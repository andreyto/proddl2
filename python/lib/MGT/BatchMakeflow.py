"""Classes for generating Makeflow inputs"""

from MGT.MakeflowArgParser import parse_makeflow_args, unparse_makeflow_args
from PRODDL.util import make_executable, is_string

from contextlib import contextmanager, closing
import os
import subprocess
from subprocess import check_call


class MakeflowWriter(object):

    _tab = " "*4
    
    def __init__(self,out,vars=None,exports=None,mode="w"):
        if is_string(out):
            self.out = open(out,mode)
            #do not repeat default exports on append -
            #this breaks the Makeflow
            if "a" in mode and exports is None:
                exports = []
        elif exports is None:
            #do not add default exports when existing
            #stream is supplied. The caller should use
            #appendInitExports() when needed.
            exports = []

        self.appendInitExports(exports=exports)
        self._write_vars(vars,0)
        #store IDs of already processed MGT jobs to avoid submitting any job twice
        self.done = set()
    
    @classmethod
    def getStandardExports(klass):
        return [
                 "BATCH_OPTIONS",
                 "MAKEFLOW_BATCH_QUEUE_TYPE",
                 "MAKEFLOW_MAX_REMOTE_JOBS"
            ]

    def appendInitExports(self,exports):
        if exports is None:
            exports = self.getStandardExports()
        w = self.out.write
        for export in exports:
            w("export {}\n".format(export))

    def _write_vars(self,vars,tab_level=0):
        if vars is None:
            vars = []
        w = self.out.write
        for var in vars:
            w(self._tab*tab_level+"{}\n".format(var))
    
    def task(self,cmd,targets=[],inputs=[],
            vars=None,
            is_local=False):
        """\\@BATCH_LOCAL=0 is vars has precedence over is_local=True
        any None values in targets or inputs will be ignored"""
        w = self.out.write
        
        targets = [ x for x in targets if x is not None ]
        inputs = [ x for x in inputs if x is not None ]
        
        w(' '.join([str(t) for t in targets])+": "+\
                ' '.join([str(t) for t in inputs])+'\n')
        
        if vars is None:
            vars = []
        
        if is_local:
            for var in vars:
                if var.strip().startswith("@BATCH_LOCAL"):
                    break
            else:
                vars = vars + ["@BATCH_LOCAL=1"]
        
        self._write_vars(vars,1)
        
        w(self._tab+"{}\n\n".format(cmd))

    appendJob = task

    def sub(self,flow,targets=[],inputs=[],
            vars=None):
        """Append a task that itself is a Makeflow.
        This will mark the sub-makeflow as local unless
        it is already marked otherwise in the 'vars'"""
        if flow not in inputs:
            inputs = inputs + [flow]
        self.task(targets=targets,inputs=inputs,
                cmd="MAKEFLOW {}".format(flow),
                vars=vars,
                is_local=True)

        appendMakeflow = sub

    def appendMgtJob(self,job):
        """Append to current makeflow one MGTAXA job object.
        No recursion to append dependencies is performed.
        @param job BatchJob object, having at least these
        attributes: BatchJob(jobId,scriptName=scriptName,cwd=cwd,outputs=(flagOk,),depend=depend)
        where jobId is globally unique.
        """
        jobId = job.jobId
        if jobId not in self.done:
            inputs = []
            for dep in job.depend:
                inputs += dep.outputs
            targets = job.outputs
            cmd = "bash " + job.scriptName
            self.appendJob(targets=targets,inputs=inputs,cmd=cmd)
            self.done.add(job.jobId)

        
    def appendMgtJobs(self,jobs):
        """Append to current makeflow from a graph of MGTAXA job objects.
        This will recurse first to append dependencies.
        @param jobs a list of BatchJob objects, @see appendMgtJob for requirements.
        """
        for job in jobs:
            self.appendMgtJobs(job.depend)
            self.appendMgtJob(job)

    def close(self):
        if self.out:
            self.out.close()
            self.out = None

def writeMakeflowRunScript(
        makeflow_bin,
        workflow,
        vars,
        args,
        out,
        env=None,
        wrapper=None,
        mode="w",
        stdout="-",
        stderr="-",
        quiet=False
        ):
    """Write a shell script that will run this makeflow.
    @param makeflow_bin Makeflow executable path
    @param workflow workflow file path
    @param env file to source as shell environment
    @param vars list of environment variable assignments "VAR=VAL"
    @param args string with all arguments to makeflow executable
    @param out file path or file object for writing the script into
    @param mode to open out if out if a file path
    @param stdout redirect Makeflow standard output to that file;
    '-' (default) means no redirection; None means $workflow.out.log
    @param stderr redirect Makeflow standard error to that file;
    '-' (default) means no redirection; None means $workflow.err.log
    @param quiet minimize chatter from script [False]"""
    out_close = False
    if is_string(out):
        out = open(out,mode)
        out_close = True
    w = out.write
    w("#!/bin/bash\n")
    if env is not None:
        w(". {}\n".format(env))
    for var in vars:
        w("export {}\n".format(var))
    if stdout is None:
        stdout = workflow+".out.log"
    if stderr is None:
        stderr = workflow+".err.log"
    redir = ""
    redir_msg = ""
    if stdout != "-":
        redir += " 1> '{}'".format(stdout)
        if not quiet:
            redir_msg += "echo Makeflow standard output is redirected to '{}'\n".format(stdout)
    if stderr != "-":
        redir += " 2> '{}'".format(stderr)
        if not quiet:
            redir_msg += "echo Makeflow standard error is redirected to '{}'\n".format(stderr)
    if redir_msg:
        w(redir_msg)
    w('echo Starting execution of Makeflow "{}"\n'.format(workflow if not quiet else ""))
    if wrapper is None:
        wr_str = ''
    else:
        wr_str = '"{}" '.format(wrapper)
    w('{wr_str}"{makeflow_bin}" {args} "{workflow}"{redir}\n'.\
            format(
                makeflow_bin = makeflow_bin,
                args = args,
                workflow = workflow,
                redir=redir,
                wr_str=wr_str))
    #use of subshell below is to reset the $? value after echo without
    #calling exit from the main shell
    w("""status=$?
    if [ $status -ne 0 ]; then
        echo "Makeflow execution failed" >&2
    else
        echo "Makeflow execution finished"
    fi
    (exit $status)
    """)
    if out_close:
        out.close()

@contextmanager
def makeflow(
        makeflow_bin,
        wrapper,
        workflow,
        makeflow_args="",
        workflow_script=None,
        web=False,
        run=True
        ):
    if os.path.isfile(workflow):
        os.remove(workflow)
    if workflow_script is None:
        workflow_script = "{}.bat".format(workflow)
    #TODO: parse Makeflow options first and take log name from
    #where if present
    makeflowLog = workflow+".makeflowlog"
    if os.path.isfile(makeflowLog):
        os.remove(makeflowLog)
    workflow_work = workflow+".tmp"
    with closing(MakeflowWriter(workflow_work)) as mkw:
        
        ## the calling code writes Makeflow tasks
        yield mkw
        
        mkf_opt, mkf_other = parse_makeflow_args(
                makeflow_args
                )
        mkf_vars= [
            'MAKEFLOW_BATCH_QUEUE_TYPE="{}"'.format(mkf_opt.batch_type),
            'MAKEFLOW_MAX_REMOTE_JOBS={}'.format(mkf_opt.max_remote)
        ]
        if mkf_opt.batch_options is not None:
            mkf_vars.insert(0,
            'BATCH_OPTIONS="{}"'.format(mkf_opt.batch_options)
            )
        mkf_args_new = " ".join( [ '"{}"'.format(arg) for arg \
                in unparse_makeflow_args(mkf_opt,mkf_other) ] )
        #In Web mode, divert Makeflow output to files - it is
        #bulky and useless for the Web user, plus, Makeflow prints
        #"nothing else to do" to stderr on success and spooks Galaxy.
        if web:
            mkf_stdout = None
            mkf_stderr = None
        else:
            mkf_stdout = "-"
            mkf_stderr = "-"
        writeMakeflowRunScript(
                makeflow_bin = makeflow_bin,
                workflow = workflow,
                vars = mkf_vars,
                args = mkf_args_new,
                out = workflow_script,
                stdout = mkf_stdout,
                stderr = mkf_stderr,
                quiet = True if web else False
                )
        make_executable(workflow_script)
    os.rename(workflow_work,workflow)
    if run:
        check_call(["bash",workflow_script])

