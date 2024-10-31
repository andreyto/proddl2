"""Classes for generating Makeflow inputs"""





from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from builtins import object
from MGT.MakeflowArgParser import parse_makeflow_args, unparse_makeflow_args
from MGT.UUID import genId
from MGT.util import make_executable, is_string

from contextlib import contextmanager, closing
import os
import subprocess
import shutil
from subprocess import check_call

def _makeflow_task_update_targets(targets,cwd=None):
    if cwd:
        targets = [ os.path.join(cwd,x) if not os.path.isabs(x) else x for x in targets ]
    return targets

class MakeflowWriter(object):

    _tab = " "*4
    
    def __init__(self,out,vars=None,exports=None,wrapper=None,env=None,mode="w"):
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
        self.wrapper = wrapper
        self.env = env
    
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
    
    def task(self,cmd=None,targets=[],inputs=[],
            vars=None,
            is_local=False,
            cwd=None,
            status_file=None,
            stdout=None,
            stderr=None,
            id_task=None,
            wrapper=None,
            env=None,
            cwd_remove_before=None,
            conf=None):
        """
        @cmd command to execute; if cmd is None, conf must be a dict that has 
        a key cmd
        @param cwd if not None, use cd to change to this directory before 
        running command.
        Note that this does not affect cwd when the job is submitted.
        @param status_file if not None, try to remove this file before
        running command and write '1' to it after finishing the command
        with zero exit code. File manipulation will be done from a directory
        that is current before changing to cwd, so if you want the file to
        be created under cwd, pass status_file=cwd/some_name.
        @note \\@BATCH_LOCAL=0 in vars has precedence over is_local=True
        any None values in targets or inputs will be ignored"""

        if conf is None:
            conf = {}

        assert cmd or (conf and "cmd" in conf), "cmd must be defined"

        if not cmd:
            cmd = conf["cmd"]
        if not targets:
            targets = conf.get("targets",[])
        if not inputs:
            inputs = conf.get("inputs",[])
        if not cwd:
            cwd = conf.get("cwd",None)
        if not status_file:
            status_file = conf.get("status_file",None)
        if not stdout:
            stdout = conf.get("stdout",None)
        if not stderr:
            stderr = conf.get("stderr",None)
        if not id_task:
            id_task = conf.get("id_task",None)
        if cwd_remove_before is None:
            cwd_remove_before = conf.get("cwd_remove_before",None)

        if cwd_remove_before is None:
            cwd_remove_before = False

        #Makeflow breaks down with cryptic messages if : is present
        #in the command line (it thinks it is another rule header)
        cmd = cmd.replace(":",r"\:")
        cmd = cmd.replace("=",r"\=")
        cmd = cmd.replace("\\\n","")
        cmd = cmd.replace("\n","; ")
        
        if wrapper is None:
            wrapper = self.wrapper

        if env is None:
            env = self.env
        
        w = self.out.write

        if targets is None or is_string(targets):
            targets = [ targets ]
        if inputs is None or is_string(inputs):
            inputs = [ inputs ]
        
        targets = [ x for x in targets if x is not None ]
        targets = _makeflow_task_update_targets(targets,cwd)
        
        inputs = [ x for x in inputs if x is not None ]
        
        if status_file is not None:
            assert status_file and not os.path.isdir(status_file)
            targets.append(status_file)

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
        
        if wrapper is not None:
            wrapper = str(wrapper).strip()
            if(wrapper):
                cmd = "{} {}".format(wrapper,cmd)
        
        cmd = cmd.strip()
        
        if env is not None:
            env = str(env).strip()
            if(env):
                cmd = "source '{}'; ({})".format(env,cmd)

        if stdout is None:
            if id_task is None:
                id_task = genId()
            if cwd is not None:
                stdout = "{}/{}.task_out".format(cwd,id_task)
            else:
                stdout = "{}.task_out".format(id_task)
        
        if stderr is None:
            if id_task is None:
                id_task = genId()
            if cwd is not None:
                stderr = "{}/{}.task_err".format(cwd,id_task)
            else:
                stderr = "{}.task_err".format(id_task)
        
        redir = ""
        
        if stdout != "-":
            redir += " 1> '{}'".format(stdout)
        if stderr != "-":
            redir += " 2> '{}'".format(stderr)
        
        if redir:
            cmd = "({cmd}){redir}".format(cmd=cmd,redir=redir)

        if cwd is not None:
            cmd = "mkdir -p '{}'; pushd '{}' && ({}) && popd".format(cwd,cwd,cmd)
            if cwd_remove_before:
                cmd = "rm -rf '{}'; {}".format(cwd,cmd)
        
        if status_file is not None:
            cmd = "rm -f '{status_file}' && ({cmd}) && (echo 1 > '{status_file}')".\
                    format(status_file=status_file,cmd=cmd)
        
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

def make_task_sandbox(top_dir,id_task=None,conf={}):
    conf = conf.copy()
    if id_task is None or id_task == "":
        id_task = conf.get("id_task",None)
    if id_task is None or id_task == "":
        id_task = genId()
    cwd = os.path.join(top_dir,id_task)
    shutil.rmtree(cwd,ignore_errors=True)
    os.makedirs(cwd)
    file_root = os.path.join(cwd,id_task)
    conf.setdefault("status_file",file_root+".workflow_status")
    conf.setdefault("stdout",file_root+".task_out")
    conf.setdefault("stderr",file_root+".task_err")
    conf.setdefault("cwd",cwd)
    conf.setdefault("id_task",id_task)
    return conf

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
        quiet=False,
        script_prolog=""
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
    @param quiet minimize chatter from script [False]
    @param script_prolog Insert this string of commands at the
    beginning of the run script (to perform cleanup, for example).
    This is inserted after the `env` component."""
    out_close = False
    if is_string(out):
        out = open(out,mode)
        out_close = True
    w = out.write
    w("#!/bin/bash\n")
    if env is not None:
        w(". {}\n".format(env))
    if script_prolog is not None:
        w("{}\n".format(script_prolog))
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

##TODO: see how to pass vars and exports to MakeflowWriter ctor
@contextmanager
def makeflow(
        workflow,
        makeflow_bin="makeflow",
        env=None,
        wrapper=None,
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
    with closing(MakeflowWriter(workflow_work,wrapper=wrapper,env=env)) as mkw:
        
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
        script_prolog = ""
        if mkf_opt.batch_type is not None:
            # In case the user copied a prior run dir,
            # the existing Makeflow's file torque.wrapper will not be updated
            # by Makeflow, and direct job status output into the previous
            # directory. Makeflow will never see completion of any job as a result.
            # We remove any such files here.
            script_prolog += "rm -f {}.wrapper".format(mkf_opt.batch_type)
        writeMakeflowRunScript(
                makeflow_bin = makeflow_bin,
                workflow = workflow,
                env = env,
                wrapper = wrapper,
                vars = mkf_vars,
                args = mkf_args_new,
                out = workflow_script,
                stdout = mkf_stdout,
                stderr = mkf_stderr,
                quiet = True if web else False,
                script_prolog = script_prolog
                )
        make_executable(workflow_script)
    os.rename(workflow_work,workflow)
    if run:
        check_call(["bash",workflow_script])

