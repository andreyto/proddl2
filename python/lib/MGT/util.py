from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MDDV package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
import os, tarfile, json, glob, shutil
import sys
import contextlib
import re
import itertools
import importlib
import fnmatch
import subprocess
import tempfile
from multiprocessing.pool import ThreadPool
import six
import hashlib

def curr_time_as_filename():
    """Return current datetime formatted as suitable for use as a file name"""
    import datetime
    return datetime.datetime.today().strftime("%y-%m-%d_%H-%M")


def make_tmp_file(*l,**kw):
    """Create and open a temporary file that will exist after this program closes it.
    @return a tuple (file object,file name).
    It does the same as tempfile.NamedTemporaryFile but the file is not automatically
    deleted after being closed. Because it works through calls to mkstemp and os.fdopen,
    the returned file object does not have a file name in its 'name' attribute.
    @param prefix If provided, the file will start with that prefix (inside the dir directory)
    @param suffix If provided, the file will have that suffix
    @param create_parents - if True (default) - create parent directories (require 'dir' option)
    @param dir - create file in this directory
    @param mode (default 'w') - open file in this mode
    @param bufsize - open with this buffer size"""
    opts1 = {}
    opts1.update(kw)
    opts1.setdefault("create_parents",True)
    if opts1.pop("create_parents"):
        try:
            dirName = opts1["dir"]
        except KeyError:
            pass
        else:
            makedir(dirName)
    l2 = []
    opts1.setdefault("mode","w")
    for k in ("mode","bufsize"):
        if k in opts1:
            l2.append(opts1[k])
            del opts1[k]
    if opts1.pop("with_time",False):
        opts1["prefix"] = opts1.get("prefix","") + currTimeAsFileName()+'.'
    (fd,name) = tempfile.mkstemp(*l,**opts1)
    return (os.fdopen(fd,*l2),name)

def make_work_file(path_base,location="original",dir=None,return_path_only=True):
    """Return a temporary file that is path_base plus some unique suffix.
    @param location [original|cwd] If original, in the same dir as original,
    if cwd, then in the current working directory
    @param dir overrides location if not None
    @param return_path_only If True, close the file and return string with name,
    otherwise return a tuple (file object, file name)
    """
    import os
    dirName,baseName = os.path.split(path_base)
    if not dirName.strip():
        dirName = os.getcwd()
    if not dir is None:
        dirName = dir
    else:
        if location == "cwd":
            dirName = os.getcwd()
    fobj,path = make_tmp_file(suffix='.tmp', prefix=baseName+"_", dir=dirName)
    if return_path_only:
        fobj.close()
        return path
    else:
        return (fobj,path)

def cat_files(inp_files,out_file):
    import shutil
    with open(out_file,'wb') as out:
        for inp_file in inp_files:
            with open(inp_file,'rb') as inp:
                shutil.copyfileobj(inp,out)

def make_paired_file_names(root="read",ext="fastq",n=2):
    if n<2:
        return "{}.{}".format(root,ext)
    else:
        return [ "{}_{}.{}".format(root,i,ext) for i in range(n) ]

class dir_helper:

    def __init__(self,dirname):
        self.dirname = dirname

    def make_paired_files(self,base="read",ext="fastq",n=2):
        return make_paired_file_names(root=os.path.join(self.dirname,base),ext=ext,n=n)

    def make_file(self,base="read",ext="fastq"):
        return make_paired_file_names(root=os.path.join(self.dirname,base),ext=ext,n=1)

    def make_paired_files_func(self,base="read",ext="fastq",n=2):
        return lambda sfx,n=n: self.make_paired_files(base="{}_{}".format(base,sfx),ext=ext,n=n)

    def make_file_func(self,base="read",ext="fastq"):
        return self.make_paired_files_func(base=base,ext=ext,n=1)

@contextlib.contextmanager
def chdir(dirname=None,create=False,exist_ok=True,to_abs=True):
    curdir = os.getcwd()
    try:
        do_chdir = False
        if dirname is not None:
            if create:
                os.makedirs(dirname,exist_ok=exist_ok)
            if to_abs:
                dirname = abspath(dirname)
            do_chdir = True
        else:
            dirname = curdir
        if do_chdir:
            os.chdir(dirname)
        yield dir_helper(dirname if to_abs else ".")
    finally:
        os.chdir(curdir)


def abspath(path):
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    return path

def is_rel_subdir(path):
    """Check that path can be used as a simple subdirectory.
    It should not be an absolute path, and should not have any .. or . components."""
    return (not (".." in path.split("/") or os.path.isabs(path))) and path==os.path.relpath(path)

def absrealpath(path):
    return abspath(os.path.realpath(path))

def is_real_subdir(lower,upper,_lower_prepared=False,_upper_prepared=False):
    """Check that path `lower` is really a subdirectory of path `upper`.
    This resolves symlinks. The main application is to prevent escapes
    through symlinks when `lower` is supposed to be restricted to be under `upper`.
    """
    if not _lower_prepared:
        lower = absrealpath(lower)
    if not _upper_prepared:
        upper = absrealpath(upper)
    return os.path.commonpath([lower,upper]) == upper

def is_real_subdir_any(lower,uppers,_upper_prepared=False):
    """Iterate is_real_subdir over multiple upper dirs"""
    lower = absrealpath(lower)
    for upper in uppers:
        if is_real_subdir(lower,upper,
            _lower_prepared=True,
            _upper_prepared=_upper_prepared):
            return True
    return False

def makedir(path):
    """Make a directory if it does not already exist, creating intermediate subdirectories."""
    if not os.path.isdir(path):
        os.makedirs(path)

def make_pardir(path):
    """Make a parent directory of a file if the parent does not already exist, creating intermediate subdirectories."""    
    pardir = os.path.dirname(path)
    if pardir:
        makedir(pardir)

def split_path_all(path):
    """Split path into all components.
    Note that for absolute path like /path it will return ['/','path'], just like os.path.split()"""
    #Copied from https://www.safaribooksonline.com/library/view/python-cookbook/0596001673/ch04s16.html
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts

_path_unix_friendly_chars_strictest = "-_.0-9a-zA-Z@"
_rx_not_path_unix_friendly_chars_strictest = re.compile("[^{}]".\
                                                    format(_path_unix_friendly_chars_strictest))

def is_path_unix_friendly(path):
    for i_part,part in enumerate(split_path_all(path)):
        if (not (i_part==0 and part=="/")) and re.search(_rx_not_path_unix_friendly_chars_strictest,part):
            return False
    return True

def path_unix_friendly_chars():
    return _path_unix_friendly_chars_strictest

def assert_path_unix_friendly(path):
    if not is_path_unix_friendly(path):
        raise ValueError("File path '{}' contains characters that can cause problems with various "
                     "tools operating from a Linux shell. Please rename all such files and subdirectories "
                     "to make sure that they only contain the following symbols: '{}' "
                     "(dash between symbols means all symbols in that range).". \
                     format(path, path_unix_friendly_chars()))

def assert_id_friendly(x,strictness="file_system_alpha"):
    assert strictness in ("file_system_alpha"), "Unimplemented strictness value: {}".format(strictness)
    if strictness == "file_system_alpha":
        if not (x[:1].isalpha() and is_path_unix_friendly(x) and not '@' in x):
            raise ValueError("ID string '{}' should start with a letter and contain only symbols in the range: '{}'".\
                             format(x, path_unix_friendly_chars()))

## from: https://gist.github.com/simonw/229186
def is_hard_link(src, dst):
    """Check that two files are hard links to the same thing.
    Hardlinks are symmetrical, so the order of dst and str does not
    matter.
    Does it work on Windows?"""
    s1 = os.stat(src)
    s2 = os.stat(dst)
    return (s1[stat.ST_INO], s1[stat.ST_DEV]) == \
            (s2[stat.ST_INO], s2[stat.ST_DEV])

def is_sym_link(src,dst):
    """Check that dst is a symlink to src.
    On Windows, always returns False"""
    ret = False
    if os.islink(dst):
        dst_targ = os.readlink(dst)
        ## now we need to compare that dst_targ and src are the same
        ## file. If we used normpath or samefile, and src was itself
        ## a symlink, we would be comparing for the final src,
        ## but this is not what we need here: we need to see if
        ## dst_targ points to src itself, not to its destination.
        ## Therefore, we use the line below.
        are_same = os.samestat(os.lstat(dst_targ),os.lstat(src))
        if are_same:
            ret = True
    return ret

def link_or_copy(src,dst,order="hsc",existing="replace"):
    """Try to create hard link, symlink or copy.
    @param order Defines the methods to try and their order
    (h - hard link, s - soft link, c - copy).
    @param existing [replace,keep,raise] What to do if the
    target already exists. Note: "replace" operation is currently
    not atomic - in case of error conditions, you can end up
    with existing target deleted and not replaced; "keep"
    operation does not check if the existing target has the same
    type as the first method in the `order` string, it just
    return the first method (it will keep the existing symlink and
    @return The symbol for the method that succeeded
    """
    actions = dict(
            h = dict(copy=lambda: os.link(src,dst),
                keep=lambda: not os.path.islink(dst) and is_hard_link(src,dst)),
            s = dict(copy=lambda: os.symlink(src,dst),
                keep=lambda: is_sym_link(src,dst)),
            c = dict(copy=lambda: shutil.copy(src,dst),
                keep=lambda: not os.path.islink(dst) and not is_hard_link(src,dst))
            )
    for o in order:
        if os.path.exists(dst):
            if existing == "replace":
                os.unlink(dst)
            elif existing == "keep":
                if actions[o]["keep"]():
                    return o
            else:
                raise OSError("Destination already exists: {}".format(dst))
            
        try:
            actions[o]["copy"]()
            return o
        except:
            pass
    raise OSError("Unable to link or copy '{}' into "\
            "'{}' with method order '{}'".format(src,dst,order))


def rmtree_one(path,ignore_errors=False,onerror=None,file_ok=True):
    """Attempts a fast removal of a potetially very large directory tree.

    Tries using several external programs (pretreat the tree) until success,
    and finish with shutil.rmtree. The extra paremeters are
    ignored in pretreatment steps,
    which always ignore any errors.
    """
    #print("DEBUG: removing {}".format(path))
    def _try_call_remove(ret,path,*l,**kw):
        ## non-existing executables will always raise OSError - wrap
        ## all calls in try-except:pass to avoid doing extra shutil.which()
        ##TODO: error messages still show up in standard error and look scary.
        ##Consider doing shutil.which once and storing result in a module-level
        ##variable.
        if ret !=0:
            if os.path.exists(path):
                try:
                    ret = subprocess.call(*l,**kw)
                except:
                    pass
        return ret
    if os.path.isdir(path):
        ret = 1
        try:
            with tempfile.TemporaryDirectory(suffix=".dir",\
                prefix="empty.") as empty_dir:
                ret = _try_call_remove(ret,path,["rsync","-rd","--force","--delete",empty_dir+"/",path])
        except:
            pass            
        ret = _try_call_remove(ret,path,["find",path,"-type","f","-delete"])
        ret = _try_call_remove(ret,path,["rm","-rf",path])
        if os.path.exists(path):
            shutil.rmtree(path,ignore_errors=ignore_errors,onerror=onerror)
    else:
        if not os.path.exists(path):
            if not ignore_errors:
                if not onerror:
                    raise OSError("No such file or directory: {}".format(path))
                else:
                    onerror(path)
        elif file_ok:
            os.unlink(path)
        else:
            raise NotADirectoryError(path)
    return path

def rmtree(path,ignore_errors=False,onerror=None,threads=1,thread_targets="*",file_ok=False):
    if threads > 1:
        it_path = glob_files(os.path.join(path,thread_targets) if thread_targets else path,return_sorted=False)
        with ThreadPool(threads) as pool:
            for x in pool.imap_unordered(lambda p: rmtree_one(p,ignore_errors=ignore_errors,
                onerror=onerror,file_ok=True),it_path):
                pass
    for p in glob_files(path,return_sorted=False):
        rmtree_one(p,ignore_errors=ignore_errors,onerror=onerror,file_ok=file_ok)

def urljoin_path(base,url):
    import urllib.parse
    #urlparse.urljoin is weird: 
    #In [6]: urljoin(urljoin("/","static/MDDV"),"jbrowse")
    #Out[6]: '/static/jbrowse'
    if not base.endswith(":"):
        if not base.endswith("/"):
            base += "/"
    return urllib.parse.urljoin(base,url)

def to_url_params(params):
    """You might have to pass OrderedDict if the order of parameters
    is important. Alternatively, urllib.quote_plus can be applied
    directly to a string."""
    import urllib.request, urllib.parse, urllib.error
    return urllib.parse.urlencode(params)

def add_to_path(dir,var="PATH",prepend=False,env=None):
    """Add a directory to the PATH environment variable"""
    dir = str(dir)
    if env is None:
        env = os.environ
    if var in env:
        if prepend:
            first = dir
            second = env[var]
        else:
            first = env[var]
            second = dir
        env[var] = os.pathsep.join((first,second))
    else:
        env[var] = dir

def tar_check_safety(tar,file_name=None):
    
    def _tar_info_str(tarinfo):
        return " ; ".join((str(_) for _ in [tarinfo.name,tarinfo.type]))
    
    def _err_msg(tarinfo,msg):
        return "Archive {} failed safety check - {} {}".format(file_name if file_name else "<Name unavailable>",msg,_tar_info_str(tarinfo))

    for tarinfo in tar:
        if os.path.isabs(tarinfo.name) or \
                os.path.isabs(os.path.normpath(tarinfo.name)):
            raise ValueError(_err_msg(tarinfo,
            "Absolute file name detected"))
        elif ".." in split_path_all(tarinfo.name) or ".." \
                in split_path_all(os.path.normpath(tarinfo.name)):
            raise ValueError(_err_msg(tarinfo,
                    "Upper directory reference is detected"))
        elif not (tarinfo.isreg() or tarinfo.isdir()):
            #e.g. if archive was artificially manipulated to contain 
            #first A/B where B is a symlink to ../../something,
            #and then A/B/C, then C might be created as ../../something/C 
            #(my guess).
            raise ValueError(_err_msg(tarinfo,
                    "Non-regular files or dirs can lead to exploits"))

def tar_extractall_safe(archive,path=None,check_path_empty=False,create_path=True,skip_safety_check=False):
    if path is None:
        path = os.getcwd()
    elif create_path:
        os.makedirs(path,exist_ok=True)
    if check_path_empty:
        assert not any(True for _ in os.listdir(path)), "Empty directory is expected"
    tar = tarfile.open(archive, "r") #will auto-detect compression
    try:
        if not skip_safety_check:
            tar_check_safety(tar,file_name=archive)
        tar.extractall(path=path)
    finally:
        tar.close()
    return path


def tar_extractall_safe_single_dir(archive,path=None,check_path_empty=False,create_path=True):
    tar_extractall_safe(archive=archive,path=path,check_path_empty=check_path_empty,create_path=create_path)
    subdirs = list(os.listdir(path))
    assert len(subdirs) == 1,\
            "Expected a single directory in archive %s" \
            % (path,)
    return os.path.join(path,subdirs[0])

def dir_to_tar(dir,archive,with_dir=False,format="gztar"):
    """
    Archive everything under a directory into a file.

    :param dir: Directory to archive
    :param archive: Name of output archive (full name, with the extension)
    :param with_dir: if True, archive hierarchy will start with dir name itself, otherwise
    archive will be a "tar bomb" - not containing the directory name itself
    :return: value of the input archive parameter
    """
    import shutil
    tmp_file = make_work_file(archive)
    base_dir = "."
    if with_dir:
        ##abspath needed so that split would return uniformly structured results
        ##in various corner cases like ".", "./" and "../"
        dir,base_dir = os.path.split(os.path.abspath(dir))
        if base_dir == "":
            ## only happens when dir=="/"
            base_dir = "."
    arch = shutil.make_archive(base_name=tmp_file, format=format,
                        root_dir=dir,base_dir=base_dir)
    os.rename(arch,archive)
    os.remove(tmp_file)

def tar_to_dir(archive,dir,skip_safety_check=False,overwrite=False):
    """Expand targz file into a directory.

    Safety checks are done - links and special files will cause an error.
    Directory will be created or can exist but be empty unless overwrite is set"""
    return tar_extractall_safe(archive=archive,
                        path=dir,
                        check_path_empty=not overwrite,
                        create_path=True,
                        skip_safety_check=skip_safety_check)

def none_from_str(s):
    if s is not None:
        if s == "None":
            return None
    return s


def re_match_any(patterns,s,flags=0,format="re"):
    """Look for the first match from a list of regexes or shell globs.

    Args:
        format (str): One of
          - re: pattern in Python regex
          - glob: pattern is a shell glob as in Python fnmatch module
    """
    assert format in ("re","glob")
    if is_string(patterns):
        patterns = (patterns,)
    for pattern in patterns:
        if format == "glob":
            pattern = fnmatch.translate(pattern)
        m = re.match(pattern,s,flags=flags)
        if m:
            return m
    return None



def is_callable(x):
    """Test that object is a function or other callable"""
    ## callable() is not available in Python 3.0-3.1
    import sys
    if sys.version_info[0] == 3 and sys.version_info[1] < 2:
        return hasattr(x, '__call__')
    else:
        return callable(x)

def is_string(x):
    import six
    return isinstance(x,six.string_types)

def is_sequence(x):
    import collections
    if is_string(x):
        return False
    return isinstance(x, collections.Sequence)

def is_mapping(x):
    import collections
    return isinstance(x, collections.Mapping)

def is_bool(x):
    return isinstance(x,bool)

def rreplace(s, old, new, count=-1):
    parts = s.rsplit(old, count)
    return new.join(parts)

def make_executable(file_name):
    """Set executable permission bit"""
    os.chmod(file_name,os.stat(file_name).st_mode|0o755)

def make_id_from_file_name(file_name,strip_sfx="."):
    base = os.path.basename(file_name)
    if strip_sfx is not None:
        base = base.rsplit(strip_sfx,1)[0]
    return base

def make_nested_dir(top,seed,levels=2):
    """
    Create a hierarchy of directories under a given top directory.

    This exists to avoid placing too many temporary directories under a single top
    in a flat structure, which can slow down metadata updates such as deletes on the 
    local file system.

    The seed parameter allows for deterministic placement of the created directory.
    The seed is hashed into hex digest and the directory structure is created from
    the initial letters of the digest.

    :param top : string, top directory for the hierarchy
    :param seed : string, the hierarchy will be generated from this seed string
    :rtype : string, path to temporary directory - will be created when necessary.
    """
    dir_source = hashlib.md5(six.b(str(seed))).hexdigest()
    dir_chain = top
    for i in range(max(min(levels,len(dir_source)),1)):
        dir_chain = os.path.join(dir_chain, dir_source[i])
        if not os.path.exists(dir_chain):
            try:
                os.mkdir(dir_chain)
            except os.error:
                if not os.path.exists(dir_chain):
                    raise
    return dir_chain

def sorted_conditional(iterable,do_sort=True):
    if do_sort:
        return sorted(iterable)
    else:
        return iterable

def glob_files(
        files_globs=None,
        files_globs_list=None,
        file_globs_sep=",",
        recurse=False,
        return_split_path=False,
        return_sorted=True,
        entry_type="fd",
        filt=None,
        filt_re_format="glob",
        filt_re_flags=0,
        filt_re_apply="basename",
        no_glob_match="ignore"
        ):
    """Iterate through both a list of globs and a list of globs from a file.

    Overall behavior of this function mimics Linux 'find files_glob -type entry_type <other predicates>'
    where 'filt' parameter is responsible for implementing the <other predicates>, and 'recurse'
    parameter restricts depth of search to zero when set to False.

    @param files_globs Either a glob string or list of glob strings
    @param files_globs_list Either a file name where each line is a glob 
    or a list of globs
    @param recurse If any matching glob is a directory, recurse into the directory tree and walk through the subtrees
    @param return_split_path For each found entry, return a tuple (resolved_glob_path,recursed_entry_under_glob_path). The
    full entry path can be reconstructed with os.path.join(). For uniformity, for non-recursed entries the returned
    value will be ("",resolved_glob_path), so that os.path.join() would have the same effect as for the recursed entries.
    @param return_sorted Return matched globs and recursed entries in lexicographically sorted order. If multiple
    glob patterns were provided, their original order will be always kept as-is - sorting will only be done within the
    pattern. Note that sorting internally buffers output in a list, so the memory requirements will grow with large
    directory trees.
    @param entry_type Any of {"f","df","d"}. If recursing, yield only these types of entries (files, dirs or both)
    @param filt If defined, shuld be one of:
        - a regular expression or shell glob (controlled by filt_re_format) that will be matched against each found path
        - a function that will be called with the found path.
    The path will be only yielded if the filt
    returns True. This is applied both to the top glob paths and to all elements found during recursion.
    @param filt_re_format {"glob","re"} How to interpret the string passed as a filter. Glob will be converted
    to a regex with fnmatch.translate, and passed to re.match.
    @param filt_re_flags flag argument to re.match when a string pattern filter is applied
    @param filt_re_apply {"basename","dirname","path"} To what component of the found path a string pattern filter
    should be applied. Note that if filt is a function, the function is always called of the full path.
    @param no_glob_match {"ignore","raise","na"} What to do if, for a given glob, no
    match that passes the filter is yielded:
        - ignore - do nothing, continue to the next glob
        - raise - raise ValueError exception
        - na - return a single record for this glob, None or (None,None) depending on
        the value of return_split_path.
    Use this option when you want to track which globs from a list of globs found
    anything, and which - did not. In particular, it is useful when you need to search for
    a fixed list of file names.
    """
    def _walk(inp_match, recurse=True, return_split_path=False,entry_type="fd", filt=None):
        if recurse and os.path.isdir(inp_match):
            for root, dirs, files in os.walk(inp_match):
                if "f" in entry_type:
                    for name in files:
                        x = os.path.join(root, name)
                        if not filt or filt(x):
                            if return_split_path:
                                yield (inp_match,os.path.relpath(x,inp_match))
                            else:
                                yield x
                if "d" in entry_type:
                    for name in dirs:
                        x = os.path.join(root, name)
                        if not filt or filt(x):
                            if return_split_path:
                                yield (inp_match, os.path.relpath(x, inp_match))
                            else:
                                yield x
        else:
            if ("f" in entry_type and os.path.isfile(inp_match)) or \
                    ("d" in entry_type and os.path.isdir(inp_match)):
                if not filt or filt(inp_match):
                    if return_split_path:
                        yield ("",inp_match)
                    else:
                        yield inp_match
    
    def _handle_no_match(any_matched,inp_match,return_split_path=False):
        if any_matched:
            return None
        if no_glob_match == "raise":
            raise ValueError("No matching path found for glob: {}".format(inp_match))
        elif no_glob_match == "na":
            ret = None
            if return_split_path:
                ret = (inp_match,ret)
            return [ret]
        elif no_glob_match == "ignore":
            return None
        else:
            raise ValueError("Unknown no_glob_match value: {}".format(no_glob_match))


    if files_globs is None:
        files_globs = []
    if is_string(files_globs):
        if not file_globs_sep:
            files_globs = [ files_globs ]
        else:
            files_globs = files_globs.split(file_globs_sep)
    if filt is not None:
        if is_string(filt):
            if filt_re_apply == "basename":
                filt_spl = os.path.basename
            elif filt_re_apply == "dirname":
                filt_spl = os.path.dirname
            elif filt_re_apply == "path":
                filt_spl = lambda path: path
            else:
                raise ValueError("Unknown value for filt_re_apply: ".format(filt_re_apply))
            patt = filt
            filt = lambda x: re_match_any(patterns=patt,s=filt_spl(x),flags=filt_re_flags,format=filt_re_format)
    for files_glob in files_globs:
        any_matched = False 
        for f in sorted_conditional(glob.iglob(files_glob),return_sorted):
            for el in sorted_conditional(_walk(f, recurse=recurse, return_split_path=return_split_path,
                            entry_type=entry_type, filt=filt),return_sorted):
                any_matched = True
                yield el
        nm = _handle_no_match(any_matched,files_glob)
        if nm:
            yield nm[0]

    if files_globs_list is not None:
        if is_string(files_globs_list):
            with open(files_globs_list,"r") as inp_globs:
                for files_glob in inp_globs:
                    any_matched = False
                    for f in sorted_conditional(glob.iglob(files_glob),return_sorted):
                        for el in sorted_conditional(_walk(f, recurse=recurse, return_split_path=return_split_path,
                                        entry_type=entry_type, filt=filt),return_sorted):
                            any_matched = True
                            yield el
                    nm = _handle_no_match(any_matched,files_glob)
                    if nm:
                        yield nm[0]
        else:
            for files_glob in files_globs_list:
                any_matched = False
                for f in sorted_conditional(glob.iglob(files_glob),return_sorted):
                    for el in sorted_conditional(_walk(f, recurse=recurse, return_split_path=return_split_path,
                                    entry_type=entry_type, filt=filt),return_sorted):
                        any_matched = True
                        yield el
                nm = _handle_no_match(any_matched,files_glob)
                if nm:
                    yield nm[0]


def load_json(file_name):
    with open_text_py23(file_name,'r') as f:
        return json.load(f)

def save_json(o,file_name,pretty=False,**kw):
    if pretty:
        kw.update(dict(sort_keys=True,
                       indent=4,
                       separators=(',', ': ')))
    with open_text_py23(file_name,'w') as f:
        json.dump(o,f,**kw)

def get_arg_as_json(x):
    if is_string(x):
        x = json.loads(x)
    return x

def file_str_format(inp_file=None, inp_str=None, out_file=None, vars={}):
    """Read file as string (or use input str), 
    apply str.format() and save (or return as str)"""

    vars = vars.copy()

    s = ""

    if inp_file is not None:
        with open(inp_file,"r") as inp:
            s += inp.read()
    if inp_str is not None:
        s += inp_str

    s = s.format(**vars)
    
    if out_file is not None:
        with open(out_file,"w") as out:
            out.write(s)
    
    return s

def encode_urlsafe(s):
    import base64
    return base64.urlsafe_b64encode(s)

def list_find(x,key,start=None,end=None,step=1):
    """Find key in list x and return the index.
    :param x: list or other sequence
    :param key: if a callable, should return True on matching elements; otherwise is compared with == operator
    :param start: start searching from this index
    :param end: search up to this index (non-inclusive)
    :param step: iterate with this step during search (must be positive or None)
    :return: Index of the first found element or a negative integer if key is not found

    start argument can be used to search repeatedly for more occurrences of the key
    """
    not_found = -1
    if start is None:
        start = 0
    ## end==None is iterpreted as len() by islice()
    ## expecting here that islice is internally optimized for random-access containers like a list,
    ## but never having tested if this is true
    en_iter = zip(itertools.count(start,step),itertools.islice(x,start,end,step))
    if not is_callable(key):
        pos = next((ind for (ind,val) in en_iter if val == key), not_found)
    else:
        pos = next((ind for (ind,val) in en_iter if key(val)), not_found)
    return pos

from itertools import groupby # for unique function.
from bisect import bisect_left, insort_left # for unique function.

def unique(seq, stable=False):
    """unique(seq, stable=False): return a list of the elements in seq in arbitrary
    order, but without duplicates.
    If stable=True it keeps the original element order (using slower algorithms).
    Unlike set(), it will also work for unhashable and even unsortable
    elements"""
    # Developed from Tim Peters version:
    #   http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52560


    # Special case of an empty s:
    if len(seq)==0: return list()

    # if it's a set:
    if isinstance(seq, set): return list(seq)

    if stable:
        # Try with a set:
        seqSet= set()
        result = []
        try:
            for e in seq:
                if e not in seqSet:
                    result.append(e)
                    seqSet.add(e)
        except TypeError:
            pass # move on to the next method
        else:
            #if uniqueDebug: print "Stable, set."
            return result

        # Since you can't hash all elements, use a bisection on sorted elements
        result = []
        sortedElem = []
        try:
            for elem in seq:
                pos = bisect_left(sortedElem, elem)
                if pos >= len(sortedElem) or sortedElem[pos] != elem:
                    insort_left(sortedElem, elem)
                    result.append(elem)
        except TypeError:
            pass  # Move on to the next method
        else:
            #if uniqueDebug: print "Stable, bisect."
            return result
    else: # Not stable
        # Try using a set first, because it's the fastest and it usually works
        try:
            u = set(seq)
        except TypeError:
            pass # move on to the next method
        else:
            #if uniqueDebug: print "Unstable, set."
            return list(u)

        # Elements can't be hashed, so bring equal items together with a sort and
        # remove them out in a single pass.
        try:
            t = sorted(seq)
        except TypeError:
            pass  # Move on to the next method
        else:
            #if uniqueDebug: print "Unstable, sorted."
            return [elem for elem,group in groupby(t)]

    # Brute force:
    result = []
    for elem in seq:
        if elem not in result:
            result.append(elem)
    #if uniqueDebug: print "Brute force (" + ("Unstable","Stable")[stable] + ")."
    return result


# Following a suggestion from Alex Martelli, sometimes this uniquePick
# is faster for more than about 300 unsortable and unhashable elements:

from pickle import dumps

def unique_pick(seq):
    result = []
    seen = set()
    for elem in seq:
        key = dumps(elem, protocol=-1)
        if key not in seen:
             seen.add(key)
             result.append(elem)
    return result

def gethostname():
    import socket
    return socket.getfqdn()

def rm_fname_extension(x,sep="."):
    return x.rsplit(sep,1)[0]


def open_text_py23(file,mode,*l,**kw):
    """Kludge to open file in binary mode for writing on Python 2 to get modules like json and csv work in Python 2 after futurize.
    Use only when strictly necessary! It will most likely result in error if you have Unicode strings on Py 2"""
    import sys
    if sys.version_info[0] < 3:
        if 'b' not in mode:
            mode = mode + 'b'
    return open(file,mode,*l,**kw)

def open_gzip_text23(file,mode,format="auto",*l,**kw):
    """Kludge to open gzip files in text mode for reading on Python 2"""
    import sys
    import gzip
    if sys.version_info[0] >= 3:
        if 'b' not in mode:
            if 't' not in mode:
                mode = mode + 't'
    else:
        if 't' in mode:
            mode = mode.replace('t','')
    if format=="auto":
        if file.endswith(".gz"):
            format = "gzip"
        else:
            format = "txt"
    if format=="gzip":
        return gzip.open(file,mode,*l,**kw)
    else:
        return open(file,mode,*l,**kw)

def import_name(name):
    res = name.rsplit(".",1)
    if len(res) > 1:
        mod,obj = res
    else:
        mod = ""
        obj = res[0]
    if mod:
        mod = importlib.import_module(mod)
        obj = getattr(mod,obj)
    else:
        obj = globals()[obj]
    return obj

def is_stdio_fn(x):
    """Return true if the argument designates stdio stream request (None or the '-' symbol)"""
    return (x is None) or (is_string(x) and x == "-")

@contextlib.contextmanager
def as_stream_cm(x,mode="r"):
    """Convert argument to a file stream object.
    
    If a string and not '-', open with a given mode.
    If None or '-', return stdio or stdout depending on the mode.
    Else retun as-is.
    
    This is a context manager. It will only auto-close streams that
    it opened."""

    ret = x
    do_close = False
    try:
        if is_stdio_fn(x):
            if "r" in mode:
                ret = sys.stdin
            elif "w" in mode or "a" in mode:
                ret = sys.stdout
            else:
                raise ValueError(f"Unknown file open mode: {mode}")
        elif is_string(x):
            ret = open(x,mode)
            do_close = True
        yield ret
    finally:
        if do_close:
            ret.close()



## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        dir_to_tar,
        tar_to_dir,
        is_path_unix_friendly,
        split_path_all,
        glob_files,
        rmtree
    ])


if __name__ == "__main__":
    _main()
