### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PRODDL package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Generate shell rc file to be sourced before installing or using 
package if non-standard install dir is used"""
import sys, os

pkg_name = "PRODDL"
conf_ext = ".json"

def is_in_path(target,path):
    if isinstance(path,str):
        path = path.split(os.pathsep)
    in_path = False
    for el in path:
        if os.path.abspath(target) == os.path.abspath(el) or \
                (os.path.exists(el) and \
                 os.path.exists(target) and \
                 os.path.samefile(target,el)):
            in_path = True
            break
    return in_path

def add_to_path(target,path=None,prepend=True):
    val = target
    if path:
        if prepend:
            val = os.pathsep.join([val,path])
        else:
            val = os.pathsep.join([path,val])
    return val

if __name__ == "__main__":

    try:
        action = sys.argv[1]
    except IndexError:
        print("Usage: %s action [args]" % (os.path.basename(sys.argv[0]),))
        sys.exit(1)

    assert action in ("gen-rc",)

    if action == "gen-rc":
        schema = sys.argv[2]
        target = sys.argv[3]
        rc = sys.argv[4]

        assert schema in ("install-dir",)

        if not os.path.isabs(target):
            target = os.path.abspath(target)

        if schema == "install-dir":
            PATH = os.environ.get("PATH","")
            PYTHONPATH = os.environ.get("PYTHONPATH","")
            with open(rc,"w") as out:
                if not is_in_path(target,sys.path):
                    print("export PYTHONPATH={0}{1}$PYTHONPATH".format(target,os.pathsep), file=out)
                if not is_in_path(target,PATH):
                    print("export PATH={0}{1}$PATH".format(target,os.pathsep), file=out)
                print("export {0}_CONF={1}".format(pkg_name,
                        os.path.join(target,pkg_name+conf_ext)), file=out)

