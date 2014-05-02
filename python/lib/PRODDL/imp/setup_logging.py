import IMP.base
import logging, sys, datetime

class _log_formatter(logging.Formatter):
    converter=datetime.datetime.fromtimestamp
    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created)
        if datefmt:
            s = ct.strftime(datefmt)
        else:
            t = ct.strftime("%Y-%m-%d %H:%M:%S")
            s = "%s.%03d" % (t, record.msecs)
        return s

def add_handler():
    log_h_console = logging.StreamHandler()
    log_h_console.setFormatter(_log_formatter(fmt='%(asctime)s %(message)s',datefmt='%Y-%m-%d,%H:%M:%S.%f'))
    logging.root.addHandler(log_h_console)

    #logging.basicConfig(stream=sys.stdout)
    #logging.root.setLevel(logging.DEBUG)
    logging.root.setLevel(logging.ERROR)

    # it gets awfully slow with internal checks
    #IMP.base.set_check_level(IMP.base.USAGE)
    IMP.base.set_log_level(IMP.base.SILENT)
    IMP.base.set_log_timer(True)
    return (log_h_console)

def remove_handler(log_h):
    logging.root.removeHandler(log_h)

