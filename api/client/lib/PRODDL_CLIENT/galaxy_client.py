import json, sys, os, logging, threading
from collections import OrderedDict
from bioblend.galaxy import GalaxyInstance

log = logging.getLogger(__name__)
        
def robust_rename(src,dst):
    try:
        os.rename(src,dst)
    except OSError:
        #will happen if dst exists on Windows
        #remove dst and try again - this is not
        #atomic anymore but there is no other way
        if os.path.exists(dst):
            os.remove(dst)
        os.rename(src,dst)

def wf_input_key_from_label(wf,label):
    found = None
    inputs = wf["inputs"]
    for key,val in inputs.items():
        if val["label"] == label:
            if found is None:
                found = key
            else:
                raise ValueError("Duplicate label: {}".format(label))
    if found is None:
        raise KeyError(label)
    return found

def dataset_to_file(gi,id,file_name):
    file_name_tmp = file_name+".tmp"
    #if use_default_filename is True, file_name_tmp must be a dir
    gi.datasets.download_dataset(
            id,
            file_name_tmp,
            use_default_filename=False,
            wait_for_completion=True
            )

    robust_rename(file_name_tmp,file_name)

def status_printer(msg,data=None):
    print str(msg)

class status_thread(threading.Thread):

    def __init__(self,kwargs={},**kw):

        kwargs = kwargs.copy()
        
        self.status_callback = kwargs.setdefault("status_callback",status_printer)
        
        threading.Thread.__init__(self,
                kwargs=kwargs,
                **kw)
    
    def status():
        return sef.status_callback


class proddl:

    dock_tool_id = "proddl_dock"

    def __init__(self,options):
        with open(options,"r") as _:
            self.opt = json.load(_)
        self.connect()

    def connect(self):
        
        log.info("Initiating Galaxy connection")

        opt_gal = self.opt["galaxy"]

        gi = GalaxyInstance(url=opt_gal["url"], 
                key=opt_gal["key"])

        self.gi = gi

    def run_dock(self,
            receptor_pdb,
            ligand_pdb,
            models_pdb,
            n_models=10,
            status_callback=status_printer):

        status_form = "proddl dock for output %s: {}" % (models_pdb,)

        opt_gal = self.opt["galaxy"]
        gi = self.gi
        
        history = gi.histories.create_history(name="PRODDL remote docking history")

        try:

            rec_uploaded = gi.tools.upload_file(
                receptor_pdb,
                history_id=history["id"],
                file_name="Receptor",
                dbkey="?",
                file_type="pdb",
            )
            assert len(rec_uploaded) == 1
            
            lig_uploaded = gi.tools.upload_file(
                ligand_pdb,
                history_id=history["id"],
                file_name="Ligand",
                dbkey="?",
                file_type="pdb",
            )
            assert len(lig_uploaded) == 1

            status_callback(status_form.format("inputs uploaded"))
            
            wf = gi.workflows.show_workflow(opt_gal["proddl_workflow_id"])
            data_map = { 
                    wf_input_key_from_label(wf,"Receptor") : 
                        {'id': rec_uploaded['outputs'][0]['id'], 'src': 'hda'},
                    wf_input_key_from_label(wf,"Ligand") : 
                        {'id': lig_uploaded['outputs'][0]['id'], 'src': 'hda'},
                        }

            params = {self.dock_tool_id: {'param': 'inp_n_models', 'value': str(n_models)}}
            #replacement_params = { "inp_n_models" : str(n_models) }
            
            res = gi.workflows.run_workflow(
                    opt_gal["proddl_workflow_id"], 
                    data_map,
                    params = params,
                    #replacement_params = replacement_params,
                    history_id=history["id"])
            
            outputs = res["outputs"]
            assert len(outputs) == 1, "Expecting a single output in workflow"

            status_callback(status_form.format("docking job submitted"))

            dataset_to_file(gi,outputs[0],models_pdb)

            status_callback(status_form.format("models saved"))

        finally:
            gi.histories.delete_history(history["id"])
        
        return status_callback

    def dock(self,
            receptor_pdb,
            ligand_pdb,
            models_pdb,
            n_models=10,
            wait=False
            ):
        kw = locals()
        del kw["self"]
        del kw["wait"]
        if wait:
            return self.run_dock(**kw)
        else:
            t = status_thread(target=self.run_dock,kwargs=kw)
            t.daemon = True
            t.start()
            return t

if __name__ == "__main__":
    pd = proddl("data/proddl_client.json")
    t = pd.dock(
            "data/test/2ptc_E.pdb",
            "data/test/2ptc_I.pdb",
            "models.pdb",
            n_models=10,
            wait=False
            )
    t.join()

