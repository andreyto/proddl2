//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/types.hpp"
#include "PRODDL/Common/options_io_json.hpp"
#include "PRODDL/Common/argparse.hpp"
#include "PRODDL/External/Pdb++/pdb++.hpp"
#include "PRODDL/IO/rigid.hpp"

#include "PRODDL/Common/g_options.hpp"
#include "PRODDL/Common/logger.hpp"


#include <cstdlib>

#include <iostream>

//#include <gzstream.h>

#include <fstream>

#include <string>

#include <iostream>
#include <iterator>

namespace PRODDL {

    typedef PRODDL_T_NUM T_num;

    void load_pdb_atoms(
            std::string file_name,
            std::vector<PDBPP::PDB>& records,
            Types<T_num>::Points& coords
            ) {

        typedef Types<T_num> TypesT;
        typedef TypesT::Point Point;
        typedef TypesT::Points Points;

        records.clear();

        std::ifstream file_in(file_name.c_str());

        ATALWAYS( file_in.good(),"Unable to open input pdb file: " + file_name);

        PDBPP::PDB record;

        bool cont = true;

        while (cont && file_in >> record) {
            switch (record.type()) {

                case PDBPP::PDB::ATOM:
                    records.push_back(record);
                    break;

                case PDBPP::PDB::END:
                    //allow for format extensions that add stuff after the END (e.g. Rosetta)
                    cont = false;
                    break;
            }
        }

        coords.resize(records.size());

        for(int i = 0; i < records.size(); i++) {

            const PDBPP::PDB& r = records[i];

            Point& c = coords(i);

            c(0) = r.atom.xyz[0];
            c(1) = r.atom.xyz[1];
            c(2) = r.atom.xyz[2];

        }

    }

    void export_pdb_nmr(
            std::string model_inp,
            int model_ind_start,
            int model_ind_end,
            std::string pdb_inp_rec,
            std::string pdb_inp_lig,
            std::string model_out
            ) {


        typedef Types<T_num> TypesT;
        typedef TypesT::Point Point;
        typedef TypesT::Points Points;
        typedef IORigid<T_num> IORigidT; 

        std::vector<PDBPP::PDB> records_rec, records_lig;

        Points coords_rec, coords_lig;

        load_pdb_atoms(pdb_inp_rec,records_rec,coords_rec);
        load_pdb_atoms(pdb_inp_lig,records_lig,coords_lig);

        std::ofstream m_out(model_out.c_str());

        ATALWAYS( m_out.good(), "Unable to open output model file: " + model_out);

        PDBPP::PDB rec_mod(PDBPP::PDB::MODEL);
        PDBPP::PDB rec_endmdl(PDBPP::PDB::ENDMDL);
        PDBPP::PDB rec_ter(PDBPP::PDB::TER);

        PDBPP::PDB rec_out;

        IORigidT inp_tr;

        int n_mod_tot = inp_tr.readCoords(model_inp);

        if (model_ind_end > n_mod_tot) {
            model_ind_end = n_mod_tot;
        }

        IORigidT::RotTranValue tr_v;

        for(int i_mod = model_ind_start; i_mod < model_ind_end; i_mod++) {

            //@todo limit max model by field width
            rec_mod.model.num = i_mod - model_ind_start + 1;

            m_out << rec_mod;

            for(int i_rec = 0; i_rec < records_rec.size(); i_rec++) {
                m_out << records_rec[i_rec];
            }


            if(records_rec.size()) {
                int rec_size = records_rec.size();
                rec_ter.ter.serialNum = records_rec[rec_size-1].atom.serialNum;
                rec_ter.ter.residue = records_rec[rec_size-1].atom.residue;
                m_out << rec_ter;
            }

            inp_tr.getCoords(i_mod,1,&tr_v);

            Points coords_out(coords_lig.copy());

            tr_v.tran(coords_out);

            for(int i_rec = 0; i_rec < records_lig.size(); i_rec++) {

                rec_out = records_lig[i_rec]; // a copy

                const Point& c = coords_out(i_rec);

                rec_out.atom.xyz[0] = c(0);
                rec_out.atom.xyz[1] = c(1);
                rec_out.atom.xyz[2] = c(2);

                m_out << rec_out;

            }

            if(records_lig.size()) {
                int rec_size = records_lig.size();
                rec_ter.ter.serialNum = records_lig[rec_size-1].atom.serialNum;
                rec_ter.ter.residue = records_lig[rec_size-1].atom.residue;
                m_out << rec_ter;
            }

            m_out << rec_endmdl;

        }

    }

} // namespace PRODDL

namespace PRODDL {

    Options gOptions;

} // namespace PRODDL

namespace {

    namespace po = boost::program_options;


    //typedef PRODDL::Docking<T_num> dk;

    void
        setGlobalOptions(const PRODDL::Options& options) {

            int runTimeLogLevel = ATLOG_LEVEL;

            options.getdefault("logLevel",runTimeLogLevel,ATLOG_LEVEL);

            PRODDL::gOptions = options;

            PRODDL::Logger::setRunTimeLevel(runTimeLogLevel);

        }


    int parse_arguments(int ac, char* av[], po::variables_map& vm)
    {
        using namespace std;
        try {

            po::options_description desc("Options for exporting the results");
            desc.add_options()
                ("help", "produce help message")
                ("options", po::value<string>(), "file with options common for this run")
                ("model-inp", po::value<string>(), "input file with models")
                ("model-ind-start", po::value<int>(), "start index in model list")
                ("model-ind-end", po::value<int>(), "end index in model list")
                ("pdb-inp-rec", po::value<string>(), "input PDB file with receptor")
                ("pdb-inp-lig", po::value<string>(), "input PDB file with ligand")
                ("format-out", po::value<string>(), "output format for models")
                ("model-out", po::value<string>(), "output file for models")
                ;

            po::store(po::parse_command_line(ac, av, desc), vm);

            if (vm.count("help")) {
                cout << desc << "\n";
                return true;
            }

            po::notify(vm);    

        }
        catch(exception& e) {
            cerr << "error: " << e.what() << "\n";
            return false;
        }
        catch(...) {
            cerr << "Exception of unknown type when parsing arguments!\n";
            return false;
        }

        return true;
    }


    void process_arguments(const po::variables_map& vm) {
        using namespace PRODDL;
        using namespace std;

        Options opt;

        if(vm.count("options")) {
            load_options_from_json_file(vm["options"].as<string>(),opt);
        }

        setGlobalOptions(opt);

        string format_out = vm["format-out"].as<string>();

        if(format_out == "pdb_nmr") {
            export_pdb_nmr(
                    vm["model-inp"].as<string>(),
                    vm["model-ind-start"].as<int>(),
                    vm["model-ind-end"].as<int>(),
                    vm["pdb-inp-rec"].as<string>(),
                    vm["pdb-inp-lig"].as<string>(),
                    vm["model-out"].as<string>()
                    );
        }
        else {
            AT_THROW(po::invalid_option_value("Option 'format-out' has invalid value: " + format_out));
        }

    }

} // namespace

int main(int ac, char* av[]) {
    PRODDL::Logger::init();
    ATLOG_STD_EXCEPTIONS_TRY();
    po::variables_map vm;
    int parse_status = parse_arguments(ac, av, vm);
    if (!parse_status) {
        return -1;
    }
    else {
        if(vm.count("help")) {
            return 0;
        }
    }
    process_arguments(vm);
    ATLOG_STD_EXCEPTIONS_CATCH();
    return 0;
}

