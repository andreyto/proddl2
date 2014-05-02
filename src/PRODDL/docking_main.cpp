//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/Common/options_io_json.hpp"
#include "PRODDL/Common/argparse.hpp"
#include "PRODDL/docking.hpp"

#include "PRODDL/Common/g_options.hpp"
#include "PRODDL/Common/logger.hpp"



#include <iostream>
#include <iterator>

namespace PRODDL {

  Options gOptions;

} // namespace PRODDL

namespace {

namespace po = boost::program_options;

typedef PRODDL_T_NUM T_num;

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

        po::options_description desc("Options for FFT scan process");
        desc.add_options()
            ("help", "produce help message")
            ("options", po::value<string>(), "file with options common for this run")
			("molforce-params", po::value<string>(), "molforce parameters file")
			("task", po::value<string>(), "type of task to perform")
			("fft-rot-grid-start", po::value<int>(), "start index in rot-grid")
			("fft-rot-grid-end", po::value<int>(), "end index in rot-grid")
			("fft-rot-scan-res", po::value<string>(), "output file of fft scan for one task")
			("fft-rot-scan-list", po::value<string>(), "file with a list of fft rot scan tasks")
			("fft-res", po::value<string>(), "output file for entire fft scan")
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
    
	string task = vm["task"].as<string>();

	if(task == "rot-scan") {
		set_option_from_arg<int>(vm,opt,"fft-rot-grid-start",true);
		set_option_from_arg<int>(vm,opt,"fft-rot-grid-end",true);
		set_option_from_arg<string>(vm,opt,"fft-rot-scan-res",true);
	}
	else if(task == "gather") {
		set_option_from_arg<string>(vm,opt,"fft-rot-scan-list",true);
		set_option_from_arg<string>(vm,opt,"fft-res",true);
	}
	else {
		AT_THROW(po::invalid_option_value("Option 'task' has invalid value: " + task));
	}

	setGlobalOptions(opt);

	string molforce_params_file = vm["molforce-params"].as<string>();

	Docking<T_num>::MolForceParams mfp;

	Docking<T_num>::load_from_hdf5(molforce_params_file,mfp);

	if(task == "gather") {
		Docking<T_num>::Master app;
		app.init(mfp);
		app.run();
		string res_file;
		opt.get("fft_res",res_file);
		app.writeResults(res_file,'b');
	}
	else if(task == "rot-scan") {
		Docking<T_num>::Slave app;
		app.init(mfp);
		app.run();
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

