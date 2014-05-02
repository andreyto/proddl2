#include "PRODDL/Common/options.hpp"
#include "PRODDL/Common/debug.hpp"
#include "PRODDL/External/jsoncpp/jsoncpp.hpp"
#include <fstream>

namespace PRODDL {

	namespace implement {
		void load_json_container_to_options(const Json::Value& node_cont,Options& opt) {
			ATALWAYS(node_cont.type() == Json::objectValue,"Option containers should be JSON objects");
			std::vector<std::string> keys = node_cont.getMemberNames();

			for (int i = 0; i < node_cont.size(); i++)
			{
				const std::string key = keys[i];
				const Json::Value val = node_cont[key];
				// switch type of value
				switch (val.type())
				{
				case Json::nullValue:
					//TODO: Maybe instead delete existing keys if present, then the user can
					//load multiple files and have an option of deleting previously defined
					//entries by providing keys with "null" value - similar to how this
					//is done in BASH
					ATDIE(("Keys with undefined values are not allowed: "+key).c_str());
					break;
				case Json::arrayValue:
					ATDIE(("Keys with array values are not supported: "+key).c_str());
					break;
				case Json::objectValue:
					//recursive call
					load_json_container_to_options(val,opt.set(key,Options()));
					break;

				case Json::stringValue:
					opt.set(key,val.asString());
					break;
				case Json::intValue:
					opt.set(key,val.asInt());
					break;
				case Json::uintValue:
					opt.set(key,val.asUInt());
					break;
				case Json::realValue:
					opt.set(key,val.asDouble());
					break;
				case Json::booleanValue:
					opt.set(key,val.asBool());
					break;
				default:
					ATDIE(("Unknown value type for key: "+key).c_str());
				}
			}

		}
	} // namespace implement
	
	void load_options_from_json_file(std::string inp_file,Options& opt)
	{
		Json::Reader reader;
		Json::Value node_cont;
		std::ifstream inp_stream(inp_file.c_str());
		ATALWAYS(reader.parse(inp_stream, node_cont, false),("Failed to load config file: " + inp_file).c_str());
		inp_stream.close();
		implement::load_json_container_to_options(node_cont,opt);
	}
	
	void load_options_from_json_string(std::string inp_str,Options& opt)
	{
		Json::Reader reader;
		Json::Value node_cont;
		ATALWAYS(reader.parse(inp_str, node_cont, false),("Failed to parse config string: " + inp_str).c_str());
		implement::load_json_container_to_options(node_cont,opt);
	}

} // namespace PRODDL