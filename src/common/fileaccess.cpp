#include "fileaccess.hpp"

using std::vector, std::string;
using std::uint32_t;
using std::flush, std::endl;

// bool code::LDPC::search_txt_to_uint(const string &path, string query, uint64_t &val){
// 	std::ifstream file(path, std::ios::in);
// 	if(!file.is_open()) throw std::filesystem::filesystem_error("cannot open file", path, std::make_error_code(std::errc::no_such_file_or_directory));

// 	string buf;
// 	while(std::getline(file, buf)){
// 		std::istringstream ss(buf);

// 		string key;
// 		ss>>key>>val;
// 		if(key == query) return true;
// 	}
// 	return false;
// }

bool util::search_txt_to_str(const string &path, string query, string &str){
	std::ifstream file(path, std::ios::in);
	if(!file.is_open()) throw std::filesystem::filesystem_error("cannot open file", path, std::make_error_code(std::errc::no_such_file_or_directory));

	string buf;
	while(std::getline(file, buf)){
		std::istringstream ss(buf);

		string key;
		ss>>key>>str;
		if(key == query) return true;
	}
	return false;
}
