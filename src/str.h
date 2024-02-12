#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <sstream>

namespace sprank {

inline std::vector<std::string> string2vector(const std::string sline) {
	std::istringstream ss(sline);
	std::string buf;
	std::vector<std::string> token;
	while(ss >> buf) token.push_back(buf);
	return token;
}

// inline bool starts_with(const std::string& str, const std::string& start) {
// 	return str.size() >= start.size() && str.substr(0, start.size()) == start;
// }

inline std::string trim(const std::string &str, const std::string &whitespace = " \t\n\r") {
	const auto str_begin = str.find_first_not_of(whitespace);
	if (str_begin == std::string::npos){
		return ""; // no content
	}
	const auto str_end = str.find_last_not_of(whitespace);
	const auto str_range = str_end - str_begin + 1;
	return str.substr(str_begin, str_range);
}

}
