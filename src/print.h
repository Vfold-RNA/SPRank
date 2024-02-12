#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <iomanip>
#include "vec3d.h"

namespace sprank {

inline void print(Float x, std::ostream& out = std::cout) {
	out << x;
}
inline void print(int x, std::ostream& out = std::cout) {
	out << x;
}
inline void print(Size_Type x, std::ostream& out = std::cout) {
	out << x;
}
template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
	out << "[";
	for(Size_Type i = 0; i< v.size(); ++i) {
		if(i != 0)
			out << " ";
		out << v[i];
	}
	out << "]";
}
template<typename T>
void print(const std::set<T>& s, std::ostream& out = std::cout) {
	out << "[";
	for(const auto sv : s) {
		out << sv << " ";
	}
	out << "]";
}

inline void print(const Vec3d& v, std::ostream& out = std::cout) {
	out << "(";
	for(Size_Type i = 0; i< v.size(); ++i) {
		if(i != 0)
			out << ", ";
		out << v[i];
	}
	out << ")";
}

template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
	print(x, out);
	out << '\n';
}

}
