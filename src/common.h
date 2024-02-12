#pragma once

#include <limits>
#include <cmath>
#include <vector>

namespace sprank {

using Float = double;
using Floats = std::vector<Float>;
using Size_Type = std::size_t;
const Size_Type k_max_size_type = (std::numeric_limits<Size_Type>::max)();
// using Size_Types = std::vector<std::size_t>;
// using Float_Pair = std::pair<Float,Float>;
// using Float_Pairs = std::vector<Float_Pair>;

const Float k_pi = Float(3.1415926535897931);
const Float k_max_float = (std::numeric_limits<Float>::max)();
const Float k_min_float = (std::numeric_limits<Float>::min)();
const unsigned k_max_unsigned = (std::numeric_limits<unsigned>::max)();
const int k_max_int = (std::numeric_limits<int>::max)();
const Float k_epsilon = std::numeric_limits<Float>::epsilon();
const Float k_tolerance = Float(0.001);
const Float k_not_a_num = std::sqrt(Float(-1)); // FIXME? check

inline const bool eq(Float a, Float b) {
	return std::abs(a - b) < k_tolerance;
}
inline const Float square(const Float a) {
	return a*a;
}

// template<unsigned n>
// inline Float int_pow(Float x) {
// 	return int_pow<n-1>(x)*x;
// }
// template<>
// inline Float int_pow<0>(Float x) {
// 	return 1;
// }

// inline const bool not_max(Float x) {
// 	return (x < 0.1 * k_max_float);
// }

}

