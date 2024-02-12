#pragma once

#include <iostream>
#include "common.h"

namespace sprank {

class Vec3d {
	Float data[3];
public:
	Vec3d();//default
	Vec3d(const Float x, const Float y, const Float z);
	const Float& operator[](Size_Type i) const;//return data[i];
	      Float& operator[](Size_Type i);//return data[i];
	const Float norm_square() const;
	const Float norm() const;
	void normalize();
	const bool is_normalized() const;
	const Float operator*(const Vec3d& v) const;//dot product
	const Vec3d& operator+=(const Vec3d& v);
	const Vec3d& operator-=(const Vec3d& v);
	const Vec3d& operator+=(Float s);
	const Vec3d& operator-=(Float s);
	const Vec3d operator+(const Vec3d& v) const;
	const Vec3d operator-(const Vec3d& v) const;
	const Vec3d& operator*=(Float s);
	const Vec3d operator*(const Float s) const;
	const bool operator==(const Vec3d& rhs) const;
	void assign(Float s);
	Size_Type size() const;//return 3;
};
const Vec3d k_zero_vec3d(0, 0, 0);
const Vec3d k_max_vec3d(k_max_float, k_max_float, k_max_float);
const Vec3d k_nan_vec3d(k_not_a_num, k_not_a_num, k_not_a_num);

const Vec3d operator*(const Float s, const Vec3d& v);
const Vec3d vec3d_cross_product(const Vec3d& a, const Vec3d& b);
const Vec3d vec3d_elementwise_product(const Vec3d& a, const Vec3d& b);

const Float vec3d_distance_square(const Vec3d& a, const Vec3d& b);
// const Float square(const Vec3d& v);

}
