#include <cassert>
#include <cmath>
#include "vec3d.h"
#include "common.h"

namespace sprank {

Vec3d::Vec3d() {
	data[0] = data[1] = data[2] = k_not_a_num;
}
Vec3d::Vec3d(const Float x, const Float y, const Float z) {
	data[0] = x;
	data[1] = y;
	data[2] = z;
}
const Float& Vec3d::operator[](Size_Type i) const {
	assert(i < 3);
	return data[i];
}
Float& Vec3d::operator[](Size_Type i) {
	assert(i < 3);
	return data[i];
}
const Float Vec3d::norm_square() const {
	return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
}
const Float Vec3d::norm() const {
	return std::sqrt(norm_square());
}

void Vec3d::normalize() {
	const Float s = this->norm();
	if(s > k_epsilon) {
		this->data[0] = this->data[0]/s;
		this->data[1] = this->data[1]/s;
		this->data[2] = this->data[2]/s;
	}
	assert(eq(this->norm(),1));
}

const bool Vec3d::is_normalized() const {
	return eq(this->norm(),1.0);
}

const Float Vec3d::operator*(const Vec3d& v) const {
	return data[0]*v[0] + data[1]*v[1] + data[2]*v[2];
}
const Vec3d& Vec3d::operator+=(const Vec3d& v) {
	data[0] += v[0];
	data[1] += v[1];
	data[2] += v[2];
	return *this;
}
const Vec3d& Vec3d::operator-=(const Vec3d& v) {
	data[0] -= v[0];
	data[1] -= v[1];
	data[2] -= v[2];
	return *this;
}
const Vec3d& Vec3d::operator+=(Float s) {
	data[0] += s;
	data[1] += s;
	data[2] += s;
	return *this;
}
const Vec3d& Vec3d::operator-=(Float s) {
	data[0] -= s;
	data[1] -= s;
	data[2] -= s;
	return *this;
}
const Vec3d Vec3d::operator+(const Vec3d& v) const {
	return Vec3d(data[0] + v[0],
				 data[1] + v[1],
				 data[2] + v[2]);
}
const Vec3d Vec3d::operator-(const Vec3d& v) const {
	return Vec3d(data[0] - v[0],
					data[1] - v[1],
					data[2] - v[2]);
}
const Vec3d& Vec3d::operator*=(Float s) {
	data[0] *= s;
	data[1] *= s;
	data[2] *= s;
	return *this;
}
const Vec3d Vec3d::operator*(const Float s) const {
	return Vec3d(s*data[0], s*data[1], s*data[2]);
}

const bool Vec3d::operator==(const Vec3d& rhs) const {
	return eq(this->data[0], rhs[0]) && eq(this->data[1], rhs[1]) && eq(this->data[2], rhs[2]);
}

void Vec3d::assign(Float s) {
	data[0] = data[1] = data[2] = s;
}

Size_Type Vec3d::size() const { return 3; }

//////////////////non-member/////////////////////

const Vec3d operator*(const Float s, const Vec3d& v) {
	return Vec3d(s * v[0], s * v[1], s * v[2]);
}

const Vec3d vec3d_cross_product(const Vec3d& a, const Vec3d& b) {
	return Vec3d(a[1]*b[2] - a[2]*b[1],
		         a[2]*b[0] - a[0]*b[2],
		         a[0]*b[1] - a[1]*b[0]);
}

const Vec3d vec3d_elementwise_product(const Vec3d& a, const Vec3d& b) {
	return Vec3d(a[0] * b[0],
		         a[1] * b[1],
			     a[2] * b[2]);
}

const Float vec3d_distance_square(const Vec3d& a, const Vec3d& b) {
	return (a[0] - b[0])*(a[0] - b[0]) + \
		   (a[1] - b[1])*(a[1] - b[1]) + \
		   (a[2] - b[2])*(a[2] - b[2]);
}

// const Float square(const Vec3d& v) {
// 	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
// }

}
