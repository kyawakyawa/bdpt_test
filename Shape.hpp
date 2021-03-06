#pragma once

#include <vector>

#include "Vec3.hpp"
#include "Material.hpp"
#include "Ray.hpp"
#include "Intersection_point.hpp"

struct Shape{
	std::vector< Material > materials;//材質
	int light_id = -1;//光源のインデックス

	inline Shape() = default;

	inline virtual Intersection_point* get_intersection(const Ray &ray) const = 0;//レイとの交差についての情報を返す

	inline virtual R get_S() const = 0;

	inline virtual void sample_one_point (Vec3 &point,Vec3 &normal,R &pdf,FColor &emission) const = 0;

	inline virtual ~Shape() = default;
};
