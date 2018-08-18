#pragma once

#include <vector>
#include <chrono>
#include "Shape.hpp"
#include "Intersection_info.hpp"
#include "Light_source.hpp"
#include "Area_light_source.hpp"
#include "Light_info.hpp"
#include "Random.hpp"
#include "Camera.hpp"

struct Scene{
	std::vector<Shape*> shapes;//物体
	std::vector<Light_source*> lights;//光源
	const int threads = 4; //スレッドの数
	Camera *camera;//カメラ

	Scene(Camera *c) : camera(c) {

	}

	inline void add(Shape *shape){//物体を追加する
		shapes.push_back(shape);
	}

	inline void add_as_light(Shape *shape){
		add(shape);

		shape->light_id = lights.size();

		lights.push_back(new Area_light_source(shape));
	}

	~Scene(){
		for(Shape *shape : shapes)
			delete shape;
		for(Light_source *light : lights)
			delete light;
	}
};
