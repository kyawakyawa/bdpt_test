#pragma once

#include "Scene.hpp"
#include "Intersection_info.hpp"
#include "Light_info.hpp"

struct Pt {
    Pt(){}

    void pt(Scene &scene,int samples) {
        scene.camera->init();

        for(int k = 0;k < samples;k++) {
            #pragma omp parallel for schedule(dynamic, 1) num_threads(4) 
            for(int i = 0;i < scene.camera->pixel_h;i++) {
                for(int j = 0;j < scene.camera->pixel_w;j++) {
                    path_tracing(scene,i,j);
                }
            }
        }

        for(int i = 0;i < scene.camera->pixel_h;i++) {
            for(int j = 0;j < scene.camera->pixel_w;j++) {
                scene.camera->img[i * scene.camera->pixel_w + j] = scene.camera->img[i * scene.camera->pixel_w + j] / (R)samples;
            }
        }
    }

    void path_tracing(Scene &scene,int i,int j) {
        const Vec3 z0 = scene.camera->sample_one_point_on_lens();
        const Vec3 &normal = scene.camera->dir;

        const R p_A_z0 = scene.camera->pdf_l;

        const R W = scene.camera->sensibility;

        //FColor alpha = FColor(W,W,W) / p_A_z0;

        const Vec3 zp = scene.camera->sample_one_point_on_image_sensor(i,j);
        const Vec3 e = (scene.camera->position - zp).normalized();

        const Vec3 z1 = scene.camera->position + e * scene.camera->dist_lens_object_plane / (e * normal);

        const Vec3 omega = (z1 - z0).normalized();

        //const R cos_ = omega * normal;
        //const R dp_square = scene.camera->dist_sensor_lens * scene.camera->dist_sensor_lens;

        //const R p_sigma = dp_square * scene.camera->pdf_i / std::pow(cos_,4);

        int depth = 0;
        const int min_depth = 4;
        const int max_depth = 8;

        //int previous_cos = cos_;
        //int previous_p_sigma = p_sigma;

        Ray ray(z0,omega);

        const R cos_ = (z0 - zp).normalized() * normal;

        // alpha *= std::pow((z0 - zp).normalized() * normal,4) * cos_ / dp_square;
        FColor alpha = FColor(W,W,W) * cos_ * cos_ / (z0 - zp).abs_square() / scene.camera->pdf_i / scene.camera->pdf_l;

        while(true) {

            Intersection_info *intersection_info = get_intersection_of_nearest(ray,scene);

            if(intersection_info == nullptr)break;

            const Intersection_point *intersection = intersection_info->intersection_point;
		    const Material material = intersection->material;
		    const Vec3 normal = ((ray.direction * intersection->normal < 0.0) ? 1.0 : -1.0) * intersection->normal;

            scene.camera->img[i * scene.camera->pixel_w + j] += alpha * material.Le;

            Vec3 omega;//次のサンプリング

            //omegaに次のサンプリングの方向のセット,previous_p_sigmaの更新,previous_cosの更新,alphaの更新を行う
            switch (material.type){

            case MT_DEFAULT: {
    			Vec3 u;
    			if (std::abs(normal.x) > 1e-6) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
    				u = (cross(Vec3(0.0, 1.0, 0.0), normal)).normalized();
    			else
    				u = (cross(Vec3(1.0, 0.0, 0.0), normal)).normalized();
    			Vec3 v = cross(normal,u);

    			R u1 = Random::rando() * 2.0 * 3.1415926535,u2 = Random::rando(),u3 = std::sqrt(u2);

    			omega = 
    				u * std::cos(u1) * u3 +
    				v * std::sin(u1) * u3 +
    				normal * std::sqrt(1 - u2);


                //sub_path_vertex.push_back(Vertex_data_for_bdpt(intersection->position,normal,material));
                //sub_path_weight.push_back(material.kd * *(sub_path_weight.end() - 1) / P);

                //previous_cos = omega * normal;

                //previous_p_sigma = previous_cos * M_1_PI;

                alpha *= material.kd;

		    } break;
            case MT_PERFECT_REF: {
				omega = -2.0 * (normal * ray.direction) * normal + ray.direction;
				alpha *= material.kd;
			} break;
			case MT_REFRACTION: {
				const R n = ((ray.direction * intersection->normal < 0.0) ? 1.0 / material.n : material.n);
				const R dn = ray.direction * normal;
				const R cos2t = 1.0 - n * n * (1.0 - dn * dn);

				//const Vec3 p = intersection->position;
				const Vec3 r = -2.0 * (normal * ray.direction) * normal + ray.direction;

				if(cos2t < 0.0){
					//ray = Ray(p,r);
                    omega = r;
					alpha *= material.kd;
					break;
				}

				Vec3 T = (n * ray.direction - (n * dn + std::sqrt(cos2t) ) * normal).normalized();
				const R a = 1.0 - material.n,b = 1.0 + material.n;
				const R F0 = (a * a) / (b * b);

				const R Fr = F0 + (1.0 - F0) * std::pow(1.0 + ((ray.direction * intersection->normal < 0.0) ? dn : T * normal),5);
				const R Tr = (1.0 - Fr) * n * n;

				const R probability = 0.25 * 0.5 * Fr;
				//if(depth > 2){
					if(Random::rando() < probability){
						//ray = Ray(p,r);
                        omega = r;
						alpha *= material.kd * (1.0 / probability) * Fr;
					}else{
						//ray = Ray(p,T);
                        omega = T;
						alpha *= material.kd * (1.0 / (1.0 - probability)) * Tr;
					}
				//}else {
					//L += material.kd * (pathtracing(Ray(p,r),depth + 1) * Fr + pathtracing(Ray(p,T),depth + 1) * Tr) * (1.0 / P);
				//}
			} break;
            default : {
                ;
            }break;
            }

            R PRR = std::min((R)1,std::max(std::max(material.kd.red,material.kd.green),material.kd.blue));

    		if(depth < min_depth)
			    PRR = 1.0;
		
		    if (depth >= max_depth)//最大値で大幅に確率を下げる
			    PRR *= pow(0.5, depth - max_depth);
		
		

    		if(PRR < Random::rando()){
			    delete intersection_info;
                break;
		    }

            alpha = alpha / PRR;

            ray = Ray(intersection->position,omega);

            depth++;

            delete intersection_info;
        }
    }


    /*static*/ inline Intersection_info* get_intersection_of_nearest(const Ray &ray,const Scene &scene) {
		R min_t = INF;
		Intersection_info *intersection_info = new Intersection_info();

		for(Shape *shape : scene.shapes){
			Intersection_point *intersection = shape->get_intersection(ray);

			if(intersection != nullptr && min_t > intersection->distance){
				intersection_info->rewrite(shape,intersection);
				min_t = intersection->distance;
			}else
				delete intersection;
		}

		if(intersection_info->shape == nullptr){
			delete intersection_info;
			return nullptr;
		}

		return intersection_info;
    }
};