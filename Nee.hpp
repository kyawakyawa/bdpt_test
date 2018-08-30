#pragma once

#include "Scene.hpp"
#include "Intersection_info.hpp"
#include "Light_info.hpp"

struct Nee {
    Nee(){}

    /*static*/ R *light_P = nullptr;

    void nee(Scene &scene,int samples) {
        scene.camera->init();

        if(light_P == nullptr)set_light(scene);

        for(int k = 0;k < samples;k++) {
            #pragma omp parallel for schedule(dynamic, 1) num_threads(8) 
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

    /*static*/ void set_light(Scene &scene) {
   		if(scene.lights.size() < 1)
			return;
		light_P = new R[scene.lights.size() + 1];

		light_P[0] = 0.0;

		for(int i = 0;i < scene.lights.size();i++) {
			light_P[i + 1] = light_P[i] + scene.lights[i]->get_weight();
		}
		const R d = 1.0 / light_P[scene.lights.size()];
		for(int i = 1;i < scene.lights.size() + 1;i++) {
			light_P[i] *= d;
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

        const R cos_ = omega * normal;
        const R dp_square = scene.camera->dist_sensor_lens * scene.camera->dist_sensor_lens;

        //const R p_sigma = dp_square * scene.camera->pdf_i / std::pow(cos_,4);

        int depth = 0;
        const int min_depth = 4;
        const int max_depth = 8;

        //int previous_cos = cos_;
        //int previous_p_sigma = p_sigma;

        Ray ray(z0,omega);
        // alpha *= std::pow((z0 - zp).normalized() * normal,4) * cos_ / dp_square;

        FColor alpha = FColor(W,W,W) * std::pow((z0 - zp).normalized() * normal,4) * cos_ / dp_square / scene.camera->pdf_i / p_A_z0;

        bool flag = true;


        while(true) {

            Intersection_info *intersection_info = get_intersection_of_nearest(ray,scene);

            if(intersection_info == nullptr)break;

            const Intersection_point *intersection = intersection_info->intersection_point;
		    const Material material = intersection->material;
		    const Vec3 normal = ((ray.direction * intersection->normal < 0.0) ? 1.0 : -1.0) * intersection->normal;

            if(flag && ray.direction * intersection->normal < 0.0){
                scene.camera->img[i * scene.camera->pixel_w + j] += alpha * material.Le;
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

            Vec3 omega;//次のサンプリング

            //omegaに次のサンプリングの方向のセット,previous_p_sigmaの更新,previous_cosの更新,alphaの更新を行う
            switch (material.type){

            case MT_DEFAULT: {
                flag = false;

                FColor lgp;Vec3 omegal;
			    nee(scene,intersection->position,normal,lgp,omegal);
			    scene.camera->img[i * scene.camera->pixel_w + j] += alpha * material.kd * M_1_PI * lgp;

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

                alpha *= material.kd;

		    } break;
            case MT_PERFECT_REF: {
                flag = true;
				omega = -2.0 * (normal * ray.direction) * normal + ray.direction;
				alpha *= material.kd;
			} break;
			case MT_REFRACTION: {
                flag = true;
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

    inline void nee(const Scene &scene,const Vec3 &position,const Vec3 &normal,FColor &lgp,Vec3 &omegal) const {

		const R u = Random::rando();

		for(int i = 0;i < scene.lights.size();i++) {
			if(u < light_P[i + 1] || i == (int)scene.lights.size() - 1){
				Light_info *light_info = scene.lights[i]->light_point();

				const Vec3 path = light_info->point - position;
				omegal = path.normalized();

				if(is_shadow(scene,position,light_info->point)) {
					delete light_info;
					break;
				}

				const R wl_dot_nl = (-omegal) * light_info->normal;
				if(wl_dot_nl < EPS){
					delete light_info;
					break;
				}

				const R w1_dot_n1 = omegal * normal;
				if(w1_dot_n1  < EPS){
					delete light_info;
					break;
				}

				lgp = light_info->emission
						* w1_dot_n1 * wl_dot_nl
						* (1.0 / (path.abs_square() * (light_P[i + 1] - light_P[i]) * light_info->pdf));
				//std::cerr << (light_P[i + 1] - light_P[i]) * light_info->pdf << std::endl;
				delete light_info;
				return;
			}
		}
		lgp = FColor(0,0,0);
	}    
    /*static*/ inline bool is_shadow(const Scene &scene,const Vec3 &y,const Vec3 &z) const {
		const R max_t = (y - z).abs() - EPS;
		for(Shape *shape : scene.shapes){
			Intersection_point *intersection = shape->get_intersection(Ray(y,z - y));
			if(intersection != nullptr && max_t > intersection->distance ){
				delete intersection;
				return true;
			}
			delete intersection;
		}
		return false;
    }


};