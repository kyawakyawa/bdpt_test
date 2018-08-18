#pragma once

#include<vector>

#include "Scene.hpp"
#include "Intersection_info.hpp"
#include "Light_info.hpp"

struct Vertex_data_for_bdpt {
    const Vec3 position;
    const Vec3 normal;
    //const Vec3 omega_in;
    const Material bsdf;
    //const FColor weight;
    //const R p;
    //const R russian_roulette;

    Vertex_data_for_bdpt(const Vec3 &position_,const Vec3 &normal_/*,Vec3 &omega_in_*/,const Material &bsdf_/*,const FColor &weight_*//*,R p_*//*,const R russian_roulette_*/) :
    position(position_),normal(normal_)/*,omega_in(omega_in_)*/,bsdf(bsdf_)/*,weight(weight_)*//*,p(p_)*//*,russian_roulette(russian_roulette_)*/ {};
};

struct Bdpt {
    /*static*/ std::vector<Vertex_data_for_bdpt> light_sub_path_vertex;
    /*static*/ std::vector<FColor> light_sub_path_weight;

    /*static*/ std::vector<Vertex_data_for_bdpt> eye_sub_path_vertex;
    /*static*/ std::vector<FColor> eye_sub_path_weight;

    /*static*/ R *light_P;

    /*static*/ R pdf_y0,pdf_z0;

    Bdpt() {
        light_P = nullptr;
    }

    /*static*/ void bdpt(Scene &scene,int samples) {
        scene.camera->init();

        if(light_P == nullptr)set_light(scene);

        for(int k = 0;k < samples;k++) {
            for(int i = 0;i < scene.camera->pixel_h;i++) {
                for(int j = 0;j < scene.camera->pixel_w;j++) {
                    light_tracing(scene,i,j);
                    path_tracing(scene,i,j);
                    merge(scene,i,j);
                }
            }
        }

        for(int i = 0;i < scene.camera->pixel_h;i++) {
            for(int j = 0;j < scene.camera->pixel_w;j++) {
                scene.camera->img[i * scene.camera->pixel_w + j] = 
                    (scene.camera->img_e[i * scene.camera->pixel_w + j] + scene.camera->img_l[i * scene.camera->pixel_w + j]) / samples;
            }
        }
    }

    /*static*/ void light_tracing(Scene &scene,const int i,const int j) {
        light_sub_path_vertex.clear();
        light_sub_path_weight.clear();

        light_sub_path_weight.push_back(FColor(1,1,1));

        const R u = Random::rando();

		for(int i = 0;i < scene.lights.size();i++) {
			if(u < light_P[i + 1] || i - 1 == scene.lights.size()){
				Light_info *light_info = scene.lights[i]->light_point();

                light_sub_path_vertex.push_back(Vertex_data_for_bdpt(
                    light_info->point
                    ,light_info->normal
                    ,Material(FColor(1,1,1))
                ));

                pdf_y0 = (light_P[i + 1] - light_P[i]) * light_info->pdf;

                light_sub_path_weight.push_back(light_info->emission / pdf_y0);


                Vec3 u;
				if (std::abs(light_info->normal.x) > 1e-6) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
					u = (cross(Vec3(0.0, 1.0, 0.0), light_info->normal)).normalized();
				else
					u = (cross(Vec3(1.0, 0.0, 0.0), light_info->normal)).normalized();
				Vec3 v = cross(light_info->normal,u);

				R u1 = Random::rando() * 2.0 * 3.1415926535,u2 = Random::rando(),u3 = std::sqrt(u2);

				Vec3 omega = 
					u * std::cos(u1) * u3 +
					v * std::sin(u1) * u3 +
					light_info->normal * std::sqrt(1 - u2);

                light_sub_path_weight.push_back(light_sub_path_weight[0]);

                get_subpath(Ray(light_info->point,omega),light_sub_path_vertex,light_sub_path_weight,scene,i,j,true);

                delete light_info;
                break;
            }
        }
    }

    /*static*/ void path_tracing(Scene &scene,const int i,const int j) {
        eye_sub_path_vertex.clear();
        eye_sub_path_weight.clear();

        eye_sub_path_weight.push_back(FColor(1,1,1));

        const Vec3 x0 = scene.camera->sample_one_point_on_lens();
        const Vec3 &normal = scene.camera->dir;

        eye_sub_path_vertex.push_back(Vertex_data_for_bdpt(
            x0
            ,normal
            ,Material(FColor(1,1,1))
        ));

        const R W = scene.camera->sensibility;

        pdf_z0 = scene.camera->pdf_l;

        eye_sub_path_weight.push_back(FColor(W,W,W) / pdf_z0);

        const Vec3 xp = scene.camera->sample_one_point_on_image_sensor(i,j);
        const Vec3 e = (scene.camera->position - xp).normalized();

        const Vec3 x1 = scene.camera->position + e * scene.camera->dist_lens_object_plane / (e * scene.camera->dir);

        const Vec3 omega = (x1 - x0).normalized();

        eye_sub_path_weight.push_back(std::pow(omega * normal,(R)4.0) / ((R)3.1415926535 * scene.camera->dist_sensor_lens * scene.camera->dist_sensor_lens * scene.camera->pdf_i) * eye_sub_path_weight[0]);
        get_subpath(Ray(x0,omega),eye_sub_path_vertex,eye_sub_path_weight,scene,i,j,false);

    }

    /*static*/ void get_subpath(Ray ray,std::vector<Vertex_data_for_bdpt> &sub_path_vertex,std::vector<FColor> &sub_path_weight,Scene &scene,const int i,const int j,const bool is_light_tracing){
        int depth = 0;
        const int min_depth = 2;
        const int max_depth = 8;
        while(true) {
            
            Intersection_info *intersection_info = get_intersection_of_nearest(ray,scene);

		    if(intersection_info == nullptr){//物体が存在しない
                break;
		    }
		
		    const Intersection_point *intersection = intersection_info->intersection_point;
		    const Material material = intersection->material;
		    const Vec3 normal = ((ray.direction * intersection->normal < 0.0) ? 1.0 : -1.0) * intersection->normal;
		
		    //FColor L = (ray.direction * intersection->normal < 0.0) ? material.Le : FColor(0,0,0);

		    R P = std::max(std::max(material.kd.red,material.kd.green),material.kd.blue);

		    if(depth < min_depth)
			    P = 1.0;
		
		    if (depth >= max_depth)//最大値で大幅に確率を下げる
			    P *= pow(0.5, depth - max_depth);
		
		

		    if(P < Random::rando()){//再帰が深すぎ
			    delete intersection_info;
                break;
		    }

            if(is_light_tracing) {
                ;
            }else {
                const Shape *shape = intersection_info->shape;

                if(shape->light_id >= 0) {//衝突したのが光源だったら
                    const int id = shape->light_id;
                    const Light_source *light = scene.lights[id];

                    const int store = pdf_y0;//weight計算の光源側のpdfが変わるので記憶しておく
                    pdf_y0 = (light_P[i + 1] - light_P[i]) / shape->get_S();
                }
            }

            switch (material.type){

            case MT_DEFAULT: {
			Vec3 u;
			if (std::abs(normal.x) > 1e-6) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
				u = (cross(Vec3(0.0, 1.0, 0.0), normal)).normalized();
			else
				u = (cross(Vec3(1.0, 0.0, 0.0), normal)).normalized();
			Vec3 v = cross(normal,u);

			R u1 = Random::rando() * 2.0 * 3.1415926535,u2 = Random::rando(),u3 = std::sqrt(u2);

			Vec3 omega = 
				u * std::cos(u1) * u3 +
				v * std::sin(u1) * u3 +
				normal * std::sqrt(1 - u2);

            sub_path_vertex.push_back(Vertex_data_for_bdpt(intersection->position,normal,material));
            sub_path_weight.push_back(material.kd * *(sub_path_weight.end() - 1) / P);

            ray = Ray(intersection->position,omega);

		    } break;
            default : {
                ;
            }break;
            }

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

    /*static*/ void merge(Scene &scene,const int i,const int j) {
        for(int s = 1;s <= light_sub_path_vertex.size();s++) {
            for(int t = 1;t <= eye_sub_path_vertex.size();t++) {
                Vertex_data_for_bdpt y = light_sub_path_vertex[s - 1];
                Vertex_data_for_bdpt z = eye_sub_path_vertex[t - 1];
                FColor fy,fz;

                switch (y.bsdf.type) {

                case MT_DEFAULT : {
                    fy = y.bsdf.kd / 3.1415926535;
                }break;
                default : {
                    ;
                }
                }

                switch (z.bsdf.type) {

                case MT_DEFAULT : {
                    fz = z.bsdf.kd / 3.1415926535;
                }break;
                default : {
                    ;
                }
                }

                const FColor c = fy * fz * G(y,z);

                std::vector<Vertex_data_for_bdpt> x;

                for(int i = 0;i < s;i++) x.push_back(light_sub_path_vertex[i]);
                for(int i = t - 1;i >= 0;i--) x.push_back(eye_sub_path_vertex[i]);

                const R w = get_weight(x,s,t);

                if(y.normal * z.normal >= 0 && is_shadow(scene,y.position,z.position)) {
                    continue;
                }

                if(t <= 1) {
                    ;
                }else {
                    scene.camera->img_e[i * scene.camera->pixel_w + j] += w * light_sub_path_weight[s] * c * eye_sub_path_weight[s];
                }
            }
        }
    }

    /*static*/ inline R G(const Vertex_data_for_bdpt &x,const Vertex_data_for_bdpt &y) {
        Vec3 omega = (x.position - y.position).normalized();

        return std::abs(omega * x.normal) * std::abs(omega * y.normal) / (x.position - y.position).abs_square();
    }

    /*static*/ inline R get_weight(const std::vector<Vertex_data_for_bdpt> &path_vertex,const int s,const int t) {
        const int k = (int)path_vertex.size() - 1;

        R *Gi = new R[k];

        for(int i = 0;i < k;i++) {
            Gi[i] = G(path_vertex[i + 1],path_vertex[i]);
        }

        R *P = new R[k + 1];//Pi = pi+1 / pi (i=0,1,...,k)

        R denominator,numerator;
        for(int i = 1;i < k;i++) {

            switch (path_vertex[i + 1].bsdf.type) {
            case MT_DEFAULT : {
                denominator = Gi[i];
            }break;
            default : {
                ;
            }
            }

            switch (path_vertex[i - 1].bsdf.type) {
            case MT_DEFAULT : {
                numerator = Gi[i - 1];
            }break;
            default : {

            }
            }

            P[i] = numerator / denominator;


        }


        switch (path_vertex[1].bsdf.type) {
        case MT_DEFAULT : {
            denominator = Gi[0];
        }break;
        default : {
            ;
        }
        }

        numerator = pdf_y0;

        P[0] = numerator / denominator;

        switch (path_vertex[k - 1].bsdf.type) {
        case MT_DEFAULT : {
            numerator = Gi[k - 1];
        }break;
        default : {
            ;
        }
        }

        denominator = pdf_z0;

        P[k] = numerator / denominator;

        R w = 1;
        R m = 1;

        for(int i = s - 1;i >= 0;i--) {
            m /= P[i];
            w += m * m;
        }
        m = 1;
        for(int i = s + 1;i <= k;i++) {
            m *= P[i];
            w += m * m;
        }

        delete[] Gi;
        delete[] P;

        return 1 / w;
    }

    /*static*/ inline bool is_shadow(const Scene &scene,const Vec3 &y,const Vec3 &z) {
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

    ~Bdpt() {
        delete[] light_P;
    }

};