#pragma once

#include<vector>

#include "Scene.hpp"
#include "Intersection_info.hpp"
#include "Light_info.hpp"

struct Vertex_data_for_bdpt {
    const Vec3 position;
    const Vec3 normal;
    const Material bsdf;
    const FColor alpha;//頂点のウェイト
    R p_light;//光源側からトレースする時一つ前の頂点からこの頂点に来る確率密度(面積測度)
    R p_eye;//視線側からトレースする時一つ前の頂点からこの頂点に来る確率密度(面積測度)

    Vertex_data_for_bdpt(const Vec3 &position_,const Vec3 &normal_,const Material &bsdf_,const FColor &alpha_,const R p_light_,const R p_eye_) :
    position(position_),normal(normal_),bsdf(bsdf_),alpha(alpha_),p_light(p_light_),p_eye(p_eye_) {};
};

struct Bdpt {
    /*static*/ std::vector<Vertex_data_for_bdpt> light_sub_path_vertex;

    /*static*/ std::vector<Vertex_data_for_bdpt> eye_sub_path_vertex;

    /*static*/ R *light_P = nullptr;

    /*static*/ //R pdf_y0,pdf_z0;

    Bdpt() {
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

        const R u = Random::rando();

		for(int i = 0;i < scene.lights.size();i++) {
			if(u < light_P[i + 1] || i == (int)scene.lights.size() - 1){
				Light_info *light_info = scene.lights[i]->light_point();

                const R pdf_y0 = (light_P[i + 1] - light_P[i]) * light_info->pdf;

                const FColor alpha = light_info->emission / pdf_y0;

                light_sub_path_vertex.push_back(Vertex_data_for_bdpt(
                    light_info->point
                    ,light_info->normal
                    ,Material(FColor(M_PI,M_PI,M_PI))
                    ,alpha
                    ,pdf_y0
                    ,1//暫定的に入力
                ));

                Vec3 u;
				if (std::abs(light_info->normal.x) > EPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
					u = (cross(Vec3(0.0, 1.0, 0.0), light_info->normal)).normalized();
				else
					u = (cross(Vec3(1.0, 0.0, 0.0), light_info->normal)).normalized();
				Vec3 v = cross(light_info->normal,u);

				R u1 = Random::rando() * 2.0 * 3.1415926535,u2 = Random::rando(),u3 = std::sqrt(u2);

				Vec3 omega = 
					u * std::cos(u1) * u3 +
					v * std::sin(u1) * u3 +
					light_info->normal * std::sqrt(1 - u2);

                //light_sub_path_weight.push_back(light_sub_path_weight[0]);

                R cos_ = light_info->normal * omega;

                get_subpath(M_PI * alpha,cos_ * M_1_PI,cos_,Ray(light_info->point,omega),light_sub_path_vertex,scene,i,j,true);

                delete light_info;
                break;
            }
        }
    }

    /*static*/ void path_tracing(Scene &scene,const int i,const int j) {
        eye_sub_path_vertex.clear();

        const Vec3 z0 = scene.camera->sample_one_point_on_lens();
        const Vec3 &normal = scene.camera->dir;

        const R pdf_z0 = scene.camera->pdf_l;

        const R W = scene.camera->sensibility;

        const FColor alpha = FColor(W,W,W) / pdf_z0;

        eye_sub_path_vertex.push_back(Vertex_data_for_bdpt(
            z0
            ,normal
            ,Material(FColor(1,1,1))
            ,alpha
            ,1
            ,pdf_z0
        ));



        const Vec3 zp = scene.camera->sample_one_point_on_image_sensor(i,j);
        const Vec3 e = (scene.camera->position - zp).normalized();

        const Vec3 z1 = scene.camera->position + e * (scene.camera->dist_lens_object_plane / (e * scene.camera->dir));

        const Vec3 omega = (z1 - z0).normalized();

        const R cos_ = omega * normal;
        const R dp_square = scene.camera->dist_sensor_lens * scene.camera->dist_sensor_lens;

        const R p_sigma = dp_square * scene.camera->pdf_i / std::pow(cos_,4);

        const FColor alpha2 = alpha * std::pow((z0 - zp).normalized() * normal,4) * cos_ / dp_square / scene.camera->pdf_i;

        get_subpath(alpha2,p_sigma,cos_,Ray(z0,omega),eye_sub_path_vertex,scene,i,j,false);

    }

    /*static*/ void get_subpath(FColor alpha,R previous_p_sigma,R previous_cos,Ray ray,std::vector<Vertex_data_for_bdpt> &sub_path_vertex,Scene &scene,const int i,const int j,const bool is_light_tracing){
        while(true) {
            
            Intersection_info *intersection_info = get_intersection_of_nearest(ray,scene);

		    const Intersection_point *intersection = intersection_info->intersection_point;
		    const Material material = intersection->material;
		    const Vec3 normal = ((ray.direction * intersection->normal < 0.0) ? 1.0 : -1.0) * intersection->normal;
            const R dist_square = intersection->distance * intersection->distance;
            const R cos_ = normal * (-ray.direction);

            /*if(is_light_tracing) {

                int i = 0,j = 0;
                const R lens_t = scene.camera->get_intersection_with_lens(ray,i,j);
                if(lens_t > EPS && (intersection_info == nullptr || intersection->distance - lens_t > EPS))
                    std::cerr << i << " " << j << std::endl;

                if(lens_t > EPS && (intersection_info == nullptr || intersection->distance - lens_t > EPS)) {
                    scene.camera->img_l[i * scene.camera->pixel_w + j] = FColor(1,0,0);//FColor(depth == 0,0,0) ;FColor(std::max(1.0 - depth * 0.1,0.0),std::min(depth * 0.1,1.0),0);
                    break;
                }

            }*/

		    if(intersection_info == nullptr){//物体が存在しない
                break;
		    }

            const R PRR = (intersection_info->shape->light_id >= 0 ) ? 1 : std::min((R)1,std::max(std::max(material.kd.red,material.kd.green),material.kd.blue));

    		if(PRR < Random::rando()){
			    delete intersection_info;
                break;
		    }

            alpha = alpha / PRR;

            Vertex_data_for_bdpt v(
                intersection->position
                ,intersection->normal
                ,material
                ,alpha
                ,1//暫定
                ,1//暫定
            );//頂点データ

            const R rev_P = p_area_measure(v,-ray.direction,cos_,previous_cos,dist_square);//逆向きの確率密度（面積測度）を計算

            if(is_light_tracing) {
                sub_path_vertex[sub_path_vertex.size() - 1].p_eye = rev_P;
            }else {
                sub_path_vertex[sub_path_vertex.size() - 1].p_light = rev_P;
            }

            const R P = previous_p_sigma * cos_ / dist_square;//p_area_measure(sub_path_vertex[sub_path_vertex.size() - 1],ray.direction,previous_cos,cos_,dist_square);
            //順方向の確率密度（面積測度)を計算

            if(is_light_tracing) {
                v.p_light = P;
            }else {
                v.p_eye = P;
            }

            sub_path_vertex.push_back(v);

            if(!is_light_tracing) {//パストレだったら
                const Shape *shape = intersection_info->shape;

                if(shape->light_id >= 0) {//衝突したのが光源だったら、s=0の処理を行う
                    const int id = shape->light_id;

                    std::vector<Vertex_data_for_bdpt> x;

                    for(int i = sub_path_vertex.size() - 1;i >= 0;i--) {
                        x.push_back(sub_path_vertex[i]);
                    }

                    x[0].p_light = (light_P[id + 1] - light_P[id]) / shape->get_S();//ポリゴンが光源の時まずい

                    scene.camera->img_e[i * scene.camera->pixel_w + j] += material.Le * x[0].alpha * get_weight(x,0,x.size());
                }
            }

            Vec3 omega;//次のサンプリング

            //omegaに次のサンプリングの方向のセット,previous_p_sigmaの更新,previous_cosの更新,alphaの更新を行う
            switch (material.type){

            case MT_DEFAULT: {
    			Vec3 u;
    			if (std::abs(normal.x) > EPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
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

                previous_cos = omega * normal;

                previous_p_sigma = previous_cos * M_1_PI;

                alpha *= material.kd;

		    } break;
            default : {
                ;
            }break;
            }

            ray = Ray(intersection->position,omega);

            delete intersection_info;
            break;
        }
    }

    //x1からx2を選ぶ面積測度の確率密度PA(x2)を求める。ただし点x2はx1からomega12方向にレイを飛ばした時最初に衝突する点
    inline R p_area_measure(const Vertex_data_for_bdpt &x1
                            ,const Vec3 &omega12//x1→x2の方向ベクトル(正規化されたもの)
                            ,const R cos1//x1の法線とx1→x2方向がなす角のcos
                            ,const R cos2//x2の法線とx2→x1方向がなす角のcos
                            ,const R dist_square) //x1とx2の距離の二乗
                            {
        R P = 1;

        switch (x1.bsdf.type) {
        case MT_DEFAULT : {
            P = M_1_PI * cos1 * cos2 / dist_square;
        }break;
        default : ;

        }

        return P;
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
            for(int t = 2;t <= eye_sub_path_vertex.size();t++) {
                Vertex_data_for_bdpt y = light_sub_path_vertex[s - 1];
                Vertex_data_for_bdpt z = eye_sub_path_vertex[t - 1];

                if(y.normal * z.normal < 0.0) continue;//裏向き

                const Vec3 omega_y_in = (s > 1) ? (light_sub_path_vertex[s - 2].position - y.position).normalized() : Vec3(0,0,0);
                const Vec3 omega_y_out = (z.position - y.position).normalized();

                const FColor fy = calc_bsdf(y,omega_y_in,omega_y_out);

                const Vec3 omega_z_in = -omega_y_out;
                const Vec3 omega_z_out = (t > 1) ? (eye_sub_path_vertex[t - 2].position - z.position).normalized() : Vec3(0,0,0);

                const FColor fz = calc_bsdf(z,omega_z_in,omega_z_out);

                const R cos_y = y.normal * omega_y_out;
                const R cos_z = z.normal * omega_z_in;

                const R dist_square = (y.position - z.position).abs_square();

                const FColor c = fy * fz * (cos_y * cos_z / dist_square);

                std::vector<Vertex_data_for_bdpt> x;

                for(int i = 0;i < s;i++) x.push_back(light_sub_path_vertex[i]);
                for(int i = t - 1;i >= 0;i--) x.push_back(eye_sub_path_vertex[i]);

                x[s - 1].p_eye = p_area_measure(z,omega_z_in,cos_z,cos_y,dist_square);
                x[s].p_light = p_area_measure(y,omega_y_out,cos_y,cos_z,dist_square);

                const R w = get_weight(x,s,t);

                Ray ray(y.position,z.position - y.position);

                Intersection_info *intersection = get_intersection_of_nearest(ray,scene);

                if(intersection != nullptr && (intersection->intersection_point->position - z.position).abs() >= EPS) {
                    delete intersection;
                    continue;
                }
                delete intersection;

                if(t <= 1) {
                    ;
                }else {
                    scene.camera->img_e[i * scene.camera->pixel_w + j] += w * light_sub_path_vertex[s - 1].alpha * c * eye_sub_path_vertex[t - 1].alpha;
                }
            }
        }
    }

    inline FColor calc_bsdf(const Vertex_data_for_bdpt &x,const Vec3 &omega_in,const Vec3 &omega_out) {
        FColor ret;
        switch (x.bsdf.type) {

        case MT_DEFAULT : {
            ret = x.bsdf.kd * M_1_PI;
        }break;
        default : {
        ;
        }
        }
        return ret;
    }

    /*static*/ inline R get_weight(const std::vector<Vertex_data_for_bdpt> &path_vertex,const int s,const int t) {
        const int k = (int)path_vertex.size() - 1;

        R *P = new R[k + 1];//Pi = pi+1 / pi (i=0,1,...,k)

        for(int i = 0;i <= k;i++) {
            P[i] = path_vertex[i].p_light / path_vertex[i].p_eye;
        }

        R w = 1;

        R m = 1;

        for(int i = s - 1;i >= 0;i--) {
            m *= P[i];
            w += m * m;
        }

        m = 1;

        for(int i = s + 1;i <= k - 2/*-2はt=0,1を除外*/;i++) {
            m *= P[i];
            w += m * m;
        }

        delete[] P;

        return 1 / w;
    }

    /*static*/ /*inline bool is_shadow(const Scene &scene,const Vec3 &y,const Vec3 &z) {
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
    }*/

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