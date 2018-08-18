#pragma once

#include "Vec3.hpp"
#include "FColor.hpp"
#include "Ray.hpp"
#include "Intersection_point.hpp"
#include "Random.hpp"

struct Camera {
    const Vec3 position = Vec3(0,0,0);//レンズの中心の位置
    const Vec3 dir = Vec3(0,0,-1);
    const Vec3 up = Vec3(0,1,0);

    const int pixel_w = 480;
    const int pixel_h = 480;

    long sample_num = 0;

    FColor *img = nullptr;// 合わせた結果
    FColor *img_l = nullptr;// t = 0,1 の時の結果の記憶
    FColor *img_e = nullptr;// t >= 2 の時の結果の記憶

    const Vec3 x,y;

    const R width_of_sensor = 40;
    const R height_of_sensor;
    const R dist_sensor_lens;
    const Vec3 sensor_center;

    const R focal_length;//レンズの焦点距離
    const R lens_radius;//レンズの半径
    const R f_number;

    const R dist_lens_object_plane;

    const R pdf_i;
    const R pdf_l;

    const R sensibility;

    Camera();
    Camera(const int w,const int h);
    Camera(const int w,const int h,const Vec3 &p,const Vec3 &d,const Vec3 &u,const R fov_,const R focus_distance,const R f_number_) : position(p),dir(d),up(u),pixel_w(w),pixel_h(h)
    ,img(new FColor[h * w]),img_l(new FColor[h * w]),img_e(new FColor[h * w])
    ,x(cross(dir,up).normalized()),y(cross(dir,x).normalized())
    ,height_of_sensor(width_of_sensor * h / w)
    ,dist_sensor_lens(0.5 * width_of_sensor / std::tan(fov_ / 360.0 * M_PI))
    ,sensor_center(-dir * dist_sensor_lens)

    ,focal_length(focus_distance * dist_sensor_lens / (focus_distance + dist_sensor_lens)),lens_radius(focal_length / f_number_ / 2.0),f_number(f_number_)

    ,dist_lens_object_plane(focus_distance)

    ,pdf_i(1.0 / (height_of_sensor * width_of_sensor / pixel_h / pixel_w)),pdf_l(1.0 / M_PI / lens_radius / lens_radius)
    ,sensibility(0.001 * pdf_i) {};

    void init() {
        for(int i = 0;i < pixel_h * pixel_w;i++)img[i] = img_l[i] = img_e[i] = FColor(0,0,0);
    }

    Intersection_point* get_intersection_with_lens(const Ray &ray) {
        const Vec3 &d = ray.direction;
		Vec3 s2;

		////平面上のもう一つの点を取得
		if(std::abs(dir.x) > 1e-9)
			s2 = Vec3(-(dir.y + dir.z) / dir.x,1,1);
		else if(std::abs(dir.y) > 1e-9)
			s2 = Vec3(1,-(dir.x + dir.z) / dir.y,1);
		else
			s2 = Vec3(1,1,-(dir.x + dir.y) / dir.z);

		const Vec3 s = ray.start - (s2 + position);

		if(d * dir == 0)//レイと平面が並行
			return nullptr;

		const R t = -(s * dir) / (d * dir);

		if(t < 0.0)//平面が視線の後側
			return nullptr;

        if(d * dir > 0)//レンズの裏側から当たる
            return nullptr;

        const Vec3 intersection_position = ray.start + t * d;

        if((position - intersection_position).abs() >= lens_radius) //レンズの外
            return nullptr;


		return new Intersection_point(t,intersection_position,dir,Material(FColor(0,0,0)));//ダミーのマテリアル
    }

    Vec3 sample_one_point_on_lens() {
        R u1,u2;
        do {
            u1 = Random::rando();
            u2 = Random::rando();
        }while(u1 * u1 + u2 * u2 > 1.0);

        return position + u1 * lens_radius * x + u2 * lens_radius * y;
    }

    Vec3 sample_one_point_on_image_sensor(const int i,const int j) {

        if(i < 0 || i >= pixel_h || j < 0 || j >= pixel_w) return Vec3(0,0,0);

        const R dx = Random::rando();
        const R dy = Random::rando();
        const Vec3 xp = position + sensor_center
                        + x * width_of_sensor * ((j + dx) / pixel_w - 0.5)
                        + y * height_of_sensor * ((i + dy) / pixel_h - 0.5);

        return xp;
    }

    void out_ppm_stdio() const {
        std::printf("P3\n%d %d\n255\n",pixel_w,pixel_h);
        //for(int i = 0;i < pixel_h * pixel_w;i++) 
            //img[i].print255();
        for(int i = pixel_h - 1;i >= 0;i--) 
            for(int j = pixel_w - 1;j >= 0;j--) 
                img[i * pixel_w + j].print255();
    }

    ~Camera() {
        delete[] img;
        delete[] img_l;
        delete[] img_e;
    }
};