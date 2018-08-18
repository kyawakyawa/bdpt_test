#include <cmath>
#include <iostream>
#include "Bdpt.hpp"
#include "Sphere.hpp"

int main() {
    Camera camera(320,240,Vec3(50.0, 40.8, 190.0),Vec3(0,0,-1),Vec3(0,1,0),2 * std::atan(2.0 / 3) * 180.0 / 3.1415926535,140,42.0/17);
    /*std::cout << camera.sensor_center + camera.position << std::endl;
    std::cout << camera.height_of_sensor << std::endl;
    std::cout << camera.dist_sensor_lens << std::endl;
    std::cout << camera.dist_lens_object_plane << std::endl;
    std::cout << camera.lens_radius << std::endl;*/

    //Sphere(7.5,Vec(50.0, 72.5, 81.6),    Color(16,16,16), Color(),              REFLECTION_TYPE_DIFFUSE), //照明
	//Sphere(1e5, Vec( 1e5+1, 40.8, 81.6), Color(),      Color(0.75, 0.25, 0.25), REFLECTION_TYPE_DIFFUSE), // 左
	//Sphere(1e5, Vec(-1e5+99, 40.8, 81.6),Color(),      Color(0.25, 0.25, 0.75), REFLECTION_TYPE_DIFFUSE), // 右
	//Sphere(1e5, Vec(50, 40.8, 1e5),      Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE), // 奥
	//Sphere(1e5, Vec(50, 40.8, -1e5+250), Color(),      Color(),                 REFLECTION_TYPE_DIFFUSE), // 手前
	//Sphere(1e5, Vec(50, 1e5, 81.6),      Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE), // 床
	//Sphere(1e5, Vec(50, -1e5+81.6, 81.6),Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE), // 天井
	//Sphere(20,Vec(50, 20, 50),           Color(),      Color(0.25, 0.75, 0.25), REFLECTION_TYPE_DIFFUSE), // 緑球

    Scene scene(&camera);

    scene.add_as_light(new Sphere(Vec3(50.0,72.5,81.6),7.5,Material(FColor(0,0,0),FColor(16,16,16))));
    scene.add(new Sphere(Vec3( 1e5+1, 40.8, 81.6),1e5,Material(FColor(0.75,0.25,0.25))));
    scene.add(new Sphere(Vec3( -1e5+99, 40.8, 81.6),1e5,Material(FColor(0.25,0.25,0.75))));
    scene.add(new Sphere(Vec3( 50, 40.8, 1e5),1e5,Material(FColor(0.75,0.75,0.75))));
    scene.add(new Sphere(Vec3( 50, 40.8, -1e5+250),1e5,Material(FColor(0,0,0))));
    scene.add(new Sphere(Vec3( 50, 1e5, 81.6),1e5,Material(FColor(0.75,0.75,0.75))));
    scene.add(new Sphere(Vec3( 50, -1e5+81.6, 81.6),1e5,Material(FColor(0.75,0.75,0.75))));
    scene.add(new Sphere(Vec3( 50, 20 ,50),20,Material(FColor(0.25,0.75,0.25))));

    Bdpt bdpt;

    bdpt.bdpt(scene,16);

    scene.camera->out_ppm_stdio();

    return 0;


}