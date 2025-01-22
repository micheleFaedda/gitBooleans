//
// Created by Michele on 20/01/25.
//
#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include <mesh_booleans/booleans.h>
#include <arrangements/code/processing.h>


using namespace cinolib;
using namespace std;
bool debug = true;
bool demo = false;
bool debug_impl = true;
bool flag_arrangement_debug = false;


int main(int argc, char **argv)
{

    GLcanvas gui;
    DrawableTrimesh<>new_mesh;

    DrawableTrimesh<> m("/Users/michele/Documents/GitHub/gitBooleans/cmake-build-debug/diff_t0_t1.obj");
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui,"Diff_test");

    gui.push(&m);
    gui.push(&controls);
    m.updateGL();

    std::vector<int> p_ids;
    std::vector<phmap::flat_hash_set<uint>> patches;
    parsePatches("/Users/michele/Documents/GitHub/gitBooleans/results/debug/patches_t0_t1.txt", patches);
    int p = -1;


    int t_deb = 0;
    double ray_x_v0 = 0.0;
    int exp_x_v0 = 0;
    double ray_y_v0 = 0.0;
    int exp_y_v0 = 0;
    double ray_z_v0 = 0.0;
    int exp_z_v0 = 0;

    double ray_x_v1 = 0.0;
    int exp_x_v1 = 0;
    double ray_y_v1 = 0.0;
    int exp_y_v1 = 0;
    double ray_z_v1 = 0.0;
    int exp_z_v1 = 0;

    double step = 0.1;        // Small step
    double step_fast = 1.0;   // Large step
    const char* format = "%.6f"; // Format with six decimal places
    bool find = false;
    DrawableSegmentSoup ray;
    ray.thickness = 20.0;
    ray.use_gl_lines = true;
    gui.push(&ray);

    gui.callback_app_controls = [&]()
    {

        if(ImGui::Button("Reset")){
            for(uint i = 0; i < m.num_polys(); ++i){
                m.poly_data(i).color = cinolib::Color(cinolib::Color::WHITE());
            }
            m.updateGL();
        }
        if(ImGui::Button("Show Diff")){
            std::vector<std::vector<std::vector<unsigned int>>> parts_to_color_debug;
            std::string filename = "/Users/michele/Documents/GitHub/gitBooleans/results/debug/parts_to_color_t0_t1.txt";

            parseFileToParts(filename, parts_to_color_debug, debug_impl);

            for(uint i = 0; i < parts_to_color_debug.size(); ++i) {
                for (uint j = 0; j < parts_to_color_debug[i].size(); ++j) {

                    int t_id =  m.poly_id({parts_to_color_debug[i][j][0], parts_to_color_debug[i][j][1], parts_to_color_debug[i][j][2]});
                    if (t_id == -1) std::cerr << "Error: triangle not found" << std::endl;

                    if(i == 0){ //WHITE
                        m.poly_data(t_id).color = Color(0.9f,0.9f,0.9f,0.8f);
                    }else if(i == 1){ //RED
                        m.poly_data(t_id).color = Color(201/255.f,79/255.f,86/255.f);
                    } else if(i == 2){  //GREEN
                        m.poly_data(t_id).color = Color(122/255.f,239/255.f,72/255.f);
                    } else if (debug_impl && i==3){ //WHITE TRANSPARENT
                        m.poly_data(t_id).color = Color(255/255.f,255/255.f,255/255.f,0.9f);
                    }
                }
            }
            m.updateGL();
        }
        if(ImGui::Button("Next Patch"))
        {
            // button clicked: do something
            for(uint i = 0; i < m.num_polys(); ++i){
                m.poly_data(i).color = cinolib::Color(cinolib::Color::WHITE());
            }
            std::cout<< "p : "<< p << std::endl;
            p++;

            //print the triangles id into patches_debug_diff_tIds
            for(auto print : p_ids){
                std::cout << "Triangle patch: " << print << std::endl;
            }

            for(uint t_id : patches.at(p)){

               /* uint v0 = bool_tris.at(t_id*3);
                uint v1 = bool_tris.at(t_id*3+1);
                uint v2 = bool_tris.at(t_id*3+2);

                std::cout << "Triangle: " << bool_tris.at(t_id*3) << " " << bool_tris.at(t_id*3+1) << " " << bool_tris.at(t_id*3+2) << std::endl;
                int p_id = m.poly_id({v0, v1, v2});
                */m.poly_data(t_id).color = cinolib::Color(cinolib::Color::YELLOW());
                //}
            }
            if(p == patches.size()-1){
                p = -1;
            }
            m.updateGL();
        }


        if(ImGui::Button("Show Triangle")){
            m.poly_data(t_deb).color = cinolib::Color(cinolib::Color::YELLOW());
            cinolib::Marker marker;
            marker.color = cinolib::Color(cinolib::Color::RED());
            marker.text = std::to_string(t_deb);
            marker.disk_radius = 0.3f;
            marker.pos_3d = m.poly_centroid(t_deb);
            gui.push(marker);
            m.updateGL();
        }

        if(ImGui::InputDouble("Ray x v0", &ray_x_v0, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp x v0", &exp_x_v0, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("Ray y v0", &ray_y_v0, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp y v0", &exp_y_v0, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("Ray z v0", &ray_z_v0, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp z v0", &exp_z_v0, 0, m.num_verts())){

        }

        if(ImGui::InputDouble("Ray x v1", &ray_x_v1, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp x v1", &exp_x_v1, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("Ray y v1", &ray_y_v1, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp y v1", &exp_y_v1, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("Ray z v1", &ray_z_v1, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp z v1", &exp_z_v1, 0, m.num_verts())){

        }

        if(ImGui::Button("Show ray")){
            vec3d v0(6.26349e+09, 5.36871e+08, 4.47392e+09);
            vec3d v1(6.26349e+09, 5.70425e+09 , 4.47392e+09);

            vec3d v2(0,0,0);
            ray.push_seg(v0, v2);
            //ray.push_seg(vec3d(ray_x_v0 * pow(10, exp_x_v0), ray_y_v0 * pow(10, exp_y_v0), ray_z_v0 * pow(10, exp_z_v0)),
                         //vec3d(ray_x_v1 * pow(10, exp_x_v1), ray_y_v1 * pow(10, exp_y_v1), ray_z_v1 * pow(10, exp_z_v1)));
           // ray.draw();
           gui.push_marker(v0,"v0", Color::RED(), 50);
           gui.push_marker(v1,"v1", Color::RED(), 50);

        }


    };
;
    //ray.push_seg(vec3d(ray_x_v0 * pow(10, exp_x_v0), ray_y_v0 * pow(10, exp_y_v0), ray_z_v0 * pow(10, exp_z_v0)),
    //vec3d(ray_x_v1 * pow(10, exp_x_v1), ray_y_v1 * pow(10, exp_y_v1), ray_z_v1 * pow(10, exp_z_v1)));
    // ray.draw();

    return gui.launch();

}