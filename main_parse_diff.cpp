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

    string name = "mostro0_1_all";

    string file_out = "diff_"+name+".obj";
    string file_parts_to_color = "parts_to_color_"+name+".txt";
    string file_patches = name+".txt";

    //string file_mesh = "/Users/michele/Documents/GitHub/gitBooleans/cmake-build-release/diff_"+name+".obj";
    string file_mesh = "/Users/michele/Documents/GitHub/gitBooleans/cmake-build-debug/input_error.obj";
    DrawableTrimesh<> m(file_mesh.c_str());
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui,"Diff_test");

    gui.push(&m);
    gui.push(&controls);
    m.updateGL();

    std::vector<int> p_ids;
    std::vector<phmap::flat_hash_set<uint>> patches;
    parsePatches("/Users/michele/Documents/GitHub/gitBooleans/results/debug/patches_"+name+".txt", patches);
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
    Marker v0, v1, v2, r1;
    ray.thickness = 20.0;
    ray.use_gl_lines = true;
    gui.push(&ray);


    static int selected_item = 0;
    const char* items[] = { "Patch Analysis", "Show Triangle", "Pick Vertex",
                            "Diff Analysis", "Ray Analysis", "Push Point",
                            "Reset", "Pick Triangle", "Pick another triangle",
                            "Draw point"};
    const int items_count = sizeof(items) / sizeof(items[0]);

    gui.callback_app_controls = [&]()
    {
        if (ImGui::BeginCombo("Seleziona un'opzione", items[selected_item]))
        {
            for (int i = 0; i < items_count; i++)
            {
                bool is_selected = (selected_item == i);
                if (ImGui::Selectable(items[i], is_selected))
                {
                    selected_item = i; // Aggiorna l'elemento selezionato
                }
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        if(selected_item == 6){
            if(ImGui::Button("Reset")){
                for(uint i = 0; i < m.num_polys(); ++i){
                    m.poly_data(i).color = cinolib::Color(cinolib::Color::WHITE());
                }
                m.updateGL();
            }
        }
        if(selected_item == 3){
            if(ImGui::Button("Show Diff")){
                std::vector<std::vector<std::vector<unsigned int>>> parts_to_color_debug;
                std::string filename = "/Users/michele/Documents/GitHub/gitBooleans/results/debug/parts_to_color_"+name+".txt";

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
        }
        if(selected_item == 0){
            if(ImGui::Button("Next Patch"))
            {
                // button clicked: do something
                for(uint i = 0; i < m.num_polys(); ++i){
                    m.poly_data(i).color = cinolib::Color(255/255.f,255/255.f,255/255.f,0.9f);
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
        }

        if(selected_item == 1) {
            if (ImGui::Button("Show Triangle")) {
                m.poly_data(t_deb).color = cinolib::Color(cinolib::Color::YELLOW());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(t_deb);
                marker.disk_radius = 0.3f;
                marker.pos_3d = m.poly_centroid(t_deb);
                gui.push(marker);
                m.updateGL();
            }
        }

        if(selected_item == 4 ){
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
        }


            gui.callback_mouse_right_click = [&](int modifiers) -> bool {
                if (selected_item == 7){
                    if (modifiers & GLFW_MOD_ALT) {
                        vec3d p;
                        vec2d click = gui.cursor_pos();
                        if (gui.unproject(click, p)) // transform click in a 3d point
                        {
                            uint pid = m.pick_poly(p);

                            if (m.poly_data(pid).color == cinolib::Color::RED()) {
                                m.poly_data(pid).color = cinolib::Color::WHITE();
                                gui.pop_all_markers();

                            } else {
                                std::cout << "p_id: " << pid << std::endl;
                                m.poly_data(pid).color = cinolib::Color::RED();

                                v0.pos_3d = m.vert(m.poly_vert_id(pid, 0));
                                v0.text = "v0";
                                v0.disk_radius = 0.3f;
                                v0.color = Color::BLUE();

                                v1.pos_3d = m.vert(m.poly_vert_id(pid, 1));
                                v1.text = "v1";
                                v1.disk_radius = 0.3f;
                                v1.color = Color::BLACK();

                                v2.pos_3d = m.vert(m.poly_vert_id(pid, 2));
                                v2.text = "v2";
                                v2.disk_radius = 0.3f;
                                v2.color = Color::GREEN();

                                gui.push(v0);
                                gui.push(v1);
                                gui.push(v2);

                                m.vert_data(m.poly_vert_id(pid, 1)).color = Color::BLACK();
                                m.vert_data(m.poly_vert_id(pid, 2)).color = Color::GREEN();

                            }
                            m.updateGL();
                        }
                    }
                }

                if(selected_item == 8){
                    if(modifiers & GLFW_MOD_ALT)
                    {
                        vec3d p;
                        vec2d click = gui.cursor_pos();
                        if(gui.unproject(click, p)) // transform click in a 3d point
                        {
                            uint pid = m.pick_poly(p);

                            if(m.poly_data(pid).color == cinolib::Color::GREEN()){
                                m.poly_data(pid).color = cinolib::Color::WHITE();
                                gui.pop_all_markers();

                            }else{

                                m.poly_data(pid).color = cinolib::Color::GREEN();

                            }
                            m.updateGL();
                        }
                    }
                }
                if(selected_item == 9){
                    if(modifiers & GLFW_MOD_ALT){
                        vec3d p;
                        vec2d click = gui.cursor_pos();
                        if(gui.unproject(click, p)) // transform click in a 3d point
                        {
                            r1.pos_3d = p;
                            r1.text = "r1";

                            std::cout << "r1 coords : " << p.x() << " " << p.y() << " " << p.z() << std::endl;
                            gui.push(r1);
                        }
                    }
                }

                return false;
            };

    };






    /*gui.callback_mouse_right_click2 = [&](int modifiers) -> bool
    {
        if(modifiers & GLFW_MOD_ALT)
        {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p)) // transform click in a 3d point
            {
                uint pid = m.pick_poly(p);

                if(m.poly_data(pid).color == cinolib::Color::GREEN()){
                    m.poly_data(pid).color = cinolib::Color::WHITE();
                    gui.pop_all_markers();

                }else{

                    m.poly_data(pid).color = cinolib::Color::GREEN();

                }
                m.updateGL();
            }
        }
        return false;
    };*/

  /*  Marker r1, r2;
    gui.callback_mouse_left_click = [&](int modifiers) -> bool {

        if(modifiers & GLFW_MOD_ALT){
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p)) // transform click in a 3d point
            {
                r1.pos_3d = p;
                r1.text = "r1";

                std::cout << "r1 coords : " << p.x() << " " << p.y() << " " << p.z() << std::endl;
                gui.push(r1);
            }
        }


        return false;
    };*/



    //ray.push_seg(vec3d(ray_x_v0 * pow(10, exp_x_v0), ray_y_v0 * pow(10, exp_y_v0), ray_z_v0 * pow(10, exp_z_v0)),
    //vec3d(ray_x_v1 * pow(10, exp_x_v1), ray_y_v1 * pow(10, exp_y_v1), ray_z_v1 * pow(10, exp_z_v1)));
    // ray.draw();

    return gui.launch();

}