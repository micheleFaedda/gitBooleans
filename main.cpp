//
// Created by Michele on 30/05/24.
//

#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <arrangements/external/Indirect_Predicates/include/implicit_point.h>
#include <mesh_booleans/booleans.h>
#include <cstdlib>
#include <arrangements/code/processing.h>
#include <io_functions.h>

using namespace cinolib;
bool debug = true;
bool demo = false;
bool debug_impl = true;
bool patch_view = false;

int main(int argc, char **argv)
{
    string file_path;
    string file_path2;

    if(!debug && argc < 3) {
        std::cout << "syntax error!" << std::endl;
        std::cout << "./gitBooleans input1.obj input2.obj" << std::endl;
        return -1;

    }else if(!debug){
        file_path = argv[1];
        file_path2 = argv[2];
    }else{
        //file_path = "../data/mostro4.obj";
        file_path = "/Users/michele/Documents/GitHub/InteractiveAndRobustMeshBooleans/folder_test/Tinghi10K/112917_sf_a.obj";
        //file_path2 = "../data/mostro5.obj";
        file_path2 = "/Users/michele/Documents/GitHub/InteractiveAndRobustMeshBooleans/folder_test/mesh_rotated/112917_sf_a.obj";
    }
    string name = "mostro0_0mod";

    string file_out = "diff_"+name+".obj";
    string file_parts_to_color = "parts_to_color_"+name+".txt";
    string file_patches = "patches_"+name+".txt";

    vector<string> files = {file_path, file_path2};

    BoolOp op = UNION;

    vector<double> in_coords, bool_coords;
    vector<uint> in_tris, bool_tris;
    vector<uint> in_labels;
    vector<bitset<NBIT>> bool_labels;

    Profiler profiler;
    profiler.push("Total time");
    //start timer
    loadMultipleFiles(files, in_coords, in_tris, in_labels);

    cinolib::write_OBJ("diff_original.obj", in_coords, in_tris, {});

    //booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);

    initFPU();

    point_arena arena;
    std::vector<genericPoint*> arr_verts; // <- it contains the original expl verts + the new_impl verts
    std::vector<uint> arr_in_tris, arr_out_tris;
    std::vector<std::bitset<NBIT>> arr_in_labels;
    std::vector<DuplTriInfo> dupl_triangles;
    Labels labels;
    std::vector<phmap::flat_hash_set<uint>> patches;
    cinolib::Octree octree; // built with arr_in_tris and arr_in_labels

    customArrangementPipeline(in_coords, in_tris, in_labels, arr_in_tris, arr_in_labels, arena, arr_verts,
                              arr_out_tris, labels, octree, dupl_triangles);



    //customBooleanPipeline(arr_verts, arr_in_tris, arr_out_tris, arr_in_labels, dupl_triangles, labels, patches, octree, op, bool_coords, bool_tris, bool_labels);
    FastTrimesh tm(arr_verts, arr_out_tris, true);


    computeAllPatches(tm, labels, patches, true);

    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);


    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);

    computeInsideOutCustom(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);


    /******************************************************************************************************/
    //DIFF CODE
    std::vector<std::vector<std::vector<uint>>> parts_to_color;
    std::vector<int> num_tri_to_color_per_part;
    uint num_parts = debug_impl ? 4 : 3 ;
    num_parts = patch_view ? (num_parts + 1): num_parts;
    std::cout << "num_parts: " << num_parts << std::endl;

    parts_to_color.resize(num_parts);
    num_tri_to_color_per_part.resize(num_parts);

    tm.resetTrianglesInfo();
    uint num_tris_in_final_result = 0;

    vector<int> p_ids;

    for(uint t_id = 0 ; t_id < tm.numTris(); ++t_id){
        if(labels.surface[t_id][1]){ //if the triangle belong to B
            if(labels.surface[t_id].count() == 2 && labels.inside[t_id].count() == 0){
                //TODO: GRAY PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[0]++;
                for(uint p_id = 0; p_id < patches.size(); ++p_id){
                    if(patches[p_id].contains(t_id)){
                        p_ids.push_back(p_id);
                        break;
                    }
                }
                continue;

            }else if(labels.inside[t_id][0]){ //parts of B inside A
                //TODO: RED PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[1]++;
                for(uint p_id = 0; p_id < patches.size(); ++p_id){
                    if(patches[p_id].contains(t_id)){
                        p_ids.push_back(p_id);
                        break;
                    }
                }
                continue;

            }else if(!labels.inside[t_id][0]){
                //TODO: GREEN PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[2]++;
                for(uint p_id = 0; p_id < patches.size(); ++p_id){
                    if(patches[p_id].contains(t_id)){
                        p_ids.push_back(p_id);
                        break;
                    }
                }
                continue;
            }
        }else if(patch_view){
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_result++;
            num_tri_to_color_per_part[3]++;
        }

        if(debug_impl && labels.surface[t_id][0]){
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_result++;
            int pos = patch_view ? 4 : 3;
            num_tri_to_color_per_part[pos]++;
            for(uint p_id = 0; p_id < patches.size(); ++p_id){
                if(patches[p_id].contains(t_id)){
                    p_ids.push_back(p_id);
                    break;
                }
            }
        }
    }
    //remove duplicates in p_ids
    std::sort(p_ids.begin(), p_ids.end());
    p_ids.erase(std::unique(p_ids.begin(), p_ids.end()), p_ids.end());

    std::vector <uint> patches_debug_diff_tIds;
    patches_debug_diff_tIds.resize(num_tris_in_final_result);
    for(uint i = 0; i < tm.numTris(); i++){
        if(tm.triInfo(i) == 1){
            patches_debug_diff_tIds.push_back(i);
        }
    }

    /******************************************************************************************************/

    //computeFinalExplicitResult(tm, labels, num_tris_in_final_solution, bool_coords, bool_tris, bool_labels, true);

    bool_tris.resize(num_tris_in_final_result * 3);
    bool_labels.resize(num_tris_in_final_result);

    std::vector <uint> vert_to_color;
    vert_to_color.resize(3);
/*

    //reserve space for the parts to color and the vertices of each part
    for(uint i = 0; i < parts_to_color.size(); ++i){
        parts_to_color[i].reserve(num_tri_to_color_per_part[i]);

        for(uint j = 0; j < parts_to_color[i].size(); ++j)
            parts_to_color[i][j].reserve(3);
    }
*/
    std::vector<int>  vertex_index(tm.numVerts(), -1);
    uint num_vertices = 0;
    uint tri_offset = 0;

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(tm.triInfo(t_id) == 0) continue; // triangle not included in final version

        const uint *triangle = tm.tri(t_id);
        for(uint i = 0; i < 3; i++)
        {
            uint old_vertex = triangle[i];
            if (vertex_index[old_vertex] == -1) {
                vertex_index[old_vertex] = num_vertices++;
            }

            bool_tris[3 * tri_offset + i] = vertex_index[old_vertex];
            vert_to_color[i] = vertex_index[old_vertex];
        }
/*
        const uint t_id_aux = t_id;
        for(uint j = 0; j < patches.size(); ++j){

            if(patches[j].contains(t_id_aux)){
                std::cout<<"entrato" << std::endl;
                parts_to_color[j].push_back(vert_to_color);
                break;
            }
        }*/
        //parts_to_color[0].push_back(vert_to_color);

        if(labels.surface[t_id][1]){ //if the triangle belong to B
            if(labels.surface[t_id].count() == 2 && labels.inside[t_id].count() == 0){
                //TODO: GRAY PARTS
                parts_to_color[0].push_back(vert_to_color);


            }else if(labels.inside[t_id][0]){ //parts of B inside A
                //TODO: RED PARTS
                parts_to_color[1].push_back(vert_to_color);


            }else if(!labels.inside[t_id][0]){
                //TODO: GREEN PARTS
                parts_to_color[2].push_back(vert_to_color);
            }
        }else if(patch_view){

            parts_to_color[3].push_back(vert_to_color);
        }

        if(labels.surface[t_id][0] && debug_impl){
            int pos_app = patch_view ? 4 : 3;
            parts_to_color[pos_app].push_back(vert_to_color);
        }

        bool_labels[tri_offset] = labels.surface[t_id];
        tri_offset++;
    }

    // loop over vertices
    bool_coords.resize(3 * num_vertices);
    for(uint v_id = 0; v_id < (uint)tm.numVerts(); v_id++) {
        if (vertex_index[v_id] == -1) continue;
        double* v = bool_coords.data() + (3 * vertex_index[v_id]);
        tm.vert(v_id)->getApproxXYZCoordinates(v[0], v[1], v[2]);
    }

    // rescale output
    double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();
    for(double &c : bool_coords) c /= multiplier;
    /**********************************************************************************************/

    profiler.pop(); //end timer



    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});
    // Percorso dell'eseguibile e del file da passare come argomento
    const char* exe = "../mesh_booleans_inputcheck";

    // Costruzione del comando da eseguire
    std::string command = std::string(exe) + " " + file_out.c_str();

    // Esecuzione del comando
    int result = system(command.c_str());

    if(result != 0){
        std::cerr << "Error in the execution of the command" << std::endl;
    }


    //save the green, gray and red parts in a file to be used in the GUI

    //create the file with name parts_to_color +
    savePartsToFile(parts_to_color,
                   "/Users/michele/Documents/GitHub/gitBooleans/results/debug/"+file_parts_to_color,debug_impl, patch_view);

    if(patch_view)
        savePatchesTriangles("/Users/michele/Documents/GitHub/gitBooleans/results/debug/" + file_patches,
                         p_ids, patches);


    ///____________ GUI _____________________________________________________________________________________________
    GLcanvas gui;
    DrawableTrimesh<> m(bool_coords, bool_tris);
    SurfaceMeshControls<DrawableTrimesh<>> mesh_controls (&m, &gui,file_out.c_str());

    gui.push(&m);
    gui.push(&mesh_controls);
    Marker point;
    int point_id = -1;

    double step = 0.1;        // Small step
    double step_fast = 1.0;   // Large step
    const char* format = "%.6f"; // Format with six decimal places

    double x_pos = 0.0;
    double y_pos = 0.0;
    double z_pos = 0.0;

    int exp_x_v0 = 0;
    int exp_y_v0 = 0;
    int exp_z_v0 = 0;

    int t_deb;
    int p = -1;
    Marker v0, v1, v2;

    static int selected_item = 0;
    const char* items[] = { "Patch Analysis", "Show Triangle", "Pick Vertex", "Diff Analysis", "Ray Analysis", "Push Point", "Reset", "Pick Triangle" };
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

        if(selected_item == 0 && patch_view){
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

                uint v0 = bool_tris.at(t_id*3);
                uint v1 = bool_tris.at(t_id*3+1);
                uint v2 = bool_tris.at(t_id*3+2);

                std::cout << "Triangle: " << bool_tris.at(t_id*3) << " " << bool_tris.at(t_id*3+1) << " " << bool_tris.at(t_id*3+2) << std::endl;;
                m.poly_data(t_id).color = cinolib::Color(cinolib::Color::YELLOW());
                //}
            }
            if(p == patches.size()-1){
                p = -1;
            }
            m.updateGL();
            }
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
            for(uint i = 0; i < parts_to_color.size(); ++i) {
                for (uint j = 0; j < parts_to_color[i].size(); ++j) {

                    int t_id =  m.poly_id({parts_to_color[i][j][0], parts_to_color[i][j][1], parts_to_color[i][j][2]});
                    if (t_id == -1) std::cerr << "Error: triangle not found" << std::endl;

                    if(i == 0){ //WHITE
                            m.poly_data(t_id).color = Color(0.9f,0.9f,0.9f,0.8f);
                    }else if(i == 1){ //RED
                            m.poly_data(t_id).color = Color(201/255.f,79/255.f,86/255.f);
                    } else if(i == 2){  //GREEN
                        m.poly_data(t_id).color = Color(122/255.f,239/255.f,72/255.f);
                    } else if (i==3){ //WHITE TRANSPARENT
                        m.poly_data(t_id).color = Color(255/255.f,255/255.f,255/255.f,0.0f);
                    } else if (debug_impl && i==4) { //WHITE TRANSPARENT
                        m.poly_data(t_id).color = Color(255 / 255.f, 255 / 255.f, 255 / 255.f, 0.9f);
                    }
                }
            }
            m.updateGL();
        }
        }
        if(selected_item == 1){
        if(ImGui::InputInt("Triangle id", &t_deb, 0, m.num_verts())){
            for(uint i = 0; i < m.num_polys(); ++i){
                //m.poly_data(i).color = cinolib::Color(cinolib::Color::WHITE());
                m.poly_data(i).color = cinolib::Color(1.0f,1.0f,1.0f,0.0f);
            }
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
        if(selected_item == 5){
        if(ImGui::InputDouble("x vert", &x_pos, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp x vert", &exp_x_v0, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("y vert", &y_pos, step, step_fast, format)){

        }
        if(ImGui::InputInt("y vert", &exp_y_v0, 0, m.num_verts())){

        }
        if(ImGui::InputDouble("z vert", &z_pos, step, step_fast, format)){

        }
        if(ImGui::InputInt("Exp z vert", &exp_z_v0, 0, m.num_verts())){

        }

        if(ImGui::Button("Push point")){

            point.pos_3d = vec3d(x_pos * pow(10, exp_x_v0), y_pos * pow(10, exp_y_v0), z_pos * pow(10, exp_z_v0));
            point.text = std::to_string(point_id);
            point.disk_radius = 0.4f;
            point.color = Color::RED();
            gui.push(point);
            m.updateGL();

        }
    }

        gui.callback_mouse_right_click = [&](int modifiers) -> bool
        {
            if(selected_item == 7) {
                if (modifiers & GLFW_MOD_SHIFT) {
                    vec3d p;
                    vec2d click = gui.cursor_pos();
                    if (gui.unproject(click, p)) // transform click in a 3d point
                    {
                        uint pid = m.pick_poly(p);

                        if (m.poly_data(pid).color == cinolib::Color::RED()) {
                            m.poly_data(pid).color = cinolib::Color::WHITE();
                            gui.pop_all_markers();

                        } else {

                            m.poly_data(pid).color = cinolib::Color::RED();
                            explicitPoint3D v0_exp = tm.vert(m.poly_vert_id(pid, 0))->toExplicit3D();
                            explicitPoint3D v1_exp = tm.vert(m.poly_vert_id(pid, 1))->toExplicit3D();
                            explicitPoint3D v2_exp = tm.vert(m.poly_vert_id(pid, 2))->toExplicit3D();

                            std::cout << "v0: " << v0_exp.X() << " " << v0_exp.Y() << " " << v0_exp.Z() << std::endl;
                            std::cout << "v1: " << v1_exp.X() << " " << v1_exp.Y() << " " << v1_exp.Z() << std::endl;
                            std::cout << "v2: " << v2_exp.X() << " " << v2_exp.Y() << " " << v2_exp.Z() << std::endl;

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

                            tm.vert(m.poly_vert_id(pid, 1))->toExplicit3D();

                            m.vert_data(m.poly_vert_id(pid, 1)).color = Color::BLACK();
                            m.vert_data(m.poly_vert_id(pid, 2)).color = Color::GREEN();

                        }
                        m.updateGL();
                    }
                }
                return false;
            }
        };
    };



    m.updateGL();

    return gui.launch();
}

