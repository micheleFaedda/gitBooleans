//
// Created by Michele on 30/05/24.
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
bool debug_impl = false;

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
        //file_path = "../data/test/Horse/Horse_conv/horse0.obj";
        //file_path = "../data/t0.obj";
        file_path = "../data/mostro0.obj";
        //file_path = "../data/test/cube.obj";
        //file_path2 = "../data/test/Horse/Horse_conv/horse1.obj";
        file_path2 = "../data/mostro1.obj";
        //file_path2 = "../data/t1.obj";
        //file_path2 = "../data/test/pyramid_transform.obj";

    }
    string name = "mostro0_1";

    string file_out = "diff_cube_"+name+".obj";
    string file_parts_to_color = "parts_to_color_"+name+".txt";
    string file_patches = name+".txt";

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

    cinolib::write_OBJ("input_error.obj", in_coords, in_tris, {});

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
    FastTrimesh tm(arr_verts, arr_out_tris, false);


    computeAllPatches(tm, labels, patches, false);

    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);
    computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);


    /******************************************************************************************************/
    //DIFF CODE
    std::vector<std::vector<std::vector<uint>>> parts_to_color;
    std::vector<int> num_tri_to_color_per_part;
    uint num_parts = debug_impl ? 4 : 3 ;
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
        }

        if(debug_impl && labels.surface[t_id][0]){
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_result++;
            num_tri_to_color_per_part[3]++;
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

    //computeFinalExplicitResult(tm, labels, num_tris_in_final_solution, bool_coords, bool_tris, bool_labels, true);c

    uint num_vertices = 0;
    std::vector<int>  vertex_index(tm.numVerts(), -1);
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
        }

        if(labels.surface[t_id][0] && debug_impl){
            parts_to_color[3].push_back(vert_to_color);
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

    //save the green, gray and red parts in a file to be used in the GUI

    //create the file with name parts_to_color +
    savePartsToFile(parts_to_color,
                    "/Users/michele/Documents/GitHub/gitBooleans/results/debug/parts_to_color_" + file_parts_to_color,debug_impl);

    savePatchesTriangles("/Users/michele/Documents/GitHub/gitBooleans/results/debug/patches_" + file_patches,
                         p_ids, patches);

    GLcanvas gui;
    DrawableTrimesh<> m(bool_coords, bool_tris);
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui,file_out.c_str());

    gui.push(&m);
    gui.push(&controls);


    int t_deb;
    int p = -1;

    gui.callback_app_controls = [&]()
    {
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

            for(uint t_id : patches.at(p_ids.at(p))){

                uint v0 = bool_tris.at(t_id*3);
                uint v1 = bool_tris.at(t_id*3+1);
                uint v2 = bool_tris.at(t_id*3+2);

                std::cout << "Triangle: " << bool_tris.at(t_id*3) << " " << bool_tris.at(t_id*3+1) << " " << bool_tris.at(t_id*3+2) << std::endl;
                int p_id = m.poly_id({v0, v1, v2});
                m.poly_data(p_id).color = cinolib::Color(cinolib::Color::YELLOW());
                //}
            }
            if(p == patches.size()-1){
                p = -1;
            }
            m.updateGL();
        }
        if(ImGui::Button("Reset")){
            for(uint i = 0; i < m.num_polys(); ++i){
                m.poly_data(i).color = cinolib::Color(cinolib::Color::WHITE());
            }
            m.updateGL();
        }
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
                    } else if (debug_impl && i==3){ //WHITE TRANSPARENT
                        m.poly_data(t_id).color = Color(255/255.f,255/255.f,255/255.f,0.9f);
                    }
                }
            }
            m.updateGL();
        }
        if(ImGui::InputInt("Text Input", &t_deb, 0, m.num_verts())){
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
    };

    m.updateGL();

    return gui.launch();
}

