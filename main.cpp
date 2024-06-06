//
// Created by Michele on 30/05/24.
//

#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include <mesh_booleans/booleans.h>


using namespace cinolib;
using namespace std;
bool debug = true;
bool flag_arrangement_debug = false;

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
        file_path = "../data/cube.obj";
        file_path2 = "../data/cube1_3.obj";
    }


    vector<string> files = {file_path, file_path2};

    BoolOp op = UNION;
    string file_out = "diff.obj";

    vector<double> in_coords, bool_coords;
    vector<uint> in_tris, bool_tris;
    vector<uint> in_labels;
    vector<bitset<NBIT>> bool_labels;


    loadMultipleFiles(files, in_coords, in_tris, in_labels);

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
    computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);

    /******************************************************************************************************/
    //DIFF CODE
    std::vector<std::vector<std::vector<uint>>> parts_to_color;
    std::vector<int> num_tri_to_color_per_part;
    uint num_parts = 3;
    parts_to_color.resize(num_parts);
    num_tri_to_color_per_part.resize(num_parts);

    tm.resetTrianglesInfo();

    uint num_tris_in_final_result = 0;

    for(uint t_id = 0 ; t_id < tm.numTris(); ++t_id){
        if(labels.surface[t_id][1]){ //if the triangle belong to B
            if(labels.surface[t_id].count() == 2 && labels.inside[t_id].count() == 0){
                //TODO: GRAY PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[0]++;

            }else if(labels.inside[t_id][0]){ //parts of B inside A
                //TODO: RED PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[1]++;

            }else if(!labels.inside[t_id][0]){
                //TODO: GREEN PARTS
                tm.setTriInfo(t_id, 1);
                num_tris_in_final_result++;
                num_tri_to_color_per_part[2]++;
            }
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

    //reserve space for the parts to color and the vertices of each part
    for(uint i = 0; i < parts_to_color.size(); ++i){
        parts_to_color[i].reserve(num_tri_to_color_per_part[i]);

        for(uint j = 0; j < parts_to_color[i].size(); ++j)
            parts_to_color[i][j].reserve(3);
    }

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
    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});

    GLcanvas gui;
    DrawableTrimesh<>new_mesh;

    DrawableTrimesh<> m(file_out.c_str());
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui);

    gui.push(&m);
    gui.push(&controls);

    for(uint i = 0; i < parts_to_color.size(); ++i) {
        for (uint j = 0; j < parts_to_color[i].size(); ++j) {

            int t_id =  m.poly_id({parts_to_color[i][j][0], parts_to_color[i][j][1], parts_to_color[i][j][2]});
            if (t_id == -1) std::cerr << "Error: triangle not found" << std::endl;

            if(i == 0){
                m.poly_data(t_id).color = Color::GRAY(2.0);
            }else if(i == 1){
                m.poly_data(t_id).color = Color::RED();
            }  else if(i == 2){
                m.poly_data(t_id).color = Color::GREEN();
            }
        }
    }

    m.updateGL();
    return gui.launch();
}