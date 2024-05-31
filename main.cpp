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

int main(int argc, char **argv)
{

    string file_path = "../data/armadillo_T0.obj";
    string file_path2 = "../data/armadillo_T1.obj";

    vector<string> files = {file_path, file_path2};




    BoolOp op = SUBTRACTION;
    string file_out = "../results/difference.obj";

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




    customBooleanPipeline(arr_verts, arr_in_tris, arr_out_tris, arr_in_labels, dupl_triangles, labels,
                          patches, octree, op, bool_coords, bool_tris, bool_labels);



    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});


    GLcanvas gui;
    DrawableTrimesh<>new_mesh;

    DrawableTrimesh<> m(file_out.c_str());
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui);

    gui.push(&m);
    gui.push(&controls);

    m.updateGL();

    return gui.launch();

    return 0;
}