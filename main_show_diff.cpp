//
// Created by Michele on 03/06/24.
//


#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;

int main(int argc, char **argv) {


    std::string mantained_parts = "../results/grey_parts.obj";
    std::string deleted_parts = "../results/red_parts.obj";
    std::string added_parts = "../results/green_parts_inside.obj";

    GLcanvas gui;
    DrawableTrimesh<>new_mesh;

    DrawableTrimesh<> m(mantained_parts.c_str());
    DrawableTrimesh<> m_rem(deleted_parts.c_str());
    m_rem.poly_set_color(Color::RED());


    gui.push(&m);
    gui.push(&m_rem);

    DrawableTrimesh<> m2(mantained_parts.c_str());
    DrawableTrimesh<> m2_add(added_parts.c_str());
    m2_add.poly_set_color(Color::GREEN());
    m2.translate(vec3d(0.7,0,0));
    m2_add.translate(vec3d(0.7,0,0));


    gui.push(&m2);
    gui.push(&m2_add);


    SurfaceMeshControls<DrawableTrimesh<>> controls_m(&m, &gui,"m_mantained");
    SurfaceMeshControls<DrawableTrimesh<>> controls_m_rem(&m_rem, &gui,"m_removed");
    SurfaceMeshControls<DrawableTrimesh<>> controls_m2(&m2, &gui,"m2_mantained");
    SurfaceMeshControls<DrawableTrimesh<>> controls_m2_add(&m2_add, &gui, "m2_added");

    gui.push(&controls_m);
    gui.push(&controls_m_rem);
    gui.push(&controls_m2);
    gui.push(&controls_m2_add);

    m.updateGL();
    m_rem.updateGL();
    m2.updateGL();
    m2_add.updateGL();

    return gui.launch();

}