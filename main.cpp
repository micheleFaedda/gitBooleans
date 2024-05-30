//
// Created by Michele on 30/05/24.
//

#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include <mesh_booleans/booleans.h>


using namespace cinolib;

int main(int argc, char **argv)
{
    std::string file_path = "../data/armadillo_T0.obj";
    GLcanvas gui;
    DrawableTrimesh<> m(file_path.c_str());
    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui);


    gui.push(&m);
    gui.push(&controls);

    /*gui.callback_app_controls = [&]()
    {
        if(ImGui::Button("show me vert")) {

        }
    }*/
    m.updateGL();
    gui.push(&m);




    return gui.launch();
}