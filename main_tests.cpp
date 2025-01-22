//
// Created by Michele on 03/06/24.
//


#include <iostream>
#include <limits>
#include <cinolib/rationals.h>
#include <implicit_point.h>
#include "code/intersect_custom.h"
#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/drawable_triangle_soup.h>
using namespace cinolib;
void savePartsToFile(const std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                     const std::string& filename,
                     bool debug_impl);
void savePartsToFile(const std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                     const std::string& filename,
                     bool debug_impl) {
    // Open a file for writing
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Define the color labels
    std::vector<std::string> labels = {"GRAY Parts", "RED Parts", "GREEN Parts", "WHITE Parts"};

    // Adjust the number of parts to save
    size_t num_parts = debug_impl ? 4 : 3;

    // Iterate through the parts
    for (size_t i = 0; i < num_parts; ++i) {
        // Write the label
        if (i < labels.size()) {
            outFile << labels[i] << "\n";
        } else {
            outFile << "UNKNOWN Parts\n";
        }

        // Write the indices of vertices
        for (const auto& vertex : parts_to_color[i]) {
            for (unsigned int index : vertex) {
                outFile << index << " ";
            }
            outFile << "\n";
        }
    }

    // Close the file
    outFile.close();
}
bool parseFileToParts(const std::string& filename,
                      std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                      bool debug_impl);
bool parseFileToParts(const std::string& filename,
                      std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                      bool debug_impl) {
    // Open the file for reading
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    // Define labels and corresponding indices
    std::unordered_map<std::string, size_t> labelToIndex = {
            {"GRAY Parts", 0},
            {"RED Parts", 1},
            {"GREEN Parts", 2},
            {"WHITE Parts", 3} // New label for debug_impl
    };

    // Adjust the parts_to_color size
    size_t num_parts = debug_impl ? 4 : 3;
    parts_to_color.resize(num_parts);

    // Parsing state
    std::string line;
    size_t currentPart = 0;

    // Parse the file line by line
    while (std::getline(inFile, line)) {
        // Check if the line is a label
        if (labelToIndex.find(line) != labelToIndex.end()) {
            currentPart = labelToIndex[line];
        } else {
            // Parse indices in the current line
            std::istringstream iss(line);
            std::vector<unsigned int> vertices;
            unsigned int vertex;
            while (iss >> vertex) {
                vertices.push_back(vertex);
            }
            // Add parsed vertices to the correct part
            if (currentPart < parts_to_color.size()) {
                parts_to_color[currentPart].push_back(vertices);
            }
        }
    }

    inFile.close();
    return true;
}


int main(int argc, char **argv) {

    std::vector<std::vector<std::vector<unsigned int>>> parts_to_color = {
            {{0, 1, 2}, {3, 4, 5}}, // GRAY Parts
            {{6, 7, 8}},            // RED Parts
            {{9, 10, 11}, {12, 13}},// GREEN Parts
            {{14, 15}}              // Additional for debug_impl
    };

    // Save to a file with debug_impl = false
    savePartsToFile(parts_to_color, "output_normal.txt", false);

    // Save to a file with debug_impl = true
    savePartsToFile(parts_to_color, "output_debug.txt", true);




    // Data structure to store the parsed parts
    std::vector<std::vector<std::vector<unsigned int>>> parts_to_color_debug;

    // Filename to parse
    std::string filename = "/Users/michele/Documents/GitHub/gitBooleans/cmake-build-debug/";

    // Parse the file with debug mode enabled
    if (parseFileToParts(filename, parts_to_color_debug, true)) {
        // Print the parsed data for verification
        for (size_t i = 0; i < parts_to_color_debug.size(); ++i) {
            std::cout << "Part " << i << ":\n";
            for (const auto& vertexGroup : parts_to_color_debug[i]) {
                for (unsigned int vertex : vertexGroup) {
                    std::cout << vertex << " ";
                }
                std::cout << "\n";
            }
        }
    } else {
        std::cerr << "Failed to parse the file.\n";
    }
    /*bigrational zero = bigrational(0,0,0)

    const uint64_t one = 412425389222715451;
    const uint64_t two = 316268131472977041;
    const bignatural one_natural = bignatural(one);
    const bignatural two_natural = bignatural(two);
    bigrational comp = bigrational(one_natural, two_natural,0);
    bool is_zero = comp > bigrational(0);
    std::cout << "Is zero: " << is_zero << std::endl;


    bigrational Xmin_s1 = std::max(zero, zero);
    std::cout << "Xmin_s1: " << Xmin_s1 << std::endl;
*/


    //std::vector <uint> faces;
    //std::vector <double> coords;

    /*
    //triangle t_id = 32
    const bigrational x0(3);
    const bigrational y0(0,0,0);
    const bigrational z0(3);

    const bigrational x1(0,0,0);
    const bigrational y1(0,0,0);
    const bigrational z1(3);

    const bigrational x2(2999631764189355,2251799813685248, 1);
    const bigrational y2(0,0,0);
    const bigrational z2(8088317212689633,4503599627370496,1);

    const std::vector<bigrational> t0 = {x0, y0, z0};
    const std::vector<bigrational> t1 = {x1, y1, z1};
    const std::vector<bigrational> t2 = {x2, y2, z2};

    const bigrational x0_ray (17826409253836285,9007199254740992, 1);
    const bigrational y0_ray (2,3,1);
    const bigrational z0_ray (36017237179440431,54043195528445952,1);

    const double x0_d = -7.70753e+09;

    std::cout << "x0_ray: " << x0_ray << std::endl;

    const bigrational x1_ray (17826409253836285,9007199254740992, 1);
    const bigrational y1_ray (0,0,0);
    const bigrational z1_ray (17,4,1);


    const std::vector<bigrational> ray0 = {x0_ray, y0_ray, z0_ray};
    const std::vector<bigrational> ray1 = {x1_ray, y1_ray, z1_ray};

    std::vector <double> coords;
    coords.push_back(x0.get_d());
    coords.push_back(y0.get_d());
    coords.push_back(z0.get_d());

    coords.push_back(x1.get_d());
    coords.push_back(y1.get_d());
    coords.push_back(z1.get_d());

    coords.push_back(x2.get_d());
    coords.push_back(y2.get_d());
    coords.push_back(z2.get_d());


    faces.push_back(0);
    faces.push_back(1);
    faces.push_back(2);

    std::vector <double> ray_coords;
    ray_coords.push_back(x0_ray.get_d());
    ray_coords.push_back(y0_ray.get_d());
    ray_coords.push_back(z0_ray.get_d());

    ray_coords.push_back(x1_ray.get_d());
    ray_coords.push_back(y1_ray.get_d());
    ray_coords.push_back(z1_ray.get_d());



    //triangle t_id = 2
    const bigrational x0_2(0,0,0);
    const bigrational y0_2(0,0,0);
    const bigrational z0_2(3);

    const bigrational x1_2(3);
    const bigrational y1_2(0,0,0);
    const bigrational z1_2(3);

    const bigrational x2_2(3);
    const bigrational y2_2(1);
    const bigrational z2_2(3);

    const std::vector<bigrational> t0_2 = {x0_2, y0_2, z0_2};
    const std::vector<bigrational> t1_2 = {x1_2, y1_2, z1_2};
    const std::vector<bigrational> t2_2 = {x2_2, y2_2, z2_2};


    coords.push_back(x0_2.get_d());
    coords.push_back(y0_2.get_d());
    coords.push_back(z0_2.get_d());

    coords.push_back(x1_2.get_d());
    coords.push_back(y1_2.get_d());
    coords.push_back(z1_2.get_d());

    coords.push_back(x2_2.get_d());
    coords.push_back(y2_2.get_d());
    coords.push_back(z2_2.get_d());

    faces.push_back(3);
    faces.push_back(4);
    faces.push_back(5);

    bool intersection_t0 = segment_triangle_intersect_3d(&ray0[0], &ray1[0], &t0[0], &t1[0], &t2[0]);
    std::cout << "Intersection t0: " << intersection_t0 << std::endl;

    bool intersection_segment_32 = segment_segment_intersect_3d(&ray0[0], &ray1[0], &t0[0], &t1[0]);
    bool intersection_segment_2 = segment_segment_intersect_3d(&ray0[0], &ray1[0], &t0_2[0], &t1_2[0]);
    std::cout << "Intersection segment 32: " << intersection_segment_32 << std::endl;
    std::cout << "Intersection segment 2: " << intersection_segment_2 << std::endl;

    bool intersection_t1 = segment_triangle_intersect_3d(&ray0[0], &ray1[0], &t0_2[0], &t1_2[0], &t2_2[0]);
    std::cout << "Intersection t1: " << intersection_t1 << std::endl;
*/
    //triangle t_id = 1863
    /*const double x0_t0 = -5.58796e+09;
    const double y0_t0 = -6.51274e+09;
    const double z0_t0 = 7.42589e+09;

    const double x1_t0 = -5.1519e+09;
    const double y1_t0 = -6.67568e+09;
    const double z1_t0 = 6.9899e+09;

    const double x2_t0 = -5.58113e+09;
    const double y2_t0 = -6.51105e+09 ;
    const double z2_t0 = 7.42747e+09;

    const std::vector<double> t0 = {x0_t0, y0_t0, z0_t0};
    const std::vector<double> t1 = {x1_t0, y1_t0, z1_t0};
    const std::vector<double> t2 = {x2_t0, y2_t0, z2_t0};

    //triangle t_id = 2148
    const double x0_t1 = -5.1519e+09;
    const double y0_t1 = -6.67568e+09;
    const double z0_t1 = 6.9899e+09;

    const double x1_t1 = -5.58113e+09;
    const double y1_t1 = -6.51105e+09;
    const double z1_t1 = 7.42747e+09;

    const double x2_t1 = -5.58796e+09;
    const double y2_t1 = -6.51274e+09 ;
    const double z2_t1 = 7.42589e+09;

    const std::vector<double> t0_2 = {x0_t1, y0_t1, z0_t1};
    const std::vector<double> t1_2 = {x1_t1, y1_t1, z1_t1};
    const std::vector<double> t2_2 = {x2_t1, y2_t1, z2_t1};

    //triangle t_id = 2144
    const double x0_t2 = -6.05799e+09 ;
    const double y0_t2 =  6.32808e+09;
    const double z0_t2 =  7.06737e+09;

    const double x1_t2 = -6.40978e+09;
    const double y1_t2 = 6.04367e+09;
    const double z1_t2 = 7.90116e+09;

    const double x2_t2 = -5.12595e+09;
    const double y2_t2 = 6.50793e+09;
    const double z2_t2 = 7.17377e+09;

    const std::vector<double> t0_3 = {x0_t2, y0_t2, z0_t2};
    const std::vector<double> t1_3 = {x1_t2, y1_t2, z1_t2};
    const std::vector<double> t2_3 = {x2_t2, y2_t2, z2_t2};


    //ray
    const double x0_ray = -5.44033e+09;
    const double y0_ray = 6.56649e+09;
    const double z0_ray = 7.28109e+09;

    const double x1_ray = -5.44033e+09;
    const double y1_ray = 1.03925e+10;
    const double z1_ray = 7.28109e+09;

    const std::vector<double> ray0 = {x0_ray, y0_ray, z0_ray};
    const std::vector<double> ray1 = {x1_ray, y1_ray, z1_ray};


    int n_triangles = 3;

    faces.push_back(0);
    faces.push_back(1);
    faces.push_back(2);

    faces.push_back(3);
    faces.push_back(4);
    faces.push_back(5);

    faces.push_back(6);
    faces.push_back(7);
    faces.push_back(8);






    coords.push_back(x0_t0);
    coords.push_back(y0_t0);
    coords.push_back(z0_t0);

    coords.push_back(x1_t0);
    coords.push_back(y1_t0);
    coords.push_back(z1_t0);

    coords.push_back(x2_t0);
    coords.push_back(y2_t0);
    coords.push_back(z2_t0);

    coords.push_back(x0_t1);
    coords.push_back(y0_t1);
    coords.push_back(z0_t1);

    coords.push_back(x1_t1);
    coords.push_back(y1_t1);
    coords.push_back(z1_t1);

    coords.push_back(x2_t1);
    coords.push_back(y2_t1);
    coords.push_back(z2_t1);

    coords.push_back(x0_t2);
    coords.push_back(y0_t2);
    coords.push_back(z0_t2);

    coords.push_back(x1_t2);
    coords.push_back(y1_t2);
    coords.push_back(z1_t2);

    coords.push_back(x2_t2);
    coords.push_back(y2_t2);
    coords.push_back(z2_t2);

    GLcanvas gui;
    DrawableTrimesh<>m(coords, faces);
    m.scale(100000);
    //cinolib::DrawableTriangleSoup soup(coords, faces, Color::RED());


    SurfaceMeshControls<DrawableTrimesh<>> controls(&m, &gui, "test");

    gui.push(&m);
    gui.push(&controls);
    //gui.push(&soup);

    m.updateGL();

    DrawableSegmentSoup segment;
    Color customColor = Color::BLUE();

    segment.push_seg(vec3d(ray0[0], ray0[1], ray0[2]), vec3d(ray1[0], ray1[1], ray1[2]), customColor);
    segment.draw(100);
    gui.push(&segment);


    return gui.launch();
    //plane
    const explicitPoint3D ir(0.0,0.0,0.0);
    const explicitPoint3D is(4.0,4.0,0.0);
    const explicitPoint3D it(0.0,4.0,0.0);

    //line to generate P0
    const explicitPoint3D iq_P0(0.5,0.5,-1.0);
    const explicitPoint3D ip_P0(0.5,0.5,1.0);

    //line to generate P1
    const explicitPoint3D iq_P1(0.5,1,-1.0);
    const explicitPoint3D ip_P1(0.5,1,1.0);

    //line to generate P2
    const explicitPoint3D iq_P2 (1.0,1.0,-1.0);
    const explicitPoint3D ip_P2(1.0,1.0,1.0);

    //line to generate P3
    const explicitPoint3D iq_P3(1.0,0.5,-1.0);
    const explicitPoint3D ip_P3(1.0,0.5,1.0);

    //create 4 implicit points
    genericPoint *P0 = new implicitPoint3D_LPI(iq_P0, ip_P0, ir, is, it);
    genericPoint *P1 = new implicitPoint3D_LPI(iq_P1, ip_P1, ir, is, it);
    genericPoint *P2 = new implicitPoint3D_LPI(iq_P2, ip_P2, ir, is, it);
    genericPoint *P3 = new implicitPoint3D_LPI(iq_P3, ip_P3, ir, is, it);
*/
    /************************* ORIENT3D IMPLICIT **********************/
   /* int T0_orient = genericPoint::orient3D(*P0, ir, it, is);
    int T1_orient = genericPoint::orient3D(*P1, ir, it, is);
    int T2_orient = genericPoint::orient3D(*P2, ir, it, is);
    int T3_orient = genericPoint::orient3D(*P3, ir, it, is);

    std::cout << "Result P0 Implicit: " << T0_orient << std::endl;
    std::cout << "Result P1 Implicit: " << T1_orient << std::endl;
    std::cout << "Result P2 Implicit: " << T2_orient << std::endl;
    std::cout << "Result P3 Implicit: " << T3_orient << std::endl;
*/
    /*Rationals*/
/*    bigrational x_P0, y_P0, z_P0;
    P0->getExactXYZCoordinates(x_P0,y_P0,z_P0);

    bigrational x_P1, y_P1, z_P1;
    P1->getExactXYZCoordinates(x_P1,y_P1,z_P1);

    bigrational x_P2, y_P2, z_P2;
    P2->getExactXYZCoordinates(x_P2,y_P2,z_P2);

    bigrational x_P3, y_P3, z_P3;
    P3->getExactXYZCoordinates(x_P3,y_P3,z_P3);

    bigrational x_ir, y_ir, z_ir;
    ir.getExactXYZCoordinates(x_ir,y_ir,z_ir);

    bigrational x_is, y_is, z_is;
    is.getExactXYZCoordinates(x_is,y_is,z_is);

    bigrational x_it, y_it, z_it;
    it.getExactXYZCoordinates(x_it,y_it,z_it);

    std::cout<<"Sign of x_it:" << y_it.sgn() << std::endl;
    std::cout<<"Print the value of y_it: " << y_it << std::endl;
    //creating the rational points
    const std::vector<bigrational> P0_rational = {x_P0, y_P0, z_P0};
    const std::vector<bigrational> P1_rational = {x_P1, y_P1, z_P1};
    const std::vector<bigrational> P2_rational = {x_P2, y_P2, z_P2};
    std::vector<bigrational> P3_rational = {x_P3, y_P3, z_P3};

    std::vector<bigrational> ir_rational = {x_ir, y_ir, z_ir};
    std::vector<bigrational> is_rational = {x_is, y_is, z_is};
    std::vector<bigrational> it_rational = {x_it, y_it, z_it};
*/
    /************************* ORIENT3D RATIONALS **********************/
   /* auto result_rational_P0 = orient3d(&P0_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P1 = orient3d(&P1_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P2 = orient3d(&P2_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P3 = orient3d(&P3_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);

    //PRINT THE TYPE OF result_rational_P0
    std::cout << "Result P0 Rationals: " << result_rational_P0 << std::endl;
    std::cout << "Result P1 Rationals: " << result_rational_P1 << std::endl;
    std::cout << "Result P2 Rationals: " << result_rational_P2 << std::endl;
    std::cout << "Result P3 Rationals: " << result_rational_P3 << std::endl;
*/
    /**********TEST IMPL TRIANGLE AND SEGMENT *********/

    //segment to generate the intersection
    const std::vector<bigrational> s0 = {bigrational(0.5), bigrational(0.6), bigrational(-1.0)};
    const std::vector<bigrational> s1 = {bigrational(0.5), bigrational(0.6), bigrational(1.0)};

    //test
   /* auto result_segment = segment_triangle_intersect_3d(&s0[0], &s1[0], &P0_rational[0], &P1_rational[0], &P2_rational[0]);

    std::cout << "Result Segment Triangle: " << result_segment << std::endl;
*/
/*
    bigrational intersect_point[3];
    plane_line_intersection(&P0_rational[0], &P1_rational[0], &P2_rational[0], &s0[0], &s1[0], intersect_point);

    std::cout << "Intersect point: " << intersect_point[0] << " " << intersect_point[1] << " " << intersect_point[2] << std::endl;

    double e0 = cinolib::orient3d(&P1_rational[0], &P0_rational[0], &s0[0], &intersect_point[0]).get_d();
    double e1 = cinolib::orient3d(&P0_rational[0], &P2_rational[0], &s0[0], &intersect_point[0]).get_d();
    double e2 = cinolib::orient3d(&P2_rational[0], &P1_rational[0], &s0[0], &intersect_point[0]).get_d();

    std::cout << "e0: " << e0 << std::endl;
    std::cout << "e1: " << e1 << std::endl;
    std::cout << "e2: " << e2 << std::endl;

    // Smallest positive normal numbers
    double smallest_normal_double = std::numeric_limits<double>::min();


    double smallest_subnormal_double = std::numeric_limits<double>::denorm_min();

    std::cout << "Smallest positive normal double: " << smallest_normal_double << std::endl;
    std::cout << "Smallest positive subnormal double: " << smallest_subnormal_double << std::endl;


    bigrational smallest_normal_rational = bigrational(std::numeric_limits<double>::min());
    bigrational smallest_subnormal_rational = bigrational(std::numeric_limits<double>::denorm_min());

    std::cout << "Smallest positive normal rational: " << smallest_normal_rational << std::endl;
    std::cout << "Smallest positive subnormal rational: " << smallest_subnormal_rational << std::endl;


    return 0;
*/
}

