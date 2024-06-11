//
// Created by Michele on 03/06/24.
//


#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/rationals.h>
#include <implicit_point.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Point_3.h>

#include <string>
#include <sstream>
using namespace cinolib;


//create a function
inline void createRationalPoint(){

}

int main(int argc, char **argv) {

    const explicitPoint3D *ir = new explicitPoint3D(0.0,0.0,0.0);
    const explicitPoint3D *is = new explicitPoint3D(4.0,4.0,0.0);
    const explicitPoint3D *it = new explicitPoint3D(0.0,4.0,0.0);

    const explicitPoint3D *PY = new explicitPoint3D(2.0, 2.0, 1);


    const explicitPoint3D *iq_P0 = new explicitPoint3D(0.5,0.5,-1.0);
    const explicitPoint3D *ip_P0 = new explicitPoint3D(0.5,0.5,1.0);

    const explicitPoint3D *iq_P1 = new explicitPoint3D(0.5,1,-1.0);
    const explicitPoint3D *ip_P1 = new explicitPoint3D(0.5,1,1.0);

    const explicitPoint3D *iq_P2 = new explicitPoint3D(1.0,1.0,-1.0);
    const explicitPoint3D *ip_P2 = new explicitPoint3D(1.0,1.0,1.0);

    const explicitPoint3D *iq_P3 = new explicitPoint3D(1.0,0.5,-1.0);
    const explicitPoint3D *ip_P3 = new explicitPoint3D(1.0,0.5,1.0);

    genericPoint *P0 = new implicitPoint3D_LPI(*iq_P0, *ip_P0, *ir, *is, *it);
    genericPoint *P1 = new implicitPoint3D_LPI(*iq_P1, *ip_P1, *ir, *is, *it);
    genericPoint *P2 = new implicitPoint3D_LPI(*iq_P2, *ip_P2, *ir, *is, *it);
    genericPoint *P3 = new implicitPoint3D_LPI(*iq_P3, *ip_P3, *ir, *is, *it);


    int T0_orient= genericPoint::orient3D(*P0, *ir, *it, *is);
    int T1_orient = genericPoint::orient3D(*P1, *ir, *it, *is );
    int T2_orient = genericPoint::orient3D(*P2, *ir, *it, *is );
    int T3_orient = genericPoint::orient3D(*P3, *ir, *it, *is );


    std::cout << "T0_orient: " << T0_orient << std::endl;
    std::cout << "T1_orient: " << T1_orient << std::endl;
    std::cout << "T2_orient: " << T2_orient << std::endl;
    std::cout << "T3_orient: " << T3_orient << std::endl;

    /*RAZIONALI*/
    bigrational x_P0, y_P0, z_P0;
    P0->getExactXYZCoordinates(x_P0,y_P0,z_P0);

    bigrational x_P1, y_P1, z_P1;
    P1->getExactXYZCoordinates(x_P1,y_P1,z_P1);

    bigrational x_P2, y_P2, z_P2;
    P2->getExactXYZCoordinates(x_P2,y_P2,z_P2);

    bigrational x_P3, y_P3, z_P3;
    P3->getExactXYZCoordinates(x_P3,y_P3,z_P3);

/***********************************************************************************************************************/
    // The sequence of numbers
    std::vector<uint> num;
    std::vector<uint> den;

    for(int i = 0; i < x_P0.get_num().size(); i++){
        num.push_back(x_P0.get_num()[i]);
    }
    for(int i = 0; i < x_P0.get_den().size(); i++){
        den.push_back(x_P0.get_den()[i]);
    }
    // Create an empty string to store the concatenated numbers
    std::string concatenated;

    // Loop through the numbers and concatenate them into a string
    for (int number : num) {
        concatenated += std::to_string(number);
    }

    // Convert the concatenated string to a long integer
    long long result = std::stoll(concatenated);

    // Output the result
    std::cout << "The concatenated number is: " << result << std::endl;



    std::vector<CGAL_Q> rationals;
    // Define the kernel
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

    // Define the types for rational coordinates
    typedef CGAL::Quotient<CGAL::MP_Float> Rational;

    // Define the 3D point type using rational coordinates
    typedef CGAL::Point_3<K, Rational> Point_3;


    // Create rational numbers using CGAL::Gmpq
    //CGAL::Gmpq x() // represents 3/4




    return 0;

}

