//
// Created by Michele on 03/06/24.
//


#include <iostream>

#include <cinolib/rationals.h>
#include <implicit_point.h>
using namespace cinolib;


int main(int argc, char **argv) {
    //plane
    const explicitPoint3D *ir = new explicitPoint3D(0.0,0.0,0.0);
    const explicitPoint3D *is = new explicitPoint3D(4.0,4.0,0.0);
    const explicitPoint3D *it = new explicitPoint3D(0.0,4.0,0.0);

    //line to generate P0
    const explicitPoint3D *iq_P0 = new explicitPoint3D(0.5,0.5,-1.0);
    const explicitPoint3D *ip_P0 = new explicitPoint3D(0.5,0.5,1.0);

    //line to generate P1
    const explicitPoint3D *iq_P1 = new explicitPoint3D(0.5,1,-1.0);
    const explicitPoint3D *ip_P1 = new explicitPoint3D(0.5,1,1.0);

    //line to generate P2
    const explicitPoint3D *iq_P2 = new explicitPoint3D(1.0,1.0,-1.0);
    const explicitPoint3D *ip_P2 = new explicitPoint3D(1.0,1.0,1.0);

    //line to generate P3
    const explicitPoint3D *iq_P3 = new explicitPoint3D(1.0,0.5,-1.0);
    const explicitPoint3D *ip_P3 = new explicitPoint3D(1.0,0.5,1.0);

    //create 4 implicit points
    genericPoint *P0 = new implicitPoint3D_LPI(*iq_P0, *ip_P0, *ir, *is, *it);
    genericPoint *P1 = new implicitPoint3D_LPI(*iq_P1, *ip_P1, *ir, *is, *it);
    genericPoint *P2 = new implicitPoint3D_LPI(*iq_P2, *ip_P2, *ir, *is, *it);
    genericPoint *P3 = new implicitPoint3D_LPI(*iq_P3, *ip_P3, *ir, *is, *it);

    /************************* ORIENT3D IMPLICIT **********************/
    int T0_orient = genericPoint::orient3D(*P0, *ir, *it, *is);
    int T1_orient = genericPoint::orient3D(*P1, *ir, *it, *is);
    int T2_orient = genericPoint::orient3D(*P2, *ir, *it, *is);
    int T3_orient = genericPoint::orient3D(*P3, *ir, *it, *is);

    std::cout << "Result P0 Implicit: " << T0_orient << std::endl;
    std::cout << "Result P1 Implicit: " << T1_orient << std::endl;
    std::cout << "Result P2 Implicit: " << T2_orient << std::endl;
    std::cout << "Result P3 Implicit: " << T3_orient << std::endl;

    /*Rationals*/
    bigrational x_P0, y_P0, z_P0;
    P0->getExactXYZCoordinates(x_P0,y_P0,z_P0);

    bigrational x_P1, y_P1, z_P1;
    P1->getExactXYZCoordinates(x_P1,y_P1,z_P1);

    bigrational x_P2, y_P2, z_P2;
    P2->getExactXYZCoordinates(x_P2,y_P2,z_P2);

    bigrational x_P3, y_P3, z_P3;
    P3->getExactXYZCoordinates(x_P3,y_P3,z_P3);

    bigrational x_ir, y_ir, z_ir;
    ir->getExactXYZCoordinates(x_ir,y_ir,z_ir);

    bigrational x_is, y_is, z_is;
    is->getExactXYZCoordinates(x_is,y_is,z_is);

    bigrational x_it, y_it, z_it;
    it->getExactXYZCoordinates(x_it,y_it,z_it);

    //creating the rational points
    std::vector<bigrational> P0_rational = {x_P0, y_P0, z_P0};
    std::vector<bigrational> P1_rational = {x_P1, y_P1, z_P1};
    std::vector<bigrational> P2_rational = {x_P2, y_P2, z_P2};
    std::vector<bigrational> P3_rational = {x_P3, y_P3, z_P3};

    std::vector<bigrational> ir_rational = {x_ir, y_ir, z_ir};
    std::vector<bigrational> is_rational = {x_is, y_is, z_is};
    std::vector<bigrational> it_rational = {x_it, y_it, z_it};

    /************************* ORIENT3D RATIONALS **********************/
    auto result_rational_P0 = orient3d(&P0_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P1 = orient3d(&P1_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P2 = orient3d(&P2_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);
    auto result_rational_P3 = orient3d(&P3_rational[0], &ir_rational[0], &it_rational[0], &is_rational[0]);

    std::cout << "Result P0 Rationals: " << result_rational_P0 << std::endl;
    std::cout << "Result P1 Rationals: " << result_rational_P1 << std::endl;
    std::cout << "Result P2 Rationals: " << result_rational_P2 << std::endl;
    std::cout << "Result P3 Rationals: " << result_rational_P3 << std::endl;

    return 0;

}

