//
// Created by jackcook on 7/11/20.
//

#ifndef CPPGFUNCTION_BOREHOLES_H
#define CPPGFUNCTION_BOREHOLES_H

#include <iostream>
#include <math.h>
#include <tuple>
#include <vector>

namespace gt {
    namespace boreholes {

        struct Borehole {
            // Destructor
            virtual ~Borehole() {
            }

            double H;    // height or length of the borehole (meters)
            double D;    // borehole burial depth (meters)
            double r_b;  // borehole radius (meters)
            double x;    // position (meters) of the center of the borehole along the x-axis
            double y;    // position (meters) of the center of the borehole along the y-axis

            Borehole(double H=0.0, double D=0.0, double r_b=0.0, double x=0.0, double y=0.0) : H(H), D(D), r_b(r_b), x(x), y(y) {
            }

            double distance(Borehole target);
            std::tuple<double, double> position();
        };

        std::vector<Borehole> rectangle_field(int N_1, int N_2, double B_1, double B_2, double H, double D, double r_b);

    } // boreholes name space

} // gt namespace

#endif //CPPGFUNCTION_BOREHOLES_H
