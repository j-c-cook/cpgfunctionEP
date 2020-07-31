#include <iostream>
#include "cppgfunction/boreholes.h"
#include <vector>
#include <cmath>
#include "cppgfunction/utilities.h"
#include "cppgfunction/gfunction.h"
#include "jcc/interpolation.h"

int main() {
    // ---------------------------------------------------------
    // Simulation parameters
    // ---------------------------------------------------------

    std::vector<double>  x {0, 200, 400, 600};
    std::vector<double>  y {373.0, 156.1, 113.6, 93.1};
    std::vector<double>  xp {90, 210, 310};
    std::vector<double>  yp (3);
    jcc::interpolation::interp1d(xp, yp, x, y);


    // Borehole dimensions
    // Borehole dimensions
    float D = 4.0;             // Borehole buried depth (m)
    float H = 100.;           // Borehole length (m)
    float r_b = 0.075;         // Borehole radius (m)
    float B = 7.5;             // Borehole spacing (m)

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Number of segments per borehole
    // Geometrically expanding time vector.
    float dt = 100*3600.;                  // Time step
    double tmax = 3000. * 8760. * 3600.;    // Maximum time
    int Nt = 30;                         // Number of time steps
    double ts = pow(H, 2)/(9.*alpha);            // Bore field characteristic time

    float dt_ = 3600.;
    double tmax_ = 0.;
    int Nt_ = 50;
//    std::vector<double> time(Nt_);
    std::vector<double> time(Nt);
    gt::utilities::time_geometric(time, dt, tmax, Nt);
//    gt::utilities::time_geometric(time, dt_, tmax_, Nt_);

    // ---------------------------------------------------------
    // Borehole fields
    // ---------------------------------------------------------

    // Field of 2x3 (n=6) boreholes
    int n1 = 2;
    int n2 = 2;
    int nbh = n1 * n2;
    std::vector<gt::boreholes::Borehole> field(nbh);
    gt::boreholes::rectangle_field(field, n1, n2, B, B,H, D, r_b);
//    field[1].x +=2;
//    field[2].x -=1;
//    field[2].y -=1;
//    field[3].y -=1;

//    std::vector<double> gfunction(Nt_);
    std::vector<double> gfunction(Nt);

    gt::gfunction::uniform_temperature(gfunction, field, time, alpha, 25, true, true);

    for (int i=0; i<gfunction.size(); i++) {
        std::cout << gfunction[i] << std::endl;
    }
//    gt::boreholes::Borehole b1 (H, D, r_b, x, y);
//    gt::boreholes::Borehole b2 (H, D, r_b, 5, y);
//
//    float dist = b2.distance(b1);
//
//    int n1 = 2;
//    int n2 = 2;
//    int nbh = n1 * n2;
//    std::vector<gt::boreholes::Borehole> field(nbh);
//    gt::boreholes::rectangle_field(field, n1, n2, 5, 5, 100, 5, r_b);

    return 0;
}
