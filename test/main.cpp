#include <iostream>
#include <cpgfunction/boreholes.h>
#include <vector>
#include <cmath>
#include <cpgfunction/utilities.h>
#include <cpgfunction/gfunction.h>

int main() {
    // ---------------------------------------------------------
    // Simulation parameters
    // ---------------------------------------------------------

    // Borehole dimensions
    float D = 4.0;             // Borehole buried depth (m)
    float H = 100.;           // Borehole length (m)
    float r_b = 0.075;         // Borehole radius (m)
    float B = 7.5;             // Borehole spacing (m)

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Number of segments per borehole
    int nSegments = 25;

    // Geometrically expanding time vector.
    float dt = 100*3600.;                  // Time step
    double tmax = 3000. * 8760. * 3600.;    // Maximum time
    int Nt = 30;                         // Number of time steps
    double ts = pow(H, 2)/(9.*alpha);            // Bore field characteristic time

    std::vector<double> time = gt::utilities::time_geometric(dt, tmax, Nt);

    // ---------------------------------------------------------
    // Borehole fields
    // ---------------------------------------------------------

    // Field of 2x3 (n=6) boreholes
    int n1 = 2;
    int n2 = 2;
    std::vector<gt::boreholes::Borehole> borefield = gt::boreholes::rectangle_field(n1, n2, B, B,H, D, r_b);
    borefield[1].x +=2;
    borefield[2].x -=1;
    borefield[2].y -=1;
    borefield[3].y -=1;

    std::vector<double> gfunction(Nt);

    gt::gfunction::uniform_temperature(gfunction, borefield, time, alpha, nSegments,true, true);

    for (double i : gfunction) {
        std::cout << i << std::endl;
    }

    return 0;
}
