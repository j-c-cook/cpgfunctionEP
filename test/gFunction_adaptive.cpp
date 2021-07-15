//
// Created by jackcook on 7/15/21.
//

#include <iostream>
#include <vector>
#include <cpgfunction/coordinates.h>
#include <cpgfunction/boreholes.h>

int main(){
    // ---------------------------------------------------------
    // Inputs
    // ---------------------------------------------------------

    // Borehole dimensions
    float D = 4.0;             // Borehole buried depth (m)
    float H = 100.;            // Borehole length (m)
    float r_b = 0.075;         // Borehole radius (m)
    float B = 7.5;             // Borehole spacing (m)

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Field of 3x4 (n=12) boreholes
    int n1 = 3;
    int n2 = 4;
    int nbh = n1 * n2;

    // Number of segments per borehole
    // TODO: Add in function for this

    // Total time duration (in sec, hour, month or year)
    std::string time_units = "year";
    double duration = 20;
    // TODO: Add in function for this

    // Coordinates
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::rectangle(n1, n2, B, B);
    // Create borehole field
    std::vector<gt::boreholes::Borehole> boreField =
            gt::boreholes::boreField(coordinates, r_b, H, D);

    return 0;
}  // main();