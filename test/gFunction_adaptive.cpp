//
// Created by jackcook on 7/15/21.
//

#include <iostream>
#include <vector>
#include <cpgfunction/coordinates.h>
#include <cpgfunction/boreholes.h>
#include <cpgfunction/utilities.h>

int main(){
    // ---------------------------------------------------------
    // Inputs
    // ---------------------------------------------------------

    // Borehole dimensions
    double D = 4.0;             // Borehole buried depth (m)
    double H = 100.;            // Borehole length (m)
    double r_b = 0.075;         // Borehole radius (m)
    double B = 7.5;             // Borehole spacing (m)

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
    double duration = 20.;
    std::vector<double> time = gt::utilities::time_vector(
            H, alpha, duration, time_units);

    // Coordinates
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::rectangle(n1, n2, B, B);
    // Create borehole field
    std::vector<gt::boreholes::Borehole> boreField =
            gt::boreholes::boreField(coordinates, r_b, H, D);

    return 0;
}  // main();