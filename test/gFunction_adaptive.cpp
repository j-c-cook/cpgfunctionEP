//
// Created by jackcook on 7/15/21.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <cpgfunction/coordinates.h>
#include <cpgfunction/boreholes.h>
#include <cpgfunction/utilities.h>
#include <cpgfunction/segments.h>
#include <cpgfunction/gfunction.h>

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
    // adaptive discretization
    double drilling_depth = nbh * H;
    gt::segments::adaptive adpt_disc;
    int nSegments = adpt_disc.discretize(H, drilling_depth);

    // Total time duration (month or year)
    std::string time_units = "year";
    double duration = 300.;
    // Geometrically growing time vector for duration
    std::vector<double> time = gt::utilities::time_geometric_auto(duration, time_units);

    // Coordinates
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::rectangle(n1, n2, B, B);
    // Create borehole field
    std::vector<gt::boreholes::Borehole> boreField =
            gt::boreholes::boreField(coordinates, r_b, H, D);
    // Detect number of threads (default uses all available threads)
    int n_Threads = int(thread::hardware_concurrency());
    // Compute (and time) the g-Function
    auto start = std::chrono::steady_clock::now();
    std::vector<double> gFunction =
            gt::gfunction::uniform_borehole_wall_temperature(boreField,
                                                             time,
                                                             alpha,
                                                             nSegments,
                                                             true,
                                                             n_Threads,
                                                             false);
    auto end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration_cast<
            std::chrono::milliseconds>(end - start).count() / 1000.;
    std::cout << "g-Function calculation duration: " << seconds << " seconds"
              << std::endl;

    for (double i : gFunction) {
        std::cout << i << std::endl;
    }

    return 0;
}  // main();