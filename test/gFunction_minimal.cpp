#include <iostream>
#include <cpgfunction/boreholes.h>
#include <vector>
#include <cmath>
#include <cpgfunction/utilities.h>
#include <cpgfunction/gfunction.h>
#include <cpgfunction/coordinates.h>
#include <chrono>


int main() {
    // ---------------------------------------------------------
    // Simulation parameters
    // ---------------------------------------------------------

    // Borehole dimensions
    float D = 4.0;             // Borehole buried depth (m)
    float H = 100.;            // Borehole length (m)
    float r_b = 0.075;         // Borehole radius (m)
    float B = 7.5;             // Borehole spacing (m)

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Number of segments per borehole
    int nSegments = 12;

    // Geometrically expanding time vector.
    float dt = 100*3600.;                   // Time step
    double tmax = 3000. * 8760. * 3600.;    // Maximum time
    int Nt = 30;                             // Number of time steps
    double ts = pow(H, 2)/(9.*alpha);    // Bore field characteristic time

    std::vector<double> time = gt::utilities::time_geometric(dt, tmax, Nt);

    // ---------------------------------------------------------
    // Borehole fields
    // ---------------------------------------------------------

    // Field of 3x4 (n=12) boreholes
    int n1 = 3;
    int n2 = 4;

    // Coordinates
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::rectangle(n1, n2, B, B);

    // Move boreholes slightly off of symmetric
    std::get<0>(coordinates[1]) += 2;  // move position 1 borehole x
    std::get<0>(coordinates[2]) -= 1;  // move position 2 borehole x
    std::get<1>(coordinates[2]) -= 1;  // move position 2 borehole y
    std::get<1>(coordinates[3]) -= 1;  // move position 3 borehole y

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
}
