//
// Created by jackcook on 8/18/21.
//

//
// Created by jackcook on 8/16/21.
//

#include <cpgfunction/coordinates.h>
#include <cpgfunction/boreholes.h>
#include <cpgfunction/utilities.h>
#include <cpgfunction/gfunction.h>


int main() {
    // -- Borehole geometry --
    double H = 100;  // height of the borehole (in meters)
    double D = 1;  // burial depth (in meters)
    double r_b = 0.057;  // borehole radius (in meters)

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Number of segments per borehole
    int nSegments = 12;

    double duration = 1;
    std::string units = "year";
    double const expansion = 0.5;

    vector<double> time = gt::utilities::time_vector_constant_expansion(
            H, alpha, duration, units, expansion);

    // -- Definitions --
    // Coordinate geometry
    double Bx = 5.;
    double By = 5;


    for (int i=2; i<17; i+=2) {
        int Nx = i;
        int Ny = i;
        int nSources = nSegments * Nx * Ny;

        // -- Configurations --
        // Get x,y coordinates for a rectangle
        std::string shape = "Rectangle";
        std::vector<std::tuple<double, double>> coordinates =
                gt::coordinates::configuration(shape, Nx, Ny, Bx, By);

        // -- boreField --
        std::vector<gt::boreholes::Borehole> boreField =
                gt::boreholes::boreField(coordinates, r_b, H, D);


        auto start = std::chrono::steady_clock::now();
        std::vector<double> gFunction =
                gt::gfunction::uniform_borehole_wall_temperature(boreField,
                                                                 time,
                                                                 alpha,
                                                                 nSegments,
                                                                 true);
        auto end = std::chrono::steady_clock::now();
        auto milli = chrono::duration_cast<chrono::milliseconds>
                (end - start).count();
        double seconds = double(milli) / 1000;
        std:: cout << nSources << "\t" << seconds << std::endl;
    }
}
