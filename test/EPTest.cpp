//
// Created by jackcook on 8/16/21.
//

#include <cpgfunction/coordinates.h>
#include <cpgfunction/boreholes.h>
#include <cpgfunction/utilities.h>
#include <cpgfunction/gfunction.h>

// This test is done in EnergyPlus for their unit tests


int main() {
    // -- Definitions --
    // Coordinate geometry
    int Nx = 2;
    int Ny = 2;
    double Bx = 5.;
    double By = 5;

    // -- Configurations --
    // Get x,y coordinates for a rectangle
    std::string shape = "Rectangle";
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::configuration(shape, Nx, Ny, Bx, By);

    // -- Borehole geometry --
    double H = 100;  // height of the borehole (in meters)
    double D = 1;  // burial depth (in meters)
    double r_b = 0.057;  // borehole radius (in meters)

    // -- boreField --
    std::vector<gt::boreholes::Borehole> boreField =
            gt::boreholes::boreField(coordinates, r_b, H, D);

    // Thermal properties
    double alpha = 1.0e-6;      // Ground thermal diffusivity (m2/s)

    // Number of segments per borehole
    int nSegments = 1;

    double duration = 1;
    std::string units = "year";
    double const expansion = 0.5;

    vector<double> time = gt::utilities::time_vector_constant_expansion(
            H, alpha, duration, units, expansion);

    std::vector<double> gFunction =
            gt::gfunction::uniform_borehole_wall_temperature(boreField,
                                                             time,
                                                             alpha,
                                                             nSegments,
                                                             true);
    std::cout.precision(16);
    for (size_t i=0; i<gFunction.size(); i++) {
        std::cout << gFunction[i] << std::endl;
    }

}
