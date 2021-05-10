//
// Created by jackcook on 2/28/21.
//

// tests for creating coordinates

#include <cpgfunction/coordinates.h>

int main() {
    // -- Definitions --
    // Coordinate geometry
    int Nx = 3;
    int Ny = 3;
    double Bx = 6.;
    double By = 4.5;

    // -- Configurations --
    // Get x,y coordinates for a rectangle
    std::string shape = "Rectangle_configuration";
    std::vector<std::tuple<double, double>> coordinates =
            gt::coordinates::configuration(shape, Nx, Ny, Bx, By);
    std::string output_path = shape + ".json";
    gt::coordinates::export_coordinates_to_file(coordinates, output_path);

    // Get x,y coordinates for an open-rectangle
    shape = "OpenRectangle_configuration";
    coordinates = gt::coordinates::configuration(shape, Nx, Ny, Bx, By);
    output_path = shape + ".json";
    gt::coordinates::export_coordinates_to_file(coordinates, output_path);

    // Get x,y coordinates for a U shape
    shape = "U_configuration";
    coordinates = gt::coordinates::configuration(shape, Nx, Ny, Bx, By);
    output_path = shape + ".json";
    gt::coordinates::export_coordinates_to_file(coordinates, output_path);

    // Get x,y coordinates for a U shape
    shape = "L_configuration";
    coordinates = gt::coordinates::configuration(shape, Nx, Ny, Bx, By);
    output_path = shape + ".json";
    gt::coordinates::export_coordinates_to_file(coordinates, output_path);

    // Read in custom configuration from (.json) file
    shape = "custom_configuration";
    std::string input_path = "U.json";
    coordinates = gt::coordinates::configuration(shape, input_path);
    output_path = shape + ".json";
    gt::coordinates::export_coordinates_to_file(coordinates, output_path);

    return 0;
}