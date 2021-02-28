//
// Created by jackcook on 2/28/21.
//

#include <cpgfunction/configurations.h>
#include <cmath>
#include <cpgfunction/boreholes.h>


void Append_Unique(std::vector<std::tuple<double, double>>& Existing,
                   const std::vector<std::tuple<double, double>>& NewCoordinates) {
    // Append the unique coordinates to the existing vector of coordinates

    int n = Existing.size();
    int m = NewCoordinates.size();

    std::vector<std::tuple<double, double>> AppendThese;

    for (int i=0; i<n; i++) {
        for (int j=i; j<m; j++) {
            std::tuple<double, double> existing_pair = Existing[i];
            std::tuple<double, double> new_pair = NewCoordinates[j];

            double x1 = std::get<0>(existing_pair);
            double y1 = std::get<1>(existing_pair);

            double x2 = std::get<0>(new_pair);
            double y2 = std::get<1>(new_pair);

            double d = gt::Distance_Formula(x1, y1, x2, y2);

            if (d > 1.0E-06) {
                AppendThese.push_back(new_pair);
            } // FI
        }  // Next j
    }  // Next i

    Existing.insert(Existing.end(), AppendThese.begin(), AppendThese.end());
}


namespace gt {
    namespace configurations {

        std::vector<std::tuple<double, double>> Uniform::generate_coordinates() {
            // make use of the methods to create a vector of tuples for x,y coordinates
            std::vector<std::tuple<double, double>> coordinates;

            // create the coordinates from the bottom left
            std::vector<std::tuple<double, double>> bottom_left_coordinates = bottom_left();
            // insert the bottom left coordinates into the master vector
            Append_Unique(coordinates, bottom_left_coordinates);

            // create the coordinates from the bottom right
            std::vector<std::tuple<double, double>> bottom_right_coordinates = bottom_right();
            // insert the bottom left coordinates into the master vector
            Append_Unique(coordinates, bottom_right_coordinates);

            // create the coordinates from the top right
            std::vector<std::tuple<double, double>> top_right_coordinates = top_right();
            // insert the bottom left coordinates into the master vector
            Append_Unique(coordinates, top_right_coordinates);

            // create the coordinates from the top left
            std::vector<std::tuple<double, double>> top_left_coordinates = top_left();
            // insert the bottom left coordinates into the master vector
            Append_Unique(coordinates, top_left_coordinates);

            return coordinates;
        }

        std::vector<std::tuple<double, double>> Uniform::bottom_left() {
            // This is the bottom left corner control.
            // This returns the (x,y) coordinates with origin as bottom left corner

            std::vector<std::tuple<double, double>> bottom_left_coordinates;

            if (SpaceX < 1.0E-06)  {  // if the x-spacing is 0, then we have no coordinates
                return bottom_left_coordinates;
            } else {
                double div = DistanceX / SpaceX;
                int Nx = std::floor(div) + 1;
                std::tuple<double, double> pair;  // container for coordinate pairs
                for (int i=0; i<Nx; i++) {  // create the x,y coordinates and push them back to the coordinates
                    for (int j=0; j<BottomY; j++) {
                        double _x = double(i) * SpaceX;
                        double _y = double(j) * SpaceY;
                        pair = std::make_tuple(_x, _y);
                        bottom_left_coordinates.push_back(pair);
                    } // next j
                }  // next i
            } // FI
            return bottom_left_coordinates;
        }  // bottom_left function

        std::vector<std::tuple<double, double>> Uniform::bottom_right() {
            // Generate coordinates stemming from the bottom right corner
            std::vector<std::tuple<double, double>> bottom_right_coordinates;

            if (SpaceY < 1.0E-06) {  // if the y-spacing is 0, then we have no coordinates here
                return bottom_right_coordinates;
            } else {
                double div = DistanceY / SpaceY;
                int Ny = std::floor(div) + 1;
                std::vector<double> xValues;
                for (int i=0; i<RightX; i++) {
                    double xv = double(i) * SpaceX;
                    xValues.push_back(xv);
                }
                std::tuple<double, double> pair;
                for (int j=0; j<Ny; j++) {
                    for (int i=0; i<xValues.size(); i++) {
                        double _x = DistanceY - xValues[i];
                        double _y = double(j) * SpaceY;
                        pair = std::make_tuple(_x, _y);
                        bottom_right_coordinates.push_back(pair);
                    } // next i
                } // next j
                return bottom_right_coordinates;
            } // FI
        }  // bottom_right function

        std::vector<std::tuple<double, double>> Uniform::top_right() {
            // The control of the top right corner
            std::vector<std::tuple<double, double>> top_right_coordinates;

            if (SpaceX < 1.0E-06) {
                return top_right_coordinates;
            } else {
                double div = DistanceX / SpaceX;
                int Nx = std::floor(div) + 1;
                std::vector<double> yValues;
                for (int j=0; j<TopY; j++) {
                    double yv = double(j) * SpaceY;
                    yValues.push_back(yv);
                }
                std::tuple<double, double> pair;
                for (int j=0; j<yValues.size(); j++) {
                    for (int i=0; i<Nx; i++) {
                        double _x = double(i) * SpaceX;
                        double _y = DistanceY - yValues[j];
                        pair = std::make_tuple(_x, _y);
                        top_right_coordinates.push_back(pair);
                    } // next i
                } // next j
                return top_right_coordinates;
            }  // FI
        }  // top_right function

        std::vector<std::tuple<double, double>> Uniform::top_left() {
            // The control of the top left corner
            std::vector<std::tuple<double, double>> top_left_coordinates;

            if (SpaceY < 1.0E-06) {
                return top_left_coordinates;
            } else {
                double div = DistanceY / SpaceY;
                std::vector<double> yValues;
                int Ny = std::floor(div) + 1;
                for (int j=0; j<Ny; j++) {
                    double yv = double(j) * SpaceX;
                    yValues.push_back(yv);
                }
                std::tuple<double, double> pair;
                for (int j=0; j<yValues.size(); j++) {
                    for (int i=0; i<LeftX; i++) {
                        double _x = double(i) * SpaceX;
                        double _y = yValues[j];
                        pair = std::make_tuple(_x, _y);
                        top_left_coordinates.push_back(pair);
                    }  // Next i
                }  // Next j
            }  // FI
            return top_left_coordinates;
        }  // top_left function
    }  // Configurations Namespace
}  // g-function name space








































