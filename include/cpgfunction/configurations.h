//
// Created by jackcook on 2/28/21.
//


#ifndef CPGFUNCTION_CONFIGURATIONS_H
#define CPGFUNCTION_CONFIGURATIONS_H

#include <iostream>
#include <vector>
#include <tuple>

namespace gt {  // g-function namespace
    namespace configurations {  // Configurations struct
        // Uniform - Takes in 8 arguments to create any L, U, Open or rectangular configuration
        struct Uniform {  // Uniform generator struct
            // Destructor
            virtual ~Uniform() = default;
            // instances
            int BottomY;    // The number of rows from the bottom of the field in the y direction
            int LeftX;      // Tthe number of boreholes in the x-direction heading left starting from the far right of the field
            int TopY;       // The number of rows from the top heading in the y-direction down
            int RightX;     // the number of rows from the right of the field heading in the x direction to the right
            double SpaceX;     // The spacing of the borehole field in the x-direction
            double SpaceY;     // The spacing of the borehole field in the y-direction
            double DistanceX;  // The total distance in the x-direction in the field
            double DistanceY;  // The total distance in the y-direction in the field

            explicit Uniform(int BottomY = 0, int LeftX = 0, int TopY = 0, int RightX = 0, double SpaceX = 5.,
                    double SpaceY = 5., double DistanceX = 30., double DistanceY = 30.) : BottomY(BottomY),
                    LeftX(LeftX), TopY(TopY), RightX(RightX), SpaceX(SpaceX), SpaceY(SpaceY), DistanceX(DistanceX),
                    DistanceY(DistanceY){
            }

            // Methods
            std::vector<std::tuple<double, double>> generate_coordinates();  // the "driver" method
            std::vector<std::tuple<double, double>> bottom_left();
            std::vector<std::tuple<double, double>> bottom_right();
            std::vector<std::tuple<double, double>> top_right();
            std::vector<std::tuple<double, double>> top_left();
            void export_coordinates(std::vector<std::tuple<double, double>>& coordinates, std::string output_path);

        };
    }
}

#endif //CPGFUNCTION_CONFIGURATIONS_H
