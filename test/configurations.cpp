//
// Created by jackcook on 2/28/21.
//

// tests for creating coordinates

#include <cpgfunction/configurations.h>

int main() {

    // create a rectangle

    int BottomY = 10;  // nrows in Y from bottom
    int LeftX = 10;    // nrows in x from left
    int TopY = 0;     // nrows in y from top
    int RightX = 0;   // nrows in x from right

    double SpaceX = 5;   // borehole spacing in the x
    double SpaceY = 5;   // borehole spacing in the y

    double DistanceX = (LeftX - 1) * SpaceX;  // total distance covered in the x
    double DistanceY = (BottomY - 1) * SpaceY;  // total distance covered in the y

    gt::configurations::Uniform uni(BottomY=BottomY, LeftX=LeftX, TopY=TopY, RightX=RightX, SpaceX=SpaceX,
                                    SpaceY=SpaceY, DistanceX=DistanceX, DistanceY=DistanceY);

    std::vector<std::tuple<double, double>> coordinates = uni.generate_coordinates();

    int a = 1;

    return 0;
}