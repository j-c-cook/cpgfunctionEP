//
// Created by jackcook on 7/11/20.
//

#ifndef CPPGFUNCTION_UTILITIES_H
#define CPPGFUNCTION_UTILITIES_H

#include <iostream>
#include <vector>
#include <cmath>

namespace gt {
    namespace utilities {
        std::vector<double> time_geometric(float dt, double tmax, int Nt);
        void cook_spitler_time (std::vector<double> &logtime);
        void convert_time(std::vector<double> &logtime, std::vector<double> &time, double H, double alpha);
    } // namespace utilities
} // namespace gt

#endif //CPPGFUNCTION_UTILITIES_H
