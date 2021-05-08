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
        std::vector<double> time_geometric(double dt, double tmax, int Nt);
        std::vector<double> time_Eskilson(const double &H, const double &alpha);
        std::vector<double> convert_time(std::vector<double> &logtime, const double &H, const double &alpha);
        void cook_spitler_time (std::vector<double> &logtime);
        void convert_time(std::vector<double> &logtime, std::vector<double> &time, double H, double alpha);
    } // namespace utilities
} // namespace gt

#endif //CPPGFUNCTION_UTILITIES_H
