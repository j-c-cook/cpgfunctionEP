//
// Created by jackcook on 7/15/20.
//

#include <iostream>
#include <vector>
#include <cpgfunction/segments.h>

using namespace std;

#ifndef CPPGFUNCTION_INTERPOLATION_H
#define CPPGFUNCTION_INTERPOLATION_H

namespace jcc::interpolation {

    double linterp(double xp, double x0, double y0, double x1, double y1);
    void interp1d(vector<double>& xp, vector<double>& yp, vector<double>& x,
                  vector<double>& y);
    void interp1d(double &xp, double &yp, vector<double>& x, vector<double>& y);
    void interp1d(double &xp, double &yp, vector<double> &time,
                  gt::segments::SegmentResponse &SegRes, int &i, int &j, int &k);
    double bilinterp1d(double q11, double q12, double q21, double q22,
                       double x1, double x2, double y1, double y2,
                       double x, double y);

} // jcc::interpolation

#endif //CPPGFUNCTION_INTERPOLATION_H
