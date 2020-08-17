//
// Created by jackcook on 7/15/20.
//

#include <iostream>
#include <vector>
#include "../cppgfunction/SegmentResponse.h"
#include "../cppgfunction/heat_transfer.h"

using namespace std;

#ifndef CPPGFUNCTION_INTERPOLATION_H
#define CPPGFUNCTION_INTERPOLATION_H

namespace jcc::interpolation {

    double linterp(double xp, double x0, double y0, double x1, double y1);
    void interp1d(vector<double>& xp, vector<double>& yp, vector<double>& x, vector<double>& y);
    void interp1d(double &xp, double &yp, vector<double>& x, vector<double>& y);
    void interp1d(double &xp, double &yp, vector<double> &time,
                  unordered_map<gt::heat_transfer::nKey, double, gt::heat_transfer::KeyHasher> &h_map,
                  vector<gt::boreholes::Borehole> &boreSegments,
                  const int i, const int j, const int hash_mode);
    void interp1d(double &xp, double &yp, vector<double> &time,
                  gt::heat_transfer::SegmentResponse &SegRes, int &i, int &j, int &k);

} // jcc::interpolation

#endif //CPPGFUNCTION_INTERPOLATION_H
