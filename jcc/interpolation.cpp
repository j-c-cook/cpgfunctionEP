//
// Created by jackcook on 7/15/20.
//

#include "interpolation.h"
#include <csignal>

using namespace std;

namespace jcc::interpolation {

    double linterp(double xp, double x0, double y0, double x1, double y1) {
        double yp;
        yp = y0 + ((y1-y0) / (x1-x0)) * (xp-x0);
        return yp;
    } // interp

    double interp1d(vector<double>& xp, vector<double>& yp, vector<double>& x, vector<double>& y) {
        int counter = 0;
        for (int i=0; i<yp.size(); i++) {
            for (int j = counter; j<y.size();j++) {
                if (xp[i] - x[j] < 10) {
                    yp[i] = y[j];
                    break;
                } else if (xp[i] >= x[j] && xp[i] <= x[j+1]) {
                    yp[i] = linterp(xp[i], x[j], y[j], x[j+1], y[j+1]);
                    break;
                } else if (xp[i] < x[0] || xp[i] > x[x.size()-1]) {
                    throw invalid_argument( "Need to add extrapolation" );
                } else {
                    counter++;
                } // fi
            } // next j
        } // next i
    } // interp1d

} // jcc::interpolation


