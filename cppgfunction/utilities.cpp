//
// Created by jackcook on 7/11/20.
//

#include "utilities.h"

namespace gt {
    namespace utilities {
        void time_geometric(std::vector<double>& time, float dt, double tmax, int Nt) {
            // TODO: place this resizing into a "general" namespace
            int len = time.size(); // check to see if there is enough space in the vector
            // if need be, resize the vector to be the same size as the number of boreholes needed
            if (len != Nt) {
                time.resize(Nt);
            } else ; // else do nothing

            double value;
            double tmax_calc = double(Nt) * double(dt);
            if (tmax > tmax_calc) {
                double dr = 1.0e99;
                double r = 2.;
                while (std::abs(dr) > 1.0e-10) {
                    dr = std::pow(1+tmax/double(dt)*(r-1), 1/double(Nt)) - r;
                    r += dr;
                } // end while
                for (int j=0; j < time.size(); j++) {
                    value = 1 - std::pow(r, double(j+1));
                    value = value * double(dt) / (1 - r);
                    time[j] = value;
                } // end for
            } else {
                for (int j=0; j < time.size(); j++) {
                    value = double(dt) * double(j + 1);
                    time[j] = value;
                } // end for
            } // end if
        } // void time_geometric

    } // namespace utilities
} // namespace gt
