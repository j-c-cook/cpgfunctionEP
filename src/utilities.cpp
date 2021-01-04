//
// Created by jackcook on 7/11/20.
//

#include <cpgfunction/utilities.h>

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

        void cook_spitler_time (std::vector<double> &logtime){
            int np = 31; // 31 total points
            if (logtime.size() != np) {
                logtime.resize(np);
            }
            for (int i=1; i<np+1; i++) {
                logtime[i-1] = 0.0005*pow(i, 3) - 0.0303 * pow(i, 2) + 0.8491*i - 9.4028;
            }
        } // cook_spitler_time

        void convert_time(std::vector<double> &logtime, std::vector<double> &time, const double H, const double alpha) {
            int nt = logtime.size();
            if (time.size() != nt) {
                time.resize(nt);
            }
            double ts = pow(H, 2) / (9 * alpha);
            for (int i=0; i<nt; i++) {
                time[i] = exp(logtime[i]) * ts;
            }
        } // convert_time

    } // namespace utilities
} // namespace gt
