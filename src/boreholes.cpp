//
// Created by jackcook on 7/11/20.
//

#include <cpgfunction/boreholes.h>


namespace gt {

    double Distance_Formula(double x1, double y1, double x2, double y2) {
        return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
    }

    namespace boreholes {

        double Borehole::distance(Borehole target) {
            double x1 = x;
            double y1 = y;
            double x2 = target.x;
            double y2 = target.y;
            double dist = Distance_Formula(x1, y1, x2, y2);
            return std::max(r_b, double(dist));  // max needs doubles
        }

        std::tuple<double, double> Borehole::position() {
            std::tuple<double, double> t (x, y);
            return t;
        };

        std::vector<Borehole> boreField(const std::vector<std::tuple<double, double>> &coordinates, const double &r_b,
                                        const double &H, const double &D){
            std::vector<Borehole> bores(coordinates.size());

            double x;
            double y;

            for (int i = 0; i < coordinates.size(); i++) {
                x = std::get<0>(coordinates[i]);
                y = std::get<1>(coordinates[i]);
                bores[i] = Borehole(H, D, r_b, x, y);
            }  // next i

            return bores;
        }  // boreField();

    // ------ Create Fields --------
    std::vector<Borehole> rectangle_field(int N_1, int N_2, double B_1, double B_2, double H, double D, double r_b) {
        int nbh = N_1 * N_2;
        std::vector<Borehole> borefield(nbh);

        double pos_x;
        double pos_y;
        int point_nb = 0;
        for( int j = 0; j < N_2; j++ ) {
            for (int i = 0; i < N_1; i++) {
                pos_x = double(i) * B_1;
                pos_y = double(j) * B_2;
                borefield[point_nb] = Borehole(H, D, r_b, pos_x, pos_y);
                point_nb++;
            } // end i
        } // end j
        return borefield;
    } // rectangular field function

    } // boreholes name space

} // gt name space
