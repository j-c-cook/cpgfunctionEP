//
// Created by jackcook on 8/4/20.
//

#include <cpgfunction/SegmentResponse.h>
using std::size_t;
using std::hash;

namespace gt { namespace heat_transfer {
    void h_ij_Key::func(double & ti, const int _i, const int _j, const int _k,
                        vector<gt::boreholes::Borehole> &boreSegments, const int _hash_mode, bool Real) {
        hash_mode = _hash_mode;  // set hash mode for when Key hash function gets called
        double precision = 10000;  // number of floating points
        double _time = ti * precision;
        time = round(_time);

        if (hash_mode==0) {
            // if the hash mode is 0 then we are only storing similarities
            gt::boreholes::Borehole b1 = boreSegments[_i];
            gt::boreholes::Borehole b2 = boreSegments[_j];
            // EXTENSION of similarities
            // Cimmino defined a distance tolerance of 0.1 if its not the same borehole
//            double disTol = 0.1;
//            double rTol;
//            if (_i==_j) {
//                rTol = 1.0e-6 * b1.r_b;
//            } else {
//                rTol = disTol;
//            }
//            double dist = b1.distance(b2);
//            double d = b1.distance(b2) * 1 / rTol;
//            int tmp_d = round(d);
//            double DELTA = abs(remainder(d ,double(tmp_d)));
//            d -= DELTA;
//            d = d * precision;
            double d = b1.distance(b2) * precision;

            double Hij = b1.H/2 - b2.H * precision;
            double Dij;
            if (Real) {
                // the Real segment responses will be stored with the difference of burial depths
                Dij = (b1.D - b2.D) * precision;
            } else {
                // the Image segment responses will be stored with the summation of burial depths
                Dij = (b1.D + b2.D) * precision;
            }
            dij = round(d);
            delta_Hij = round(Hij);
            delta_D_ij = round(Dij);
        } else if (hash_mode==1) {
            // if the hash mode is 1 then we are storing all i->j interactions (one way)
            i = _i;
            j = _j;
            k = _k;
        } else {
            throw invalid_argument("Currently the hash_mode you input is not accepted.");
        }
    }

    bool check_key (unordered_map<nKey, double, KeyHasher> &h_map, nKey &sim_key) {
        // Key is not present
        if (h_map.find(sim_key) == h_map.end())
            return false;
        return true;
    } // check_key

    double hash_table_lookup(double &h, unordered_map<nKey, double, KeyHasher> &h_map,
                           vector<double> &time, vector<gt::boreholes::Borehole> &boreSegments,
                           const int i, const int j, const int k, const int hash_mode) {
        gt::boreholes::Borehole b1;
        gt::boreholes::Borehole b2;
        if (hash_mode==0) {
            // real and image stored seperate
            ;
//            h_ij_Key _h_map_key_calc_real(time[k], i, j, k, boreSegments, hash_mode, true);
//            Key hij_map_key_real(_h_map_key_calc_real);
//            h_ij_Key _h_map_key_calc_image(time[k], i, j, k, boreSegments, hash_mode, false);
//            Key hij_map_key_image(_h_map_key_calc_image);
//            double R;
//            double I;
//            b1 = boreSegments[i];
//            b2 = boreSegments[j];
//            if (check_key(h_map, hij_map_key_real)) {
//                // (i, j) exists for real
//                R = h_map[hij_map_key_real];
//            } else {
//                // (j, i) exists for real
//                _h_map_key_calc_real.func(time[k], j, i, k, boreSegments, hash_mode, true);
//                Key hij_map_key_real_ji(_h_map_key_calc_real);
//                R = b2.H / b1.H * h_map[hij_map_key_real_ji];
//            }
//            if (check_key(h_map, hij_map_key_image)) {
//                // (i, j) exists for image
//                I = h_map[hij_map_key_image];
//            } else {
//                // (j, i) exists for image
//                _h_map_key_calc_image.func(time[k], j, i, k, boreSegments, hash_mode, false);
//                Key hij_map_key_image_ji(_h_map_key_calc_image);
//                I = b2.H / b1.H * h_map[hij_map_key_image_ji];
//            }
////            HIJ[i][j][k+1] = R + I;
////            h_ij[i][j][k+1] = R + I;
//            h = R + I;
        } else if (hash_mode==1) {
            double RandI_ij;
            if (i <= j) {
                h_ij_Key _h_map_key_calc_realandimage(time[k], i, j, k, boreSegments, hash_mode);
                RandI_ij = h_map[make_tuple(i, j, k)];
//                Key hij_map_key_realandimage(_h_map_key_calc_realandimage);
//                RandI_ij = h_map[hij_map_key_realandimage];

            } else {
                h_ij_Key _h_map_key_calc_realandimage(time[k], j, i, k, boreSegments, hash_mode);
                Key hij_map_key_realandimage_ji(_h_map_key_calc_realandimage);
                b1 = boreSegments[i];
                b2 = boreSegments[j];
                RandI_ij = b2.H / b1.H * h_map[make_tuple(i, j, k)];
//                RandI_ij = b2.H / b1.H * h_map[hij_map_key_realandimage_ji];
            }


//            if (check_key(h_map, hij_map_key_realandimage)) {
//                // (i, j) exists for real and image
//                RandI_ij = h_map[hij_map_key_realandimage];
//            } else {
//                // (j, i) exists for real and image
//                _h_map_key_calc_realandimage.func(time[k], j, i, k, boreSegments, hash_mode, true);
//                Key hij_map_key_realandimage_ji(_h_map_key_calc_realandimage);
//                b1 = boreSegments[i];
//                b2 = boreSegments[j];
//                RandI_ij = b2.H / b1.H * h_map[hij_map_key_realandimage_ji];
//            }
            h =  RandI_ij;
//            HIJ[i][j][k+1] = RandI_ij;
//            h_ij[i][j][k+1] = RandI_ij;
        }
    }

} } // gt::heat_transfer
