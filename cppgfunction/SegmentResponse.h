//
// Created by jackcook on 8/4/20.
//

#include <iostream>
#include "boreholes.h"
#include "unordered_map"
#include <boost/functional/hash.hpp>

#ifndef CPPGFUNCTION_SEGMENTRESPONSE_H
#define CPPGFUNCTION_SEGMENTRESPONSE_H

using namespace std;

namespace gt::heat_transfer {

    using nKey = tuple<int, int, int>;

    struct Key
    {
        unsigned long int time;
        int dij;
        int delta_D_ij;
        int delta_Hij;
        int i;
        int j;
        int k;
        int hash_mode;

        bool operator==(const Key &other) const
        {
            if (hash_mode==0) {
                return (time == other.time
                && dij == other.dij
                && delta_Hij == other.delta_Hij
                && delta_D_ij == other.delta_D_ij);
            } else if (hash_mode==1) {
                return (time == other.time
                        && i == other.i
                        && j == other.j);
            } else {
                throw invalid_argument("Currently the hash_mode you input is not accepted.");
            }
        }
    };

    struct h_ij_Key : Key {
        ~h_ij_Key(){};

        void func(double & ti, const int i, const int j, int k, vector<gt::boreholes::Borehole> &boreSegments, const int hash_mode, bool Real);

        h_ij_Key(double &ti, const int i, const int j, const int k, vector<gt::boreholes::Borehole> &boreSegments, const int hash_mode, bool Real=true){
            func(ti, i, j, k, boreSegments, hash_mode, Real);
        };
    };

    struct KeyHasher
    {
        std::size_t operator()(const nKey& k) const
        {
            using std::size_t;
            using std::hash;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            return boost::hash_value(k);
//            if (k.hash_mode==0) {
//                return hash<unsigned long int>()(k.time) << 1 ^ hash<int>()(k.dij)
//                        << 1 ^ hash<int>()(k.delta_Hij) << 1 ^ hash<int>()(k.delta_D_ij);
//            } else if (k.hash_mode==1) {
//                tuple<int, int, int> t{k.i, k.j, k.k};
//                return boost::hash_value(t);
////                return hash<unsigned long int>()(k.time) << 1 ^ hash<int>()(k.i) << 1 ^ hash<int>()(k.j);
//            } else {
//                throw invalid_argument("Currently the hash_mode you input is not accepted.");
//            }
        }
    };

    bool check_key(unordered_map<nKey, double, KeyHasher> &h_map, nKey &sim_key);
    double hash_table_lookup(double &h, unordered_map<nKey, double, KeyHasher> &h_map,
                           vector<double> &time, vector<gt::boreholes::Borehole> &boreSegments,
                           int i, int j, int k, int hash_mode);
}

#endif //CPPGFUNCTION_SEGMENTRESPONSE_H
