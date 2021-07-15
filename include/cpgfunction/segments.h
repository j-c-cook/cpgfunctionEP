//
// Created by jackcook on 7/15/21.
//

#include <iostream>
#include <vector>
#include <cpgfunction/boreholes.h>

using namespace std;

#ifndef CPGFUNCTIONEP_SEGMENTS_H
#define CPGFUNCTIONEP_SEGMENTS_H

namespace gt::segments {

    struct SegmentResponse {
        ~SegmentResponse() {} // destructor

        int nSources;
        int nSum;
        vector<gt::boreholes::Borehole> boreSegments;
        vector < vector < double > > h_ij;

        SegmentResponse(int nSources,
                        int nSum,
                        int nt) :
                nSources(nSources),
                boreSegments(nSources),
                h_ij(nSum, vector<double>(nt, 0)),
                nSum(nSum)
        {} // constructor

        // storage_mode = 1 is the reduced segment response vector
        int storage_mode = 1;

//        void ReSizeContainers(int n, int nt);
        void get_h_value(double &h, int i, int j, int k);
        void get_index_value(int &index, int i, int j);
    };  // struct SegmentResponse();


}

#endif //CPGFUNCTIONEP_SEGMENTS_H
