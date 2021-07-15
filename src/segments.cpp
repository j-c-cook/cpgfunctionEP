//
// Created by jackcook on 7/15/21.
//

#include <iostream>
#include <cpgfunction/segments.h>

using namespace std;

namespace gt::segments {

    void SegmentResponse::get_h_value(double &h, const int i, const int j,
                                      const int k) {
        int index;
        switch (storage_mode) {
            case 0 :
                cout << "Case 0 not written yet" << endl;
                break;
            case 1 :
                if (i <= j) {
                    get_index_value(index, i, j);
                    h = h_ij[index][k];
                } else {
                    get_index_value(index, j, i);
                    h = boreSegments[j].H/boreSegments[i].H * h_ij[index][k];
                }
                break;
            default:
                throw invalid_argument("The case selected is not currently "
                                       "implemented.");
        }  // switch();
    }  // SegmentResponse::get_h_value();

    void SegmentResponse::get_index_value(int &index, const int i, const int j) {
        index = i * (2*nSources - i - 1) / 2 + j;
    }  // SegmentResponse::get_index_value();

}