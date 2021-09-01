//
// Created by jackcook on 8/31/21.
//

#include <cpgfunction/utilities.h>
#include <cpgfunction/boreholes.h>
#include <cpgfunction/heat_transfer.h>

#include <vector>

double absolute_error(std::vector<double> actual, std::vector<double> predicted) {
    double summation = 0.;
    for (int i=0; i<actual.size(); i++) {
        summation += std::abs(actual[i] - predicted[i]);
    }
    return summation / double(actual.size());
}


int main() {
    double r_b = 0.075;  // Borehole radius (meters)
    double H = 150.;  // Height of boreholes (meters)
    double alpha = 1.0e-06;  // Thermal diffusivity

    std::vector<double> time = gt::utilities::time_Eskilson(H, alpha);
    std::vector<double> h_ij_real(time.size(), 0.);
    std::vector<double> h_ij_mirror(time.size(), 0.);
    std::vector<double> h_ij_combined(time.size(), 0.);

    bool realSource;
    bool imagSource;

    // Case A
    double H_i = 10.; // Height of borehole i (meters)
    double D_i = 4.; // Burial depth of borehole i (meters)
    double x_i = 0.;  // x-coordinate of borehole i (meters)
    double y_i = 0.;  // y-coordinate of borehole i (meters)

    gt::boreholes::Borehole segment_i =
            gt::boreholes::Borehole(H_i, D_i, r_b, x_i, y_i);

    double H_j = 10; // Height of borehole j (meters)
    double D_j = 144.; // Burial depth of borehole j (meters)
    double x_j = 67.5; // x-coordinate of borehole j (meters)
    double y_j = 67.5; // y-coordinate of borehole j (meters)

    gt::boreholes::Borehole segment_j =
            gt::boreholes::Borehole(H_j, D_j, r_b, x_j, y_j);

    int N=10;
    double h;
    gt::heat_transfer::FLSApproximation FLSApprox = gt::heat_transfer::FLSApproximation(N);
    std::vector<double> d_ = FLSApprox.construct_dm(segment_i, segment_j);
    std::vector<double> h_ij_combined_approx;
    std::vector<double> h_ij_real_approx;
    std::vector<double> h_ij_image_approx;

    for (int k=0; k<time.size(); k++){
        h_ij_combined_approx.push_back(FLSApprox.finite_line_source(time[k], alpha, segment_i, segment_j, d_, true, true));
        h_ij_real_approx.push_back(FLSApprox.finite_line_source(time[k], alpha, segment_i, segment_j, d_, true, false));
        h_ij_image_approx.push_back(FLSApprox.finite_line_source(time[k], alpha, segment_i, segment_j, d_, false, true));
    }

    // Compute only the real source
    realSource = true;
    imagSource = false;
    for (int k=0; k<time.size(); k++)
        h_ij_real[k] = gt::heat_transfer::finite_line_source(time[k], alpha,
                                                             segment_i,
                                                             segment_j,
                                                             realSource,
                                                             imagSource);

    // Compute only the mirror source
    realSource = false;
    imagSource = true;
    for (int k=0; k<time.size(); k++)
        h_ij_mirror[k] = gt::heat_transfer::finite_line_source(time[k], alpha,
                                                             segment_i,
                                                             segment_j,
                                                             realSource,
                                                             imagSource);

    // Compute the combined segment response
    realSource = true;
    imagSource = true;
    for (int k=0; k<time.size(); k++)
        h_ij_combined[k] = gt::heat_transfer::finite_line_source(time[k], alpha,
                                                             segment_i,
                                                             segment_j,
                                                             realSource,
                                                             imagSource);

    double rmse_val;

    rmse_val = absolute_error(h_ij_real, h_ij_real_approx);
    std::cout << "Absolute error: " <<rmse_val * 100. << std::endl;
    rmse_val = absolute_error(h_ij_mirror, h_ij_image_approx);
    std::cout << "Absolute error: " <<rmse_val * 100. << std::endl;
    rmse_val = absolute_error(h_ij_combined, h_ij_combined_approx);
    std::cout << "Absolute error: " <<rmse_val * 100. << std::endl;


    return 0;
}