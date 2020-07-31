// -*- lsst-c++ -*-

//
// Created by jackcook on 7/11/20.
//

#include "gfunction.h"
#include <chrono>
#include "spline.h"
//#include <cstdio>
//#include <cstdlib>
//#include <numeric>
#include "gauss_jacobi.h"
#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>
#include "../jcc/interpolation.h"
#include <thread>
#include <boost/asio.hpp>
//#include <omp.h>

extern "C" void dgesv_( int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );

using namespace std;  // lots of vectors, only namespace to be used

namespace gt::gfunction {
    void uniform_temperature(vector<double>& gfunction, vector<gt::boreholes::Borehole> boreholes,
                             vector<double>& time, const double alpha, const int nSegments,
                             const bool use_similarities, const bool disp) {
        // TODO: place this resizing into a "general" namespace
        int len_actual = time.size(); // check to see if there is enough space in the vector
        int len_g = gfunction.size();
        // if need be, resize the vector to be the same size as the number of boreholes needed
        if (len_actual != len_g) {
            gfunction.resize(len_actual);
        } else ; // else do nothing

        if (disp) {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "Calculating g-function for uniform borehole wall temperature" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        auto startall = std::chrono::steady_clock::now();
        // Open up processes here
        // Create a vector of threads
        //may return 0 when not able to detect
        const auto processor_count = thread::hardware_concurrency();
//        // Launch the pool with n threads.
//        boost::asio::thread_pool pool(processor_count);
//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

        // Number of boreholes
        int nbh = boreholes.size();
        // Total number of line sources
        int nSources = nSegments * nbh;
        // Number of time values
        int nt = time.size();

        // Split boreholes into segments
        vector<gt::boreholes::Borehole> boreSegments(nSources);
        _borehole_segments(boreSegments, boreholes, nSegments);

        // Initialize segment-to-segment response factors (https://slaystudy.com/initialize-3d-vector-in-c/)
        // NOTE: (nt + 1), the first row will be full of zeros for later interpolation
        vector< vector< vector<double> > > h_ij(nSources ,
                vector< vector<double> > (nSources, vector<double> (nt+1, 0.0)) );
        // Calculate segment to segment thermal response factors
        auto start = std::chrono::steady_clock::now();
        gt::heat_transfer::thermal_response_factors(h_ij, boreSegments, time, alpha, use_similarities, disp);
        auto end = std::chrono::steady_clock::now();

        if (disp) {
            std::cout << "Building and solving system of equations ..." << std::endl;
        }
        // -------------------------------------------------------------------------
        // Build a system of equation [A]*[X] = [B] for the evaluation of the
        // g-function. [A] is a coefficient matrix, [X] = [Qb,Tb] is a state
        // space vector of the borehole heat extraction rates and borehole wall
        // temperature (equal for all segments), [B] is a coefficient vector.
        // -------------------------------------------------------------------------

        // -------- timings for debug
        double milli = 0;
        double segment_length_time = 0;
        double time_vector_time = 0;
        double segment_h_values_time = 0;
        double fill_A_time = 0;
        double load_history_reconstruction_time = 0;
        double temporal_superposition_time = 0;
        double fill_gsl_matrices_time = 0;
        double LU_decomposition_time = 0;

        auto start2 = std::chrono::steady_clock::now();

        // ------ Segment lengths -------
        start = std::chrono::steady_clock::now();
        std::vector<float> Hb(nSources);
        auto _segmentlengths = [&boreSegments, &Hb](const int nSources) {
            for (int b=0; b<nSources; b++) {
                Hb[b] = boreSegments[b].H;
            } // next b
        }; // auto _segmentlengths
//        boost::asio::post(pool, [nSources, &boreSegments, &Hb, &_segmentlengths]{ _segmentlengths(nSources); });
        _segmentlengths(nSources);
        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        segment_length_time += milli;

        // ------ time vectors ---------
        start = std::chrono::steady_clock::now();
        // create new time vector that starts at 0
        std::vector<double> _time(time.size()+1);
        std::vector<double> dt(time.size());

        auto _fill_time = [&_time, &time, &dt]() {
            for (int i=0; i<_time.size(); i++) {
                if (i==0) {
                    _time[0] = 0;
                    dt[i] = time[i];
                } else {
                    _time[i] = time[i-1];
                    dt[i] = time[i] - time[i-1];
                } // fi
            } // next i
        }; // auto _fill_time
//        boost::asio::post(pool, [&_fill_time, &time, &_time]{ _fill_time() ;});
        _fill_time();

        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        time_vector_time += milli;

//         TODO: come back to thinking about how to thread this or another (depracated)

//        auto _fill_dt = [&dt, &time]() {
//            for (int i=0; i<time.size(); i++) {
//                if (i==0) {
//                    dt[i] = time[i];
//                } else {
//                    dt[i] = time[i] - time[i-1];
//                } // fi
//            } // end i
//        };
//        _fill_dt();

//        pool.join(); // starting up a new idea after this, pool will close here


//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

        // ---------- segment h values -------------
        /** Starting up pool2 here **/
        // Launch the pool with n threads.
        auto tic = std::chrono::steady_clock::now();
        boost::asio::thread_pool pool2(processor_count);
        auto toc = std::chrono::steady_clock::now();
        if (disp) {
            double milli = std::chrono::duration_cast<std::chrono::milliseconds>(tic - toc).count();
            double seconds = milli;
            std::cout << "Time to open a pool : "
                      << seconds
                      << " sec" << std::endl;
        }

        start = std::chrono::steady_clock::now();
        // h_dt : an interpolation through the k direction of the cube
//        vector< vector< vector<double> > > h_dt(nSources ,
//                vector< vector<double> > (nSources, vector<double> (nt)) );
        vector< vector< vector<double> > > h_dt(nSources ,vector< vector<double> > (nSources, vector<double> (nt)) );
        // dh_ij is the difference of h_ij, ie. dh_[i][j][k] = h_ij[i][j][k]-h_ij[i][j][k-1], ie. \delta h_ij
        std::vector< std::vector< std::vector<double> > > dh_ij(nSources ,
                std::vector< std::vector<double> > (nSources, std::vector<double> (nt)) );

        // Thermal response factors evaluated at t=dt (h_dt)
        auto _interpolate = [&h_ij, &h_dt, &_time, &dt, &dh_ij](const int i) {
            std::vector<double> y(_time.size());
            for (int j=0; j <h_ij[i].size(); j++) {
                for (int k=0; k<h_ij[i][j].size(); k++) {
                    if (k==1) {
                        dh_ij[i][j][k-1] = h_ij[i][j][k];  // start at time "1"
                    } else if (k>1) {
                        dh_ij[i][j][k-1] = h_ij[i][j][k] - h_ij[i][j][k-1];
                    } // fi
                    y[k] = h_ij[i][j][k];
                } // end k

                std::vector<double> yp(_time.size());
                jcc::interpolation::interp1d(dt, yp, _time, y );

                for (int k=0; k<h_dt[i][j].size(); k++) {
                    h_dt[i][j][k] = yp[k];
                } // end k
            } // next j

        }; // _interpolate
        // h_dt for loop
        for (int i=0; i<h_ij.size(); i++) {
//            for (int j=0; j<h_ij[i].size(); j++) {
                boost::asio::post(pool2, [&_interpolate, i]{ _interpolate(i) ;});
//                 _interpolate(i, j);
//            } // end j
        } // end i

        // if _interpolate threaded, join() pools here.
        pool2.join();  // need interpolated values moving forward
        /** Starting up pool3 here **/
        // Launch the pool with n threads.
//        boost::asio::thread_pool pool3(processor_count);
//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

//        auto _dh_ij = [&h_ij, &dh_ij](const int i, const int j) {
//            for (int k=1; k<h_ij[i][j].size(); k++) {
//                if (k==1) {
//                    dh_ij[i][j][k-1] = h_ij[i][j][k];  // start at time "1"
//                } else {
//                    dh_ij[i][j][k-1] = h_ij[i][j][k] - h_ij[i][j][k-1];
//                } // fi
//            } // next k
//        }; // +dh_ij
//        // Thermal response factor increments
//        for (int i=0; i<h_ij.size(); i++) {
//            for (int j=0; j<h_ij[0].size(); j++) {
////                boost::asio::post(pool3, [&_dh_ij, i, j]{ _dh_ij(i, j) ;});
//                _dh_ij(i, j);
//            } // next j
//        } // next i

        // after threaded, join() pool here
//        pool3.join();

        h_ij = vector<vector<vector<double> >> ();  // destruct h_ij; no longer needed

        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        segment_h_values_time += milli;

        // after interpolation scheme, get rid of h_ij first
        // Initialize segment heat extraction rates
        vector<vector<double> > Q(nSources, vector<double> (nt));

        // Define A and b for utitilizing Ax=b
        /**
         * A = [    [   ],
         *          [hb[0:len(hb), 0]
         *     ]
         * b = [    [  ],
         *          [sum(hb)]
         *     ]
         * **/

        int SIZE = nSources + 1;
        // _gesv initializiation
        int nrhs = 1; // number of columns in the b Matrix
        int lda = SIZE;
        int ldb = SIZE;
        std::vector<int> _ipiv(SIZE);
        int info;
        vector<double> A_ (SIZE * SIZE);
        vector<double> b_ (SIZE);
        // LU decomposition gsl initialization
//        gsl_matrix * _A = gsl_matrix_alloc(SIZE, SIZE);
//        gsl_vector * _b = gsl_vector_alloc(SIZE);
//        gsl_vector *_x = gsl_vector_alloc (SIZE);
//        int s;
//        gsl_permutation * _p = gsl_permutation_alloc (SIZE);

        std::vector<std::vector <double>> A(nSources +1, std::vector<double> (nSources +1));
        std::vector<double> b(nSources + 1);

        // Fill A
        int n = SIZE - 1;
//        for (int j=0; j<SIZE; j++) {
//            if (j==n) {
//                gsl_matrix_set (_A, n, n, 0);
//                A[n][n] = 0;
//            } else {
//                gsl_matrix_set (_A, n, j, Hb[j]);
//                A[n][j] = Hb[j];
//            } // fi
//        } // next j

        // Energy conservation: sum([Qb*Hb)] = sum([Hb])
//        int n = A.size() - 1;
//        int n = SIZE -1;
//        _fill_A(n, SIZE);
        n = b.size() - 1;
        double Hb_sum=0;
        for (auto & _hb : Hb) {
            Hb_sum += _hb;
        }
        b[n] = Hb_sum;

        // Build and solve the system of equations at all times

        // the loop p=n depends on what occured at p=n-1, so this will be be in series
        // however threading will be interspersed throughout to make use of as many threads as possible
        // TODO: since this is in series, move the variable declarations here
//        vector<vector<double>> h_ij_dt (nSources, vector<double> (nSources));
        std::vector<double> Tb_0 (nSources);
        // Restructured load history
        // create interpolation object for accumulated heat extraction
        std::vector<std::vector<double>> q_reconstructed (nSources, std::vector<double> (nt));

        for (int p=0; p<nt; p++) {
            // current thermal response factor matrix
            auto _fill_h_ij_dt = [&h_dt, &A] (const int i, const int p) {
                int m = h_dt[0].size();
                for (int j=0; j<A[i].size(); j++) {
                    if (j==A[i].size()-1) {
                        A[i][j] = -1;
                    } else {
                        A[j][i] = h_dt[i][j][p];
                    }
//                    h_ij_dt[j][i] = h_dt[i][j][p];

                } // next j
                ;
            }; // _fill_h_ij_dt
            // _fill_A
            //        auto _fill_A = [&A, &Hb, &_A](const int i, const int SIZE) { // TODO: keep in mind this function can make use of threading

            // ------------- fill A ------------
            start = std::chrono::steady_clock::now();
            auto _fillA = [&Hb, &h_dt, &A, &A_](int i, int p, int SIZE) {
                int n = SIZE - 1;
                for (int j=0; j<SIZE; j++) {
                    if (i == n) { // then we are referring to Hb
                        if (j==n) {
//                            gsl_matrix_set (_A, i, n, 0);
                            A_[i+j*SIZE] = 0;
//                            A[i][n] = 0;
                        } else {
//                            gsl_matrix_set (_A, i, j, Hb[j]);
                            A_[i+j*SIZE] = Hb[j];
//                            A[i][j] = Hb[j];
                        } // fi
                    } else {
                        if (j==A[i].size()-1) {
//                            gsl_matrix_set (_A, i, j, -1);
                            A_[i+j*SIZE] = -1;
//                            A[i][j] = -1;
                        } else {
//                            gsl_matrix_set (_A, j, i, h_dt[i][j][p]);
                            A_[j+i*SIZE] = h_dt[i][j][p];
//                            A[j][i] = h_dt[i][j][p];
                        } // fi
                    } // fi
                } // next k
            };
            boost::asio::thread_pool pool3(processor_count);
            for (int i=0; i<SIZE; i++) {
                boost::asio::post(pool3, [&_fillA, i, p, SIZE]{ _fillA(i, p, SIZE) ;});
//                _fillA(i, p, SIZE);
            }
            pool3.join();
//            n = SIZE - 1;
//            for (int i=0; i<SIZE; i++) {
//                for (int j=0; j<SIZE; j++) {
//                    if (i == n) { // then we are referring to Hb
//                        if (j==n) {
//                        gsl_matrix_set (_A, i, n, 0);
////                            A[i][n] = 0;
//                        } else {
//                        gsl_matrix_set (_A, i, j, Hb[j]);
////                            A[i][j] = Hb[j];
//                        } // fi
//                    } else {
//                        if (j==A[i].size()-1) {
//                        gsl_matrix_set (_A, i, j, -1);
////                            A[i][j] = -1;
//                        } else {
//                        gsl_matrix_set (_A, j, i, h_dt[i][j][p]);
////                            A[j][i] = h_dt[i][j][p];
//                        } // fi
//                    } // fi
//                } // next k
//            }
             // _fill_A
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            fill_A_time += milli;

//            for (int i=0; i<h_dt.size(); i++) {
//                for (int j=0; j<A[i].size(); j++) {
//                    if (j==A[i].size()-1) {
////                        gsl_matrix_set (_A, i, j, -1);
//                        A[i][j] = -1;
//                    } else {
////                        gsl_matrix_set (_A, j, i, h_dt[i][j][p]);
//                        A[j][i] = h_dt[i][j][p];
//                    }
////                    h_ij_dt[j][i] = h_dt[i][j][p];
//                } // next j
//            }
//
//            for (int j=0; j<SIZE; j++) {
//                if (j==n) {
//                    gsl_matrix_set (_A, n, n, 0);
//                    A[n][n] = 0;
//                } else {
//                    gsl_matrix_set (_A, n, j, Hb[j]);
//                    A[n][j] = Hb[j];
//                } // fi
//            } // next j

//            for (int i=0; i<h_dt.size(); i++) {
//                _fill_h_ij_dt(i, p); // TODO: thread this
//            } // next i

//            // fill A with h_ij, and the last row with -1
//            for (int i=0; i<h_ij_dt.size(); i++) {
//                for (int j=0; j<A[i].size(); j++) {
//                    if (j==A[i].size()-1) {
//                        A[i][j] = -1;
//                    } else {
//                        A[i][j] = h_ij_dt[i][j];
//                    } // fi
//                } // next j
//            } // next i

            // TODO: join() threads

            // ----- load history reconstruction -------
            start = std::chrono::steady_clock::now();
            load_history_reconstruction(q_reconstructed,time, _time, Q, dt, p);
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            load_history_reconstruction_time += milli;

            // ----- temporal superposition
            start = std::chrono::steady_clock::now();
            _temporal_superposition(Tb_0, dh_ij, q_reconstructed, p);
            // fill b with -Tb
            for (int i=0; i<Tb_0.size(); i++) {
                b[i] = -Tb_0[i];
            }
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            temporal_superposition_time += milli;

            int m = A.size();
            int n = A[0].size();
            vector<double> x(b.size());
//            _solve_eqn(x, A, b);
            /** was _solve_eqn **/

            // ---- fill gsl matrix A and b -----
            start = std::chrono::steady_clock::now();
//            for (int i=0; i<A.size(); i++) {
//                for (int j=0; j<A[i].size(); j++) {
//                    gsl_matrix_set (_A, i, j, A[i][j]);
//                } // next j
//            } // next i

            for (int i=0; i<b.size(); i++) {
//                gsl_vector_set(_b, i, b[i]);
                b_[i] = b[i];
            } // next i

//            for (int i=0; i<SIZE; i++) {
//                std::cout << i << ' ';
//                for (int j=0; j<SIZE; j++) {
//                    std::cout << gsl_matrix_get(_A, i, j) << ' ';
//                }
//                std::cout << '\n';
//            }

            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            fill_gsl_matrices_time += milli;

            // ----- LU decomposition -----
            start = std::chrono::steady_clock::now();
            dgesv_( &n, &nrhs, &*A_.begin(), &lda, &*_ipiv.begin(), &*b_.begin(), &ldb, &info );

//            gsl_linalg_LU_decomp (_A, _p, &s);
//
//            gsl_linalg_LU_solve (_A, _p, _b, _x);

            for (int i=0; i<SIZE; i++) {
//                x[i] = gsl_vector_get(_x, i);
                x[i] = b_[i];
            } // next i
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            LU_decomposition_time += milli;

            // ---- Save Q's for next p ---
            for (int j=0; j<Q.size(); j++) {
                Q[j][p] = x[j];
            } // next j
            // the borehole wall temperatures are equal for all segments
            double Tb = x[x.size()-1];
            gfunction[p] = Tb;
            int a = 1;
        } // next p
        segment_length_time /= 1000;
        time_vector_time /= 1000;
        segment_h_values_time /= 1000;
        segment_length_time /= 1000;
        fill_A_time /= 1000;
        load_history_reconstruction_time /= 1000;
        temporal_superposition_time /= 1000;
        fill_gsl_matrices_time /= 1000;
        LU_decomposition_time /= 1000;

        cout << "------ timings report -------" << endl;
        cout << " t\t " << " t/p\t" << "name" << endl;
        cout << segment_length_time << "\t" << segment_length_time << "\t" << "segment length time" << endl;
        cout << time_vector_time << "\t" << time_vector_time << "\t" << "time vector time" << endl;
        cout << segment_h_values_time << "\t" << segment_h_values_time << "\t" << "segment h values time" << endl;
        cout << fill_A_time << "\t" << fill_A_time / double(nt) << "\t" << "time to fill vector A" << endl;
        cout << load_history_reconstruction_time << "\t" << load_history_reconstruction_time / double(nt)
        << "\t" << "load hist reconstruction" << endl;
        cout << temporal_superposition_time << "\t" << temporal_superposition_time / double(nt)
             << "\t" << "temporal superposition time:" << endl;
        cout << fill_gsl_matrices_time << "\t" << fill_gsl_matrices_time / double(nt)
             << "\t" << "gsl fill matrices time" << endl;
        cout << LU_decomposition_time << "\t" << LU_decomposition_time/double(nt)
             << "\t" << "LU decomp time" << endl;

//        gsl_permutation_free (_p);
//        gsl_vector_free (_x);
//        gsl_vector_free (_b);
//        gsl_matrix_free (_A);

        auto end2 = std::chrono::steady_clock::now();
        if (disp) {
            double milli1 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
            double seconds1 = milli1 / 1000;
            double milli2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - startall).count();
            double seconds2 = milli2 / 1000;
            std::cout << "Elapsed time in seconds : "
                      << seconds1
                      << " sec" << std::endl;
            std::cout << "Total time for g-function evaluation : "
                    << seconds2
                    << " sec" << std::endl;
        }
        int a = 1;
    } // void uniform_temperature

    void _borehole_segments(std::vector<gt::boreholes::Borehole>& boreSegments,
            std::vector<gt::boreholes::Borehole>& boreholes, const int nSegments) {
        double H;
        double D;
        int count = 0;
        // Split boreholes into segments
        for(auto& b : boreholes) {
            // TODO: maybe thread this later on
            for (int i=0; i<nSegments; i++) {
                H = b.H / double(nSegments);
                D = b.D + double(i) * b.H / double(nSegments);
                boreSegments[count] = gt::boreholes::Borehole(H, D, b.r_b, b.x, b.y);
                count++;
            }  // end for
        } // end for
    } // void _borehole_segments

    void load_history_reconstruction(std::vector<std::vector<double>>& q_reconstructed,
            vector<double>& time, vector<double>& _time, vector<vector<double> >& Q,
            vector<double>& dt, const int p) {
        // for s in range p+1
        int nSources = Q.size();

        // Inverted time steps
        std::vector<double> dt_reconstructed (p+1);
        for (int i=p; i>=0; i--) {
            dt_reconstructed[p-i] = dt[i];  // reverse the dt
//            if (i==0) {
//                dt_reconstructed[i] = time[0];
//            } else {
//                dt_reconstructed[i] = dt[p-i+1];  // reverse the dt
//            }
        }
        // t_restructured is [0, cumsum(dt_reversed)]
        std::vector<double> t_reconstructed(p+2);  // will start at 0
        for (int i=1; i<=p+1; i++) {
            t_reconstructed[i] = dt_reconstructed[i-1];
        }
        for (int i=1; i<=p+1; i++) {
            t_reconstructed[i] = t_reconstructed[i] + t_reconstructed[i-1];
        }
         // local time vector
        std::vector<double> t(p+3);
        for (int i=0; i<t.size(); i++) {
            if (i==t.size()-1) {
                t[i] = _time[i-1] + _time[1];
            } else {
                t[i] = _time[i];
            }
        }
        int _tsize = t.size();
        // Q*dt
        std::vector<std::vector <double> > Q_dt (nSources, std::vector<double> (t.size()));
        auto _Q_dot_dt = [&Q_dt, &Q, &dt, &p, &_tsize, &t](const int i) {
            for (int j = 1; j<_tsize; j++) {
                if (j>=p+1) {
                    Q_dt[i][j] = Q_dt[i][j-1];
                } else {
                    double aa =Q[i][j-1];
                    double ab = dt[j-1];
                    Q_dt[i][j] = Q[i][j-1] * dt[j-1] + Q_dt[i][j-1];
                }  // fi
            } // next j
        };
        for (int i=0; i<nSources; i++) {
            _Q_dot_dt(i);  // TODO: thread this
        } // next i
        int a =1;
        // TODO: join()

        auto _interpolate = [&Q_dt, &q_reconstructed, &t, &t_reconstructed, &dt_reconstructed, &p](const int i) {
            int n = t.size();
            std::vector<double> y(n);
            for (int j=0; j<n; j++) {
                y[j] = Q_dt[i][j];
            }

            std::vector<double> yp(n);
//            std::cout << p << std::endl;
            if (p==2) {
                int a = 1;
            }
            jcc::interpolation::interp1d(t_reconstructed, yp, t, y);

            n = q_reconstructed[0].size();
            for (int j=0; j<p; j++) {
                double c = yp[j];
                double d = yp[j+1];
                double e = dt_reconstructed[j];
//                q_reconstructed[i][j] = (s(t_reconstructed[j+1]) - s(t_reconstructed[j])) / dt_reconstructed[j];
                q_reconstructed[i][j] = (d - c) / e;
            }
        }; // _interpolate
        for (int i=0; i<nSources; i++) {
            _interpolate(i); // TODO: thread this
        }
        // TODO: join()
    } // load_history_reconstruction

    void _temporal_superposition(vector<double>& Tb_0, vector<vector<vector<double> > >& dh_ij,
            std::vector<std::vector<double>>& q_reconstructed, int p)
            {
        const auto processor_count = thread::hardware_concurrency();
        // Launch the pool with n threads.
        boost::asio::thread_pool pool(processor_count);
//

        std::fill(Tb_0.begin(), Tb_0.end(), 0);
//        for (double & i : Tb_0) { // set all values in Tb_0 = 0
//            i = 0;
//        }
        // Number of heat sources
        int nSources = q_reconstructed.size();
        // Number of time steps
//        int nt = q_reconstructed[0].size();
        int nt = p + 1;

//        int k, j;
//#pragma omp parallel for private(j)
//        for (int k=0; k<nt; k++) {
//            for (int j=0; j<nSources; j++) {
//                for (int i=0; i<nSources; i++) {
////#pragma omp critical (TEMPORAL_CRITICAL)
//                    {
//                        Tb_0[i] += dh_ij[i][j][k] * q_reconstructed[j][nt-k-1] ;
//                    }
//
//                } // next i
//            } // next j
//        } // next k
//#pragma omp parallel for collapse(3)
//        for (int i=0; i<nSources; i++) {
//            for (int j=0; j<nSources; j++) {
//                for (int k=0; k<nt; k++) {
//                    Tb_0[i] += dh_ij[i][j][k] * q_reconstructed[j][nt-k-1] ;
//                }
//            }
//        }
        auto _borehole_wall_temp = [&dh_ij, &q_reconstructed, &Tb_0](const int i, const int nSources, const int nt){
            for (int j =0; j<nSources; j++) {
                for (int k=0; k<nt; k++) {
                    Tb_0[i] += dh_ij[i][j][k] * q_reconstructed[j][nt-k-1] ;
                }
            }
        };
        for (int i=0; i<nSources; i++) {
            boost::asio::post(pool, [&_borehole_wall_temp, i, nSources, nt]
            { _borehole_wall_temp(i, nSources, nt); });
//            _borehole_wall_temp(i, nSources, nt);
        }

//        for (int i=0; i<nSources; i++) {
//            for (int k=0; k<nt; k++) {
//                for (int j=0; j<nSources; j++) {
//                    double a = dh_ij[i][j][k];
//                    double ab = q_reconstructed[i][nt-j-1];
//                    Tb_0[i] += dh_ij[i][j][k] * q_reconstructed[i][nt-k-1];
//                } // next j
//            } // next k
//        } // next i

        pool.join();

        int a = 1;
    }  // _temporal_superposition

//    void _solve_eqn(std::vector<double>& x, std::vector<std::vector<double>>& A, std::vector<double>& b) {
//        int SIZE = x.size();
//        gsl_matrix * _A = gsl_matrix_alloc(SIZE, SIZE);
//        for (int i=0; i<A.size(); i++) {
//            for (int j=0; j<A[i].size(); j++) {
//                gsl_matrix_set (_A, i, j, A[i][j]);
//            } // next j
//        } // next i
//
//        gsl_vector * _b = gsl_vector_alloc(SIZE);
//        for (int i=0; i<b.size(); i++) {
//            gsl_vector_set(_b, i, b[i]);
//        } // next i
//
//        gsl_vector *_x = gsl_vector_alloc (SIZE);
//
//        int s;
//
//        gsl_permutation * p = gsl_permutation_alloc (SIZE);
//
//        gsl_linalg_LU_decomp (_A, p, &s);
//
//        gsl_linalg_LU_solve (_A, p, _b, _x);
//
//        for (int i=0; i<SIZE; i++) {
//            x[i] = gsl_vector_get(_x, i);
//        } // next i
//
//        gsl_permutation_free (p);
//        gsl_vector_free (_x);
//    } // _solve_eq

} // namespace gt::gfunction
