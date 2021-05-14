//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/gesv.h>

extern void print_matrix( std::string desc, int m, int n, double* a, int lda );


int main(){
    int n = 3;  // number of elements in each row (square matrix)
    int nrhs = 1;  // number of columns in b matrix (right hand side)
    int lda = n;
    int ldb = n;
    vector<int> i_piv(n, 0);  // pivot column vector
    int info;
    vector<double> a = {5, 2, 8, 9, 7, 2, 10, 3, 4};
    vector<double> b = {22, 13, 17};

    la::_gesv::gesv(n, nrhs, a, lda, i_piv, b, ldb, info);

    print_matrix("Solution: ", n, nrhs, &*b.begin(), ldb);

    return 0;
}


/* Auxiliary routine: printing a matrix */
void print_matrix( std::string desc, int m, int n, double* a, int lda ) {
    int i, j;
    std::cout << desc << std::endl;
//    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.4f", a[i+j*lda] );
        printf( "\n" );
    }
}