#include <iostream>
#include <stdio.h>
#include "PARA.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
using namespace std;

/* to compile this file: 
g++ -I/usr/local/include -c gsl_Jinv.cpp && g++ -L/usr/local/lib gsl_Jinv.o -lgsl -lgslcblas -lm -o jinv && ./jinv
*/

int main (void){
    double pFpy[2][2] = {{0,0},{0,1}};
    double x = 5.0;
    gsl_matrix *m = gsl_matrix_alloc(N_sp+1, N_sp+1);
    gsl_matrix *m_inv = gsl_matrix_alloc(N_sp+1, N_sp+1);
//*****************************************************************************

// set m:
    for (int ire=0; ire<N_sp+1; ire++){
        for (int isp=0; isp<N_sp+1; isp++){
            if (ire==0 && isp==0){
              gsl_matrix_set(m, ire, isp, 1.);
            }
            else if (ire==0 || isp==0){
              gsl_matrix_set(m, ire, isp, 0.);
            }
            else{
              gsl_matrix_set(m, ire, isp, pFpy[ire][isp]);
            }
        cout <<"ire "<<ire<<"isp "<<isp<< " "<<gsl_matrix_get(m,ire,isp) <<endl;
    }
}
//*****************************************************************************



//*****************************************************************************
// inverse of m:
    gsl_vector *v = gsl_vector_alloc(N_sp+1);

//*****************************************************************************
    gsl_matrix_free(m); gsl_matrix_free(m_inv);
    gsl_vector_free(v);
    return 0;
}