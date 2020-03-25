// Newton method
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include "Newton.h"
#include "gsl_inverse.h"
#include "my_linalg.h"
#include "reaction.h"
#include "PARA.h"

# include <cmath>
#include <gsl/gsl_matrix.h>
using namespace std;

int isp, ire, jsp;
void get_inverse (double *mA_inv, double *mA, int const N);
int const N = N_sp +1;

//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
//  1)   H+    +   e     ->   H     +   ph.
//  2)   H     +   e     ->   H-    +   ph.
//  3)   H-    +   H     ->   H2    +   e              
//  4)   H2  +   H    <-> 3 H  
//  5)   3 H       ->     H2     +      H
//  6)    H2    +     J21        ->         2 H  
//  7)    H-    +     ph       ->   H    +     e
//  8)    H    +     e-       ->   H+    +     2e-

void SOL_IMPLICIT(double* dy, double *y0, double* y1, double dt, double nH, double T_K, double* xk, double J_LW, double Tb){
    double* j = new double [N*N]; double* j_inv = new double [N*N];
    double F[N];
    double r_f[N], r_f_fw[N], r_f_bw[N];
    // double xk[N_react1];
    double delta_y;
    double y_tmp[N], err_0, r_f_big, dr_f,
           y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp;
    double dr_fdy[N][N];
    
    react_coef(xk, nH, y1[1], y1[2], T_K, J_LW, Tb);
    react_rat(r_f, y1, xk, nH, T_K);
    //printf("NEWTON: k_pdH2=%3.2e, k_pdHm=%3.2e, k_pdH2p=%3.2e\n",xk[21],xk[22],xk[23]);

    // numerically calculate dr_fdy 
    for (jsp=1;jsp<N;jsp++){
        if(abs(y1[jsp]) < 1.e-100) delta_y = epE10;
        else delta_y = epE8 * y1[jsp];
        for (isp=1; isp<N; isp++){
            if(isp == jsp) y_tmp[isp] = y1[isp]+ .5*delta_y;
            else y_tmp[isp] = y1[isp];
        }

        react_rat(r_f_fw, y_tmp, xk, nH, T_K);
        y_tmp[jsp] = y1[jsp] - 0.5 *delta_y;
        react_rat(r_f_bw, y_tmp, xk, nH, T_K);

        for (isp=1; isp<N; isp++){
            dr_fdy[isp][jsp] = (r_f_fw[isp] - r_f_bw[isp]) /delta_y;
            //printf("dr_fdy = %3.2e\n",dr_fdy[isp][jsp]);
            r_f_big = max( abs(r_f_fw[isp]) , abs(r_f_bw[isp]) );
            if (r_f_big != 0.) {
                dr_f = abs(r_f_fw[isp] -r_f_bw[isp] );
                if (dr_f/r_f_big < 1.e-15) dr_fdy[isp][jsp] = 0.;
            }
        }
    }

    // set Jij = partial dyi/dt partial yj
    j[0] = 1.;
    for (jsp=1; jsp<N; jsp++) j[jsp] = 0.;
    for (isp=1; isp<N; isp++) j[isp*N] = 0.;

    for (isp = 1; isp<N; isp++){
        //printf("[");
        for (jsp = 1; jsp<N; jsp++){
            if(isp != jsp) j[isp*N+jsp] = -dt*dr_fdy[isp][jsp];
            else j[isp*N+jsp] = 1. - dt*dr_fdy[isp][jsp];
            //printf("%3.2e, ",j[isp*N+jsp]);
        }
        //printf("]\n");
    }
    //printf("\n");
    get_inverse (j_inv, j, N);

    F[0] = 0;
    for (isp=1;isp<N;isp++) F[isp] = y1[isp] - dt*r_f[isp] -y0[isp];

    dot(dy, N, j_inv, F);
    for (isp=1;isp<N;isp++) y1[isp] -= dy[isp];

    delete [] j; delete [] j_inv;
}