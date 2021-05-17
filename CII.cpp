#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include "gsl_inverse.h"
#include "my_linalg.h"
#include "PARA.h"
#include <iomanip>
#include <time.h>
#include "CII.h"

# include <cmath>
#include <gsl/gsl_matrix.h>
using namespace std;

static double g_0 = 2.0, g_1 = 4.0, A_10 = 2.4e-6, Q_10, C_10, C_01, tau0, tau_cnt;
static double err_eps = 1.e-4;
static int i,j;
static int const N = 1;
static double esc[N], error[N], desc[N], esc_f[N], error_f[N], A[N*N], A_inv[N*N];
static double f_1, f_0;

double CIIcool(double nH, double T_K, double Nc_CII, double y_m, double y_a, double y_e,
               double esc_10, double tau_c) {
//! 1. n_gas, 2. T_K, 3. column density of CII?, 4. y_m=yH2, 5. y_a=yH, 6. y_e 7. 
//  common /line2/Q_10,,C_10,C_01, tau0
   tau_cnt = tau_c;
//     set coefficients

   double DT_10 = 92.;
   double DE_10 = DT_10*k_B;
   double nu_10 = DE_10/h_p;

   Q_10=Q_bg(DT_10); // wli: = 0

   double gamma_e = 2.8e-7/sqrt(T_K*1.e-2);
   double gamma_H = 8.e-10*pow(T_K*1.e-2, 7.e-2);
   double gamma_H2 = 5.e-1*gamma_H;
   C_10 = nH * (y_e*gamma_e+y_a*gamma_H+y_m*gamma_H2);
   C_01 = (g_1/g_0)*C_10*exp(-DT_10/T_K);

   double m_C = 12.*m_H;
   double v_th = sqrt(2.*k_B*T_K/m_C);
   tau0 = (A_10/8./pi)*pow(3.e10/nu_10,3)*Nc_CII/v_th;

   CIIpop(f_0,f_1,esc_10);
   // printf(" final from CII pop : f_0=%15.10e\t,f_1=%15.10e\tesc_10=%15.10e\n",f_0,f_1,esc_10);
   double S_10 = 1.0/(g_1*f_0/(g_0*f_1)-1.0);
   double Ld_CII = DE_10*A_10*f_1*esc_10*(1.-Q_10/S_10)/nH;
   return Ld_CII;
}

void CIIpop(double &f_0, double& f_1, double &esc_10){
//       dimension esc(1),error(1),desc(1),esc_f(1),error_f(1),A(1,1)
   esc[0] = esc_10;
   for (int itr=0;itr<1000;itr++){
      twolevel(f_0,f_1,esc,error);
      // printf("%15.10e\t%15.10e\t%15.10e\t%15.10e\n",f_0,f_1,esc[0],error[0]);

// c     evaluate error
      double err_max = 0.; 
      for (i=0;i<N;i++){
         err_max = max(abs(error[i]),err_max);
      }
      if (err_max < err_eps) break;

      double delta_esc;
      for (i=0;i<N;i++){
         if (esc[i] == 0.) delta_esc = 1.e-10;
         else delta_esc = 1.e-5*esc[i];

         for (int jj=0; jj<N; jj++){
            if (jj == i) esc_f[jj] = esc[jj]+delta_esc;
            else esc_f[jj] = esc[jj];
         }

         twolevel(f_0,f_1,esc_f,error_f);

   // C ------ set the matrix A
         for (j=0;j<N;j++){
            A[i*N+j] = (error_f[i]-error[i]) / (esc_f[j]-esc[j]);
         }
       
      }
      
      // printf("A=%15.10e\n",A[0]);
      // exit(0);

   // C ------- set the vector desc
      for (i=0;i<N;i++){
         desc[i] = -error[i];
      }
    //   printf("desc=%15.10e\n",desc[0]);

// C ------- solve linear equations
//          call gaussj(A,1,1,desc,1,1,0)
      get_inverse(A_inv, A, N);
      dot(desc,N,A_inv,error);

      // printf("A_inv=%15.10e\n",A[0]);
      // printf("desc=%15.10e\t err_max=%15.10e\n",desc[0],err_max);
      // exit(0);

// c     !wli in solving reaction: call gaussj(a=A,n=N_sp,np=N_sp,b=ddy,m=1,mp=1,ind=1)

      double fact = 1.;
      if (itr > 20){
         if (err_max > 1.0) fact = 1.e-2;
      }
      for (i=0;i<N;i++){
         if(esc[i]*desc[i] != 0.) fact = min(fact,0.4*abs(esc[i]/desc[i]));
      }
// c     improve esc
      for (i=0;i<N;i++) esc[i]=esc[i] - fact*desc[i];
   }
      
   esc_10 = esc[0];

   // printf(" final in CII pop : f_0=%15.10e\t,f_1=%15.10e\tesc_10=%15.10e\n",f_0,f_1,esc_10);
   // printf("desc=%15.10e\t err_max=%15.10e\n",desc[0],err_max);
   // exit(0);
}



void twolevel(double &f_0, double &f_1, double *esc,double *error){
//     solves level population (f_0, f_1) of two-level system
//     for given escape probability (esc)
//     returns the difference of old and new esc as (error). 

//     common /line2/Q_10,A_10,C_10,C_01,g_0,g_1,tau0,tau_cnt
    double esc_10 = esc[0];
//     population of level 0 & 1
    double R_01=(g_1/g_0)*A_10*esc_10*Q_10+C_01;
    double R_10=esc_10*A_10*(1.0+Q_10)+C_10;
    f_1=R_01/(R_10+R_01);
    f_0=R_10/(R_10+R_01);
    double tau_10=tau0*(f_0*g_1/g_0-f_1);
    esc[0]=beta_esc(tau_10,tau_cnt);
    error[0] = esc_10-esc[0];
}


double Q_bg(double T_nu){
//      parameter(T_rad=27.3d0)
//      x=T_nu/T_rad
//      if(x.gt.1.d2) then
//         Q_bg_CMB=0.d0
//      else
//         Q_bg_CMB=1.d0/(dexp(x)-1.d0)
//      endif
//      Q_bg=Q_bg_CMB
//      Q_bg=0.d0
    return 0.;
}

double beta_esc(double tau_L, double tau_C){
    if (tau_L < 0) return 1.;
    else if (tau_L < 1.e-5) return exp(-tau_C);
    else return exp(-tau_C)*(1.-exp(-tau_L))/tau_L;
}

/*
g++ -L/usr/local/lib -lgsl -lgslcblas -lm CII.cpp OI.cpp gsl_inverse.o my_linalg.o  -o cooling_cpp.out && ./cooling_cpp.out
*/

// int main(){

//     double xnH = 1.0e3, T_K = 1.0e3, xNc_CII = 1.e10, y_m = 1.e-3;
//     double y_a = 1.0, y_e = 1.e-3, esc_10 = 0.5, tau_c = 0.5;
//     printf("%15.10e\n",CIIcool(xnH, T_K,  xNc_CII,  y_m,  y_a, y_e, esc_10,  tau_c));
//    //  printf("%5.4e\n",beta_esc(.5,.5));
//     return 0;
// }