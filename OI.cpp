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
#include "OI.h"

# include <cmath>
#include <gsl/gsl_matrix.h>
using namespace std;

static double g_0 = 5.0, g_1 = 3.0, g_2 = 1.0;
static double nu_10,nu_20,nu_21, Q_10,Q_20,Q_21, A_10,A_20,A_21,C_10,C_20,C_21,C_01,C_02,C_12;
static double err_eps = 1.e-5;
static double Nclmn, v_th, tau_cnt;
static int i,j;
static int const N = 3;
static double esc[N], error[N], desc[N], esc_f[N], error_f[N], A[N*N], A_inv[N*N];
static double f_0, f_1, f_2;

double OIcool (double nH,double T_K,double Nc_OI,double y_m,double y_a,double y_e,double *esc,double tau_c){
   tau_cnt = tau_c;
   if (T_K < 10.) return 0.;
//      set coefficients

   double DT_10 = 2.3e2;
   double DT_20 = 3.28e2;
   double DT_21 = 9.8e1;

   double DE_10 = DT_10*k_B;
   double DE_20 = DT_20*k_B;
   double DE_21 = DT_21*k_B;
    
   nu_10 = DE_10/h_p;
   nu_20 = DE_20/h_p;
   nu_21 = DE_21/h_p;

   A_10 = 9.0e-5;
   A_20 = 1.0e-10;
   A_21 = 1.7e-5;

   Q_10 = Q_bg(DT_10);
   Q_20 = Q_bg(DT_20);
   Q_21 = Q_bg(DT_21);

   double gamma10_e = 1.4e-8;
   double gamma20_e = 1.4e-8;
   double gamma21_e = 5.0e-9;

   double gamma10_H = 9.2e-11*pow(T_K*1.e-2, 0.67);
   double gamma20_H = 4.3e-11*pow(T_K*1.e-2, 0.8);
   double gamma21_H = 1.1e-10*pow(T_K*1.e-2,0.44);

   double gamma10_H2 = 5.e-1*gamma10_H;
   double gamma20_H2 = 5.e-1*gamma20_H;
   double gamma21_H2 = 5.e-2*gamma21_H;

   C_10 = nH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2);
   C_20 = nH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2);
   C_21 = nH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2);
   C_01=(g_1/g_0)*C_10*exp(-DT_10/T_K);
   C_02=(g_2/g_0)*C_20*exp(-DT_20/T_K);
   C_12=(g_2/g_1)*C_21*exp(-DT_21/T_K);

   Nclmn = Nc_OI;
   double m_O = 16.*m_H;
   v_th = sqrt(2.*k_B*T_K/m_O);

   pop3lev(f_0,f_1,f_2,esc);

   double S_10 = 1./(g_1*f_0/(g_0*f_1)-1.);
   double S_20 = 1./(g_2*f_0/(g_0*f_2)-1.);
   double S_21 = 1./(g_2*f_1/(g_1*f_2)-1.);

   double x_10 = esc[0]*(1.-Q_10/S_10);
   double x_20 = esc[1]*(1.-Q_20/S_20);
   double x_21 = esc[2]*(1.-Q_21/S_21);

   double Ld_OI = (DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20+DE_21*A_21*f_2*x_21)/nH;
   return Ld_OI;

}
      

void pop3lev(double &f_0, double &f_1, double &f_2, double *esc){
   //     determin population under given esc
   for (int itr=0;itr<100;itr++){
      thrlev(f_0,f_1,f_2,esc,error);

      //     evaluate error
      double err_max = 0.;
      for (i=0;i<N;i++) err_max = max(abs(error[i]),err_max);
   //     error small enough ?
      if (err_max < err_eps) break;
   //     if not, improve guess for esc
   //     by Newton-Raphson scheme
      double delta_esc;
      for (j=0;j<N;j++){
         if(esc[j] == 0.){
            delta_esc = 1.e-10;
         }
         else {
            delta_esc=1.e-5*esc[j];
         }

         for (int jj=0;jj<N;jj++){
            if (jj==j){
               esc_f[jj] = esc[jj]+delta_esc;
            }
            else esc_f[jj] = esc[jj];
         }
       
         thrlev(f_0,f_1,f_2,esc_f,error_f);
   // ------ set the matrix A
         for (i=0;i<N;i++) A[i*N+j] = (error_f[i]-error[i])/(esc_f[j]-esc[j]);
      }

   // ------- set the vector desc
      for (i=0;i<N;i++) desc[i] = -error[i];
      //   printf("desc=%15.10e\n",desc[0]);

   // C ------- solve linear equations
   //          call gaussj(A,1,1,desc,1,1,0)
      get_inverse(A_inv, A, N);
      dot(desc,N,A_inv,error);

         // printf("A_inv=%15
   // ------- solve linear equations
      // !!!!!!!!!! gaussj(A,3,3,desc,1,1,0)
   //     !wli in solving reaction: call gaussj(a=A,n=N_sp,np=N_sp,b=ddy,m=1,mp=1,ind=1)

       double fact = 1.;
  
      if (itr > 20){
         if (err_max > 1.) fact = 1.e-2;
      }

      for (i=0;i<N;i++){
         if (esc[i]*desc[i] != 0.){
            fact = min(fact,0.4*abs(esc[i]/desc[i]));
         }
      }

   //     next guess for esc
      for (i=0;i<N;i++) esc[i] = esc[i] - fact*desc[i];

   }

}


void thrlev(double &f_0, double &f_1, double &f_2, double *esc, double *error){
   double esc_10 = esc[0];
   double esc_20 = esc[1];
   double esc_21 = esc[2];
   
   
   double R_10 = esc_10*A_10*(1.+Q_10)+C_10;
   double R_20 = esc_20*A_20*(1.+Q_20)+C_20;
   double R_21 = esc_21*A_21*(1.+Q_21)+C_21;
   double R_01 = (g_1/g_0)*esc_10*A_10*Q_10+C_01;
   double R_02 = (g_2/g_0)*esc_20*A_20*Q_20+C_02;
   double R_12 = (g_2/g_1)*esc_21*A_21*Q_21+C_12;
   
   f_0 = ( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
         /( (R_01+R_02+R_20)*(R_10+R_12+R_21) 
         -(R_01-R_21)*(R_10-R_20) );
   f_1 = (f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21);
   f_2 = (f_0*R_02+f_1*R_12)/(R_21+R_20);

   double   Nc_0=Nclmn*f_0;
   double   Nc_1=Nclmn*f_1;
   double   Nc_2=Nclmn*f_2;

   double tau_10=(A_10/8./pi)*pow(3.e10/nu_10,3)
         *(Nc_0*g_1/g_0-Nc_1)/v_th ;
   double tau_20=(A_20/8./pi)*pow(3.e10/nu_20,3)
         *(Nc_0*g_2/g_0-Nc_2)/v_th ;
   double tau_21=(A_21/8./pi)*pow(3.e10/nu_21,3)
         *(Nc_1*g_2/g_1-Nc_2)/v_th ;

   double   err_10 = esc_10-beta_esc(tau_10,tau_cnt);
   double   err_20 = esc_20-beta_esc(tau_20,tau_cnt);
   double   err_21 = esc_21-beta_esc(tau_21,tau_cnt);
      
   error[0] = err_10;
   error[1] = err_20;
   error[2] = err_21;

}

/*
g++ -L/usr/local/lib -lgsl -lgslcblas -lm CII.cpp OI.cpp gsl_inverse.o my_linalg.o  -o OI.out && ./OI.out
*/

// int main(){

//     double xnH = 1.0e3, T_K = 1.0e3, Nc_OI = 1.e10, y_m = 1.e-3;
//     double y_a = 1.0, y_e = 1.e-3, esc_10 = 0.5, tau_c = 0.;
//     for (i=0;i<N;i++) esc[i] = 0.5;
//     printf("%15.10e\n", OIcool(xnH,T_K,Nc_OI,y_m,y_a,y_e,esc,tau_c));
//     printf("%5.4e  %5.4e  %5.4e\n", esc[0],esc[1],esc[2]);
//     return 0;
// }