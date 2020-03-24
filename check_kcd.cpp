#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include<algorithm>

#include "reaction.h"
#include "my_linalg.h"
#include "Newton.h"
#include "PARA.h"
using namespace std;

/* g++ -Wall -I/usr/local/include -c check_kcd.cpp 
g++  -L/usr/local/lib check_kcd.o Newton.o PARA.o my_linalg.o gsl_inverse.o reaction.o -lgsl -lgslcblas -lm -o check_kcd
./check_kcd
*/

//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
int N = 10000, i;
double y_H2 = 1.e-2;
double y_e = 1.e-4;
double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2p = 1.0e-12, y_Hm = 1.0e-12;
double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
double ys[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};
double gamma_h1, gamma_l1, gamma_h2, gamma_l2, nc_1, nc_2, 
        pp, gamma_sum, xk_L, xk_H, n_cr, a, n_cr_H, n_cr_H2, 
        zH, zH2, chi_H2, z_n;
double k;
int n_max=5, n;


ofstream myfile;
double xk_Martin(double y_H2, double y_H, double nH, double T_K){
    double lgT = log10(T_K);
    gamma_h1=-1.784239e2-6.842243e1*log10(T_K)
    +4.320243e1*pow(lgT,2)-4.633167*pow(lgT,3)
    +6.970086e1*log10(1.0+4.087038e4/T_K);
    gamma_h2=-2.370570e4/T_K;
    gamma_l1=1.288953e2-5.391334e1*lgT
        +5.315517*pow(lgT,2)
        -1.973427e1*log10(1.0+1.678095e4/T_K);
    gamma_l2=-2.578611e4/T_K;
    nc_1=1.482123e1-4.890915*lgT
        +4.749030e-1*pow(lgT,2)-1.338283e2/T_K;
    nc_2=-1.164408+nc_1;
    pp=8.227443e-1+5.864073e-1*exp(-T_K/1850)
        -2.056313*exp(-T_K/440);
    nc_1=pow(1.e1,nc_1);
    nc_2=pow(1.0e1,nc_2);
    gamma_sum=gamma_h1
        -(gamma_h1-gamma_l1)/(1.0+pow(nH/nc_1,pp))
        +gamma_h2
        -(gamma_h2-gamma_l2)/(1.0+pow(nH/nc_2,pp));
    k=pow(1.0e1,gamma_sum);
    return k;
}
double xk_GA(double y_H2, double y_H, double nH, double T_K){
    chi_H2  =  5.1965e4; //binding energy of H2
    double lgT4 = log10(T_K/1.e4);
    n_cr_H  = pow( 10., 3.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) );
    n_cr_H2 = pow( 10., 4.845 - 1.3*lgT4 + 1.62*pow(lgT4,2) );
    n_cr = 1.0/(y_H/n_cr_H + 2.0*y_H2/n_cr_H2); // Eq(14) GA08
    a=1./(1.+nH/n_cr);

    xk_L  = 6.67e-12*sqrt(T_K)/exp(1 + 63593./T_K);
    xk_H = 3.52e-9*exp(-43900./T_K);
    xk_H = 6.5e-7/sqrt(T_K)/exp(chi_H2/T_K)*(1.-exp(-6000./T_K));//Omukai 2001
//     connect between v=0 and LTE  
    if(a==1.) k = xk_L;
    else if (a==0.) k = xk_H;
    else k = pow(xk_H, 1.-a) * pow(xk_L,a);
    return k;
}

double xk_GA_O(double y_H2, double y_H, double nH, double T_K){
    chi_H2  =  5.1965e4; //binding energy of H2
    double lgT4 = log10(T_K/1.e4);
    n_cr_H  = pow( 10., 3.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) );
    n_cr_H2 = pow( 10., 4.845 - 1.3*lgT4 + 1.62*pow(lgT4,2) );
    n_cr = 1.0/(y_H/n_cr_H + 2.0*y_H2/n_cr_H2); // Eq(14) GA08
    n_cr = pow( 10., 4.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) ); //Omukai 2001
    a=1./(1.+nH/n_cr);

    xk_L  = 1.12e-10*exp(-7.035e4/T_K);//Omukai 2001
    xk_H = 6.5e-7/sqrt(T_K)/exp(chi_H2/T_K)*(1.-exp(-6000./T_K));//Omukai 2001
//     connect between v=0 and LTE  
    if(a==1.) k = xk_L;
    else if (a==0.) k = xk_H;
    else k = pow(xk_H, 1.-a) * pow(xk_L,a);
    return k;
}

double xk_back(double kcd, double T_K){
    double lgT = log10(T_K); 
    chi_H2  =  5.1965e4; //binding energy of H2
    zH2 = pow( 10., 2.20859 -1.8089*lgT + 0.451858* pow(lgT,2.0) );
    zH=0.0;
    for (n=1; n<=n_max; n++){
        z_n = 2.*pow(n,2)*exp(-1.57798e5*(1.-1./pow(n,2))/T_K);
        zH = zH+z_n;
    }
    k = kcd *zH2/pow(zH,2)*1.493e-20 /pow(T_K,1.5)*exp(chi_H2/T_K);
    printf("%3.2e, %3.2e, %3.2e\n",lgT, zH2, k);
    return k;
}

int main(){
    double nH0 = .1; double T_K0 = 200;
    double nH = nH0;
    double T_K = T_K0;
//kcd v.s. nH
    myfile.open("../data/kcd_nH.txt", ios::out | ios::trunc );
    myfile<<"nH\tGA\tOmukai\tMartin\t3bGA\t3bOmukai\t3bMartin"<<endl;
    T_K = 1.e4;
    for (i=0; i<N; i++){
        nH *= 1.003;
        myfile <<"\t"<<nH;
        myfile  <<"\t"<<xk_GA(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_GA_O(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_Martin(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_back(xk_GA(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                <<"\t"<<xk_back(xk_GA_O(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                <<"\t"<<xk_back(xk_Martin(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                << endl;
    }
    myfile.close();

//kcd v.s. T
    myfile.open("../data/kcd_T.txt", ios::out | ios::trunc );
    myfile<<"nH\tGA\tOmukai\tMartin\t3bGA\t3bOmukai\t3bMartin"<<endl;
    nH = 3000;
    T_K = 10.;
    for (i=0; i<N/2; i++){
        T_K *= 1.003;
        myfile <<"\t"<<T_K;
        myfile  <<"\t"<<xk_GA_O(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_GA_O(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_Martin(y_H2, y_H, nH,T_K)*y_H2*y_H*pow(nH,2)
                <<"\t"<<xk_back(xk_GA(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                <<"\t"<<xk_back(xk_GA_O(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                <<"\t"<<xk_back(xk_Martin(y_H2, y_H, nH,T_K), T_K)*pow(y_H*nH,3)
                << endl;
    }
    myfile.close();
    return 0;
}