//compare reaction coefficients original v.s. Glover & Abel 2008
// also check n_cr=Lambda_LTE/Lambda_0 with in reaction.cpp

// g++ check_cool.cpp PARA.cpp thermo.cpp -o check_cool && ./check_cool

# include <iostream>
# include <stdio.h>
#include <fstream>
#include <string.h>
# include <cmath>
#include<algorithm>

#include "PARA.h"
#include "thermo.h"
using namespace std;

double LambdaGA_H2(double nH, double T_K, double* y){
    double y_H = y[1], y_H2 = y[2], y_He = y[7], y_Hp = y[4], y_e = y[3];
    double T3 = T_K/1.e3;
    //double logT3 = log10(T3);
    double logT3 = log(T3);
    double LH=0, LH2=0, LHe=0, LHp=0, Le=0; // cooling by collisions with H, H2, H+, e-
    double L0, L_LTE, L, ftau;
    int const n_a = 6; //must be const.
    double a[n_a] = {0,0,0,0,0,0};
    int i, j;
//  1) H
    if (T_K<=10.) for (i=0;i<n_a;i++) a[i]=0.;
    else if (10.<T_K && T_K<=100.){
        a[0] = -16.818342;
        a[1] =  37.383713;
        a[2] = 58.145166;
        a[3] = 48.656103;
        a[4] = 20.159831;
        a[5] =  3.8479610;
    }
    else if (100.<T_K && T_K<=1000.){
        a[0] = -24.311209;
        a[1] = 3.5692468;
        a[2] = -11.332860;
        a[3] = -27.850082;
        a[4] = -21.328264;
        a[5] = 4.2519023;
    }
    else if (1000.<T_K && T_K<=6000.){  // wli relaxed the criteria in original Table
        a[0] = -24.311209;
        a[1] = 4.6450521;
        a[2] = -3.7209846;
        a[3] = 5.9369081;
        a[4] = -5.5108047;
        a[5] = 1.5538288;
    }
    else for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    //for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    
    for (i=0; i<n_a; i++) LH += a[i]*pow(logT3,i);
//  2) H2 // wli criteria cancelled
    if (100.<T_K && T_K<=6000.){
        a[0] = -23.962112;
        a[1] = 2.09433740;
        a[2] = -0.77151436;
        a[3] =  0.43693353;
        a[4] = -0.14913216;
        a[5] = -0.033638326;
    }
    else for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    //for (i=0;i<n_a;i++) a[i]=0.; // set to 0

    a[0] = -23.962112;
    a[1] = 2.09433740;
    a[2] = -0.77151436;
    a[3] =  0.43693353;
    a[4] = -0.14913216;
    a[5] = -0.033638326;

    for (i=0; i<n_a; i++) LH2 += a[i]*pow(logT3,i);

//  3) He // wli criteria cancelled 
    if (10.<T_K && T_K<=6000.){
        a[0] = -23.689237;
        a[1] = 2.1892372;
        a[2] = -0.81520438;
        a[3] = 0.29036281;
        a[4] = -0.16596184;
        a[5] = 0.19191375;
    }
    else for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    //for (i=0;i<n_a;i++) a[i]=0.; // set to 0

    for (i=0; i<n_a; i++) LHe += a[i]*pow(logT3,i);
//  4) H+ // wli criteria cancelled 
    if (10.<T_K && T_K<=10000.){
        a[0] = -21.716699;
        a[1] = 1.3865783;
        a[2] = -0.37915285;
        a[3] = 0.11453688;
        a[4] = -0.23214154;
        a[5] = 0.058538864;
    }
    else for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    //for (i=0;i<n_a;i++) a[i]=0.; // set to 0

    for (i=0; i<n_a; i++) LHp += a[i]*pow(logT3,i);

//  5) e
    if (T_K<=10) for (i=0;i<n_a;i++) a[i] = 0.;
    else if (10.<T_K && T_K<=200.){
        a[0] = -34.286155;
        a[1] = -48.537163;
        a[2] = -77.121176;
        a[3] = -51.352459;
        a[4] = -15.169160;
        a[5] = -0.98120322;
    }
    else if (200.<T_K && T_K<=10000.){ //in table 200<T<10000 K; criteria cancelled 
        a[0] = -22.190316;
        a[1] = 1.5728955;
        a[2] = -0.21335100;
        a[3] = 0.96149759;
        a[4] = -0.91023195;
        a[5] = 0.13749749;
    }
    else for (i=0;i<n_a;i++) a[i]=0.; // set to 0
    //for (i=0;i<n_a;i++) a[i]=0.; // set to 0

    for (i=0; i<n_a; i++) Le += a[i]*pow(logT3,i);
//L0 in erg/s (per H2 molecule)
    //LH = pow(10.,LH); LH2 = pow(10.,LH2); LHe = pow(10.,LHe); LHp = pow(10.,LHp); Le = pow(10.,Le);
    LH = exp(LH); LH2 = exp(LH2); LHe = exp(LHe); LHp = exp(LHp); Le = exp(Le);
    // LH2*y_H2 +
    L0 = (LH*y_H +  LHe*y_He + LHp*y_Hp + Le*y_e)*nH; //L0 in erg/s (per H2 molecule)  GA08-eq37
    printf("LH=%3.2e, LH2=%3.2e, LHe=%3.2e, LHp=%3.2e, Le=%3.2e\n",LH,LH2,LHe,LHp,Le);
// LTE condition, Hollenbach & McKee (1979) (adopted by Yoshida et al. 2006)
// Λ_rot
    double La = 9.5e-22*pow(T3,3.76)/exp(pow(0.13/T3,3.0));
    double Lb = 1.0+0.12*pow(T3,2.1);
    double Lc = 3.0e-24*exp(-(0.51/T3));
    double L_rot = (La/Lb)+Lc;
// Λ_vib
    La = 6.7e-19/exp((5.86/T3));
    Lb = 1.6e-18/exp((11.7/T3));
    double L_vib = La+Lb;
// Λ_LTE in erg/s (per H2 molecule)
    L_LTE = (L_rot+L_vib);
// Λ
    L = L_LTE/(1.+L_LTE/L0);
    printf("GA: L0=%3.2e, L_LTE=%3.2e, L=%3.2e\n",L0/nH, L_LTE/nH, L/nH);
    ftau = min(1., pow(nH/8.e9,-0.45));
    return L*ftau*(y_H2*nH);  // in erg/cm^3/s
}

// Yoshida et al.(2006) -> Λ in erg /cm^3/s
double LambdaY_H2(double nH, double T_K, double* y){
    double y_H2 = y[2];
    double T3 = T_K/1.e3;
//     Λ_rotを決める
    double La = 9.5e-22*pow(T3,3.76)/exp(pow(0.13/T3,3.0));
    double Lb = 1.0+0.12*pow(T3,2.1);
    double Lc = 3.0e-24*exp(-(0.51/T3));
    double L_rot = (La/Lb)+Lc;
//     Λ_vibを決める
    La = 6.7e-19/exp((5.86/T3));
    Lb = 1.6e-18/exp((11.7/T3));
    double L_vib = La+Lb;
//     Λ_LTEを決める
    double L_LTE = (L_rot+L_vib)/nH; // number density of H atoms
//     Λ_H2(n→0)を決める
    La = -103.0+97.59*log10(T_K)-48.05*pow(log10(T_K),2.0); 
    Lb = 10.80* pow(log10(T_K),3.0) - 0.9032*pow(log10(T_K),4.0); 
    double L0 = pow(10.0,(La+Lb)); 
//     n_cr/n_Hを決める
    La=(L_LTE/L0); 
    //printf("n_cr = %3.2e\n",nH*La);
//     Λ_H2(H2の冷却関数)in erg cm^3 / s
    double L_H2=L_LTE/(1.0+La); //cout<<"\nL_H2= "<<L_H2<<endl;
    printf("Y : L0=%3.2e, L_LTE=%3.2e, L=%3.2e\n",L0, L_LTE, L_H2);
    double ftau = (1<pow(nH/8.e9,-0.45))?1:pow(nH/8.e9,-0.45);
    return L_H2*ftau *(y_H2*nH)*nH; //in erg /cm^3 / s
}
//**************************************************************************************
int main(){    
    int N = 10000, i;
    double y_H2 = 1.e-2;
    double y_e = 1.e-4;
    double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2p = 1.0e-12, y_Hm = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double ys[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};

    ofstream myfile;
    myfile.open("data/LambdaH2_nH.txt", ios::out | ios::trunc );
    myfile<<"nH\tGA\tYoshida\t Yoshida_Hm \t reaction"<<endl;
    double nH0 = .1; double T_K0 = 40;
    double nH = nH0;
    double T_K = T_K0, lgT4;
    double n_cr_H, n_cr_H2, n_cr;
//Λ v.s. nH
    for (i=0; i<N/2; i++){
        nH *= 1.003;
        T_K = T_K0*pow(nH/nH0, 2./3.);
        T_K = 5000;
        lgT4 = log10(T_K/1.e4);
        n_cr_H  = pow( 10., 3.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) );
        n_cr_H2 = pow( 10., 4.845 - 1.3*lgT4 + 1.62*pow(lgT4,2) );
        n_cr = 1.0/(y_H/n_cr_H + 2.0*y_H2/n_cr_H2); // Eq(14) GA08
        
        myfile <<"\t"<<nH;
        myfile  <<"\t"<<LambdaGA_H2(nH,T_K,ys)
                <<"\t"<<LambdaY_H2(nH,T_K,ys)
                <<"\t"<<1.e6/sqrt(T_K)*(1.6*exp(-pow(400./T_K,2))+1.4*y_H2*exp(-12000./(T_K+1200.)))
                <<"\t"<<n_cr
                << endl;
    }
    myfile.close();    
//Λ v.s. T
    nH = 1.e7; T_K = 100;
    myfile.open("data/LambdaH2_T.txt", ios::out | ios::trunc );
    myfile<<"# nH="<<nH;
    myfile<<"\nT_K\tGA\tYoshida\t Yoshida_Hm \t reaction"<<endl;
    for (i=0; i<N; i++){
        lgT4 = log10(T_K/1.e4);
        n_cr_H  = pow( 10., 3.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) );
        n_cr_H2 = pow( 10., 4.845 - 1.3*lgT4 + 1.62*pow(lgT4,2) );
        n_cr = 1.0/(y_H/n_cr_H + 2.0*y_H2/n_cr_H2); // Eq(14) GA08

        T_K += 10.; 
        myfile<<"\t"<<T_K;
        myfile  <<"\t"<<LambdaGA_H2(nH,T_K,ys)
                <<"\t"<<LambdaY_H2(nH,T_K,ys)
                <<"\t"<<1.e6/sqrt(T_K)*(1.6*exp(-pow(400./T_K,2))+1.4*y_H2*exp(-12000./(T_K+1200.)))
                <<"\t"<<n_cr
                << endl;
    }
    myfile.close();

//Λ v.s. species
    nH = 0.045; T_K = 100;
    ys[2] =1.e-3; ys[8] = 8.e-2;
    y_H2 = 1.e-3; y_Hep = 8.e-2;
    myfile.open("data/Lambda_T.txt", ios::out | ios::trunc );
    myfile<<"# nH="<<nH;
    // myfile<<"\nT_K\tLambda_{H2}\tLambda_H\tLambda_{He+}\tLambda_{H,He+}\tHe1s"<<endl;
    for (i=0; i<N; i++){
        T_K *=1.001; 
        myfile<<"\t"<<T_K;
        myfile  <<"\t"<<LambdaGA_H2(nH, T_K, ys)/pow(nH,2)
                <<"\t"<<LambdaY_H2(nH,T_K,ys)/pow(nH,2)
                <<"\t"<<Lambda_H(nH, T_K, y_H, y_e)/pow(nH,2)
                <<"\t"<<Lambda_Hep(nH, T_K, y_Hep, y_e)/pow(nH,2)
                <<"\t"<<(Lambda_H(nH, T_K, y_H, y_e)+Lambda_Hep(nH, T_K, y_Hep, y_e) )/pow(nH,2)
                << endl;
    }
    myfile.close();
    return 0;
}