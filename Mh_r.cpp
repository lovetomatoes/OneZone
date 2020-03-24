#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include "PARA.h"
using namespace std;

// g++ Mh_r.cpp PARA.cpp -o Mr && ./Mr

double f_NFW(double c, double x){
    return -c*x/(1+c*x) + log(1+c*x);
}

double R200(double Mh, double z){
    double M_r,rho_m,c,rs,rvir,delta0;
    rho_m = rho_m0*pow(1+z,3);
    rvir = pow( Mh/(4./3*pi*200*rho_m),1./3. );
    return rvir;
}

double M_enc(double r, double Mh, double z){
    double M_r,rho_m,c,rs,r200,delta0;
    c = 10;
    rho_m = rho_m0*pow(1+z,3);
    r200 = R200(Mh, z);
    rs = r200/c;
    delta0 = 200./3.*pow(c,3)/f_NFW(c,1);
    M_r = 4*pi*rho_m*delta0*pow(rs,3)*f_NFW(c,r/r200);
    return M_r;
}

double Vc(double Mr,double r){
    return sqrt(G*Mr/r);
}

double Phi(double r, double Mh, double z){
    double rho_m,c,rs,r200,delta0;
    c = 10;
    rho_m = rho_m0*pow(1+z,3);
    r200 = R200(Mh, z);
    rs = r200/c;
    delta0 = 200./3.*pow(c,3)/f_NFW(c,1);

    return 4.*pi*G*rho_m*delta0*pow(rs,3) /r *log(1+r/rs);
}

int main(){
    int a[10] = {1,2,3,4,5,6,7,8,9,10};
    int* p = a;
    for (int i=0; i<9; i++) printf("%d\t",*(++p));

    for (int i=0;i<10;i++) {cout<<"haha\t"; cout<<"right?\n";}
    double z = 10, Mh = 1.e8*Ms, rvir, Vvir;
    rvir = R200(Mh,z);
    Vvir = Vc(Mh,rvir);
    int i,j,k, N = 1000;
    double r[N], M[N], V[N], P[N];
    ofstream file;
    file.open("enclosed_Mass.txt", ios::out | ios::trunc );
    for (i=0; i<N; i++){
        r[i] = (double)(i+1)/N*rvir;
        M[i] = M_enc(r[i],Mh,z);
        V[i] = Vc(M[i],r[i]);
        P[i] = Phi(r[i],Mh,z);
        if (i==0) file<<"\t x"<<"\t Mx "<<"\t Vcx \t GM/r \t Vrx\n";
        file<<" "<<r[i]/rvir<<" "<<M[i]/Mh<<" "<<V[i]/Vvir<<" "<<G*M[i]/r[i]<<" "<<sqrt(2*P[i])/Vvir<<endl;
        if (i==N-1) file<<"\t x"<<"\t Mx "<<"\t Vcx \t GM/r \t Vrx\n";
    }
    printf("M_enc = %10.9e Mh\n", M_enc(rvir,Mh,z) /Mh);
    printf("rvir = %3.2e pc\n", rvir);
    printf("Mh = %3.2e Ms\n", Mh/Ms);
    printf("Vc = %3.2e km/s \n", Vvir/km);

    return 0;
}