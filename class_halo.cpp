// parent halo properties;
#include <iostream>
#include <stdio.h>
#include "PARA.h"
#include "dyn.h"
#include "class_halo.h"
#include <cmath>
#include <fstream>
using namespace std;

/* to compile this file: 
g++ class_halo.cpp PARA.cpp dyn.cpp -o halo.out && ./halo.out
*/

HALO:: HALO(double Mh0, double z0){
    Mh = Mh0;
    z = z0;
    // concentration parameter c from Dekel & Birnboim 2006 Eq(22)
    c = 18*pow(Mh0/(1.e11*Ms), -0.13)/(1+z);
    //c=4; //wli: trying fixed concentration factor!!!!
    double d = Omega_mz(z) - 1;
    Delta_crit = 18.0*pi*pi + 82*d - 39*d*d;

    delta0 = Delta_crit/3.*pow(c,3)/F_NFW(1);

    rho_crit = RHO_crit(z); // mean density of DM at z
    rho_c = rho_crit * delta0;

    Rvir = pow( Mh/(4./3*pi*Delta_crit*rho_crit),1./3. );
    Rs = Rvir/c;
    Vc = sqrt(G*Mh/Rvir);
    t_dyn = Rvir/Vc;
    Tvir = G*Mh*(mu*m_H)/(2.*k_B*Rvir);
    gc = 2*c/(log(1+c) - c/(1+c));
    alpha = Tvir/pow(Mh , 2./3);

    i = 0;
    id = 0;
}

HALO:: ~HALO(void){
}

// x = r/Rvir; c = Rvir/Rs
double HALO:: F_NFW(double x){
    return -c*x/(1+c*x) + log(1+c*x);
}

double HALO:: Rho_r(double r){ //
    return rho_crit*delta0/(r/Rs) / pow(1+r/Rs,2);
}

double HALO:: M_enc(double r){
    double M_r = 4*pi*rho_crit*delta0*pow(Rs,3)*F_NFW(r/Rvir);
    return M_r;
}

double HALO:: Phi(double r){
    // lim r -> 0
    //return -4*pi*G*rho_crit*delta0*Rs*Rs;
    return -4*pi*G*rho_crit*delta0*pow(Rs,3)/r*log(1+r/Rs);
}

double HALO:: V_inf2(double rout, double rin){
    return 2*(Phi(rout)-Phi(rin));
}

/* //origin, checking...(forget goal)
int main(){
    int N = 1000;
    int i,j;
    double Mh0, z0;
    double x, r, v;
    Mh0 = 1.e8*Ms;
    z0 = 10;
    HALO halo(Mh0, z0);
    ofstream file;
    file.open("data/Vs.txt", ios::out | ios::trunc );
    for (i=0;i<N+1;i++){
        if (i==0) file<<" x"<<" Vr/Vc"<<endl;
        x = (double)i/N;
        r = x*halo.Rvir;
        v = sqrt(-2*halo.Phi(r));
        file <<" "<<x<<" "<<v/halo.Vc<<endl;
    }
    file.close();

    ofstream file1;
    file1.open("data/Vs1.txt", ios::out | ios::trunc );
    double z1,z2,zi;
    z1 = 0; z2 = 100;
    for (i=0;i<N+1;i++){
        if (i==0) file1 <<" z"<<" Vr0/Vc"<<" g(c)"<<endl;
        zi = (1-(double)i)/N*z1 + (double)i/N*z2;
        HALO halo(Mh0, zi);
        v = sqrt(8*pi*G*halo.rho_crit*halo.delta0*pow(halo.Rs,2));
        file1 <<" "<<zi<<" "<<v/halo.Vc<<" "<<sqrt(2*halo.c/(log(1+halo.c) - halo.c/(1+halo.c) ))<<endl;
    }
    file1.close();
    HALO halo1(1.e8*Ms,9);
    printf("z=9, M=1.e8Ms, Rvir = %5.2f kpc, Vc = %5.2f km/s, Tvir = %5.2f K\n", halo1.Rvir/kpc,halo1.Vc/km,halo1.Tvir);
    printf("z=9, age of Universe = %3.2f Myr\n", t_from_z(9)/Myr);    
    return 0;
} */

double Mh_Tz(double Tvir, double z){
    double alpha_T = 2.324e4; //Tvir = 2.324e4 M8^(2./3.) (1+z10)^-1
    return 1.e8*Ms*pow(Tvir/alpha_T*11/(1+z),1.5);
}

double Mh_Vc(double Vc, double z){
    double alpha_Vc = 1.7945e6; // Vc = 1.7945e6 M8^(1/3) z11^(1/2) cm/s
    return 1.e8*Ms*pow( Vc/alpha_Vc/sqrt((1+z)/11.), 3);
}

/* 
// halo Vc, Rvir, n_crit at z...
int main(){
    double z = 17.2;
    double Mh = 2.36e7*Ms;
    HALO halo(Mh,z);
    printf("n_crit:%6.4e, Rvir=%6.4e kpc, Vc=%6.4e km/s\n",halo.rho_crit/(mu*m_H), halo.Rvir/kpc, halo.Vc/km);
} */