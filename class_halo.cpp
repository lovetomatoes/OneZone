// parent halo properties;
#include <iostream>
#include <iomanip>
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
    //  fit of Bullock et al. (2001)
    c = 18*pow(Mh0/(1.e11*Ms), -0.13)/(1+z);
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
    Kvir = k_B*Tvir/(mu*m_H)/pow(fb*Delta_crit*RHO_DM(z),2./3.); // wli: K normalization not precise, see Voit 05 fig1, Visbal 14
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
    return rho_c/(r/Rs) / pow(1+r/Rs,2);
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

double Mh_Tz(double Tvir, double z){
    double alpha_T = 2.324e4; //Tvir = 2.324e4 M8^(2./3.) * (1+z10)
    return 1.e8*Ms*pow(Tvir/alpha_T*11/(1+z),1.5);
}

double Mh_Vc(double Vc, double z){
    double alpha_Vc = 1.7945e6; // Vc = 1.7945e6 M8^(1/3) z11^(1/2) cm/s
    return 1.e8*Ms*pow( Vc/alpha_Vc/sqrt((1+z)/11.), 3);
}

double Tv_Vc(double Vc){
    return 0.5*(mu*m_H)*pow(Vc,2)/k_B;
}


/* to compile this file: 
g++ class_halo.cpp PARA.cpp dyn.cpp -o halo.out && ./halo.out
*/

// halo Vc, Rvir, n_crit at z...
// int main(){
//     double z = 20;
//     double Mh = Mh_Vc(4*km, z); Mh = Mh_Tz(1000,z);
// // checking Omukai 2001 rho_dm evolution
//     // printf("t from bigbang:%3.2e s\n",t_from_z(0));
//     // printf("Omega_0=%3.2e\n", pow(t_from_z(0)/3.1e17*h0, -2.) );
//     printf("Mh=%3.2e\n",Mh/Ms);
//     HALO halo(Mh,z);
//     printf("Vc=%3.2e, ",halo.Vc/km);
//     printf("n_crit:%6.4e, Rvir=%6.4e kpc, Tvir=%6.4e K\n",halo.rho_crit/(mu*m_H), halo.Rvir/kpc, halo.Tvir);
//     printf("K1=%3.2e,K2=%3.2e,Kvir=%3.2e\n",K_Tn(100,.1), K_Tn(1000,.01), halo.Kvir);
    
//     Mh = 1.e4*Ms;
//     z = 0; Mh = 1.e15/h0*Ms;
//     HALO halo1(Mh,z);
//     printf("K200:%5.3e\n",halo1.Kvir);

//     z = 10;
//     Mh = 1.e5*Ms;
//     double Mh1 = 1.e8*Ms, Mhrat = exp(log(Mh1/Mh)/20.);
//     ofstream f1;
//     f1.open("a.txt", ios::out | ios::trunc );
//     f1<<setiosflags(ios::scientific)<<setprecision(5);
//     f1<<setw(12)<<"Mh_Ms"<<setw(12)<<"Kcore"<<setw(12)<<"K_ISM"<<endl;
//     while (Mh<=Mh1){
//         HALO halo2(Mh,z);
//         f1<<setw(12)<<Mh/Ms<<setw(12)<<halo2.Kvir/10.<<setw(12)<<K_ISM(z)<<endl;
//         Mh *= Mhrat;
//     }
//     f1.close();
// }
