// calculate the dynamical processes & values therein

# include <iostream>
# include <stdio.h>
# include "dyn.h"
# include "PARA.h"
# include "class_halo.h"
# include <cmath>
using namespace std;

// Cp = sqrt( 32*G / (3*pi) );
double n_ff(double z, double nH0, double rho_DM, double Dt){
// explicit method
    double dlnNdt = Cp * sqrt(rho_DM + (mu*m_H)*nH0);
    return nH0 * exp(dlnNdt*Dt);
/* 
//implicit method one step, parameter: Dt
    return nH0/( 1 - Dt * C*sqrt(nH0)); // t_ff = 1/C/sqrt(nH0)
 */
/* 
//analytical formula, parameter: t_ff0
    double t_ff0 = 1./sqrt(nH0)/C;
    return 1./pow( C*(t_ff0-t_act/2), 2 );
 */
}

double N_CORE(double z)
{
    return 7*pow((1+z)/11.,3);
}

double N_ADB(double S, double T){
    // adiabatic, entropy K = k_B*T*n^(-2/3)
    return pow(k_B*T / S, 1.5); // n \propto T^1.5
}

double Omega_mz(double z){
    return Omega_m0*pow(1+z,3.) /(Omega_m0*pow(1+z,3.) + Omega_L0);
}

double RHO_DM(double z){
    return 3*pow(H0,2)/(8*pi*G)*pow(1+z,3)*Omega_m0;
}

double Hz(double z){
    //cout<<" "<<H0<<" "<<z<<" "<<Omega_L0<<endl;
    //cout<<H0*sqrt( Omega_m0*pow(1+z,3) + Omega_L0 )<<endl;
    return H0*sqrt( Omega_m0*pow(1+z,3) + Omega_L0 );
}

double RHO_crit(double z){
    //z_dependence (1+z)^3
    return 3*pow(H0,2)/(8*pi*G)*pow(1+z,3)*Omega_m0/Omega_mz(z);
}

double t_freefall(double nH){
    return 1./C/sqrt(nH);
}

double z_ana(double z0, double t_from_z0){
    return pow(3*H0*sqrt(Omega_m0)/2*t_from_z0 + pow(1+z0,-1.5), -2./3) - 1;
}

double t_from_z(double z){ // age of universe at redshift z
    return 2./(3*H0*sqrt(Omega_m0)) * pow(1+z, -1.5);
}

