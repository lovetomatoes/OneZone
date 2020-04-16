// must declare & define here, or defalt initialization to 0;
//must const to avoid duplicate in other .o files
#include<cmath>

const double G=6.67408e-8, c=2.9979245e10, k_B=1.38064852e-16, m_H=1.66053904e-24;
const double gamma_adb=5./3, pi=3.141593;
const double mu = 1.2;
const double C = sqrt( 32*G*m_H/ (3*pi) );// t_ff = C^-1 * n^-.5
const double Cp = sqrt( 32*G / (3*pi) );// t_ff = C^-1 * n^-.5
const double erg2eV = 6.2415e11;
const double km = 1.e5;
const double pc = 3.e18, kpc = 3.e21, Mpc = 3.e24;
const double Myr = 1.e6*(365*24*3600);
const double yr = 365*24*3600;
const double Ms = 2.e33;
const double h0 = .677;
const double H0 = h0*100*km/Mpc;

// Planck constant
const double h_p = 6.63e-27, eV = 1.60217657e-12, Angstron = 1.e-8;
const double Mb = 1.e-18;

const double Omega_m0 = 0.311; // match calculation in code_tree/analysis_trees.f
const double Omega_L0 = 1 - Omega_m0;

const double fb = 0.16;
const double rho_m0 = 3*pow(H0,2)/(8*pi*G);// not used;
const double nb200_0 = 18*pow(pi,2)*rho_m0/m_H*fb;

//const int N_sp = 5, N_react = 8;
const int N_sp = 9, N_react = 42;
const int N_sp1 = N_sp+1, N_react1 = N_react+1;

const double ep10 = 0.1, ep5 = .2, ep2 = .5;
const double epE3 = .001, epE2 = .01, epE8 = 1.e-8, epE6 = 1.e-6, epE10 = 1.e-10;
const double epE1 = 0.1, epE4 = 1.e-4, epE5 = 1.e-5, epE7 = 1.e-7, epE9 = 1.e-9;

const double z_rcb = 1089, sigma1 = 30*km;