#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "class_halo.h"
#include "PARA.h"
#include "RK4.h"
#include "LE_iso.h"
#include "dyn.h"
//#include "class_gas.h"
using namespace std;

/* 
rm Mg*.txt
g++ -c LE_iso.cpp
g++ -c RK4.cpp
g++ LE_iso.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_iso.out
./le_iso.out
*/

//input para
//output: a profile of density distributions, Phi(x), x=r/a, Phi = log(rho/rho_g0)

double Tvir = 1.6e4;
//double z = 20; double Mh = 1.e7*Ms;
double z1 = 20; double Mh1 = 1.e5*Ms;
//double z = 30; double Mh = Mh_Tz(Tvir,z);
//double z = 25; double Mh = Mh_Tz(Tvir,z);
//double z = 10; double Mh = Mh_Tz(Tvir,z);
//double z = 0; double Mh = Mh_Tz(Tvir,z);
double R_EQ(double Tg, double rhoc, double rs){
    return 9.*pi*G*mu*m_H*rhoc*(rs*rs)/(k_B*Tg);
}
double A_TR(double Tg, double rhoc, double R){
    return sqrt(k_B*Tg/ (4*pi*G*mu*m_H*rhoc*R) );
}
extern "C" void profile(char* filename, double Tg, double R, double z=z1, double Mh=Mh1){
    HALO halo1(Mh,z);
    int const N = 1000000;
    int const n = 2;
    int i;
    double rho_g0 = halo1.rho_c * R;
    double M_intg = 0;
    double a = sqrt(k_B*Tg/(4*pi*G*mu*m_H*rho_g0));
    double alpha = a/halo1.Rs;
    // set x1, largest xi 
    double x1 = halo1.Rvir /a;
    //x1 = 100; self-gravity
    printf("scaling factor a = %3.2e cm, alpha = %3.2e, delta_rhoc=%3.2e, Rs=%3.2e cm Rvir=%3.2e\n", a,alpha,halo1.rho_c,halo1.Rs,halo1.Rvir);
    double* x = new double [N];
    
    double** y = new double* [N];
    for (i=0;i<N;i++){
        y[i] = new double [n];
    }
    int c = 2;
    double* v = new double [c];
    v[0] = alpha; v[1] = R;

    //B for checking DM gravity w/ Suto, Sasaki, Makino 1998
    double B = 4*pi*G*(mu*m_H)*halo1.delta0*halo1.rho_crit*pow(halo1.Rs,2)/(k_B*Tg);
    printf("B=%3.2e",B);
    M_intg = 0;
    double* dydx0 = new double [n];
    fstream file;
    file.open(filename, ios::out | ios::trunc );

    file<<"r r_Rvir r_Rs r_r1 phi psi rhog_over_rhog0 n_gas rhoDM_over_rhog0 n_DM M_intg M_fit Bfx\n";
    // dx
    double dx;
    double err=.01;
    // boundary conditions
    i = 0;
    dx = x1/N;
    // boundary x[0]=0, but 0 causes sigularity, using x[0]=dx/1000<<dx
    //x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    x[i] = dx/1.e8; y[i][0] = 0; y[i][1] = 0; // isothermal case
    printf("center boundary: r0=%3.2e (Rvir)\t dr=%3.2e (Rvir)\n",a*x[0]/halo1.Rvir,a*dx/halo1.Rvir);

    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        check_conv(err,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        x[i]=x[i-1]+dx;
        //integrate within R_vir
        if (x[i]<=halo1.Rvir/a) M_intg += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);
        double r1 = A_TR(Tg,halo1.rho_c,R_EQ(Tg,halo1.rho_c,halo1.Rs));
    ////file<<"r r_Rvir r_Rs r_r1 phi psi rhog_over_rhog0 n_gas rhoDM_over_rhog0 n_DM M_intg Bfx\n";
        file<<a*x[i]<<" "<<a*x[i]/halo1.Rvir<<" "<<alpha*x[i]<<" "<<a*x[i]/r1;
        file<<" "<<y[i][0]<<" "<<y[i][1];
        file<<" "<<exp(-y[i][0]); //rhog_over_rhog0
        file<<" "<<rho_g0*exp(-y[i][0])/(mu*m_H); // n_gas
        file<<" "<< 1./(R* alpha*x[i]* pow((1+alpha*x[i]),2) );// rhoDM_over_rhog0
        file<<" "<<halo1.rho_crit*halo1.delta0/( alpha*x[i] * pow((1+alpha*x[i]),2) )/(mu*m_H);// n_DM
        file<<" "<<M_intg;
        file<<" "<<4.*pi/3.*rho_g0*pow(a*x[i],3);
        file<<" "<<B*(1- log(1+alpha*x[i])/ (alpha*x[i]))<<endl;
    }
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    file.close();
    printf("z = %3.2f\tMh = %3.2e Ms\t rs%3.2e cm\trhoc=%3.2e/cc\n",halo1.z,halo1.Mh/Ms,halo1.Rs,halo1.rho_c);
    printf("rho_g0=%3.2e, rho_crit=%3.2e, delta_0=%3.2e\n",rho_g0, halo1.rho_crit, halo1.delta0);
}

double Mg(char* filename, double Tg, double R, double z=z1, double Mh=Mh1){
    printf("R in solving Mg is: %3.2e\n",R);
    HALO halo1(Mh,z);
    int const N = 1000000;
    int const n = 2;
    int i;
    double rho_g0 = halo1.rho_c * R;
    double M_intg = 0;
    double a = sqrt(k_B*Tg/(4*pi*G*mu*m_H*rho_g0));
    double alpha = a/halo1.Rs;
    // set x1, largest xi 
    double x1 = halo1.Rvir *10 /a;

    double* x = new double [N];
    
    double** y = new double* [N];
    for (i=0;i<N;i++){
        y[i] = new double [n];
    }
    int c = 2;
    double* v = new double [c];
    v[0] = alpha; v[1] = R; 

    //B for checking DM gravity w/ Suto, Sasaki, Makino 1998
    double B = 4*pi*G*(mu*m_H)*halo1.delta0*halo1.rho_crit*pow(halo1.Rs,2)/(k_B*Tg);
    double* dydx0 = new double [n];
    fstream file;
  
    //________________________________________________________
    // commented out by wli 20200219; now use file to print M_intg v.s. xi
    //________________________________________________________
    if (FILE *f = fopen(filename, "r")) fclose(f);
    else {
        file.open(filename, ios::out | ios::trunc);
        file<<"R ng_0 M_intg ng_vir\n";
        file.close();
    }
    file.open(filename, ios::out | ios::app);

    // dx
    double dx;
    double err=.01;
    // boundary conditions
    i = 0;
    dx = x1/N;
    x[i] = dx/1.e8; y[i][0] = 0; y[i][1] = 0; // isothermal case
    int count_in = 0;
    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        check_conv(err,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        x[i]=x[i-1]+dx;

        //integrate within R_vir
        if (x[i]<=halo1.Rvir/a) M_intg += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);
        else break;
    } 
    ////file<<"R alpha M_intg n_vir\n";    
    file<<" "<<R<<" "<<rho_g0/(mu*m_H)<<" "<<M_intg/Ms<<" "<<rho_g0*exp(-y[i][0])/(mu*m_H)<<endl;
    file.close();
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
   
    //printf("z = %3.2f\tMh = %3.2e Ms\t rs=%3.2e cm\trhoc=%3.2e/cc\trvir=%3.2e cm\n",halo1.z,halo1.Mh/Ms,halo1.Rs,halo1.rho_c,halo1.Rvir);
    printf("z = %3.2f\tMh = %3.2e Ms\t Mg_vir = %3.2e Ms\n",halo1.z,halo1.Mh/Ms,M_intg/Ms);
    printf("in Mg func: R = %3.2e\t ng_0 = %3.2e\n", R, halo1.rho_c*R/(mu*m_H));   

    return M_intg;
}

double N_VIR(double Tg, double R, double z=z1, double Mh=Mh1){
    //printf("in N_VIR: R=%3.2e, Tg=%3.2e, z=%3.2e, Mh=%3.2e\n",R,Tg,z,Mh/Ms);
    HALO halo1(Mh,z);
    int const N = 100000;
    int const n = 2;
    int i;
    double rho_g0 = halo1.rho_c * R;
    double M_intg = 0;
    double a = sqrt(k_B*Tg/(4*pi*G*mu*m_H*rho_g0));
    double alpha = a/halo1.Rs;
    // set x1, largest xi from 1.01 *Rvir
    double x1 = halo1.Rvir *1.01 /a;

    double* x = new double [N];
    double** y = new double* [N];
    for (i=0;i<N;i++){
        y[i] = new double [n];
    }
    int c = 2;
    double* v = new double [c];
    v[0] = alpha; v[1] = R; 

    //B for checking DM gravity w/ Suto, Sasaki, Makino 1998
    double B = 4*pi*G*(mu*m_H)*halo1.delta0*halo1.rho_crit*pow(halo1.Rs,2)/(k_B*Tg);
    double* dydx0 = new double [n];

    // dx & boundary conditions
    double dx = x1/N;
    i = 0;
    x[i] = dx/1.e8; y[i][0] = 0; y[i][1] = 0; // isothermal case
    int count_in = 0;
    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        check_conv(epE2,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx);
        x[i]=x[i-1]+dx;

        //integrate within R_vir
        if (x[i]<=halo1.Rvir/a) M_intg += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);
        else break;
    }
    double n_vir = rho_g0*exp(-y[i][0])/(mu*m_H);
    printf("in N_VIR: n_vir=%3.2e\n",n_vir);
    delete [] x; delete [] y; delete [] v; delete [] dydx0;

    // printf("z = %3.2f\tMh = %3.2e Ms\t Mg_vir = %3.2e Ms\n",halo1.z,halo1.Mh/Ms,M_intg/Ms);
    // printf("in Mg func: R = %3.2e\t ng_0 = %3.2e\n", R, halo1.rho_c*R/(mu*m_H));   

    return n_vir;
}

void Nvir2N0(double& n_sol, double& nvir_max, double ni, double Tg, double z, double Mh){
    double R0 = .1, delta_R = epE2*R0, R1;
    double R_max, n0_max;
    double nvir_fw, nvir_bw ;
    double nvir_0, nvir_1;
    double dndR, d2ndR;
    double dR = R0;
    int it = 0;
    HALO halo(Mh,z);

// local maximum of n_vir
    while ( pow(dR,2)>=epE2*pow(R0,2) && it<6){ //dR < 0.1 R0 or it=6, 结束计算
        delta_R = epE2*R0;
        nvir_fw = N_VIR(Tg, R0+.5*delta_R, z, Mh); nvir_bw = N_VIR(Tg, R0-.5*delta_R, z, Mh); 
        nvir_0 = N_VIR(Tg,R0,z,Mh);
        dndR = ( nvir_fw - nvir_bw )/ delta_R;
        d2ndR = 4.*( nvir_fw + nvir_bw - 2.*nvir_0)/ pow(delta_R, 2);

        dR = - dndR/d2ndR; R1 = R0 + dR;
    // update R0
        R0 = R1;
        it++;
        //printf("\nin LOOP for MAXIMUM: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
    }
    // R_max = R0; n0_max = R0*halo.rho_c/(mu*m_H);
    nvir_max = N_VIR(Tg,R0,z,Mh);
    printf("\nR_max=%3.2e, nvir_max=%3.2e\n",R0,nvir_max);
    printf("z=%3.2e, Mh=%3.2e\n",z,Mh/Ms);
// unstable criterion: nvir_max < n_mean (cosmic mean density at z) 
    double n_mean = RHO_crit(z)/(mu*m_H);
    printf("nmean=%3.2e, ni=%3.2e:  UNSTABLE\n",n_mean,ni);
    if (n_mean>nvir_max) n_sol = 0; //unstable
    // if (true) {} // wli: 不解n0 只要nvir_max
    else{ // get solution of given outer boundary n_mean
        it = 0;
        R0 = ni*(mu*m_H)/halo.rho_c;
        dR = R0;
        while ( pow(dR,2)>=epE2*pow(R0,2) && it<5){ //dR < 0.1 R0 or it=6, 结束计算
            delta_R = epE2*R0;
            nvir_fw = N_VIR(Tg, R0+.5*delta_R, z, Mh); nvir_bw = N_VIR(Tg, R0-.5*delta_R, z, Mh); 
            nvir_0 = N_VIR(Tg,R0,z,Mh);
            dndR = ( nvir_fw - nvir_bw )/ delta_R;
            dR = - (nvir_0-n_mean)/dndR; R1 = R0 + dR;
        // update R0
            if (R1<0) R0 = R0/10.;
            else R0 = R1;
            it++;
            // printf("\nin LOOP: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
        }
        n_sol = R0*halo.rho_c/(mu*m_H);
        printf("solution: of n0: %3.2e\n",n_sol);
    }
}

/* 
rm Mg*.txt
g++ -c LE_iso.cpp
g++ -c RK4.cpp
g++ LE_iso.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_iso.out
./le_iso.out
*/

/* int main(){
    double R = 1.64;
    char* fname = "z_nmax.txt";
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<" z Tv_crit\n";
    double Tg = 1.5e4, Tvir=7.e3, Tv1,Tv0, Tv_sol;
    double z0 = 35, z1 = 10, Mh;
    double n_sol = 0, ni = 1, n_mean, n_mid;
    double N = 5;
    double z_rat = exp( log(z1/z0)/N ); cout<<z_rat<<endl;
    double nvir_max, nvir_m1, nvir_m0;
    for (int i=0; i<N; i++){
        Tv1 = 2.e4,Tv0 =5.e3;
        Tv1 = 3.e4,Tv0 =1.5e4;
        Tv1 = 3.e4,Tv0 =2.2e4;
        Nvir2N0(n_sol, nvir_m0, ni, Tg, z0, Mh_Tz(Tv0, z0));
        Nvir2N0(n_sol, nvir_m1, ni, Tg, z0, Mh_Tz(Tv1, z0));
        n_mean = RHO_crit(z0)/(mu*m_H);
        if (nvir_m0<n_mean or nvir_m1>n_mean) printf("\n!!!!!initial value not correct\n");
        while (pow(Tv1/Tv0-1, 2)>epE2){
            Nvir2N0(n_sol, n_mid, ni, Tg, z0, Mh_Tz((Tv0+Tv1)/2., z0));
            if (n_mid>=n_mean) {
                Tv0 = (Tv0+Tv1)/2.;
                Nvir2N0(n_sol, nvir_m0, ni, Tg, z0, Mh_Tz(Tv0, z0));
            }
            else {
                Tv1 = (Tv0+Tv1)/2.;
                Nvir2N0(n_sol, nvir_m1, ni, Tg, z0, Mh_Tz(Tv1, z0));
            }
        }
        Tv_sol = (Tv0+Tv1)/2.;
        double cs2 = gamma_adb*k_B*8000/(mu*m_H);
        f<<" "<<z0<<" "<<Tv_sol<<endl;
        z0 *= z_rat;
    }
    f.close();

    return 0;
} */