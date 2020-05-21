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
void profile(char* filename, double Tg, double R, double z=z1, double Mh=Mh1){
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
        DyDx_iso(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
        check_conv(err,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
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
        file<<" "<<M_intg/Ms; //M(r) in solar mass
        file<<" "<<4.*pi/3.*rho_g0*pow(a*x[i],3);
        file<<" "<<B*(1- log(1+alpha*x[i])/ (alpha*x[i]))<<endl;
    }
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    file.close();
    printf("z = %3.2f\tMh = %3.2e Ms\t rs%3.2e cm\trhoc=%3.2e/cc\n",halo1.z,halo1.Mh/Ms,halo1.Rs,halo1.rho_c);
    printf("rho_g0=%3.2e, rho_crit=%3.2e, delta_0=%3.2e\n",rho_g0, halo1.rho_crit, halo1.delta0);
}

void BOUNDARY(double& N_VIR, double& MG_VIR, double Tg, double R, double z=z1, double Mh=Mh1){
    MG_VIR = 0;
    HALO halo1(Mh,z);
    int const N = 100000;
    int const n = 2;
    int i;
    double rho_g0 = halo1.rho_c * R;
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
    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx_iso(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
        check_conv(epE2,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
        x[i]=x[i-1]+dx;

        //integrate within R_vir
        if (x[i]<=halo1.Rvir/a) MG_VIR += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);
        else break;
    }
    N_VIR = rho_g0*exp(-y[i][0])/(mu*m_H);
    delete [] x; delete [] y; delete [] v; delete [] dydx0;

    // printf("in BOUNDARY: R=%3.2e, Tg=%3.2e, z=%3.2e, Mh=%3.2e, Tvir=%3.2e\n",R,Tg,z,Mh/Ms,halo1.Tvir);
    // printf("in BOUNDARY: N0=%3.2e, N_VIR=%3.2e, MG_VIR=%3.2E\n", halo1.rho_c*R, N_VIR, MG_VIR/Ms);
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
    double Mg_vir=0;
// local maximum of n_vir
    while ( pow(dR,2)>=epE2*pow(R0,2) && it<6){ //dR < 0.1 R0 or it=6, 结束计算
        delta_R = epE2*R0;
        // nvir_fw & nvir_bw
        BOUNDARY(nvir_fw, Mg_vir, Tg, R0+.5*delta_R, z, Mh); BOUNDARY(nvir_bw, Mg_vir, Tg, R0-.5*delta_R, z, Mh);
        BOUNDARY(nvir_0, Mg_vir, Tg, R0, z, Mh);
        dndR = ( nvir_fw - nvir_bw )/ delta_R;
        d2ndR = 4.*( nvir_fw + nvir_bw - 2.*nvir_0)/ pow(delta_R, 2);

        dR = - dndR/d2ndR; R1 = R0 + dR;
    // update R0
        R0 = R1;
        it++;
        //printf("\nin LOOP for MAXIMUM: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
    }
    // R_max = R0; n0_max = R0*halo.rho_c/(mu*m_H);
    BOUNDARY(nvir_max, Mg_vir, Tg, R0, z, Mh);
    // printf("\nR_max=%3.2e, nvir_max=%3.2e\n",R0,nvir_max);
    // printf("z=%3.2e, Mh=%3.2e\n",z,Mh/Ms);
// unstable criterion: nvir_max < n_mean (cosmic mean density at z) 
    double n_mean = RHO_crit(z)/(mu*m_H);
    double n_nfw = fb*halo.Rho_r(halo.Rvir)/(mu*m_H);
    printf("nmean=%3.2e, n_nfw=%3.2e, nvir_max=%3.2e:\n",n_mean,n_nfw, nvir_max);
    n_mean = n_nfw;
    if (n_mean>nvir_max) {
        printf("UNSTABLE!!!!!!!!!!!!!\n");
        n_sol = 0; //unstable
    }
    // if (true) {} // wli: 不解n0 只要nvir_max
    else{ // get solution of given outer boundary n_mean
        it = 0;
        R0 = ni*(mu*m_H)/halo.rho_c;
        dR = R0;
        while ( pow(dR,2)>=epE2*pow(R0,2) && it<5){ //dR < 0.1 R0 or it=6, 结束计算
            delta_R = epE2*R0;
            BOUNDARY(nvir_fw, Mg_vir, Tg, R0+.5*delta_R, z, Mh); BOUNDARY(nvir_bw, Mg_vir, Tg, R0-.5*delta_R, z, Mh);
            BOUNDARY(nvir_0, Mg_vir, Tg, R0, z, Mh);
            dndR = ( nvir_fw - nvir_bw )/ delta_R;
            dR = - (nvir_0-n_mean)/dndR; R1 = R0 + dR;
        // update R0
            if (R1<0) R0 = R0/10.;
            else R0 = R1;
            it++;
            // printf("\nin LOOP: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
        }
        n_sol = R0*halo.rho_c/(mu*m_H);
        printf("Nvir2N0: solution: of n0: %3.2e\n",n_sol);
    }
}

// 改自Nvir2N0; essencially the same. 暂未验证: Mg_max 对应的n0, 任一Mg对应的 n_sol
void Mg2N0(double& n_sol, double& Mg_max, double ni, double Tg, double z, double Mh){
    double R0 = .1, delta_R = epE2*R0, R1;
    double R_max, n0_max;
    double Mg_fw, Mg_bw ;
    double Mg_0;
    double dMdR, d2MdR;
    double dR = R0;
    int it = 0;
    HALO halo(Mh,z);
    double nvir=0;
  // local maximum of Mg
    while ( pow(dR,2)>=epE2*pow(R0,2) && it<6){ //dR < 0.1 R0 or it=6, 结束计算
        delta_R = epE2*R0;
        // nvir_fw & nvir_bw
        BOUNDARY(nvir, Mg_fw, Tg, R0+.5*delta_R, z, Mh); BOUNDARY(nvir, Mg_bw, Tg, R0-.5*delta_R, z, Mh);
        BOUNDARY(nvir, Mg_0, Tg, R0, z, Mh);
        dMdR = ( Mg_fw - Mg_bw )/ delta_R;
        d2MdR = 4.*( Mg_fw + Mg_bw - 2.*Mg_0)/ pow(delta_R, 2);

        dR = - dMdR/d2MdR; R1 = R0 + dR;
    // update R0
        R0 = R1;
        it++;
        //printf("\nin LOOP for MAXIMUM: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
    }
    // R_max = R0; n0_max = R0*halo.rho_c/(mu*m_H);
    BOUNDARY(nvir, Mg_max, Tg, R0, z, Mh);
    printf("\nR_max=%3.2e, nvir_max=%3.2e\n",R0,Mg_max);
    printf("z=%3.2e, Mh=%3.2e\n",z,Mh/Ms);
  // unstable criterion: baryonic mass v.s. Mg_max
    if (fb*Mh>Mg_max) n_sol = 0; //unstable
    // if (true) {} // wli: 不解n0 只要Mg_max
    else{ // get solution of given Mg=fb*Mh
        it = 0;
        R0 = ni*(mu*m_H)/halo.rho_c;
        dR = R0;
        while ( pow(dR,2)>=epE2*pow(R0,2) && it<5){ //dR < 0.1 R0 or it=6, 结束计算
            delta_R = epE2*R0;
            BOUNDARY(nvir, Mg_fw, Tg, R0+.5*delta_R, z, Mh); BOUNDARY(nvir, Mg_bw, Tg, R0-.5*delta_R, z, Mh);
            BOUNDARY(nvir, Mg_0, Tg, R0, z, Mh);
            dMdR = ( Mg_fw - Mg_bw )/ delta_R;
            dR = - (Mg_0-fb*Mh)/dMdR; R1 = R0 + dR;
        // update R0
            if (R1<0) R0 = R0/10.;
            else R0 = R1;
            it++;
            // printf("\nin LOOP: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
        }
        n_sol = R0*halo.rho_c/(mu*m_H);
        printf("Mg2N0: solution: of n0: %3.2e\n",n_sol);
    }
}



/* 
g++ -c LE_iso.cpp
g++ -c RK4.cpp
g++ LE_iso.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_iso.out
./le_iso.out
*/


/* int main(){// *****************  Tvir critical value v.s. z (Tg_eff changes w/ bsm&z as input parameter) **************
    char* fname = "H2_3s_z_Tvcrit.txt";
    // fname = "threshold_3s.txt";
    int n_sigma = 3;
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<"z Tv_crit Mh_crit\n";
    // f<<" z Tv_min Mh_min Tv0_min Mh0_min\n";
    double Tv1,Tv0, Tv_sol;
    double z0 = 35, z1 = 25, Mh;
    z0 = 25, z1 = 15;
    double v_bsm;
    // v_bsm = sigma1*(1+z0)/(1+z_rcb); // z in [35:25]: 1sigma v_bsm [0.99, 0.71]
    // cout<<2*sigma1*(1+z1)/(1+z_rcb)/km<<endl;//1.7
    double n_sol = 0, ni = 1, n_mean, n_mid;
    double alpha = 4.7;
    double Tg, cs, f_Ma;
    double nvir_max, nvir_m1, nvir_m0;
  // calculation
    double N = 5;
    double z_rat = exp( log(z1/z0)/N );
    for (int i=0; i<N; i++){
        v_bsm = n_sigma*sigma1*(1+z0)/(1+z_rcb);
        cs = sqrt(gamma_adb* k_B*2000./(mu*m_H) ); //sound speed of 2000K gas
        f_Ma = 1+pow(alpha*v_bsm/cs,2);
        Tg = 2000*f_Ma;
        // Tv0 = 0.5*(mu*m_H)*pow(5*km,2)/k_B;
        // Tv1 = 0.5*(mu*m_H)*pow(10*km,2)/k_B;
        Tv0 = 5000; Tv1 = 20000;
  // 二分法求Tvir_crit; print
        Nvir2N0(n_sol, nvir_m0, ni, Tg, z0, Mh_Tz(Tv0, z0));
        Nvir2N0(n_sol, nvir_m1, ni, Tg, z0, Mh_Tz(Tv1, z0));
        n_mean = RHO_crit(z0)/(mu*m_H);
        if (nvir_m0<n_mean or nvir_m1>n_mean) {
            printf("\n!!!!!initial value not correct\n");
            return 0;
        }
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
        f<<" "<<z0<<" "<<Tv_sol<<" "<<Mh_Tz(Tv_sol,z0)/Ms<<endl;
  // 算minumum cooling mass & Tvir; print
        // double Vc0 = 3.7*km;
        // double Vc_eff = sqrt( pow(Vc0,2)+pow(alpha*v_bsm,2) );
        // f<<" "<<z0<<" "<<Tv_Vc(Vc0)<<" "<<Mh_Vc(Vc0,z0)/Ms<<" "<<Tv_Vc(Vc_eff)<<" "<<Mh_Vc(Vc_eff,z0)/Ms<<endl;
  // update z0
        z0 *= z_rat;
    }
    f.close();

    cout<<"f_Ma = "<<f_Ma<<endl;
    return 0;
} */

/* int main(){ // *****************  Tvir critical value v.s. Tg (z as input parameter) **************
    char* fname = "Tg_Tvcrit.txt";
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<"Tg Tv_crit\n";
    double Tg = 7.e3, Tg1 = 2.e4, Tv1,Tv0, Tv_sol;
    double z = 25, Mh;
    double n_sol = 0, ni = 1, n_mean, n_mid;
    double N = 10;
    double Tgrat = exp( log(Tg1/Tg)/N ); cout<<Tgrat<<endl;
    double nvir_max, nvir_m1, nvir_m0;
    for (int i=0; i<N; i++){
        Tv1 = 3.e4,Tv0 = 5.e3;
        Nvir2N0(n_sol, nvir_m0, ni, Tg, z, Mh_Tz(Tv0, z));
        Nvir2N0(n_sol, nvir_m1, ni, Tg, z, Mh_Tz(Tv1, z));
        n_mean = RHO_crit(z)/(mu*m_H);
        n_mean *= 10;
        if (nvir_m0<n_mean or nvir_m1>n_mean) {
            printf("\n!!!!!initial value not correct\n");
            break;
        }
        while (pow(Tv1/Tv0-1, 2)>epE2){
            Nvir2N0(n_sol, n_mid, ni, Tg, z, Mh_Tz((Tv0+Tv1)/2., z));
            if (n_mid>=n_mean) {
                Tv0 = (Tv0+Tv1)/2.;
                Nvir2N0(n_sol, nvir_m0, ni, Tg, z, Mh_Tz(Tv0, z));
            }
            else {
                Tv1 = (Tv0+Tv1)/2.;
                Nvir2N0(n_sol, nvir_m1, ni, Tg, z, Mh_Tz(Tv1, z));
            }
        }
        Tv_sol = (Tv0+Tv1)/2.;
        double cs2 = gamma_adb*k_B*8000/(mu*m_H);
        f<<" "<<Tg<<" "<<Tv_sol<<endl;
        Tg *= Tgrat;
    }
    f.close();

    return 0;
}
 */

// 检查Nvir2N0, Mg2N0是否给出正确的nvir_max, Mg_max 
/* int main(){
    double R = 1.e-2;
    double Tg = 1.5e4;
    double R1 = 1.e2, Rrat = exp(log(R1/R)/20.);
    double n0_sol_M, Mg_max, n0_sol_N, nvir_max, n0_sol;
    double ni = 1;
    double z = 20, Mh = 2.e7*Ms;
    HALO halo(Mh,z);
    double n_mean = RHO_crit(z)/(mu*m_H);
    double n_vir, Mg_vir;

  // Tg v.s. n0(Mg) 没用
    // Tg = 8.e3;
    // double Tg1 = 2.e4, Tgrat = exp(log(Tg1/Tg)/20.);
    // z = 20, Mh = 2.e7*Ms;
    // char* fname = "MJ_Tg.txt";
    // fstream f;
    // f.open(fname, ios::out | ios:: trunc);
    // f<<"Tg n0 Mg\n";
    // Mh = 1.e5*Ms;
    // double Mh1 = 5.e7*Ms, Mhrat = exp(log(Mh1/Mh)/20.);
    // while (Tg<Tg1){
    //     // Nvir2N0(n0_sol,nvir_max,ni,Tg,z,Mh);
    //     Mg2N0(n0_sol, Mg_max, ni, Tg, z, Mh);
    //     BOUNDARY(n_vir, Mg_vir,Tg,n0_sol*(mu*m_H)/halo.rho_c,z,Mh);
    //     f<<Tg<<" "<<n0_sol<<" "<<Mg_vir/Ms<<endl;
    //     Tg *= Tgrat;
    // }
    // f.close();

  // marginally unstable case;
    Tg = 2.3e4;
    z = 20; Mh = 4.e7*Ms;
    HALO halo1(Mh,z);
    char* f_Nsol = "profile_Mh.txt";
    Nvir2N0(n0_sol_N,nvir_max,ni,Tg,z,Mh);
    profile(f_Nsol,Tg,n0_sol_N*(mu*m_H)/halo.rho_c,z,Mh);
    printf("Rvir = %3.2e pc, Tvir = %3.2e K\n",halo1.Rvir/pc, halo1.Tvir);

  // R(n0) v.s. M/n_vir
    // char* fname = "boundary_iso.txt";
    // fstream f;
    // f.open(fname, ios::out | ios::trunc );
    // f<<" R ng0 Mg nvir Pvir\n";    
    // while (R<R1){
    //     BOUNDARY(n_vir, Mg_vir, Tg, R, z, Mh);
    //     f<<" "<<R<<" "<<R*halo.rho_c/(mu*m_H)<<" "<<Mg_vir/Ms<<" "<<n_vir<<" "<<n_vir*k_B*Tg<<endl;
    //     R *= Rrat;
    // }
    // f.close();

  // 检查求解函数
    // char* f_Nsol = "cc25profile_Nsol.txt";
    // char* f_Msol = "cc25profile_Msol.txt";
    // Nvir2N0(n0_sol_N,nvir_max,ni,Tg,z,Mh);
    // profile(f_Nsol,Tg,n0_sol_N*(mu*m_H)/halo.rho_c,z,Mh);

    // Mg2N0(n0_sol_M, Mg_max, ni, Tg, z, Mh);
    // profile(f_Msol,Tg,n0_sol_M*(mu*m_H)/halo.rho_c,z,Mh);

    // printf("Rfrom Nsol:%3.2e\n",n0_sol_N*(mu*m_H)/halo.rho_c);
    // printf("Rfrom Msol:%3.2e\n",n0_sol_M*(mu*m_H)/halo.rho_c);
    // printf("z=%3.2e, Mh=%3.2e\tn_mean = %3.2e, fb*Mh = %3.2e Ms, halo.c=%3.2e\n",z,Mh/Ms,n_mean, fb*Mh/Ms,halo.c);
    // printf("local max: from Mg:%3.2e Ms, from n@vir:%3.2e\n", Mg_max/Ms, nvir_max);
    // printf("N0 solution: from Mg:%3.2e, from n@vir:%3.2e\n", n0_sol_M, n0_sol_N);
    // printf("fb=%3.2e\n",Mg_vir/Mh);

  // 算由Nvir2N0得到的n0-->Mg -->fb
    // char* fname = "fb_Mh.txt";
    // fstream f;
    // f.open(fname, ios::out | ios:: trunc);
    // f<<" Mh fb\n";
    // Mh = 1.e5*Ms;
    // double Mh1 = 5.e7*Ms, Mhrat = exp(log(Mh1/Mh)/20.);
    // while (Mh<Mh1){
    //     Nvir2N0(n0_sol_N,nvir_max,ni,Tg,z,Mh);
    //     BOUNDARY(n_vir, Mg_vir,Tg,n0_sol_N*(mu*m_H)/halo.rho_c,z,Mh);
    //     f<<Mh/Ms<<" "<<Mg_vir/Mh<<endl;
    //     Mh *= Mhrat;
    // }
    // f.close();

  //print halo properties
    // double z1 = 20; double Mh1 = 1.e5*Ms;
    // HALO halo(Mh1,z1);
    // printf("1.e5Ms: Rmax =%3.2e\n",R_EQ(Tg,halo.rho_c,halo.Rs));
    // Mh1 = 2.e6*Ms;
    // HALO halo1(Mh1,z1);
    // printf("2.e6Ms: Rmax =%3.2e\n",R_EQ(Tg,halo1.rho_c,halo1.Rs));
    // Mh1 = 1.e7*Ms;
    // HALO halo2(Mh1,z1);
    // printf("1.e7Ms: Rmax =%3.2e\n",R_EQ(Tg,halo2.rho_c,halo2.Rs));

    return 0;
}

 */

// 1. density profile in 3 halos; 
/* int main(){
    double R = .01;
    double Tg = 1.e4;
    double Mg, nvir;
    double R1 = 1.e2, Rrat = exp(log(R1/R)/20.);

    double z = 20, Mh = 1.e5*Ms;
    // HALO halo1(Mh,z);
    // printf("1.e5Ms: Tvir=%3.2e\t",halo1.Tvir);
    // char* pfile1 = "profile1e5.txt";
    // profile(pfile1, Tg, R, z, Mh);

    // Mh = 2.e6*Ms;
    // HALO halo2(Mh, z);
    // printf("2.e6Ms: Tvir=%3.2e\t",halo2.Tvir);
    // char* pfile2 = "profile2e6.txt";
    // profile(pfile2, Tg, R, z, Mh);

    // Mh = 1.e7*Ms;
    // HALO halo3(Mh, z);
    // printf("1.e7Ms: Tvir=%3.2e\t",halo3.Tvir);
    // char* pfile3 = "profile1e7.txt";
    // profile(pfile3, Tg, R, z, Mh);

// 2. 画Mg_max v.s. R
    z = 20, Mh = 2.e6*Ms;
    HALO halo(Mh,z);
    fstream f;
    char* fname = "cc25RM2e6.txt";
    f.open(fname, ios::out | ios::trunc );
    f<<" R n_0 M_g n_vir\n";
    while (R<R1){
        BOUNDARY(nvir,Mg,Tg,R,z,Mh);
        f<<" "<<R<<" "<<R*halo.rho_c/(mu*m_H)<<" "<<Mg/Ms<<" "<<nvir<<endl;
        R *= Rrat;
    }
    f.close();
    cout<<endl;
    return 0;
}
 */

// 3. 算Tvir_crit v.s. redshift; Tg is input parameter.
/* int main(){
    double R = 1.64;
    char* fname = "z_Tvcrit.txt";
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<" z Tv_crit Mh_crit\n";
    double Tg = 2.3e4, Tvir=7.e3, Tv1,Tv0, Tv_sol;
    double z0 = 20, z1 = 10, Mh;
    double n_sol = 0, ni = 1, n_mean, n_mid;
    double N = 1;
    // HALO halo(4.e7*Ms,z0);
    // printf("4.e7halo Tvir:%3.2e\n",halo.Tvir);
    double z_rat = exp( log(z1/z0)/N ); cout<<z_rat<<endl;
    double nvir_max, nvir_m1, nvir_m0;
    for (int i=0; i<N; i++){
        Tv1 = 5.e4,Tv0 =1.e4;
        // Tv1 = 3.e4,Tv0 =1.5e4;
        // Tv1 = 3.e4,Tv0 =2.2e4;
        Nvir2N0(n_sol, nvir_m0, ni, Tg, z0, Mh_Tz(Tv0, z0));
        Nvir2N0(n_sol, nvir_m1, ni, Tg, z0, Mh_Tz(Tv1, z0));
        n_mean = 10* RHO_crit(z0)/(mu*m_H);
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
        f<<" "<<z0<<" "<<Tv_sol<<" "<<Mh_Tz(Tv_sol,z0)/Ms<<endl;
        z0 *= z_rat;
    }
    f.close();

    return 0;
} */


//---------------------------------------------------------------------------------------------------------
//                                  FILE treatment...
//---------------------------------------------------------------------------------------------------------
/* // detect if FILE exist 
    if (FILE *f = fopen(filename, "r")) fclose(f);
    else {
        file.open(filename, ios::out | ios::trunc);
        file<<"R ng_0 M_intg ng_vir\n";
        file.close();
    }
    file.open(filename, ios::out | ios::app);
 */