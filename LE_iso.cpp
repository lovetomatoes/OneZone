//: log scale (= LE_iso_log.cpp, proved better than LE_iso_lin.cpp) 
//: computation time much better than lin scale
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "class_halo.h"
#include "PARA.h"
#include "RK4.h"
#include "LE_iso.h"
#include "dyn.h"
using namespace std;

/* 
g++ -c LE_iso.cpp RK4.cpp && g++ LE_iso.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_iso.out && ./le_iso.out
*/

static int N=100, i, hit;
static int const n=2;
static double a, alpha, rho_g0, M_intg;
static double x1, dx, x_vir, xrat;

double Tvir = 1.6e4;
double z1 = 20; double Mh1 = 1.e5*Ms;
//double z = 30; double Mh = Mh_Tz(Tvir,z);

// a profile of density distributions, Phi(x), x=r/a, Phi = -log(rho/rho_g0)
void profile(string filename, double Tg, double R, double z=z1, double Mh=Mh1){
    fstream file;
    file.open(filename, ios::out | ios::trunc );
    file<<setiosflags(ios::scientific)<<setprecision(3);
    file<<setw(12)<<"r_pc"<<setw(12)<<"r_Rvir"<<setw(12)<<"r_Rs"<<setw(12)<<"phi"<<setw(12)<<"psi";
    file<<setw(12)<<"ng_ng0"<<setw(12)<<"ng"<<setw(12)<<"nDM_ng0"<<setw(12)<<"nDM";
    file<<setw(12)<<"M_intg"<<setw(12)<<"Bfx"<<setw(12)<<"gx";
    file<<endl;

    HALO halo1(Mh,z);
    rho_g0 = halo1.rho_c * R;
    a = sqrt(k_B*Tg/(4*pi*G*mu*m_H*rho_g0));
    alpha = a/halo1.Rs;
// B for checking DM gravity w/ Suto, Sasaki, Makino 1998
    double B = 4*pi*G*(mu*m_H)*halo1.delta0*halo1.rho_crit*pow(halo1.Rs,2)/(k_B*Tg);

    // printf("scaling factor a = %3.2e pc, alpha = %3.2e pc, delta_rhoc=%3.2e g/cc, Rs=%3.2e pc Rvir=%3.2e pc\n", 
    // a/pc,alpha/pc,halo1.rho_c,halo1.Rs/pc,halo1.Rvir/pc);

    double *x=NULL; double **y=NULL; double *dydx0=NULL;
    double *v=NULL;

    x = new double [N];
    y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    dydx0 = new double [n];

    int c = 2;
    v = new double [c];
    v[0] = alpha; v[1] = R;

// x1, dx & boundary conditions; log scale grids
    x_vir = halo1.Rvir/a;
    x1 = 2.*x_vir;
    // x1 = 100; self-gravity
    // x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    i = 0; x[i] = x_vir/1.e8; y[i][0] = 0; y[i][1] = 0; // isothermal case
// boundary x[0]=0, but 0 causes sigularity, using x[0]=x_vir/1e8<<dx
    xrat = exp( log(x1/x[0])/(N-1) );

    M_intg = 0;
    hit = 0;
    for (i=1;i<N;i++){
        x[i]=x[i-1]*xrat;
        dx = x[i] - x[i-1];
    //integrate within R_vir
        if (x[i]>x_vir && hit==0) {
            x[i] = x_vir;
            dx = x_vir - x[i-1];
            hit = 1;
        }

        DyDx_iso(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
        // if (x[i] == x_vir) printf("ng_vir = %5.3e\n",rho_g0*exp(-y[i][0])/(mu*m_H));
        if (x[i]<=x_vir) M_intg += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);

    // "r_pc r_Rvir r_Rs phi psi ng_ng0 ng nDM_ng0 nDM M_intg Bfx gx";
        file<<setw(12)<<a*x[i]/pc<<setw(12)<<a*x[i]/halo1.Rvir<<setw(12)<<alpha*x[i];
        file<<" "<<y[i][0]<<" "<<y[i][1];
        file<<" "<<exp(-y[i][0]); //rhog_over_rhog0
        file<<" "<<rho_g0*exp(-y[i][0])/(mu*m_H); // n_gas
        file<<" "<< 1./(R* alpha*x[i]* pow((1+alpha*x[i]),2) );// rhoDM_over_rhog0
        file<<" "<<halo1.rho_c/( alpha*x[i] * pow((1+alpha*x[i]),2) )/(mu*m_H);// n_DM
        file<<" "<<M_intg/Ms; //M(r) in solar mass
        file<<" "<<B*(1- log(1+alpha*x[i])/ (alpha*x[i])); // B*fx Eq(11) Suto98; gas self-gravity 
        file<<" "<<B/2.*(alpha*x[i]); // Eq(36) Suto98; g(x)= -Φ 
        file<<endl; // B*fx Eq(11) Suto98; gas self-gravity

        // printf("x[%d]/x_vir=%3.2f\n",i,x[i]/x_vir);
    }
    file.close();

    // printf("Mg_vir= %5.3e Ms\n",M_intg/Ms);

    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] v; delete [] dydx0;

    // printf("z = %3.2f\tMh = %3.2e Ms\t rs%3.2e cm\trhoc=%3.2e/cc\n",halo1.z,halo1.Mh/Ms,halo1.Rs,halo1.rho_c);
    // printf("rho_g0=%3.2e, rho_crit=%3.2e, delta_0=%3.2e\n",rho_g0, halo1.rho_crit, halo1.delta0);
}

void BOUNDARY(double& N_VIR, double& MG_VIR, double Tg, double R, double z=z1, double Mh=Mh1){
    MG_VIR = 0;
    HALO halo1(Mh,z);
    rho_g0 = halo1.rho_c * R;
    a = sqrt(k_B*Tg/(4*pi*G*mu*m_H*rho_g0));
    alpha = a/halo1.Rs;

    double *x=NULL; double **y=NULL; double *dydx0=NULL;
    double *v=NULL;

    x = new double [N];
    y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    dydx0 = new double [n];

    int c = 2;
    v = new double [c];
    v[0] = alpha; v[1] = R;

// x1, dx & boundary conditions; log scale grids
    x_vir = halo1.Rvir/a;
// truncated at Rvir exactly
    // *********** //
    x1 = x_vir;    // 
    // *********** //
    i = 0; x[i] = x_vir/1.e8; y[i][0] = 0; y[i][1] = 0; // isothermal case

    xrat = exp( log(x1/x[0])/(N-1) );

    for (i=1;i<N;i++){
        x[i]=x[i-1]*xrat;
        dx = x[i] - x[i-1];  
        DyDx_iso(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_iso);
        MG_VIR += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 *exp(-y[i][0]);
    }
    N_VIR = rho_g0*exp(-y[N-1][0])/(mu*m_H);

    // printf("outside loop: i=%d, N=%d, n_vir=%5.3e/cc, MG_VIR=%5.3e Ms,\n",
    // i,N,N_VIR,MG_VIR/Ms);

    for (i=0;i<N;i++) delete [] y[i];
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
    while ( pow(dR,2)>=epE4*pow(R0,2)){ //dR < 0.1 R0 or it=6, 结束计算
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
    }
    // printf("\nMAXIMUM: it=%d, dR/R0 = %3.2e\n",it,dR/R0);

    BOUNDARY(nvir_max, Mg_vir, Tg, R0, z, Mh); //求解nvir_max
    // printf("for nvir_max=%3.2e, n0=%3.2e\n",nvir_max,R0*halo.rho_c/(mu*m_H));
    // printf("z=%.1f, Mh=%3.2e Ms\n",z,Mh/Ms);
// unstable criterion: nvir_max < n_nfw (cosmic mean density at z) 
    double n_mean = RHO_crit(z)/(mu*m_H);
    double n_nfw = fb*halo.Rho_r(halo.Rvir)/(mu*m_H);
    // printf("nmean=%3.2e, n_nfw=%3.2e, nvir_max=%3.2e:\n",n_mean,n_nfw, nvir_max);
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
        while ( pow(dR,2)>=epE4*pow(R0,2)){ //dR < 0.1 R0 or it=6, 结束计算 wli! 取消it限制??
            delta_R = epE2*R0;
            BOUNDARY(nvir_fw, Mg_vir, Tg, R0+.5*delta_R, z, Mh); BOUNDARY(nvir_bw, Mg_vir, Tg, R0-.5*delta_R, z, Mh);
            BOUNDARY(nvir_0, Mg_vir, Tg, R0, z, Mh);
            dndR = ( nvir_fw - nvir_bw )/ delta_R;
            dR = - (nvir_0-n_mean)/dndR; R1 = R0 + dR;
        // update R0
            if (R1<0) R0 = R0/10.;
            else R0 = R1;
            it++;
        }
        // printf("\nin LOOP: it=%d, dR/R0 = %3.2e\n",it,dR/R0);
        n_sol = R0*halo.rho_c/(mu*m_H);
        // printf("Nvir2N0: n_nfw=%3.2e, solution: of n0: %3.2e\n",n_nfw,n_sol);
    }
}


/* 
g++ -c LE_iso.cpp RK4.cpp && g++ LE_iso.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_iso.out && ./le_iso.out
*/


/* int main(){
    double R = .01;
    double Tg = 1.e4;
    double Mg, nvir;
    double R1 = 1.e2, Rrat = exp(log(R1/R)/20.);
    double z = 20, Mh = 1.e5*Ms;
    fstream f;
    string fname;
    clock_t t0 = clock();

// 0. profile & boundary
    // double z0 = 25, Mh=1.e5*Ms;
    // double nsol, Mgas;
    // double Tg = 8000;
    // fname = to_string(N)+"_lin.txt";
    // profile(fname, Tg, 1., z0, Mh);
    // // cout<<" \n\nBOUNDARY: \n";
    // BOUNDARY(nsol,Mgas,Tg,1.,z0,Mh);

// 1. density profile in 3 halos; 
    // HALO halo1(Mh,z);
    // printf("1.e5Ms: Tvir=%3.2e\t",halo1.Tvir);
    // string pfile1 = "profile1e5.txt";
    // profile(pfile1, Tg, R, z, Mh);

    // Mh = 2.e6*Ms;
    // HALO halo2(Mh, z);
    // printf("2.e6Ms: Tvir=%3.2e\t",halo2.Tvir);
    // string pfile2 = "profile2e6.txt";
    // profile(pfile2, Tg, R, z, Mh);

    // Mh = 1.e7*Ms;
    // HALO halo3(Mh, z);
    // printf("1.e7Ms: Tvir=%3.2e\t",halo3.Tvir);
    // string pfile3 = "profile1e7.txt";
    // profile(pfile3, Tg, R, z, Mh);

// 2.1 画Mg (n_vir) v.s. R
    z = 20, Mh = 1.e7*Ms;
    HALO halo(Mh,z);
    fname = "boundary_Tg1e4z20Mh1e7.txt";
    f.open(fname, ios::out | ios::trunc );
    f<<" R n_0 M_g n_vir\n";
    while (R<R1){
        BOUNDARY(nvir,Mg,Tg,R,z,Mh);
        f<<" "<<R<<" "<<R*halo.rho_c/(mu*m_H)<<" "<<Mg/Ms<<" "<<nvir<<endl;
        R *= Rrat;
    }
    f.close();

// 2.2 检查Nvir2N0是否给出正确的nvir_max, n0 solution
    double n0_sol_N,nvir_max,ni=1;
    string f_Nsol = "profile_Nsol.txt";
    Nvir2N0(n0_sol_N,nvir_max,ni,Tg,z,Mh);
    profile(f_Nsol,Tg,n0_sol_N*(mu*m_H)/halo.rho_c,z,Mh);

// 3. Tvir_crit value v.s. Tg (z as input parameter)
    //:   gravity v.s. pressure, hope linear 
    //:  halo changes, using n_nfw=10*n_mean as fiducial boundary condition
    // fname = "Tg_Tvcrit.txt";
    // f.open(fname, ios::out | ios::trunc );
    // f<<"Tg Tv_crit\n";
    // Tg = 7.e3;
    // double Tg1 = 2.e4, Tv1,Tv0, Tv_sol;
    // double n_sol = 0, ni = 1, n_mean, n_nfw, n_mid;
    // double Nspan = 10;
    // double Tgrat = exp( log(Tg1/Tg)/Nspan ); cout<<Tgrat<<endl;
    // double nvir_max, nvir_m1, nvir_m0;
    // for (int i=0; i<N; i++){
    //     Tv1 = 5.e4,Tv0 = 5.e3;
    //     Nvir2N0(n_sol, nvir_m0, ni, Tg, z, Mh_Tz(Tv0, z));
    //     Nvir2N0(n_sol, nvir_m1, ni, Tg, z, Mh_Tz(Tv1, z));
    //     n_mean = RHO_crit(z)/(mu*m_H);
    //     n_nfw = 10*n_mean;
    //     if (nvir_m0<n_nfw or nvir_m1>n_nfw) {
    //         printf("\n!!!!!initial value not correct\n");
    //         break;
    //     }
    //     while (pow(Tv1/Tv0-1, 2)>epE2){
    //         Nvir2N0(n_sol, n_mid, ni, Tg, z, Mh_Tz((Tv0+Tv1)/2., z));
    //         if (n_mid>=n_nfw) {
    //             Tv0 = (Tv0+Tv1)/2.;
    //             Nvir2N0(n_sol, nvir_m0, ni, Tg, z, Mh_Tz(Tv0, z));
    //         }
    //         else {
    //             Tv1 = (Tv0+Tv1)/2.;
    //             Nvir2N0(n_sol, nvir_m1, ni, Tg, z, Mh_Tz(Tv1, z));
    //         }
    //     }
    //     Tv_sol = (Tv0+Tv1)/2.;
    //     double cs2 = gamma_adb*k_B*8000/(mu*m_H);
    //     f<<" "<<Tg<<" "<<Tv_sol<<endl;
    //     Tg *= Tgrat;
    // }
    // f.close();
    clock_t t1 = clock();
    printf("executing time: %.2fs s\n", (double)(t1-10)/CLOCKS_PER_SEC);
    return 0;
}
 */

// 4. 算由Nvir2N0得到的n0-->Mg -->fb
/* int main(){
    double R = 1.e-2;
    double Tg = 1.e4;
    // double R1 = 1.e2, Rrat = exp(log(R1/R)/20.);
    double n0_sol_M, Mg_max, n0_sol_N, nvir_max, n0_sol;
    double ni = 1;
    double z = 20, Mh = 1.e7*Ms;
    HALO halo(Mh,z);
    double n_mean = RHO_crit(z)/(mu*m_H);
    double n_vir, Mg_vir;

    string fname = "fb_Mh.txt";
    fstream f;
    f.open(fname, ios::out | ios:: trunc);
    f<<" Mh fb\n";
    Mh = 1.e5*Ms;
    double Mh1 = Mh_Tz(1.1e4,z), Mhrat = exp(log(Mh1/Mh)/20.);
    while (Mh<Mh1){
        Nvir2N0(n0_sol_N,nvir_max,ni,Tg,z,Mh);
        BOUNDARY(n_vir, Mg_vir,Tg,n0_sol_N*(mu*m_H)/halo.rho_c,z,Mh);
        f<<Mh/Ms<<" "<<Mg_vir/Mh<<endl;
        Mh *= Mhrat;
    }
    f.close();
    return 0;
}
 */

// 5. 算Tvir_crit v.s. redshift; Tg is input parameter.
/* int main(){
    double R = 1.64;
    string fname = "z_Tvcrit.txt";
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<" z Tv_crit Mh_crit\n";
    double Tg = 2.3e4, Tvir=7.e3, Tv1,Tv0, Tv_sol;
    double z0 = 20, z1 = 10, Mh;
    double n_sol = 0, ni = 1, n_mean, n_mid;
    double Nspan = 10;
    double z_rat = exp( log(z1/z0)/Nspan ); cout<<z_rat<<endl;
    double nvir_max, nvir_m1, nvir_m0;
    for (int i=0; i<Nspan; i++){
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
}
 */