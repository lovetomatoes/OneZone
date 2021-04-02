/// log scale (= LE_iso_log.cpp, similar with LE_adb_lin.cpp) 

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
#include "dyn.h"
using namespace std;

/* 
g++ -c LE_adb.cpp && g++ -c RK4.cpp && g++ LE_adb.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_adb.out && ./le_adb.out
*/

static double z1 = 20, Mh1 = 1.e5*Ms, Tvir = 1.6e4;
static int const N=2000, n=2;
static int i, hit, it;
static double rho_g0, ng0, M_intg, alpha, beta;
static double x1, dx, x_vir, xrat;
static double K0, eta;
static double n_adb = 1.5;


void BOUNDARY_adb(double& N_VIR, double& MG_VIR, double& r_out, double& T_ave, double R, 
                  bool write, string filename, double cs_eff_2, double z=z1, double Mh=Mh1){
// with a core + r^1 K profile;
    HALO halo1(Mh,z);
    rho_g0 = halo1.rho_c * R;
    K0 = max(K_ISM(z),0.1*halo1.Kvir);

    alpha = halo1.Kvir*pow(rho_g0,gamma_adb-1.) / (k_B*halo1.Tvir/(mu*m_H));
    beta = 4*pi*G*rho_g0*pow(halo1.Rvir,2)/ (k_B*halo1.Tvir/(mu*m_H));
    eta = cs_eff_2/(k_B*halo1.Tvir/(mu*m_H));

    double *x=NULL; double **y=NULL; double *dydx0=NULL;
    double *v=NULL;
    x = new double [N];
    y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    dydx0 = new double [n];

    int c = 7;
    v = new double [c];
    v[0] = n_adb;
    v[1] = R;
    v[2] = alpha;
    v[3] = beta;
    v[4] = halo1.c;
    v[5] = K0/halo1.Kvir;
    v[6] = eta;

  // boundary conditions
    x_vir = 1.;
    i = 0;
    x[i] = x_vir/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    xrat = pow(x_vir/x[0],1./double(N-1) );

    MG_VIR = 0; N_VIR = 0; T_ave = 0;
    fstream file;
    if (write){
        file.open(filename, ios::out | ios::trunc );
        file<<setiosflags(ios::scientific)<<setprecision(3);
        file<<setw(12)<<"xi"<<setw(12)<<"y0"<<setw(12)<<"y1";
        file<<setw(12)<<"n_g"<<setw(12)<<"n_DM"<<setw(12)<<"T_g";
        file<<endl;
    }

    for (i=1;i<N;i++){
        x[i]=x[i-1]*xrat;
        dx = x[i] - x[i-1];

        DyDx_adb(x[i-1], y[i-1], dydx0, c, v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_adb);

        if (y[i][0]<=1.e-6 or isnan(y[i][0])) {
            r_out = x[i];
            break;
        }

        if (x[i]<=x_vir) {
            MG_VIR += pow(halo1.Rvir,3)* 4*pi*pow(x[i],2)*dx * rho_g0 * y[i][0];  //integrate within R_vir
            // T_ave  += Kcore*(1+a*x[i]/r_core)* (mu*m_H) *pow(rho_g0*y[i][0], gamma_adb-1.) /k_B; // T_g
        }
        if (abs(x[i]-x_vir)/x_vir<epE4) {
            N_VIR = rho_g0 * y[i][0]/(mu*m_H);
            r_out = 1;
        }
        if (write) {
            file<<setw(12)<<x[i];
            file<<setw(12)<<y[i][0]<<setw(12)<<y[i][1];
            file<<setw(12)<<rho_g0*y[i][0]/(mu*m_H); // n_g
            file<<setw(12)<<halo1.rho_c/( halo1.c*x[i] * pow((1+halo1.c*x[i]),2) )/(mu*m_H); // n_DM
            // file<<setw(12)<<Kcore*(1+a*x[i]/r_core)* (mu*m_H) *pow(rho_g0*y[i][0], gamma_adb-1.) /k_B; // T_g
            file<<endl;
        }
    }

    if (write) file.close();

    T_ave /= i;

    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    // printf("in BOUNDARY: MGVIR=%3.2e\n",MG_VIR/Ms);
}

void Mg2N0_adb(double& n_sol, double cs_eff_2, double z, double Mh){
    double R0, delta_R, R1, ni=100.; //初始尝试值不能太小 原本ni=1 牛顿迭代求解出错
    double Mg_fw, Mg_bw ;
    double Mg_0;
    double dMdR;
    double dR;
    int it = 0;
    HALO halo(Mh,z);
    double nvir=0, r_out=0, Tg_ave=0;
  // get solution of given Mg=fb*Mh
    it = 0;
    R0 = ni*(mu*m_H)/halo.rho_c; // printf("R0=%3.2e\n",R0);
    dR = R0;
    // printf("\nin Mg2N0_ADB LOOP: it=%d, n0=%3.2e, dR/R0 = %3.2e\n",it, R0*halo.rho_c/(mu*m_H), dR/R0);
    while ( pow(dR,2)>=epE4*pow(R0,2) ){ //dR < 0.1 R0 or it=6, 结束计算 && it<5
        delta_R = epE2*R0;
        // printf("R0=%7.5e,delta_R=%7.5e, R0-=%7.5e, Mg_bw=%7.5e, R0+=%7.5e, Mg_fw%7.5e\n",R0,delta_R, R0-0.5*delta_R, Mg_bw/Ms, R0+0.5*delta_R, Mg_fw/Ms);

        BOUNDARY_adb(nvir,Mg_fw,r_out,Tg_ave,R0+.5*delta_R,false,"",cs_eff_2,z,Mh);
        BOUNDARY_adb(nvir,Mg_bw,r_out,Tg_ave,R0-.5*delta_R,false,"",cs_eff_2,z,Mh);
        BOUNDARY_adb(nvir,Mg_0,r_out,Tg_ave,R0,false,"",cs_eff_2,z,Mh);

        dMdR = ( Mg_fw - Mg_bw )/ delta_R;
        dR = - (Mg_0-fb*Mh)/dMdR;
        R1 = R0 + dR;
        // printf("dMdR=%7.5e, delta_R=%7.5e, Mg_0=%7.5e\n",dMdR/Ms, delta_R, Mg_0/Ms);
        // printf("fbMh=%7.5e, dR=%7.5e\n", fb*Mh/Ms, dR);

    // update R0
        if (R1<0) R0 /= 10.; //暂且取负数R1; 原则上R1有正数下限, below which 解出的Mg=0, 通过取更好的初值暂时解决了
        else R0 = R1;
        it++;
        // if (! isfinite(R0)) { //Mg_fw & Mg_bw = 0; dMdR=0, dR=inf 
        //     printf("#\t#\t#\t#\t#\t#\t#\t#\t#\t\t#\t#\t#\t infinity!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        //     printf("it=%d\n",it);
        //     // printf("Mh=%3.2eMs,z=%3.2e\n", halo.Mh/Ms, z );
        //     break;
        // }
    }

    n_sol = R0*halo.rho_c/(mu*m_H);
    // printf("Mg2N0_ADB: it=%d, R0=%7.5e, dR/R0 = %7.5e, n0=%7.5e\n",it, R0, dR/R0, R0*halo.rho_c/(mu*m_H));
    // printf("Mg2N0: solution: of n0: %3.2e, corresponding Mg:%3.2e, fbMh=%3.2e\n",n_sol, Mg_0/Ms, fb*Mh/Ms );
}

/* 
g++ -c LE_adb.cpp && g++ -c RK4.cpp && g++ LE_adb.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_adb.out && ./le_adb.out
*/

// int main(){
//     clock_t t0 = clock();

//     double R=1, n_adb = 1.5;
//     double n_sol, ni = 100., r_out, Tg_ave, n_vir, Mg_vir;
//     double z = 40, Mh = 8.e5*Ms;
//     z = 30; Mh = 2.e7*Ms;
//     HALO halo(Mh,z);
//     double cs_eff_2 = k_B*halo.Tvir/(mu*m_H);
//     cs_eff_2 = 0;
//     BOUNDARY_adb(n_vir,Mg_vir,r_out,Tg_ave,R,true,"z40M8e5R1.txt",cs_eff_2,z,Mh);
//     Mg2N0_adb(n_sol, cs_eff_2, z, Mh);

//     R = n_sol*(mu*m_H)/halo.rho_c;
//     clock_t t1 = clock();
//     printf("log n_sol=%.3e, dt=%.2f\n",n_sol, (double)(t1-t0)/CLOCKS_PER_SEC);

//   // R(n0) v.s. Mg, n_vir
//     R = 1.e-2;
//     int Nspan = 100;
//     double R1 = 1.e2, Rrat = exp(log(R1/R)/double(Nspan));
//     string fname = "RM.txt";
//     fstream f;
//     f.open(fname, ios::out | ios::trunc );
//     f<<setiosflags(ios::scientific)<<setprecision(3);
//     f<<setw(12)<<"R"<<setw(12)<<"ng0"<<setw(12)<<"Mg"<<setw(12)<<"nvir"<<setw(12)<<"Tg_ave"<<setw(12)<<"Rout"<<endl;    
//     while (R<R1){
//         BOUNDARY_adb(n_vir,Mg_vir,r_out,Tg_ave,R,false,"",cs_eff_2,z,Mh);
//         f<<setw(12)<<R<<setw(12)<<R*halo.rho_c/(mu*m_H)<<setw(12)<<Mg_vir/Ms<<setw(12)<<n_vir<<setw(12)<<Tg_ave<<setw(12)<<r_out<<endl;
//         R *= Rrat;
//     }
//     f.close();
//     return 0;
// }