/// linear scale

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
g++ -c LE_adb_lin.cpp && g++ -c RK4.cpp && g++ LE_adb_lin.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_adb_lin.out && ./le_adb_lin.out
*/

static int N=2000, i, hit;
static int const n=2;
static double a, alpha, beta, rho_g0, ng0, M_intg;
static double x1, dx, x_vir, xrat;
static double r_core, Kcore, Kc;
static double z1 = 20, Mh1 = 1.e5*Ms, Tvir = 1.6e4;

void profile_adb_Kc(string filename, double R, double n_adb, double z=z1, double Mh=Mh1){
    HALO halo1(Mh,z);
    rho_g0 = halo1.rho_c * R;
    ng0 = rho_g0/(mu*m_H);
    double gamma = 1 + 1./n_adb;

    Kc = K_ISM(z);
    a = sqrt( (n_adb+1)*Kc *pow(rho_g0,gamma-2.)/(4*pi*G) ); // adiabatic n_adb
    alpha = a/halo1.Rs;

    // set x1, largest xi
    x_vir = halo1.Rvir/a;
    x1 = 2*x_vir;
    // printf("in PROFILE:nH0=%3.2e, R=%3.2e, a=%3.2e, alpha=%3.2e, x1=%3.2e\n", nH0, R, a, alpha, x1);
    // printf("halo: delta_rhoc=%3.2e, Rs=%3.2e cm Rvir=%3.2e\n", halo1.rho_c,halo1.Rs,halo1.Rvir);
    double* x = new double [N];
    double** y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    int c = 3;
    double* v = new double [c];
    v[0] = n_adb;
    v[1] = R;
    v[2] = alpha;

    double* dydx0 = new double [n];
    fstream file;
    file.open(filename, ios::out | ios::trunc );

    file<<setw(12)<<"xi"<<setw(12)<<"r"<<setw(12)<<"y0"<<setw(12)<<"y1";
    file<<setw(12)<<"theta_n"<<setw(12)<<"n_g"<<setw(12)<<"n_DM"<<setw(12)<<"T_g";
    file<<endl;
    // dx & boundary conditions
    i = 0;
    dx = x1/N;
    // boundary x[0]=0, but 0 causes sigularity, using x[0]=dx/1000<<dx
    x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    // printf("center boundary: r0=%3.2e (Rvir)\t dr=%3.2e (Rvir)\n",a*x[0]/halo1.Rvir,a*dx/halo1.Rvir);
    // printf("IN PROFILE_ADB: center boundary: r0=%3.2e (Rvir)\t dr=%3.2e (Rvir)\n",a*x[0]/halo1.Rvir,a*dx/halo1.Rvir);
    // printf("n_adb=%3.2f\n",n_adb);
    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx_adb_Kc(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_adb_Kc);
        x[i]=x[i-1]+dx;
        if (pow(y[i][0],n_adb)<=1.e-5) break;
        if (i%10 ==0){
        file<<x[i]<<setw(12)<<a*x[i]/halo1.Rvir;
        file<<setw(12)<<y[i][0]<<setw(12)<<y[i][1];
        file<<setw(12)<<pow(y[i][0],n_adb)<<setw(12)<<ng0*pow(y[i][0],n_adb);
        file<<setw(12)<<halo1.rho_c/( alpha*x[i] * pow((1+alpha*x[i]),2) )/(mu*m_H);// n_DM
        file<<setw(12)<<Kc* (mu*m_H) *pow( rho_g0*pow(y[i][0],n_adb) , gamma-1.) /k_B; //T_g
        file<<endl;
        }
    }
    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    file.close();
}

void profile_adb_Kfit(string filename,double& N_VIR, double& MG_VIR, double R, double n_adb, double z=z1, double Mh=Mh1){
// with a core + r^1 K profile;
    HALO halo1(Mh,z);
    // int const N = 2000000;
    // int const n = 2;
    // int i;
    double rho_g0, ng0, Tgc, M_intg = 0, a, alpha, beta;
  // adiabatic, entropy K = k_B*T*n^(-2/3), n=pow(k_B*T / S, 1.5); // n \propto T^1.5
    double Kcore = 0.1*halo1.Kvir;
    double gamma = 1. + 1./n_adb;
    printf("K_ISM=%3.2e, K_vir=%3.2e, K_core=0.1Kvir=%3.2e\n",K_ISM(z),halo1.Kvir,0.1*halo1.Kvir);
    rho_g0 = halo1.rho_c * R;
    ng0 = rho_g0/(mu*m_H);
    r_core = 0.1*halo1.Rvir;
    a = r_core; alpha = a/halo1.Rs;
    beta = 4*pi*G*pow(r_core,2)/Kcore/pow(rho_g0,gamma-2);

  // set r1, largest physical r
    x_vir = halo1.Rvir/a;
    x1 = 2.*x_vir;
    double* x = new double [N];
    double** y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    int c = 4;
    double* v = new double [c];
    v[0] = n_adb;
    v[1] = R;
    v[2] = alpha;
    v[3] = beta;
    double* dydx0 = new double [n];
    fstream file;
    file.open(filename, ios::out | ios::trunc );
  // file<<"xi r y0 y1 n_g\n";
    file<<"xi r y0 y1 n_g n_DM T_g\n";
    // boundary conditions
    i = 0;
    dx = x1/N;
    // boundary x[0]=0, but 0 causes sigularity
    x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    printf("IN PROFILE_ADB_KFIT: center boundary: r0=%3.2e (Rvir)\t dr=%3.2e (Rvir)\n",a*x[0]/halo1.Rvir,dx/x_vir);
    printf("n_adb=%3.2f, beta=%3.2e\n",n_adb, beta);
      
    MG_VIR = 0; N_VIR = 0;
    for (i=1;i<N;i++){
        x[i]=x[i-1]+dx;
        DyDx_adb_fit(x[i-1], y[i-1], dydx0, c, v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_adb_fit);
        if (x[i] < halo1.Rvir/a) {
            MG_VIR += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 * y[i][0];  //integrate within R_vir
        }
        else {
            N_VIR = rho_g0 * y[i][0]/(mu*m_H);
            break;
        }

        if (y[i][0]<=1.e-6) break;
        if (i%10 ==0){
            file<<x[i]<<setw(12)<<a*x[i]/halo1.Rvir;
            file<<setw(12)<<y[i][0]<<setw(12)<<y[i][1];
            file<<setw(12)<<ng0*y[i][0];
            file<<setw(12)<<halo1.rho_c/( alpha*x[i] * pow((1+alpha*x[i]),2) )/(mu*m_H);// n_DM
            file<<setw(12)<<Kcore*(1+a*x[i]/r_core)* (mu*m_H) *pow(rho_g0*y[i][0], gamma-1.) /k_B; // T_g
            file<<endl;
        }
    }
    printf("OUTSIDE: MG_VIR=%3.2E, N_VIR=%3.2E\n",MG_VIR/Ms, N_VIR);
    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    file.close();
    printf("z = %3.2f\tMh = %3.2e Ms\t rs%3.2e cm\trhoc=%3.2e/cc\n",halo1.z,halo1.Mh/Ms,halo1.Rs,halo1.rho_c);
    printf("rho_g0=%3.2e, rho_crit=%3.2e, delta_0=%3.2e\n",rho_g0, halo1.rho_crit, halo1.delta0);
}

void BOUNDARY_adb(double& N_VIR, double& MG_VIR, double& r_out, double& T_ave, double R, double z=z1, double Mh=Mh1){
// with a core + r^1 K profile;
    HALO halo1(Mh,z);
    int i, c;
    double* v;
    rho_g0 = halo1.rho_c * R;
    double n_adb = 1.5, gamma = 1. + 1./n_adb;;
    ng0 = rho_g0/(mu*m_H);

    // set r1, largest physical r
    double* x = new double [N];
    double** y = new double* [N];
    for (i=0;i<N;i++) y[i] = new double [n];
    double* dydx0 = new double [n];


  // boundary conditions
    i = 0;
    // printf("IN BOUNDARY: center boundary: r0=%3.2e (Rvir)\t dr=%3.2e (Rvir)\n",a*x[0]/halo1.Rvir,a*dx/halo1.Rvir);
    bool fit = true;
    // fit = false;
    fit = (0.1*halo1.Kvir > K_ISM(z)); // printf("fit?%s\n",(fit)?"fit core":"const");
    if (fit){
        Kcore = 0.1*halo1.Kvir;
        r_core = 0.1*halo1.Rvir;
        a = r_core; alpha = a/halo1.Rs;
        beta = 4*pi*G*pow(r_core,2)/Kcore/pow(rho_g0,gamma-2);
        c = 4;
        v = new double [c];
        v[0] = n_adb;
        v[1] = R;
        v[2] = alpha;
        v[3] = beta;

        x_vir = halo1.Rvir/a;
        x1 = 2.*x_vir;
        dx = x1/N;
        // boundary x[0]=0, but 0 causes sigularity
        x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    }
    else{
        // printf("\n\nNOW USING KC\n\n");
        Kc = K_ISM(z);
        a = sqrt( (n_adb+1)*Kc *pow(rho_g0,gamma-2.)/(4*pi*G) ); // adiabatic n_adb
        alpha = a/halo1.Rs;
        c = 3;
        v = new double [c];
        v[0] = n_adb;
        v[1] = R;
        v[2] = alpha;

        x_vir = halo1.Rvir/a;
        x1 = 2.*x_vir;
        dx = x1/N;
        x[i] = dx/1.e8; y[i][0] = 1; y[i][1] = 0; // adiabatic case
    }

    MG_VIR = 0; N_VIR = 0; T_ave = 0;
    for (i=1;i<N;i++){
        x[i]=x[i-1]+dx;
        if (fit){
            DyDx_adb_fit(x[i-1], y[i-1], dydx0, c, v);
            rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_adb_fit);
            if ( y[i][0]<=1.e-5 or isnan(y[i][0]) ) { // r < r_vir and gas truncated
                N_VIR = 0;
                // printf("Kfit truncated i=%d\n",i); // if truncated at i=1, MG_VIR = 0; N_VIR = 0; T_ave = 0;
                break;
            }
            if (x[i] < halo1.Rvir/a) {
                MG_VIR += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 * y[i][0] ;  //integrate within R_vir
                T_ave  += Kcore*(1+a*x[i]/r_core)* (mu*m_H) *pow(rho_g0*y[i][0], gamma-1.) /k_B; // T_g
            }
            else {
                N_VIR = rho_g0 * y[i][0]/(mu*m_H);
                break;
            }
        }
        else{
            DyDx_adb_Kc(x[i-1], y[i-1], dydx0, c, v);
            rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_adb_Kc);
            if ( pow(y[i][0],n_adb)<=1.e-5 or isnan(y[i][0]) or isnan(pow(y[i][0],n_adb)) ) { // r < r_vir and gas truncated
                // printf("\nTELL4, theta_n=%3.2e\n",pow(y[i][0],n_adb));
                N_VIR = 0;
                // printf("Kc truncated i=%d\n",i); // if truncated at i=1, MG_VIR = 0; N_VIR = 0; T_ave = 0;
                break;
            }
            if (x[i] < halo1.Rvir/a) {
                MG_VIR += pow(a,3)* 4*pi*pow(x[i],2)*dx * rho_g0 * pow(y[i][0], n_adb);  //integrate within R_vir
                T_ave += Kc* (mu*m_H) *pow( rho_g0*pow(y[i][0],n_adb) , gamma-1.) /k_B;
            }
            else {
                N_VIR = rho_g0 * pow(y[i][0],n_adb)/(mu*m_H);
                break;
            }
        }
    }

    T_ave /= i; 

    if (x[i]*a < halo1.Rvir) r_out = x[i]*a / halo1.Rvir;
    else r_out = 1.;
    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] v; delete [] dydx0;
    // printf("in BOUNDARY: MGVIR=%3.2e\n",MG_VIR/Ms);
}

void Mg2N0_adb(double& n_sol, double ni, double z, double Mh){
    double R0, delta_R, R1;
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

        BOUNDARY_adb(nvir, Mg_fw, r_out,Tg_ave, R0+.5*delta_R, z, Mh); BOUNDARY_adb(nvir, Mg_bw, r_out,Tg_ave, R0-.5*delta_R, z, Mh);
        BOUNDARY_adb(nvir, Mg_0, r_out,Tg_ave, R0, z, Mh);
        dMdR = ( Mg_fw - Mg_bw )/ delta_R;
        dR = - (Mg_0-fb*Mh)/dMdR;
        R1 = R0 + dR;
        // printf("dMdR=%7.5e, delta_R=%7.5e, Mg_0=%7.5e\n",dMdR/Ms, delta_R, Mg_0/Ms);
        // printf("fbMh=%7.5e, dR=%7.5e\n", fb*Mh/Ms, dR);

    // update R0
        if (R1<0) R0 /= 10.; //不能取负数R1; 原则上R1有正数下限, below which 解出的Mg=0, 通过取更好的初值暂时解决了
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
    printf("\nin Mg2N0_ADB LOOP: it=%d, R0=%7.5e, dR/R0 = %7.5e, n0=%7.5e\n",it, R0, dR/R0, R0*halo.rho_c/(mu*m_H));
    // printf("Mg2N0: solution: of n0: %3.2e, corresponding Mg:%3.2e, fbMh=%3.2e\n",n_sol, Mg_0/Ms, fb*Mh/Ms );
}

/* 
g++ -c LE_adb_lin.cpp && g++ -c RK4.cpp && g++ LE_adb_lin.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o le_adb_lin.out && ./le_adb_lin.out
*/

// 检查 Mg2N0
int main(){
    double R, n_adb = 1.5;
    double n_sol, ni = 100., r_out, Tg_ave, n_vir, Mg_vir;
    double z = 20, Mh = 1.e6*Ms; //z20 M6 --> Kfit
    HALO halo(Mh,z);
    double n_nfw = fb*halo.Rho_r(halo.Rvir)/(mu*m_H);
    string f_Nsol; 
  // tests --especially Tg profile
    clock_t t0 = clock();
    // printf("0.1*halo1.Kvir=%3.2e, K_ISM=%3.2e, Tvir=%3.2e\n",0.1*halo.Kvir,K_ISM(z),halo.Tvir );
    Mg2N0_adb(n_sol, ni, z, Mh);

    clock_t t1 = clock();
    printf("lin n_sol=%.3e, dt=%.2f\n",n_sol, (double)(t1-t0)/CLOCKS_PER_SEC);

    R = n_sol*(mu*m_H)/halo.rho_c;
    // printf("R=%3.2e\n",R);
    // profile_adb_Kc("R1e-2KcMh6.txt",R,n_adb,z,Mh);
    // profile_adb_Kfit("pKfit.txt",n_vir, Mg_vir,R,n_adb,z,Mh);
    // BOUNDARY_adb(n_vir, Mg_vir, r_out,Tg_ave, R, z, Mh);
    // printf("r_out=%3.2e\n",r_out);
    // printf("halo Tvir=%3.2e, Tg_ave=%3.2e\n",halo.Tvir,Tg_ave);

  // R(n0) v.s. n_vir
    R = 1.e-5;
    int Nspan = 100;
    double R1 = 1.e4, Rrat = exp(log(R1/R)/double(Nspan));
    string fname = "RM.txt";
    fstream f;
    f.open(fname, ios::out | ios::trunc );
    f<<setiosflags(ios::scientific)<<setprecision(3);
    f<<setw(12)<<"R"<<setw(12)<<"ng0"<<setw(12)<<"Mg"<<setw(12)<<"nvir"<<setw(12)<<"Tg_ave"<<setw(12)<<"Rout"<<endl;    
    while (R<R1){
        BOUNDARY_adb(n_vir, Mg_vir, r_out,Tg_ave, R, z, Mh);
        f<<setw(12)<<R<<setw(12)<<R*halo.rho_c/(mu*m_H)<<setw(12)<<Mg_vir/Ms<<setw(12)<<n_vir<<setw(12)<<Tg_ave<<setw(12)<<r_out<<endl;
        R *= Rrat;
    }
    f.close();

  // 算由Mg2N0得到的n0-->n_vir? v.s. n_nfw
  // 1. Tvir v.s. n_sol
    // string fname = "Tv_n0Ktran_Nvir.txt";
    // fstream f;
    // f.open(fname, ios::out | ios:: trunc);
    // f<<"Tv n0_sol n_vir n_nfw Mh Mg_vir n_KISM n_core Kcore KISM Tgave r_out\n";
    // // double Tv = 1.e3, Tv1 = 2.e4; // for Kcore
    // // double Tv = 1.e2, Tv1 = 7.e8; // for Kc
    // double Tv = 1.e2, Tv1 = 2.e4; // for transition
    // double Tvrat = exp(log(Tv1/Tv)/40.);
    // while (Tv<Tv1){
    //     Mh = Mh_Tz(Tv,z);
    //     Mg2N0_adb(n_sol, ni, z, Mh);
    //     HALO halo1 (Mh, z);
    //     R = n_sol*(mu*m_H)/halo1.rho_c;
    //     n_nfw = fb*halo1.Rho_r(halo1.Rvir)/(mu*m_H);
    //     BOUNDARY_adb(n_vir, Mg_vir, r_out,Tg_ave, R, z, Mh);
    //     f<<Tv<<setw(12)<<n_sol<<setw(12)<<n_vir<<setw(12)<<n_nfw<<setw(12)<<Mh/Ms<<setw(12)<<Mg_vir/Ms;
    //     f<<setw(12)<<6* pow(halo1.Tvir/1000., 1.5)<<setw(12)<<N_CORE(z)<<setw(12)<<0.1*halo1.Kvir<<setw(12)<<K_ISM(z)<<setw(12)<<Tg_ave<<setw(12)<<r_out<<endl;
    //     Tv *= Tvrat;
    // }
    //  f.close();
  // 2. z v.s. n_sol
    // string fname = "z_n0Ktran_Nvir.txt";
    // fstream f;
    // f.open(fname, ios::out | ios:: trunc);
    // f<<"z n0_sol n_vir n_nfw Mh Mg_vir n_KISM n_core Kcore KISM Tvir Tgave r_out\n";
    // z = 35;  Mh = 1.e7*Ms;
    // double z1 = 10, zrat = exp(log(z1/z)/10.);
    // while (z>=z1){
    //     Mg2N0_adb(n_sol, ni, z, Mh);
    //     HALO halo1 (Mh, z);
    //     R = n_sol*(mu*m_H)/halo1.rho_c;
    //     n_nfw = fb*halo1.Rho_r(halo1.Rvir)/(mu*m_H);
    //     BOUNDARY_adb(n_vir, Mg_vir, r_out,Tg_ave, R, z, Mh);
    //     f<<z<<setw(12)<<n_sol<<setw(12)<<n_vir<<setw(12)<<n_nfw<<setw(12)<<Mh/Ms<<setw(12)<<Mg_vir/Ms;
    //     f<<setw(12)<<6* pow(halo1.Tvir/1000., 1.5)<<setw(12)<<N_CORE(z)<<setw(12)<<0.1*halo1.Kvir<<setw(12)<<K_ISM(z);
    //     f<<setw(12)<<halo1.Tvir<<setw(12)<<Tg_ave<<setw(12)<<r_out<<endl;
    //     z *= zrat;
    // }
    // f.close();


    return 0;
}
