#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include<algorithm>

#include "kpd.h"
#include "reaction.h"
#include "class_gas.h"
#include "my_linalg.h"
#include "Newton.h"
#include "thermo.h"
#include "dyn.h"
#include "LE_iso.h"
#include "LE_adb.h"

#include "class_halo.h"
#include "PARA.h"
#include "RK4.h"

using namespace std;
/* to compile this file: 
g++ -c class_gas.cpp
g++ class_gas.o -o gas
./gas
*/

// constructor; initializes
GAS:: GAS(double *frac0, int MergerModel, double J21, double Tbb, string treefile, bool spec, bool Ma_turn, int bsm){
    //N_sp = 5; N_react = 6; 
    MerMod = MergerModel; printf("MerMod=%d\n",MerMod);
    nMer = 0; iMP = 1;
    MPs = NULL;
    MPs = new MainProgenitor [200];
    aTree(nMer,treefile,MPs); printf("read tree done\n");
    z0 = MPs[iMP].z;
    z = z0;
    z_col = -1.;
    i_bsm = bsm;
    v_bsm = i_bsm*sigma1*(1+z)/(1+z_rcb);
    Mh = MPs[iMP].mhalo;
    HALO halo(Mh,z0);
    rhoc_DM = halo.rho_c;
    // printf("z0 = %3.2f\thalo concentration = %3.2f,\n", z0,halo.c);
    // printf("halo center number density = %3.2f,\n", halo.rho_c/(mu*m_H));
    // printf("halo virial velocity = %3.2f,\n", halo.Vc/km);

    inMer = false;
    evol_stage = (MerMod)?5:0;
    Gamma_mer = 0;
    M_major = 0; M_minor = 0;

    // initial density & T setting.
    T_K0 = halo.Tvir;
    nH0 = MPs[iMP].ng_adb;

    /* if (MerMod==0){
        nH0 = 4.5e-3;
        T_K0 = 20;
        nH0 = 6.*pow(T_K0/1000.,1.5); // Visbal 2014 Eq(2)
    } */
    rho0 = (mu*m_H) * nH0;
    e0 = k_B*T_K0/(gamma_adb-1)/(mu*m_H); // in erg/g
    P0 = nH0*k_B*T_K0;
    S0 = k_B*T_K0/pow(nH0,2./3.); 

    t_ff = 1./(Cp * sqrt(rhoc_DM + (mu*m_H)*nH0));
    // if (MerMod==0)  t_ff = 1./C/sqrt(nH0);
    cs = sqrt( gamma_adb*k_B*T_K0/(mu*m_H) );
    Mgas = Mh*fb;
    Mcore = 4.*pi/3.*pow(0.1*halo.Rvir,3)* nH0;
    M_BE = 1.18*sqrt(fb)*pow(cs,4)/(sqrt(P0*pow(G,3)));
    Mg_intg = 0;

    not_adb = false;
    Ma_on = Ma_turn;
    f_Ma = 1.;
    
    t_ff0 = 1./C/sqrt(nH0); // roughly, not including DM density
    t1 = 1.9999*t_ff0; //maximum time of evolution
    // t1 = 1.999999*t_ff0; //to n = 1.e10/cm^3 (when n > 1.e9 /cm^3, 3body starts to act prominently
    printf("Mh0 = %3.2e Ms, z0 = %2.2f, \n", Mh/Ms, z0);
    printf("nH0 = %3.2e, T_K0 = %2.2f, t_ff0 = %3.2eMyr \n", nH0, T_K0, t_ff0/Myr);

    i_m = 0;
    Nt = 5;
    t0 = 0;
    t_act = 0;
    Dt = 0;

    rho0 = (mu*m_H) * nH0;
    e0 = k_B*T_K0/(gamma_adb-1)/(mu*m_H); // in erg/g

    P0 = (gamma_adb-1) * rho0 * e0;
    v_tur2 = 2./3.*e0; //initialize turbulent energy, following ratio Eth/Ek = 3:1

    J_LW = J21; Tb = Tbb;
    y0 = NULL; y1 = NULL; ys = NULL;
    k = NULL; rf = NULL;
    N = N_sp + 1; 
    y0 = new double [N]; y1 = new double [N];
    ys = new double[(Nt+1)*N];
    k = new double [N_react+1]; rf = new double[N_react+1];

    ycool = 0.; yequi = 0.; ycool_crit = 1.e100;
    Jc_pd = 0.; Jc_cd = 0.; Jc_pred = 0.; Jc_pred_max = 0;
    delta_H2_compr_min = 1.e100;
    a=0; b=0; c=0; d=0; e=0;
    n_H2crit = 0;

    intoequi = false;

    Ta = new double [n_ra]; ka = new double [n_ra];
    read_k(n_ra, Ta, ka);

    kra(k[3], T_K0, n_ra, Ta, ka);
    kpd_Hm_H2p(Tb, k[22],k[23], spec); //返回的是kappa 需乘J21 k_pdHm k_pdH2p 没有self-shielding
    k[22] *= J21; k[23] *= J21;
    //printf("CONSTRUCTOR: k_pdH2=%3.2e, k_pdHm=%3.2e, k_pdH2p=%3.2e\n",k[21],k[22],k[23]);

    RJ = cs*t_ff0;//RJ = sqrt( pi*k_B*T_K0/ (G*pow(mu*m_H,2)*nH0) ); !wli: RJ not precise
    printf("cs_0 is %3.2e km/s; R_vir is %3.2e pc, RJ_0 is %3.2e pc \n",cs/1.e5,halo.Rvir/pc, RJ/pc);
    printf("sqrt(2)*cs_0 is %3.2e km/s; halo Vc is %3.2e km/s \n",cs/1.e5*sqrt(2), halo.Vc);
    MJ0 = 4.*pi/3.*rho0*pow(RJ/2.,3);

    for (int i=0; i<N; i++){
        if (i==0){
            y0[i] = 0; y1[i] = 0;
            ys[i] = 0; 
        }
        y0[i] = *(frac0++); y1[i] = y0[i];
        ys[i] = y0[i]; 
    }
    react_coef(k,nH0,y0[1],y0[2],T_K0,J21,Tb);
    react_rat(rf, y0, k, nH0, T_K0);
    printf("J_LW=%5.2f, Tb=%3.2e\n ****INITIALIZE DONE****\n\n",J_LW,Tb);
}


void GAS:: a_react_sol(bool write){
// len(y_it) = 6 --> to match the Ax=b
    double t_base = 1.e-3;
    double dlnt = log(Dt/t_base)/Nt;
    double* y_i0, *y_i1;
    double ts[Nt+1];
// calculate a series of fractions: add to list Nt+1 times
    double dy[N];
    int i=0, isp, jsp;
    ts[0] = 0;

    for (i=0; i<Nt; i++){
        ts[i+1] = t_base*exp((i+1)*dlnt);
        y_i0 = ys + i*N;
        y_i1 = ys + (i+1)*N;
        for (isp=0; isp<N; isp++){
            dy[isp] = y_i0[isp]; y_i1[isp] = y_i0[isp];// checked right ~
            //printf("y_i0[%d]=%3.2e, y_i1=%3.2e, dy=%3.2e\n",isp,y_i0[isp],y_i1[isp],dy[isp]);
        }
        int iter0 = 0;
        //while ( len_v(N, dy) > epE8*len_v(N, y_i1) ){ //original, 对于initial y_H+<=1.e-5 可能太精细 循环出不来
        //for(int i=1;i<3;i++){ //会有较大的起伏偏差
        // y_i1[0]=1;
        while ( len_v(N, dy) > epE5*len_v(N, y_i1) ){ //epE6check过convergence
            SOL_IMPLICIT(dy, y_i0, y_i1, ts[i+1]-ts[i], nH0, T_K0, k,rf, J_LW, Tb); // y_i0 passed but UNCHANGED.
            iter0++;
            //printf("LOOP TIME %d IN REACT_SOL\n",iter0);
        }
    }
    // set  y1, the reaction result at time t1 by initial set of Nt grids (a react sol)
    // y1[0] = 0 already set in constructor
    for (isp=1;isp<N;isp++) y1[isp] = y_i1[isp];
}

void GAS:: add_Nt(int N){
    Nt *= N;
    delete [] ys;
    ys = new double[(Nt+1)*this->N];
    for (int i=0; i<this->N; i++) ys[i] = y0[i];
}

void GAS:: react_sol(bool write){
    double y1_prev[N], delta_y[N];
    for (int i=0; i<N; i++) y1_prev[i] = 0.;

    kra(k[3],T_K0,n_ra,Ta,ka);
    react_coef(k,nH0,y0[1],y0[2],T_K0,J_LW,Tb);
    react_rat(rf,y0,k,nH0,T_K0);
    do{ 
        //add_Nt(10); 全程*10 essencially the same with *2
        add_Nt(2);
        //printf("in REACT_SOL: Nt=%d\n",Nt);
        a_react_sol(0);  
        for (int i=0; i<N; i++){
            delta_y[i] = y1[i] - y1_prev[i];
            y1_prev[i] = y1[i];
        }
    }while ( len_v(N, delta_y) > epE5*len_v(N, y1_prev) ); //epE6check过convergence

//charge neutrality & H neuclei conservation
    double y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp;
    y_H    = y1[1];
    y_H2   = y1[2];
    y_e    = y1[3];
    y_Hp   = y1[4];
    y_H2p  = y1[5];
    y_Hm   = y1[6];
    y_He   = y1[7];
    y_Hep  = y1[8];
    y_Hepp = y1[9];
    
    if(y_H<y_Hp) y_Hp = 1.0-y_H -2.*y_H2-2.*y_H2p-y_Hm;
    if(y_H>y_Hp) y_H  = 1.0-y_Hp-2.*y_H2-2.*y_H2p-y_Hm;

    double yHe = 8.3333333e-2; //number fraction of He neuclei;
    if(y_He>y_Hep) {
        if(y_He>y_Hepp) y_He = yHe - y_Hep - y_Hepp;
    }
    if(y_Hep>y_He) {
        if(y_Hep>y_Hepp) y_Hep = yHe - y_He - y_Hepp;
    }
    if(y_Hepp>y_He) {
        if(y_Hepp>y_Hep) y_Hepp = yHe - y_He - y_Hep;
    }

    y_e = y_Hp + y_Hep + 2.0*y_Hepp + y_H2p - y_Hm;

// update y0; chemical reaction for Dt done.
    y0[1] = y_H;
    y0[2] = y_H2;
    y0[3] = y_e;
    y0[4] = y_Hp;
    y0[5] = y_H2p;
    y0[6] = y_Hm;
    y0[7] = y_He;
    y0[8] = y_Hep;
    y0[9] = y_Hepp;
    
    /* y0[1] = y1[1];
    y0[2] = y1[2];
    y0[3] = y1[3];
    y0[4] = y1[4];
    y0[5] = y1[5];
    y0[6] = y1[6];
    y0[7] = y1[7];
    y0[8] = y1[8];
    y0[9] = y1[9]; */

    /* for(int isp=1;isp<N_sp1;isp++) printf("in REACT_SOL, IMPLICIT: y1[%d]=%3.2e\n",isp,y0[isp]);
    printf("\n"); */

//reset Nt
    Nt = 5;

    double kHm_e, kH2_Hm,   kH2p_Hp, kH2_H2p,   kH2_3b;
    double kHm_pd, kH2p_pd, kH2_pd;
    double kH2_cd_tot, kH2_cd_H, kH2_cd_e, kH2_cd_H2, kH2_cd_Hp, kH2_cd_He; 
    double kform_Hm, kform_H2p;
    double kHm_rcb, kHm_cd;
    double rform;

//
//  1)   H     +   e     ->   H+    +  2 e    col ionization
//  2)   H+    +   e     ->   H     +   ph.   recombination

//  3)   H     +   e     ->   H-    +   ph.   kHm_e
//  4)   H-    +   H     ->   H2    +   e     kH2_Hm
//  11)  H-    +   H+    -> 2 H               kHm_rcb

//  5)   H     +   H+    ->   H2+   +   ph.   kH2p_Hp
//  6)   H2+   +   H     ->   H2    +   H+    kH2_H2p

//  15) 3 H              ->   H2    +   H     kH2_3b

//  21)  H2   +   ph   -> 2 H                 kH2_pd
//  22)  H-   +   ph  ->  H   +  e            kHm_pd
//  23)  H2+  +   ph  ->  H   +  H+           kH2p_pd

//  7)   H2    +   H     -> 3 H               kH2_cd_H
//  8)   H2    +   H+    ->   H2+   +   H     kH2_cd_Hp
//  9)   H2    +   e     -> 2 H     +   e     kH2_cd_e
//  20)  H2    +   e     ->   H-    +   H
//  16)  2 H2            -> 2 H     +   H2    kH2_cd_H2
//  35)  H2    +   He    ->   He    +  2 H    kH2_cd_He

    kHm_e = k[3]; kH2_Hm = k[4];
    kH2p_Hp = k[5]; kH2_H2p = k[6];

    kH2_pd = k[21]; kHm_pd = k[22]; kH2p_pd = k[23];
    kH2_3b = k[15];
    kH2_cd_H = k[7]; kH2_cd_Hp = k[8]; kH2_cd_e = k[9] + k[20]; kH2_cd_H2 = k[16];
    kH2_cd_He = k[35];
    kH2_cd_tot = nH0* (kH2_cd_H*y0[1] + kH2_cd_Hp*y0[4] + kH2_cd_e*y0[3] + kH2_cd_H2*y0[2] + kH2_cd_He*y0[7]); 
    // kH2_cd_tot = nH0* (kH2_cd_H*y0[1]); // collisional dissociation dominated by H

    // regualte H2 at high n, together w/ H2 cd 
    kHm_rcb = k[11];
    kHm_cd = k[18];

    kform_Hm = kHm_e* kH2_Hm*y0[1]/( kH2_Hm*y0[1] + kHm_pd/nH0 + kHm_rcb*y0[4] + kHm_cd*y0[1]);
    kform_H2p = kH2p_Hp* kH2_H2p/(kH2_H2p + kH2p_pd/nH0);
    // H2 forming rate; 
    rform = max(kform_Hm*nH0*y0[1]*y0[3]+kform_H2p*nH0*y0[1]*y0[4], kH2_3b*pow(nH0,2)*pow(y0[1],3));
                // H-, small n              // H2+                            3b,  n>~10^9/cm^3
    rform = kform_Hm*nH0*y0[1]*y0[3];// + kform_H2p*nH0*y0[1]*y0[4]; // H2 formation dominated by H- over H2+
    // yequi = min(rform/kH2_pd, rform/kH2_cd_tot);
        // pd,  n<~100/cm^3   cd, large n 
    yequi = rform/(kH2_pd+kH2_cd_tot);
    // yequi = rform/(kH2_pd);
    // yequi = rform/(kH2_cd_tot);


    // sufficient cooling fraction. Λ_H2 (erg/cm^3/s) v.s. Γ_compr(g_Ma) (erg/g/s) 
    ycool = g_Ma* rho0*Gamma_compr(cs,1.,t_ff) * y0[2]/Lambda_H2(nH0,T_K0,y0);

    // ycool = yequi;
    Jc_pd = rform / (ycool * kH2_pd/J_LW);

    double beta_Hm = kHm_pd/J_LW; double beta_H2 = kH2_pd/J_LW;

    a = y0[1]*y0[3]*nH0/kH2_cd_tot*kHm_e*kH2_Hm*y0[1];
    b = kH2_Hm*y0[1] + kHm_rcb*y0[4] + kHm_cd*y0[1];
    c = beta_Hm/nH0;
    Jc_cd = (a/ycool - b)/c; 

    delta_H2_compr = abs(r_cH2 - Gamma_compr(cs,f_Ma,t_ff))/r_cH2 ;

    a = nH0*y0[1]*y0[3]*kHm_e*kH2_Hm*y0[1];
    b = kH2_Hm*y0[1] + kHm_rcb*y0[4] + kHm_cd*y0[1];
    c = beta_Hm/nH0;
    d = kH2_cd_tot*y0[1];
    e = beta_H2;
    ycool_crit = ycool;
    Jc_pred = (-(c*d+b*e)+sqrt( pow(c*d+b*e,2)-4.*(c*e)*(b*d-a/ycool_crit)))/(2.*c*e);

    if (Jc_pred> Jc_pred_max){
        n_H2crit = nH0;
        z_H2crit = z0;
        gMa_H2crit = g_Ma;
        fMa_H2crit = f_Ma;
        ycool_crit = ycool;
        Jc_pred_max = Jc_pred;
    }

}

void GAS:: setMerger(){
    //printf("iMP = %d, nMer = %d", iMP, nMer);
    if (MerMod !=  0){
        if (iMP< nMer){ // final merger not included, since nMer halos only nMer-1 intervals, //wli
        // from iMP=1 to nMer-1, MPs[iMP] has dt, dm, mratio, etc.
            inMer = true;
            if ( t_act >= MPs[iMP+1].t){
                iMP ++;
                if (MPs[iMP].major) M_major += MPs[iMP].dm;
            }
            if (MPs[iMP].t <= t_act and t_act < MPs[iMP+1].t){ 
                // interpolation btw mergers
                Mh = (MPs[iMP].mhalo*(MPs[iMP+1].t-t_act) + MPs[iMP+1].mhalo*(t_act-MPs[iMP].t)) / MPs[iMP].dt; 
            }
            else cout<<"\n!\n!\n!\n! wrong in setMer\n";
            dMdt = MPs[iMP].dm / MPs[iMP].dt;
            Mgas = Mh*fb;
            
            HALO halo(Mh, z); //更真实的Mh 和 z
            dEdt = k_B*halo.Tvir /(mu*m_H) * (dMdt*fb); // estimated by Tvir
            // from Virial theorm, thermal energy change according to Phi change
            /* -> thermal from Virial theorm */ // ( 3/4 of Gamma_mer to thermal get Tg~Tvir)
            Gamma_mer =  -1./2.*halo.Phi(halo.Rvir) * (2./3.* dMdt/Mh - Hz(z0) ); 
            if (Gamma_mer<0) cout<<"##COOLING##"<<2./3.* dMdt/Mh - Hz(z0)<<endl;
            //printf("Tg/Tvir = %3.2f\n",T_K0/halo.Tvir);
            
        }
        else {inMer = false; Gamma_mer = 0;}
    }
}

void GAS:: timescales(){
    //不行 因为t_chem太小无法演化下去
    // chemical reaction timescale; use minimum to avoid too absurd fowarding
    t_chem = 1.e100;
    int isp, ifast;
    for (isp=1; isp<7;isp++) { //不管He
    //for (isp=1; isp<N_sp1;isp++) {
        if (rf[isp]!=0. and y0[isp]!=0) {
            t_chem = min( abs(y0[isp]/rf[isp]), t_chem );
            ifast = isp;
        }
    }
    // printf("ifast=%d\n",ifast);
    //if (nH0<1.5e3 and nH0>900) printf("t_ch[%d]=%3.2e t_ff0\n",t_chem, y0[ifast]/rf[ifast]/t_ff0);

    //Hm + H -> H2+e  &  3H -> H2+H  
    double y_H = y0[1], y_H2 = y0[2], y_e = y0[3], y_Hep = y0[8], y_He = y0[7];
    double k_Hion = k[1], k_Heion = k[31];
 
    t_ion = 1/(k_Hion*y_H*nH0);
    double alpha_B = k[2]; ////  2)   H     +   e     ->   H-    +   ph.
    t_rcb = 1./(y_e*alpha_B*nH0); // recombination timescale = (ye_0*alpha_B*nH)^-1

    r_cH = Lambda_H(nH0,T_K0,y_H,y_e, k_Hion)/rho0;
    r_cH2 = Lambda_H2(nH0,T_K0,y0)/rho0;
    r_c =  r_cH2 + r_cH + Lambda_Hep(nH0, T_K0, y_Hep, y_e, y_He, k_Heion)/rho0;
    //开关 turn_off H2 cooling wli!
    // r_c -= r_cH2;
//// Merger heating -> turbulent & thermal
    Gamma_mer_th = fraction * Gamma_mer;
    Gamma_mer_k = (1-fraction) * Gamma_mer;
// merger case (evol_stage=4 is freefall) !wli!!!
    //if (evol_stage ==4) r_h = Gamma_compr(cs,1,t_ff) + Gamma_chem(nH0, T_K0, y0, k)/rho0;
    if (evol_stage ==4) r_h = Gamma_compr(cs,f_Ma,t_ff) + Gamma_mer_th + Gamma_chem(nH0, T_K0, y0, k)/rho0;
    else r_h = Gamma_mer_th + Gamma_chem(nH0, T_K0, y0, k)/rho0; //no compressional

    //printf("f_Ma=%3.2e, cs=%3.2e, t_ff=%3.2e, compr=%3.2e mer=%3.2e mer_th=%3.2e\n", f_Ma, cs, t_ff, Gamma_compr(cs,f_Ma,t_ff),Gamma_mer,Gamma_mer_th);

    t_c = e0/r_c;
    t_h = abs(e0/r_h);
    HALO halo(Mh,z);
    rhoc_DM = halo.rho_c;
    t_ff = 1./(Cp * sqrt(rhoc_DM + (mu*m_H)*nH0));

    // printf("t_ff1=%3.2e\t t_ff2=%3.2e\n",t_ff/Myr, 1./C/sqrt(nH0)/Myr);

    // Dt firstly by cooling/heating/free-fall timescales; also all for NO merger case
    Dt = 0.1* min( min(t_c,t_h), t_ff ); //不行不够细
    Dt = 0.01*min( min(t_ff,100*t_chem),min(t_c,t_h));

// merger case: Dt << merger intervals dt
    if (inMer and evol_stage !=4) Dt = min( Dt, .1*MPs[iMP].dt ); // wli : check 能否放宽
    
// no merger case, just free fall; one-zone case of former work
    if (MerMod==0) { 
        r_h = Gamma_compr(cs,f_Ma,t_ff) + Gamma_chem(nH0, T_K0, y0, k)/rho0;
        t_c = e0/r_c;
        t_h = abs(e0/r_h); //可能为负 potential well dilution
        //不需要t_rcb //如果加细Dt 用0.001仍然converge而且似乎更smooth但是慢得多 //不能放宽了 0.1明显不行
        Dt = 0.01*min( min(t_ff,100*t_chem),min(t_c,t_h));

        //printf("in TIMESCALES:t_c=%3.2e,t_h=%3.2e,t_ff=%3.2e,t_chem=%3.2e\n",t_c,t_h,t_ff,t_chem);
    }
    // time++ by Dt
    t_act += Dt;
    z0 = z;
    z = z_ana(z,Dt); //推进 redshift
    v_bsm = i_bsm*sigma1*(1.+z)/(1.+z_rcb);
}

// cloud evolution: freefall() -> react_sol() -> T_sol

void GAS:: freefall(){  //module of explicit integration over Dt
    HALO halo1(Mh,z0);
    bool adjust_iso = false;
    double nvir_max, nvir;
    double Vc0 = 3.7*km;
    double alpha = 4.7;
    switch(evol_stage){
        case 0: // not combined w/ mergers
            nH0 = n_ff(z0,nH0,0.,Dt); //设成0 为了使其不受rhoc_DM的redshift dependence 影响
            break;
        case 1: // adiabatic heating, n ~ T^1.5, constant entropy
            nH0 = N_ADB(S0,T_K0);
        // compare with maximum core density by core entropy: Visbal et al. 2014a
            if (N_ADB(S0,T_K0)>N_CORE(z0)) evol_stage = 2;
        // 新加: H2 cooling dominant--> collapse
            if (r_cH2 >= abs(r_h) ){ // not adiabatic
                //if (pow(halo1.Vc,2) >= pow(Vc0,2) + pow(alpha*v_bsm,2)){
                if (pow(halo1.Vc,2) >= pow(cs,2)*f_Ma + pow(alpha*v_bsm,2)){
                    printf("!from 1 to 4!\nz = %3.2f, Tvir = %3.2eK, Mh = %3.2eMs\n nH0=%3.2e Tg=%3.2eK\n", z0, halo1.Tvir, Mh/Ms,fb*Mh/Ms, T_K0);
                    evol_stage = 4; //wli:暂且用作H2 collapse
                }
            }
            break;
        case 2: // saturation following universe expansion; also adiabatic?
            nH0 = N_CORE(z0); 
            //printf("saturation");

            if (r_cH >= abs(r_h) or r_cH2>= abs(r_h)) not_adb = true;
            // 别的判据? <8000K 曾经加上以免 density drop; but it drops anyway
            if ( r_cH >= abs(r_h) and T_K0<8000){ //if Lya cooling dominated and just cool down < 8000K
            // if ( max(r_cH,r_cH2) >= abs(r_h) and T_K0<8000){ //if Lya cooling dominated and just cool down < 8000K
                cout<<"/////////////\t t_c<t_h & T_K<8000 \tFIRST ENTER ISO STAGE////////////////\n";
                printf("z = %3.2f\t Tvir = %3.2eK\t Mh = %3.2eMs\t Mgas = %3.2e\t Tg=%3.2eK\n", z0, halo1.Tvir, Mh/Ms,fb*Mh/Ms, T_K0);
                cout<<"ng_adb= "<<nH0<<endl;
                evol_stage = 3;
                //inside solving n_iso, n_vir comparing w/ cosmic mean determine if unstable
                Nvir2N0(n_iso, nvir_max, nH0, f_Ma*T_K0, z0, Mh);
                printf("//////////\tSOLVED FOR THE 1ST TIME\n");
                printf("n_iso=%3.2e, f_Ma=%3.2e, v_bsm=%3.2e, Vc=%3.2e\n", n_iso,f_Ma,v_bsm,halo1.Vc);
                if (!n_iso) evol_stage = 4; // unstable case Mg2ng return 0
                else nH0 = n_iso;
                Mh_prev = Mh; t_prev = t_act;
            }
        // 新加: H2 cooling dominant--> collapse
            if (r_cH2 >= abs(r_h) ){
                // if (pow(halo1.Vc,2) >= pow(Vc0,2) + pow(alpha*v_bsm,2)){
                if (pow(halo1.Vc,2) >= pow(cs,2)*f_Ma + pow(alpha*v_bsm,2)){
                    printf("!from 2 to 4!\nz = %3.2f, Tvir = %3.2eK, Mh = %3.2eMs\n nH0=%3.2e Tg=%3.2eK\n", z0, halo1.Tvir, Mh/Ms,fb*Mh/Ms, T_K0);
                    evol_stage = 4; //wli:暂且用作H2 collapse
                }
            }
            break;
        case 3: // iso T=8000K (H Lya cooling)
            adjust_iso = (t_act - t_prev >= min( 0.1*MPs[iMP].dt, t_ff) );
            if (adjust_iso) {
                dt_iso = t_act - t_prev;
                Mh_prev = Mh; t_prev = t_act;
                Nvir2N0(n_iso, nvir_max, nH0, f_Ma*T_K0, z0, Mh);
                if (!n_iso) evol_stage = 4; // unstable case Mg2ng return 0
                else {
                    nH0 = n_iso;
                    BOUNDARY(nvir,Mg_intg,f_Ma*T_K0,nH0*(mu*m_H)/halo1.rho_c,z0,Mh);
                }
            }
            break;
        case 4:
            Mg_intg = 0; // gas mass within Rvir
            if (z_col<0.) z_col = z0; // mark collapse redshift
            nH0 = n_ff(z0,nH0,rhoc_DM,Dt);
            break;
        case 5:
            nH0 = MPs[iMP].ng_adb;

            if (not intoequi and y0[2]/yequi<2.) {
                intoequi = true; 
                printf("into equi at %3.2e t_ff0, nH=%3.2e /cc\n",t_act/t_ff0,nH0);
            }

            if (intoequi and r_cH2>abs(r_h) or r_cH>abs(r_h)) { //factor 2 for H2 cooling deleted
                printf(" in CASE5: adb not hold, go to 3\t");
                evol_stage = 3;
                Nvir2N0(n_iso, nvir_max, nH0, f_Ma*T_K0, z0, Mh);
                printf("SOLVED FOR THE 1ST TIME\n");
                // printf("n_iso=%3.2e, f_Ma=%3.2e, v_bsm=%3.2e, Vc=%3.2e\n", n_iso,f_Ma,v_bsm,halo1.Vc);
                if (!n_iso) evol_stage = 4; // unstable case Mg2ng return 0
                else nH0 = n_iso;
                Mh_prev = Mh; t_prev = t_act;
            }
            break;
    }


    v_tur2 += Dt * Gamma_mer_k*2; // 2 coz e=1/2v^2

    double b = 1./3.; // [1/3, 0.5)
    f_Ma = 1;
    // turbulence
    if (Ma_on) f_Ma += v_tur2/pow(cs,2) * gamma_adb*b; // corrected f_Ma, using P = rho_g v_tur^2/3 from Chandrasekhar 1951
    // bsm velocity: alpha = 4.7; // [4, 8)
    if (i_bsm) f_Ma += pow(alpha*v_bsm/cs,2); //Eq(3) in Hirano+2018. from Fialkov2012
    // f_Ma = 3.; //调试用 1 和 3

    Ma = sqrt(v_tur2)/cs;
    // if (evol_stage==3) {
    //     printf("f_Ma=%3.2e v_bsm=%3.2e, Vc=%3.2e\n",f_Ma, v_bsm/km, cs/km, halo1.Vc/km);
    //     printf("Vc^2/ (cs^2 + v_tur^2 + v_bsm^2)=%3.2e\n",pow(halo1.Vc,2)/(pow(cs,2)+v_tur2+pow(v_bsm,2)));
    //     printf("concentration c=%3.2e\n",halo1.Rs/halo1.Rvir);
    // }

    if (evol_stage==4) g_Ma = f_Ma + Gamma_mer_th/Gamma_compr(cs,1.,t_ff);
    else g_Ma = 1.;

 //for MerMod=0 case
    if (MerMod==0) {
        f_Ma = 1;
        g_Ma = 1;
    }
    // update rho
    rho0 = (mu*m_H) * nH0;
}

void GAS:: T_sol(){
    //开关 turn_off cooling
    //if (MerMod != 0) r_c = 0;

    e0 += (r_h-r_c)*Dt;
    P0 = (gamma_adb-1) * rho0 * e0;
    T_K0 = e0*(gamma_adb-1)*(mu*m_H)/k_B;
}

void GAS:: get_para(){
    cs = sqrt(gamma_adb* k_B*T_K0/(mu*m_H) );
    RJ = sqrt( pi*k_B*T_K0/ (G*pow(mu*m_H,2)*nH0) );
    M_BE = 1.18*sqrt(fb)*pow(cs,4)/(sqrt(P0*pow(G,3)));
    MJ = 4.*pi/3.*rho0*pow(RJ/2.,3);
    MJ_eff = MJ*pow(f_Ma,1.5);
}


GAS:: ~GAS(void){
    delete [] y0; delete [] y1; delete [] ys;
    delete [] k; delete [] rf;
    delete [] Ta; delete [] ka;
    delete [] MPs;
    // file_ingas.close();
    // cout<<"releasing gas\n";
}