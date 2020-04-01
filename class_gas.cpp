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
GAS:: GAS(double *frac0, int MergerModel, double J21, double Tbb, char* treefile, bool Ma_turn){
    //N_sp = 5; N_react = 6; 
    nMer = 0;
    MPs = NULL;
    MPs = new MainProgenitor [100];
    aTree(nMer,MPs,treefile);
    z0 = MPs[0].z;
    z = z0;
    Mh = MPs[0].mhalo;

    HALO halo(Mh,z0);

    rhoc_DM = halo.rho_c;
    printf("z0 = %3.2f\thalo concentration = %3.2f,\n", z0,halo.c);
    printf("halo center number density = %3.2f,\n", halo.rho_c/(mu*m_H));
    printf("halo virial velocity = %3.2f,\n", halo.Vc/km);

    // initial density & T setting.
    nH0 = halo.rho_c*fb/(mu*m_H);
    T_K0 = halo.Tvir;
    nH0 = 6.*pow(T_K0/1000.,1.5); // Visbal 2014 Eq(2)
    printf("MerMod=%d\n",MerMod);
    /* if (MergerModel==0){
        nH0 = 4.5e-3;
        T_K0 = 20;
    } */
    rho0 = (mu*m_H) * nH0;
    e0 = k_B*T_K0/(gamma_adb-1)/(mu*m_H); // in erg/g

    P0 = nH0*k_B*T_K0;
    S0 = k_B*T_K0/pow(nH0,2./3.);

    t_ff = 1./(Cp * sqrt(rhoc_DM + (mu*m_H)*nH0));
    t_delay = t_ff; // ok for no merger

    cs = sqrt( gamma_adb*k_B*T_K0/(mu*m_H) );
    Mgas = Mh*fb;
    Mcore = 4.*pi/3.*pow(0.1*halo.Rvir,3)* nH0;
    M_BE = 1.18*sqrt(fb)*pow(cs,4)/(sqrt(P0*pow(G,3)));
    
    not_adb = false;
    Ma_on = Ma_turn;
    de_tot = 0;
    dM_tot = 0;
    f_Ma = 1.;
    printf("Mh0 = %3.2e Ms, z0 = %2.2f, \n", Mh/Ms, z0);
    printf("nH0 = %3.2e, T_K0 = %2.2f, \n", nH0, T_K0);
    t_ff0 = 1./C/sqrt(nH0); // roughly, not including DM density

    t1 = 1.9999*t_ff0; //maximum time of evolution
    // t1 = 1.999999*t_ff0; //to n = 1.e10/cm^3 (when n > 1.e9 /cm^3, 3body starts to act prominently

    i_m = 0;
    Nt = 5;
    t0 = 0;
    t_act = 0;

    rho0 = (mu*m_H) * nH0;
    e0 = k_B*T_K0/(gamma_adb-1)/(mu*m_H); // in erg/g

    P0 = (gamma_adb-1) * rho0 * e0;
    v_tur2 = 1./2.*e0; //initialize turbulent energy, following ratio Eth/Ek = 3:1
    //v_tur2 = 2./3.*e0; //initialize turbulent energy, following ratio Eth/Ek = 3:1
    //v_tur2 = 0; //initialize turbulent energy, following ratio Eth/Ek = 3:1

    J_LW = J21; Tb = Tbb;
    y0 = NULL; y1 = NULL; ys = NULL;
    k = NULL; rf = NULL;
    N = N_sp + 1; 
    y0 = new double [N]; y1 = new double [N];
    ys = new double[(Nt+1)*N];
    k = new double [N_react+1]; rf = new double[N_react+1];

    kpd_Hm_H2p(Tb, k[22],k[23]); //返回的是kappa 需乘J21 k_pdHm k_pdH2p 没有self-shielding
    k[22] *= J21; k[23] *= J21;
    printf("CONSTRUCTOR: k_pdH2=%3.2e, k_pdHm=%3.2e, k_pdH2p=%3.2e\n",k[21],k[22],k[23]);
    cs = sqrt( gamma_adb*k_B*T_K0/(mu*m_H) );
    RJ = cs*t_ff0;//RJ = sqrt( pi*k_B*T_K0/ (G*pow(mu*m_H,2)*nH0) );
    printf("cs_0 is %3.2e km/s; R_vir is %3.2e pc, RJ_0 is %3.2e pc \n",cs/1.e5,halo.Rvir/pc, RJ/pc);
    printf("sqrt(2)*cs_0 is %3.2e km/s; halo Vc is %3.2e km/s \n",cs/1.e5*sqrt(2), halo.Vc);

    MJ0 = 4.*pi/3.*rho0*pow(RJ/2.,3);
    M_major = 0;
    M_minor = 0;

    MerMod = MergerModel;
    inMer = false;
    inDelay = false;
    evol_stage = (MerMod)?1:0;
    Gamma_mer = 0;
    iMer = 0;
    file_ingas.open("data/mer_record.txt", ios::out | ios::trunc); //
    if (MerMod != 0) file_ingas<<"z t_halo t_act Mh Tvir Rvir dt t_dyn n_gas Mcore nc_DM M_BE MJ"<<endl;

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
    printf("J_LW=%3.2e, Tb=%3.2e\n",J_LW,Tb);
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
        while ( len_v(N, dy) > epE5*len_v(N, y_i1) ){
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
    }while ( len_v(N, delta_y) > epE6*len_v(N, y1_prev) );

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
}

void GAS:: setMerger(){
    //printf("iMer = %d, nMer = %d", iMer, nMer);
    ofstream fMer;
    fMer.open("data/mergers.txt", ios::out | ios::app); // append mode, to use it, mc then mk
    if (iMer == 0) fMer<<"j t t_act mhalo/Ms dm major Tvir"<<endl;

    if (MerMod !=  0){
        if (iMer< nMer-1){ // final merger not included, since nMer halos only nMer-1 intervals, //WLI
        // from iMer=0 to nMer-2, MPs[iMer] has dt, dm, mratio, etc.
            inMer = true;
            //printf("MPs[iMer].t = %3.2e\nt_act = %3.2e\nMPs[iMer+1].t = %3.2e\n",MPs[iMer].t,t_act,MPs[iMer+1].t);     
            HALO halo(Mh, z); //更真实的Mh 和 z
            // in a merger
            if (MPs[iMer].t <= t_act && t_act < MPs[iMer+1].t){
                dMdt = MPs[iMer].dm / MPs[iMer].dt;
                Mh = MPs[iMer].mhalo + dM_tot; //和merger tree 的Mh v.s. z check过, interpolation btw mergers
                Mgas = Mh*fb;
                Mcore = 4.*pi/3.*pow(0.1*halo.Rvir,3)* nH0*(mu*m_H);

                //dEdt = pow(halo.Vc,2) * (dMdt*fb); // estimated by Vc
                dEdt = k_B*halo.Tvir /(mu*m_H) * (dMdt*fb); // estimated by Tvir

            // from Virial theorm, thermal energy change according to Phi change
                /* -> thermal from Virial theorm*/ // ( 3/4 of Gamma_mer to thermal get Tg~Tvir)
                Gamma_mer =  -1./2.*halo.Phi(halo.Rvir) * (2./3.* dMdt/Mh - Hz(z0) ); 

                v_inf2 = halo.V_inf2(2*halo.Rvir, halo.Rvir/3);
                v_inf2 = (MPs[iMer].major)?halo.V_inf2(2*halo.Rvir, halo.Rs):halo.V_inf2(2*halo.Rvir,halo.Rvir);
                //v_inf2 = (MPs[iMer].major)? -2*halo.Phi(halo.Rvir/2):-2*halo.Phi(halo.Rvir);

                double alpha = v_inf2/pow(halo.Vc,2);double beta = T_K0/halo.Tvir;
                //printf("3/8*alpha^2 = %3.2f, beta=%3.2f, 3/8*alpha^2 - beta = %3.2f\n",3./8.*alpha, beta, 3./8.*alpha - beta);
                if (Gamma_mer<0) cout<<"##COOLING##"<<2./3.* dMdt/Mh - Hz(z0)<<endl;
                //printf("Tg/Tvir = %3.2f\n",T_K0/halo.Tvir);
            }
            // time evolve to next merger, update iMer, Mh, Mgas
            if ( t_act >= MPs[iMer+1].t){

                //cout<<"i_Mer = "<<iMer<<"\t f_Ma = "<<f_Ma<<"\t Mh = "<<Mh<<endl;
                iMer ++;
                de_tot = 0;
                dM_tot = 0;

                Mh = MPs[iMer].mhalo;
                HALO halo2(Mh,z0);
                if (MPs[iMer].major) M_major += MPs[iMer].dm;

                //printf("in Ms: Mgas = %3.2e \tMJ = %3.2e \tMJ_eff = %3.2e\n", Mgas/Ms, MJ/Ms, MJ_eff/Ms);
                fMer<<MPs[iMer].j<<" "<<MPs[iMer].t<<" "<<t_act<<" "<<MPs[iMer].mhalo/Ms<<" "<<MPs[iMer].dm/Ms<<" "<<MPs[iMer].major<<" "<<MPs[iMer].Tvir<<endl;
                //if (iMer == 0) file_ingas<<"z t_halo t_act Mh Tvir Rvir dt t_dyn n_gas Mcore nc_DM M_BE MJ"<<endl; //halo2.rho_c/(mu*m_H)RHO_crit(z0)
                file_ingas<<z<<" "<<MPs[iMer].t<<" "<<t_act<<" "<<Mh<<" "<<MPs[iMer].Tvir<<" "<<halo2.Rvir<<" "<<MPs[iMer].dt<<" "<<halo2.t_dyn<<" "<<nH0<<" "<<Mcore<<" "<<halo2.rho_crit/(mu*m_H)<<" "<<M_BE<<" "<<MJ<<endl;
                // dilution(destruction/reduction/destroyment) of H2 fraction
                // 1) by 1/2
                /* for (int i=2; i<N; i++) y0[i]/=2.;
                y0[1] = 1 - 2*y0[2] - y0[4] - y0[5];
                 */
                // 2) all destructed
                /* y0[1] += 2*y0[2]; 
                y0[2] = 0; */

            }
        }
        else {inMer = false; Gamma_mer = 0;}
    }
    fMer.close();
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
    //if (nH0<1.5e3 and nH0>900) printf("t_ch[%d]=%3.2e t_ff0\n",t_chem, y0[ifast]/rf[ifast]/t_ff0);

    //Hm + H -> H2+e  &  3H -> H2+H  
    double y_H = y0[1], y_H2 = y0[2], y_e = y0[3], y_Hep = y0[8], y_He = y0[7];
    double k_Hion = k[1], k_Heion = k[31];
 
    t_ion = 1/(k_Hion*y_H*nH0);
    double alpha_B = k[2]; ////  2)   H     +   e     ->   H-    +   ph.
    t_rcb = 1./(y_e*alpha_B*nH0); // recombination timescale = (ye_0*alpha_B*nH)^-1
    /* if (nH0<1.5e3 and nH0>900) printf("t_chem=%3.2e,t_rcb=%3.2e,t_ion=%3.2e t_ff0\n",
       t_chem/t_ff0, t_rcb/t_ff0, t_ion/t_ff0);   */
    /* if (t_rcb>t_ion) printf("t_rcb>t_ion: %3.2e,%3.2e t_chem=%3.2e,nH0=%3.2e\n",
        t_rcb/t_ff0, t_ion/t_ff0,t_chem/t_ff0,nH0); */  

    r_cH = Lambda_H(nH0,T_K0,y_H,y_e, k_Hion)/rho0;
    r_cH2 = Lambda_H2(nH0,T_K0,y0)/rho0;
    r_cH2 = 0;
    r_c =  r_cH2 + r_cH + Lambda_Hep(nH0, T_K0, y_Hep, y_e, y_He, k_Heion)/rho0;
    //开关 turn_off cooling
    //r_c = 0;
//// Merger heating -> turbulent & thermal
    Gamma_mer_th = fraction * Gamma_mer;
    Gamma_mer_k = (1-fraction) * Gamma_mer;

    if (evol_stage < 4) r_h = Gamma_mer_th + Gamma_chem(nH0, T_K0, y0, k)/rho0;
    else r_h = Gamma_compr(cs,f_Ma,t_ff) + Gamma_mer_th + Gamma_chem(nH0, T_K0, y0, k)/rho0;

    //printf("f_Ma=%3.2e, cs=%3.2e, t_ff=%3.2e, compr=%3.2e mer=%3.2e mer_th=%3.2e\n", f_Ma, cs, t_ff, Gamma_compr(cs,f_Ma,t_ff),Gamma_mer,Gamma_mer_th);

    t_c = e0/r_c;
    t_h = e0/r_h;
    HALO halo(Mh,z);
    rhoc_DM = halo.rho_c;
    t_ff = 1./(Cp * sqrt(rhoc_DM + (mu*m_H)*nH0));

    t_delay = t_ff; //delay of evolving time after merger NOT INCLUDED NOW.
    t_delay = 0;
    t_h = (t_h>0)?t_h:-t_h; //
    // Dt firstly by cooling/heating/free-fall timescales; also all for NO merger case
    Dt = 0.1* min( min(t_c,t_h), t_ff ); 
    //Dt = 0.1*min(t_h,t_ff); //!!!!! WLI appended here when freefall off, if including t_c can be too short and never evolve...

// merger case
    if (MerMod !=0){
        inDelay = false;
        if (iMer<nMer-1){
            // Dt << merger intervals dt
            Dt = min( Dt, .001*MPs[iMer].dt ); // wli : check 能否放宽
            // during delay, no compressional heating, only for MAJOR MERGERS...
            //if ( MPs[iMer].t<=t_act && t_act<MPs[iMer].t+t_delay ){
            if ( MPs[iMer].major && MPs[iMer].t<=t_act && t_act<MPs[iMer].t+t_delay ){
                inDelay = true;
            }
            // after delay time (after a MAJOR MERGER...), gas evolve time t0+=Dt 
            //if ( t_act>MPs[iMer].t+t_delay && t_act<MPs[iMer+1].t){
            if ( MPs[iMer].major && t_act>MPs[iMer].t+t_delay && t_act<MPs[iMer+1].t){
                t0 += Dt;
            }
        }
    }
// no merger case, just free fall; one-zone case of former work
    if (MerMod==0) { 
        r_h = Gamma_compr(cs,f_Ma,t_ff) + Gamma_chem(nH0, T_K0, y0, k)/rho0;
        t_c = e0/r_c;
        t_h = abs(e0/r_h); //可能为负 potential well dilution
        Dt = 0.1* min( min(t_c,t_h), t_ff ); //sufficiently same with Dt = 0.01*... 
                                             //WRONG for y_e<1.e-5, could use 0.01, 计算耗时多于带着t_chem判断
        //加入t_chem, t_rcb决定Dt 权重1, 0.0001经过调试 
        Dt = 0.1*min( min(t_ff,100*t_chem),min(t_c,t_h));
        //if (nH0>60 and nH0<1.e4) Dt = min(0.01*t_rcb,Dt);
        /* if (t_ion<100*t_ff and nH0<1.e4) {
            Dt = min(Dt, 0.001*t_rcb);
            printf("nH0=%3.2e\n",nH0);
        } */
        //if (nH0>60 and nH0<1.e4) Dt = min(0.0001*t_ion,Dt);//不行 -4结果错误rcb被延迟

        //if (nH0>150. and nH0<400.) Dt = min(0.0001*t_rcb,Dt);//这样放宽加细范围不行
        //if (nH0>60 and nH0<1.e4) Dt = min(0.001*t_rcb,Dt);// 对于-4不行 不够细 后续大的y_e

        // 10*t_chem, 0.001*t_rcb: 下面的对于y_e=1.e-3, -4, -5都可以 结果-4 和3 5不同 也影响Jc
        Dt = 0.1*min( min(t_ff,100*t_chem),min(t_c,t_h));
        if (nH0>60 and nH0<1.e4) Dt = min(.001*t_rcb,Dt); //改成0.0001也改不了-4 改成0.01则-5也不对

        //********************************************************************************
        // add_Nt没关系 设为2
        //以下3行 可以给出比较一致的-3 -4 -5 但是500->1000这个判断不好 是判断的t_rcb比较小的时候
        Dt = 0.1*min( min(t_ff,100*t_chem),min(t_c,t_h));
        Dt = min(.001*t_rcb,Dt);
        if (nH0<1000. and nH0>500.) Dt = min(.0001*t_rcb,Dt);
        //*********************************************************************************

        //以下得到-3下降更快 reaction算得更细了 按理应该对 试了-1 和-3重合 都快于-4和-5
        Dt = 0.1*min( min(t_ff,100*t_chem),min(t_c,t_h));
        Dt = min(.0001*t_rcb,Dt);
        
        //关掉H2 cooling 放宽了t_rcb尝试
        Dt = 0.01*min( min(t_ff,100*t_chem),min(t_c,t_h));
        //Dt = min(.001*t_rcb,Dt);
        
        //不在一条线上 -4下降的更慢 -3 -5重合 add_Nt = 5/2
        /* Dt = 0.1*min( min(t_ff,100*t_chem),min(t_c,t_h));
        if (nH0>60 and nH0<1.e4) Dt = min(.001*t_rcb,Dt); */

        /* Dt = 0.1*min( min(t_ff,100.*t_chem),min(t_c,t_h));
        if (nH0>60 and nH0<1.e4) Dt = min(0.0001*t_rcb,Dt); */

        //printf("in TIMESCALES:t_c=%3.2e,t_h=%3.2e,t_ff=%3.2e,t_chem=%3.2e\n",t_c,t_h,t_ff,t_chem);
    }
    // time++ by Dt
    t_act += Dt;
    z0 = z;
    z = z_ana(z,Dt); //推进 redshift
}

// cloud evolution: freefall() -> react_sol() -> T_sol

void GAS:: freefall(){  //module of explicit integration over Dt
    HALO halo1(Mh,z0);
    bool adjust_iso = false;
    switch(evol_stage){
        case 0:
            nH0 = n_ff(z0,nH0,rhoc_DM,Dt);
            break;
        case 1: 
            nH0 = N_ADB(S0,T_K0);
            //printf("adb");
            //compare with maximum core density by core entropy: Visbal et al. 2014a
            if (N_ADB(S0,T_K0)>N_CORE(z0)) evol_stage = 2;
            break;
        case 2: 
            nH0 = N_CORE(z0); 
            //printf("saturation");

            if (r_cH >= abs(r_h)) not_adb = true;
            // 对4个trees, 下一行两个判断一样的结果: gas在同一个z算n_iso, 算出来都是3tree unstable.
            //if (not_adb && T_K0<8000){ //
            if (r_cH >= abs(r_h) && T_K0<8000){
                cout<<"/////////////\t t_c<t_h/2 \t//////////////FIRST ENTER ISO STAGE\n";
                printf("z = %3.2f\t Tvir = %3.2eK\t Mh = %3.2eMs\t Mgas = %3.2e\t Tg=%3.2eK\n", z0, halo1.Tvir, Mh/Ms,fb*Mh/Ms, T_K0);
                cout<<"ng_max= "<<R_EQ(T_K0,halo1.rho_c,halo1.Rs)*halo1.rho_c/(mu*m_H)<<"\t";
                cout<<"ng_adb= "<<nH0<<endl;
                cout<<"R_adb ="<< nH0*(mu*m_H)/halo1.rho_c;
                evol_stage = 3;
                n_iso = Mg2ng(Mh*fb,nH0,f_Ma*T_K0,z0,Mh);
                printf("//////////\tSOLVED FOR THE 1ST TIME\n");
                if (!n_iso) evol_stage = 4; // unstable case Mg2ng return 0
                else nH0 = n_iso;
                Mh_prev = Mh; t_prev = t_act;
            }
            break;
        case 3:
            adjust_iso = (Mh>2*Mh_prev);
            adjust_iso = (t_act - t_prev >= t_freefall(nH0));
            //adjust_iso = false; 
            if (adjust_iso) {
                printf("###\n");
                printf("IN STAGE 3 f_Ma=%3.2e\n",f_Ma);
                dt_iso = t_act - t_prev;
                Mh_prev = Mh; t_prev = t_act;
                n_iso = Mg2ng(Mh*fb,nH0,f_Ma*T_K0,z0,Mh);
                if (!n_iso) evol_stage = 4; // unstable case Mg2ng return 0
                else nH0 = n_iso;
            }
            break;
        case 4:
            nH0 = n_ff(z0,nH0,rhoc_DM,Dt);
            break;
    }

    /* updating e and Mh */
    de_tot += Gamma_mer*Dt;
    dM_tot += dMdt * Dt;

    reduction = 1;// WL ADDED //cout<<nH0<<"\t"<<ncore<<endl;
    v_tur2 += Dt * Gamma_mer_k*2; // 2 coz e=1/2v^2
    f_Ma = (Ma_on)? 1 + v_tur2/pow(cs,2) * reduction :1; //wrong, didn't consider cs^2/gammas

    f_Ma = (Ma_on)? 1 + v_tur2/pow(cs,2) * gamma_adb :1; // corrected f_Ma, using P = rho_g v_tur^2
    f_Ma = (Ma_on)? 1 + v_tur2/pow(cs,2) * gamma_adb/3. :1; // corrected f_Ma, using P = rho_g v_tur^2/3 from Chandrasekhar 1951
    Ma = sqrt(v_tur2* reduction )/cs;
 
    if (MerMod==0) f_Ma = 1; //for MerMod=0 case

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
    delete [] MPs;
    file_ingas.close();
}