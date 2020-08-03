// include 1 merger tree
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <fstream>
#include <algorithm>

#include "evol.h"
#include "class_gas.h"
#include "class_halo.h"
#include "my_linalg.h"
#include "thermo.h"
#include "dyn.h"
#include "PARA.h"
using namespace std;

//* previous    1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
/****************************************************
*     SPECIES                                      *
*     1 : H      2 : H2     3 : e      4 : H+      *
*     5 : H2+    6 : H-                            *
*     7 : He     8 : He+    9 : He++               *
****************************************************/
/* double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-6, y_Hm = 1.0e-12, y_H2p = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;
    
    tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-7, y_Hm = 1.0e-13, y_H2p = 1.0e-13;
    y_Hp = 1.0e-5;
    y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    y_He = yHe - 2.*tiny; y_Hep = tiny; y_Hepp = tiny;
    y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;

    double frac0[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};
 */
double evol_tiny = 1.0e-20, evol_yHe = 8.33333e-2, evol_y_H2 = 1.0e-5, evol_y_Hm = 1.0e-12, evol_y_H2p = 1.0e-12;
double evol_y_Hp = 1.0e-4;
double evol_y_H = 1.0 - 2.*evol_y_H2 - 2.*evol_y_H2p - evol_y_Hm - evol_y_Hp;
double evol_y_He = evol_yHe - 2.*evol_tiny, evol_y_Hep = evol_tiny, evol_y_Hepp = evol_tiny;
double evol_y_e = evol_y_Hp + evol_y_H2p - evol_y_Hm + evol_y_Hep + 2.*evol_y_Hepp;

double frac0[] = {0., evol_y_H, evol_y_H2, evol_y_e, evol_y_Hp, evol_y_H2p, evol_y_Hm, evol_y_He, evol_y_Hep, evol_y_Hepp};
static int i;

void evol(string treename, string fout, int MerMod, double Tbb, double J21, bool spec, bool Ma_on, int i_bsm){
    printf("################################################################################\n");
    printf("################################################################################\n");
    printf("################################################################################\n");
    i =0;

    GAS gas(frac0,MerMod,J21,Tbb,treename,spec,Ma_on,i_bsm);
    printf("z0 = %3.2f, initial nH = %3.2e /cc & T = %3.2e K\n",gas.z0,gas.nH0,gas.T_K0);
    printf("t_ff,0 = %e Myr\n",gas.t_ff0/Myr);
    printf("t_mergers in main is :  %3.2e Myr\n", (gas.MPs[gas.nMer-1].t - gas.MPs[0].t)/Myr );
    printf("M_J is %3.2e Ms for n = %3.2e and T = %3.e K\n",gas.MJ0/Ms,gas.nH0,gas.T_K0);
    printf("at z=10, cosmic n_crit = %3.2e cm^-3\n",RHO_crit(10)/(mu*m_H));

//merger tree information in GAS::gas
    // ofstream f1;
    // f1.open("../data/atree.txt", ios::out | ios::trunc );
    // f1<<setiosflags(ios::scientific)<<setprecision(5);
    // for (i=0;i<gas.nMer;i++){
    //     HALO halo(gas.MPs[i].mhalo,gas.MPs[i].z);
    //     if (i==0) f1<<"\tj\tz\tt(Myr)\tdt(Myr)\tt_dyn(Myr)\tM8\tdM8\tq\tRvir\tRs\tc\tTvir\trho_c\trho_crit\td0\n";
    //     f1<<setw(16)<<gas.MPs[i].j<<setw(16)<<gas.MPs[i].z<<setw(16)<<gas.MPs[i].t/Myr<<setw(16)<<gas.MPs[i].dt/Myr<<setw(16)<<halo.t_dyn/Myr;
    //     f1<<setw(16)<<gas.MPs[i].mhalo/(1.e8*Ms)<<setw(16)<<gas.MPs[i].dm/(1.e8*Ms)<<setw(16)<<gas.MPs[i].mratio;
    //     f1<<setw(16)<<halo.Rvir<<setw(16)<<halo.Rs<<setw(16)<<halo.c<<setw(16)<<halo.Tvir;
    //     f1<<setw(16)<<halo.rho_c<<setw(16)<<halo.rho_crit<<setw(16)<<halo.delta0<<endl;
    // }
    // f1.close();
    /* for (i=1; i<gas.nMer-1; i++){
        printf("time = %6.3e,\tdt = %6.3e\n", gas.MPs[i].t/Myr, (gas.MPs[i].t-gas.MPs[i-1].t)/Myr);
    } */

    ofstream file;
    file.open(fout, ios::out | ios::trunc);
    bool py = true;
    bool DM = true;
    bool fract = true;
    bool react = false;
    bool tscales = true;
    bool haloinfo = true;
    bool heatingcooling = true;
    bool mer = true;
    double kform, rform, yequi, ycool, ypred;
    int itime=0;

    i = 0;
    //while (gas.t_act<2*gas.t_ff0){
    while (gas.z>6. && gas.nH0<1.e6){
        if (i==0) file<<setw(16)<<"t"<<setw(16)<<"Dt"<<setw(16)<<"z"<<setw(16)<<"nH"<<setw(16)<<"T";
        else file<<setw(16)<<gas.t_act/gas.t_ff0<<setw(16)<<gas.Dt/gas.t_ff0<<setw(16)<<gas.z<<setw(16)<<gas.nH0<<setw(16)<<gas.T_K0;
            
        //printf("nH=%3.2e, T_K=%3.2e\t k[15,+]=%3.2e, k[7,-]=%3.2e, y_H2=%3.2e\n",gas.nH0, gas.T_K0,gas.k[15],gas.k[7],gas.y0[2]);

        //cout<<z_ana(35,gas.t_act)<<gas.z<<endl; //not exactly the same...

        if (tscales) {            
            if (i==0) {
                file<<setw(16)<<"t_ff"<<setw(16)<<"t_c"<<setw(16)<<"t_h"<<setw(16)<<"t_rcb"<<setw(16)<<"t_chem";
                file<<setw(16)<<"t_ion"<<setw(16)<<"tc_H2"<<setw(16)<<"tc_H"<<setw(16)<<"MPs_dt"<<setw(16)<<"i2";
            }
            else {
                file<<setw(16)<<gas.t_ff/gas.t_ff0;
                file<<setw(16)<<gas.t_c/gas.t_ff0;
                file<<setw(16)<<gas.t_h/gas.t_ff0;
                file<<setw(16)<<gas.t_rcb/gas.t_ff0;
                file<<setw(16)<<gas.t_chem/gas.t_ff0;

                file<<setw(16)<<gas.t_ion/gas.t_ff0;
                file<<setw(16)<<gas.e0/gas.r_cH2 /gas.t_ff0;
                file<<setw(16)<<gas.e0/gas.r_cH / gas.t_ff0;
                file<<setw(16)<<gas.MPs[gas.iMP].dt/gas.t_ff0;
                file<<setw(16)<<0;
                /* 
                file<<setw(16)<<gas.e0*gas.rho0/Lambda_H2(gas.nH0,gas.T_K0,gas.y0); //t_cH2
                file<<setw(16)<<gas.e0*gas.rho0/Lambda_H(gas.nH0,gas.T_K0,gas.y0[1],gas.y0[3],gas.k[1]); //t_cH
                 */
            }
        }

        if (fract){
            if (i==0){
                file<<setw(16)<<"yH"<<setw(16)<<"yH2"<<setw(16)<<"ye"<<setw(16)<<"yH+"<<setw(16)<<"yH-";
                file<<setw(16)<<"y_euqi"<<setw(16)<<"y_cool"<<setw(16)<<"y_pd"<<setw(16)<<"y_cd"<<setw(16)<<"Jc_pred";
            }
            else {
                for (int i=1; i<=4; i++) file<<setw(16)<<gas.y0[i];
                file<<setw(16)<<gas.y0[6]; // Hm
                double kHm_e, kH2_Hm,   kH2p_Hp, kH2_H2p,   kH2_3b;
                double kHm_pd, kH2p_pd, kH2_pd;
                double kH2_cd_tot, kH2_cd_H, kH2_cd_e, kH2_cd_H2, kH2_cd_Hp, kH2_cd_He; 
                double kform_Hm, kform_H2p;
                double kHm_rcb, kHm_cd;
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

                kHm_e = gas.k[3]; kH2_Hm = gas.k[4];
                kH2p_Hp = gas.k[5]; kH2_H2p = gas.k[6];

                kH2_pd = gas.k[21]; kHm_pd = gas.k[22]; kH2p_pd = gas.k[23];
                kH2_3b = gas.k[15];
                kH2_cd_H = gas.k[7]; kH2_cd_Hp = gas.k[8]; kH2_cd_e = gas.k[9] + gas.k[20]; kH2_cd_H2 = gas.k[16];
                kH2_cd_He = gas.k[35];
                kH2_cd_tot = gas.nH0* (kH2_cd_H*gas.y0[1] + kH2_cd_Hp*gas.y0[4] + kH2_cd_e*gas.y0[3] + kH2_cd_H2*gas.y0[2] + kH2_cd_He*gas.y0[7]); 
                // kH2_cd_tot = gas.nH0* (kH2_cd_H*gas.y0[1]); collisional dissociation dominated by H

                // regualte H2 at high n, together w/ H2 cd 
                kHm_rcb = gas.k[11];
                kHm_cd = gas.k[18];

                kform_Hm = kHm_e* kH2_Hm/( kH2_Hm*gas.y0[1] + kHm_pd/gas.nH0 + kHm_rcb*gas.y0[4] + kHm_cd*gas.y0[1]);
                kform_H2p = kH2p_Hp* kH2_H2p/(kH2_H2p + kH2p_pd/gas.nH0);
                // H2 forming rate; 
                rform = max(kform_Hm*gas.nH0*gas.y0[1]*gas.y0[3]+kform_H2p*gas.nH0*gas.y0[1]*gas.y0[4], kH2_3b*pow(gas.nH0,2)*pow(gas.y0[1],3));
                            // H-, small n                        H2+                                     3b,  n>~10^9/cm^3
                rform = kform_Hm*gas.nH0*gas.y0[1]*gas.y0[3]; // H2 formation dominated by H- over H2+
                yequi = min(rform/kH2_pd, rform/kH2_cd_tot);
                yequi = rform/(kH2_pd+kH2_cd_tot);
                            // pd,  n<~100/cm^3      cd, large n 

                file<<setw(16)<<gas.yequi; //column 21
                file<<setw(16)<<gas.ycool;
                file<<setw(16)<<gas.ypd;
                file<<setw(16)<<gas.ycd; 
                file<<setw(16)<<gas.Jc_pred;
            }
        }
        /* if (itime <=2){
            if (gas.nH0>1.e0 && gas.y0[2]>=ycool){
            printf("-------------------------->\nonce y_H2 > ycool\n");
            printf("nH = %e\n",gas.nH0);
            itime += 1;
            }
            if (yequi>=1){
                printf("yequi >= 1, follow normal H2 formation track \n");
                itime ++;
            }
        } */

        if (DM){
            HALO halo(gas.Mh,gas.z);
            double d0 = 200./3.*pow(halo.c,3)/(log(1+halo.c)-halo.c/(1+halo.c));
            if (i==0){
                file<<setw(16)<<"z_Mh"<<setw(16)<<"t_Mh"<<setw(16)<<"Mh"<<setw(16)<<"Tvir"<<setw(16)<<"Vc";
                file<<setw(16)<<"cs"<<setw(16)<<"v_tur"<<setw(16)<<"v_bsm"<<setw(16)<<"f_Ma"<<setw(16)<<"g_Ma";
            } //rho_c = d0*RHO_DM(gas.z)/m_H checked right. 
            else{
                file<<setw(16)<<gas.MPs[gas.iMP].z<<setw(16)<<gas.MPs[gas.iMP].t/gas.t_ff0;
                file<<setw(16)<<gas.Mh/Ms<<setw(16)<<gas.MPs[gas.iMP].Tvir<<setw(16)<<halo.Vc/km;

                file<<setw(16)<<gas.cs/km<<setw(16)<<sqrt(gas.v_tur2)/km<<setw(16)<<gas.v_bsm/km;
                file<<setw(16)<<gas.f_Ma;
                file<<setw(16)<<gas.g_Ma;
            }
        }

        if (haloinfo){
            if (i==0) file<<setw(16)<<"ievol"<<setw(16)<<"M_J"<<setw(16)<<"M_Jeff"<<setw(16)<<"Mg_intg"<<setw(16)<<"q";
            else file<<setw(16)<<gas.evol_stage<<setw(16)<<gas.MJ/Ms<<setw(16)<<gas.MJ_eff/Ms<<setw(16)<<gas.Mg_intg/Ms<<setw(16)<<gas.MPs[gas.iMP].mratio;
        }


        if (i!=0){
            gas.setMerger();
            gas.timescales(); 
            gas.freefall(); 
            //cout<<"freefall done. "<<"n= "<<gas.nH0<<" /cc"<<endl;
            gas.react_sol(1); 
            //cout<<"react_sol done"<<endl;
            gas.T_sol();
            //cout<<"T_sol done. "<<"T= "<<gas.T_K0<<" K"<<endl;
            gas.get_para();
        }
        if (heatingcooling) {
            if (i==0) file<<setw(16)<<"r_c"<<setw(16)<<"r_cH2"<<setw(16)<<"r_cH"<<setw(16)<<"r_hcompr"<<setw(16)<<"r_hmer";
            else{
                file<<setw(16)<<gas.r_c;
                file<<setw(16)<<gas.r_cH2;
                file<<setw(16)<<gas.r_cH;
                if (gas.evol_stage==4 or gas.evol_stage==0) file<<setw(16)<<Gamma_compr(gas.cs,gas.f_Ma,gas.t_ff);
                else file<<setw(16)<<Gamma_compr(gas.cs,gas.f_Ma,gas.t_ff);
                file<<setw(16)<<gas.Gamma_mer_th; //merger heating; kinetic energy input: *1/3
            }
        }

        if (mer) {
            HALO halo(gas.Mh, gas.z);
            if (i==0) file<<setw(16)<<"iMP"<<setw(16)<<"nc_DM"<<setw(16)<<"dMdt"<<setw(16)<<"Gamma_mer";
            else file<<setw(16)<<gas.iMP<<setw(16)<<halo.rho_c/(mu*m_H)<<setw(16)<<gas.dMdt<<setw(16)<<gas.Gamma_mer;
        }
        file<<endl;
        i++;
    //printf("%3.2e\t %3.2e\n", Gamma_compr(gas.nH0,gas.T_K0,gas.f_Ma), Lambda_H2(gas.nH0,gas.T_K0,gas.y0));
    //cout<<gas.inDelay<<endl;
    }
    file.close();

// cout<<"\n\nff timescale nH=1.0 is : "<<t_freefall(1.)/Myr<<" Myr"<<endl;
    // printf("R_J is %3.2e pc for n = %3.2e and T = %3.e K\n",gas.RJ/pc,gas.nH0,gas.T_K0);
    // printf("evolve time is %3.2e Myr\n",gas.t_act/Myr);
    // printf("final Ma in main is :  %f \n", gas.Ma);
    // printf("v_tur in main is :  %f km/s\n", sqrt(gas.v_tur2/1.e10));
    // printf("final de_tot by merger in main is :  %f erg\n", gas.de_tot);
    // printf("final dT_tot by merger in main is :  %f K\n", gas.de_tot*(gamma_adb-1)*(mu*m_H)/k_B);
    // printf("final dM_tot by merger in main is :  %5.3e Ms\n", gas.dM_tot/Ms);
    // HALO halo (gas.Mh,gas.z);
    // printf(" z from 35 to 10 %f Myr\n", (t_from_z(10) - t_from_z(35))/Myr );
    // printf(" freefall timescale %f Myr\n", gas.t_ff0/Myr );
    // printf("1.93/2.2*dTvir = %f K\n",1.89/2.2*(gas.MPs[gas.iMP].Tvir - gas.MPs[0].Tvir));
    printf("J21= %5.2f, z_col=%3.2f\n\n\n\n",J21,gas.z_col);
}

// resembling main_Jc1.cpp
double getT(int argc, double* argv, bool write,int MerMod, double J, double Tb, string treename, bool spec, bool Ma_on, int i_bsm, double nH_tell){
    // printf("\n**************************\nSTART in getT:\n");
    GAS gas(frac0,MerMod,J,Tb,treename,spec,Ma_on,i_bsm);

    int index=treename.find("_");
    string tree = treename.substr(index+1); // tree_id 输出

    int i=0;
    fstream f1;
    if (write){
        string fout = "tr"+tree+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+"J"+to_string(int(J))+".txt";
        f1.open(fout, ios::out | ios::trunc );
        f1<<setiosflags(ios::scientific)<<setprecision(5);
    }
    while (gas.nH0<nH_tell){
        gas.setMerger();
        gas.timescales(); 
        gas.freefall(); 
        gas.react_sol(1); 
        gas.T_sol();
        gas.get_para();
       // if( fmod(i,1000)==0 ) printf("%3.2f\n",gas.z_col);
       if (write){
           if (i==0) {
                // cout<<"getT: writing file \" "<<fout<<"\"\n";
                f1<<setw(16)<<"t"<<setw(16)<<"ievol"<<setw(16)<<"z"<<setw(16)<<"nH"<<setw(16)<<"T";
                f1<<setw(16)<<"cs"<<setw(16)<<"v_bsm"<<setw(16)<<"v_tur"<<setw(16)<<"Vc"<<setw(16)<<"f_Ma";
                f1<<setw(16)<<"Mh"<<setw(16)<<"Tv"<<endl;
                f1<<setw(16)<<"y_H2"<<setw(16)<<"y_e"<<endl;
            }
            else if(fmod(i,100)==0) {
                HALO halo(gas.Mh,gas.z0);
                f1<<setw(16)<<gas.t_act/gas.t_ff0<<setw(16)<<gas.evol_stage<<setw(16)<<gas.z;
                f1<<setw(16)<<gas.nH0<<setw(16)<<gas.T_K0<<setw(16)<<gas.cs/km<<setw(16)<<gas.v_bsm/km;
                f1<<setw(16)<<sqrt(gas.v_tur2)/km<<setw(16)<<halo.Vc/km<<setw(16)<<gas.f_Ma;
                f1<<setw(16)<<gas.Mh/Ms<<setw(16)<<halo.Tvir;
                f1<<setw(16)<<gas.y1[2]<<setw(16)<<gas.y1[3]<<endl;
            }
       }

        i++;
    }
    if (write) f1.close();

    argv[0] = gas.n_H2crit;
    argv[1] = gas.z_H2crit;
    argv[2] = gas.Jc_pred_max;
    argv[3] = gas.fMa_H2crit;
    argv[4] = gas.gMa_H2crit;
    argv[5] = gas.z_col;

    // printf("*******IN GET_T***********\n");
    // cout<<"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0);
    // printf("\tJ=%3.2f, z_col=%3.2f, nH_tell=%3.2e, T=%3.2e\n\n\n",J, gas.z_col, nH_tell, gas.T_K0);
    return gas.T_K0;
}


void evol_Jc(string treename, string fout, double Tb, int MerMod, bool spec, bool Ma_on, int i_bsm){
    printf("############################\t JC_SOL\t ########################\n");
    printf("i_bsm=%d Ma_on=%d\n",i_bsm,(Ma_on)?1:0);
    double T_tell = 4000, nH_tell=1.e4;
    // boundary of bisection J21 
    double J0 = epE2, J1 = 1.e4; //包括所有的Tb所需range 没必要
    J0 = 500, J1 = 1.5e3;
    double T0, T1, T;
    int const c = 6;
    double y0[c];
    double y1[c];
    double y[c];

    ofstream file;
    ifstream checkf_exist(fout.c_str());
    if (checkf_exist.good()) cout<<fout.c_str()<<" exist\n";
    else {
        file.open(fout, ios::out | ios::trunc);
        file<<setw(16)<<"tree"<<setw(16)<<"Tb"<<setw(16)<<"i_bsm"<<setw(16)<<"tur"; 
        // file<<setw(16)<<"n_H2crit"<<setw(16)<<"z_H2crit";
        file<<setw(16)<<"Jc";
        file<<setw(16)<<"z_col"<<endl;
        file<<endl;
        file.close();
    }
    file.open(fout, ios::out | ios::app);

    T0 = getT(c,y0, false,MerMod, J0, Tb, treename, spec, Ma_on, i_bsm,nH_tell); 
    T1 = getT(c,y1, false,MerMod, J1, Tb, treename, spec, Ma_on, i_bsm,nH_tell);
    cout<<"***********log**************\n";
    if ( T0-T_tell>0 ) {
        while(T0-T_tell>0){
            printf("J0=%3.2f,T0=%3.2e,T_tell=%3.2e, wrong INITIAL LEFT boundary\n",J0,T0,T_tell);
            J1 = J0; T1 = T0;
            J0 *= 0.5;
            T0 = getT(c,y0,false,MerMod, J0, Tb, treename, spec, Ma_on, i_bsm,nH_tell);
        }
    }
    else if (T1-T_tell<0) {
        while (T1-T_tell<0){
            printf("J1=%3.2f,T1=%3.2e,T_tell=%3.2e, wrong INITIAL RIGHT boundary\n",J1,T1,T_tell);
            J0 = J1; T0 = T1;
            J1 *= 2.;
            T1 = getT(c,y1,false,MerMod, J1, Tb, treename, spec, Ma_on, i_bsm,nH_tell); 
        }
    }

    while (J1-J0 > 0.05*J0){
        bool write = false;
        // if (J1-J0 > 0.05*J0) write = true;
        T = getT(c,y,write,MerMod, (J0+J1)/2., Tb, treename, spec, Ma_on, i_bsm,nH_tell);
        printf("#######\t########\t########\t#########");
        printf("J0=%3.2f, T0=%3.2f,J1=%3.2f,T1=%3.2f, Jmid=%4.3f\tTmid=%3.2e\n",J0,T0,J1,T1,(J0+J1)/2.,T);
        if (T-T_tell>=0) {
            J1 = (J0+J1)/2.; T1 = T;
            for (i=0;i<c;i++) y1[i] = y[i];
        }
        else {
            J0 = (J0+J1)/2.; T0 = T;
            for (i=0;i<c;i++) y0[i] = y[i];
        }
    }


    printf("J0=%3.2f\tT0=%3.2e\tJ1=%3.2f\tT1=%3.2e\n --> Jc_sol=%3.2f z_col=%3.2f\n",J0,T0,J1,T1, J1, y1[5]);
    // int index=treename.find("fort.");
    // string tree = treename.substr(index+5); // tree_id 输出
    
    int index=treename.find("_");
    string tree = treename.substr(index+1); // tree_id 输出

    file<<setw(16)<<stoi(tree)<<setw(16)<<Tb<<setw(16)<<i_bsm<<setw(16)<<((Ma_on)?1:0);
    file<<setw(16)<<J1;
    file<<setw(16)<<y1[5];
    file<<endl;
    file.close();

        // file<<"tree"<<"Tb"<<"i_bsm"<<"tur"; 
        // file<<"n_H2crit"<<"z_H2crit";
        // file<<"Jc"<<"Jc_pred";
        // file<<"fMa_H2crit"<<"gMa_H2crit";
        // file<<"z_col"<<endl;
//     argv[0] = gas.n_H2crit;
//     argv[1] = gas.z_H2crit;
//     argv[2] = gas.Jc_pred;
//     argv[3] = gas.fMa_H2crit;
//     argv[4] = gas.gMa_H2crit;
//     argv[5] = gas.z_col;

    string fevol = "tr"+tree+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+"J"+to_string(int(J0))+".txt";
    // evol(treename, fevol, MerMod, Tb, J0, spec, Ma_on, i_bsm);
    cout<<fevol<<endl;
    fevol = "tr"+tree+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+"J"+to_string(int(J1))+".txt";
    // evol(treename, fevol, MerMod, Tb, J1, spec, Ma_on, i_bsm);
    cout<<fevol<<endl;
}