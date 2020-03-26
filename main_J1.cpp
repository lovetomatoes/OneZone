// include 1 merger tree
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include<algorithm>

#include "reaction.h"
#include "kpd.h"
#include "class_gas.h"
#include "LE_iso.h"
#include "class_halo.h"
#include "my_linalg.h"
#include "Newton.h"
#include "thermo.h"
#include "dyn.h"
#include "PARA.h"
using namespace std;

//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
int i =0;

int main(){
    printf("################################################################################\n");
    printf("################################################################################\n");

    int MerMod = 0;
    double J21 = 100.; //without mergers, WG11 shielding, Jcrit~50 for Tb=1.e4K, Jcrit~3000 for Tb=1.e5K
    double Tbb = 8.e3;

    double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-6, y_Hm = 1.0e-12, y_H2p = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;
    double nH = 1.e8;
    double T_K = 1.e3;
    //* previous    1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
    /****************************************************
    *     SPECIES                                      *
    *     1 : H      2 : H2     3 : e      4 : H+      *
    *     5 : H2+    6 : H-                            *
    *     7 : He     8 : He+    9 : He++               *
    ****************************************************/
    double frac0[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};

    char* fname = "../code_tree/fort.217";
    GAS gas(frac0,MerMod,J21,Tbb,fname,true);
    double t_ff0 = 1./C/sqrt(gas.nH0);
    double t1 = 1.9999*t_ff0;
    printf("z0 = %3.2f, initial nH = %3.2e /cc & T = %3.2e K\n",gas.z0,gas.nH0,gas.T_K0);
    printf("t_ff,0 = %e Myr\n",t_ff0/Myr);
    printf("t1 in main is :  %f t_ff0\n", t1/t_ff0);
    printf("t_delay in main is :  %f t_ff0\n", gas.t_delay/t_ff0);
    printf("t_mergers in main is :  %3.2e Myr\n", (gas.MPs[gas.nMer-1].t - gas.MPs[0].t)/Myr );
    printf("M_J is %3.2e Ms for n = %3.2e and T = %3.e K\n",gas.MJ0/Ms,gas.nH0,gas.T_K0);
    printf("at z=10, n_mean = %3.2e cm^-3\n",RHO_crit(10)/m_H);

    //merger tree information in GAS::gas
    ofstream f1;
    f1.open("../data/atree.txt", ios::out | ios::trunc );
    for (i=0;i<gas.nMer;i++){
        HALO halo(gas.MPs[i].mhalo,gas.MPs[i].z);
        if (i==0) f1<<"\tj\tz\tt(Myr)\tdt(Myr)\tt_dyn(Myr)\tM8\tdM8\tq\tRvir\tRs\tc\tTvir\trho_c\trho_crit\td0\n";
        f1<<" "<<gas.MPs[i].j<<" "<<gas.MPs[i].z<<" "<<gas.MPs[i].t/Myr<<" "<<gas.MPs[i].dt/Myr<<" "<<halo.t_dyn/Myr;
        f1<<" "<<gas.MPs[i].mhalo/(1.e8*Ms)<<" "<<gas.MPs[i].dm/(1.e8*Ms)<<" "<<gas.MPs[i].mratio;
        f1<<" "<<halo.Rvir<<" "<<halo.Rs<<" "<<halo.c<<" "<<halo.Tvir;
        f1<<" "<<halo.rho_c<<" "<<halo.rho_crit<<" "<<halo.delta0<<endl;
    }
    f1.close();

    /* for (i=1; i<gas.nMer-1; i++){
        printf("time = %6.3e,\tdt = %6.3e\n", gas.MPs[i].t/Myr, (gas.MPs[i].t-gas.MPs[i-1].t)/Myr);
    } */

    ofstream file;
    file.open("../data/evolve.txt", ios::out | ios::trunc );
    bool py = true;
    bool DM = true;
    bool fract = true;
    bool react = false;
    bool tscales = true;
    bool haloinfo = false;
    bool heatingcooling = false;
    bool mer = false;
    double kform, rform, yequi, ycool;
    int itime=0;

    i = 0;
    //while (gas.t0<gas.t1 && gas.nH0<1.e15){
    while (gas.z>10 && gas.nH0<1.e11){
        if (i==0){
            if (py) file <<" t Dt z nH T";
            else file <<" t Dt z nH T";
        }
        else file<<" "<<gas.t_act/t_ff0<<" "<<gas.Dt<<" "<<gas.z<<" "<<gas.nH0<<" "<<gas.T_K0;
            
        //printf("nH=%3.2e, T_K=%3.2e\t k[15,+]=%3.2e, k[7,-]=%3.2e, y_H2=%3.2e\n",gas.nH0, gas.T_K0,gas.k[15],gas.k[7],gas.y0[2]);

        //cout<<z_ana(35,gas.t_act)<<gas.z<<endl; //not exactly the same...
        if (fract){
            if (i==0){
                if (py) file<<" yH yH2 ye yH+ yH2+ yH- yHe yHe+ yHe++ Ma v_tur MJ_eff"; //y_equi y_cool
                else    file<<" y_H y_{H2} y_e y_{H+} y_{H2+} y_{H-} y_{He} y_{y_He+} y_{He++} Ma v_{tur} MJ_{eff}";
            }
            else {
                for (int i=1; i<N_sp+1; i++) file<<" "<<gas.y0[i];

                kform = gas.k[2]*gas.k[3]/(gas.k[3]+gas.k[7]/gas.nH0);
                rform = max(kform*gas.nH0*gas.y0[1]*gas.y0[3], gas.k[5]*pow(gas.nH0,2)*pow(gas.y0[1],3));
                            //H-, n小                              3b,  n>~10^9/cm^3
                yequi = min(rform/gas.k[6], rform/(gas.k[4]*gas.y0[1]*gas.nH0) );
                                // pd,  n<~100/cm^3      cd n大 
                //file<<" "<<yequi;
                //file<<" "<<gas.Dt*gas.rf[2]+gas.y0[2]; prediciton of yH2
                //ycool = 3.*k_B*gas.T_K0*gas.y0[2]/(2.*gas.t_ff*Lambda_H2(gas.nH0,gas.T_K0,gas.y0[2])*(mu*m_H));
                //file<<" "<<ycool;

                file<<" "<<gas.Ma<<" "<<sqrt(gas.v_tur2);
                file<<" "<<gas.MJ_eff/Ms;

            }
        }

        if (DM){
            HALO halo(gas.Mh,gas.z);
            double d0 = 200./3.*pow(halo.c,3)/(log(1+halo.c)-halo.c/(1+halo.c));
            if (i==0){
                if (py) file <<" nc_DM z_Mh t_Mh Vc inMer";
                else file <<" nc_{DM} z_{Mh} t_{Mh} Vc inMer";
            } //rho_c = d0*RHO_DM(gas.z)/m_H checked right. 
            else file<<" "<<gas.rhoc_DM/m_H<<" "<<gas.MPs[gas.iMer].z<<" "<<gas.MPs[gas.iMer].t/gas.t_ff0<<" "<<halo.Vc<<" "<<((gas.inMer)?1:0);
        }
        if (react){
            if (i==0){
                if (py) file<<" k1 k2 k3 k4 k5 k6 k7 k8";
                else file<<" k1 k2 k3 k4 k5 k6 k7 k8";
            }
            else for (int i=1; i<N_react+1; i++) file<<" "<<gas.k[i];
            //if (i==0) file<<"\trf1 \trf2 \trf3 \trf4 \trf5 ";
            //else for (int i=1; i<N_sp+1; i++) file<<" "<<gas.rf[i];
        }

        /* if (itime <=2){
            if (gas.nH0>1.e0 && gas.y0[2]>=ycool){
            printf("-------------------------->\nonce y_H2 > ycool\n");
            printf("nH = %e\n",gas.nH0);
            itime += 1;
            }
            if (yequi>=1){
                printf("yequi >= 1, follow normal H2 formation track \n");
                itime += 1;
            }
        } */

        if (haloinfo){
            if (i==0) {
                if (py) file<<" iMer Mh Tvir dMh q";
                else file<<" iMer Mh T_{vir} dMh q";
            }
            else file<<" "<<gas.iMer<<" "<<gas.Mh/Ms<<" "<<gas.MPs[gas.iMer].Tvir<<" "<<gas.MPs[gas.iMer].dm/Ms<<" "<<gas.MPs[gas.iMer].mratio;
        }

        if (tscales) {            
            if (i==0) {
                if (py) file<<" t_ff t_cH2 t_cH t_h t_rcb t_chem";
                else file<<" t_{ff} t_{c,H2} t_{c,H} t_h t_{rcb} t_{chem}";
            }
            else {
                file<<" "<<gas.t_ff;
                file<<" "<<gas.e0*gas.rho0/Lambda_H2(gas.nH0,gas.T_K0,gas.y0);
                file<<" "<<gas.e0*gas.rho0/Lambda_H(gas.nH0,gas.T_K0,gas.y0[1],gas.y0[3],gas.k[1]);
                file<<" "<<gas.t_h<<" "<<gas.t_rcb;
                file<<" "<<0;
            }
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
            if (i==0){
                if (py) file<<" r_c r_cH2 r_cH r_hcompr r_hmer";
                else file<<" r_c r_{c,H2} r_{c,H} r_{h,compr} r_{h,mer}";

            }
            else{
                file<<" "<<gas.r_c;
                file<<" "<<Lambda_H2(gas.nH0,gas.T_K0,gas.y0);
                file<<" "<<Lambda_H(gas.nH0,gas.T_K0,gas.y0[1],gas.y0[3],gas.k[1]);
                file<<" "<<Gamma_compr(gas.cs,gas.f_Ma,gas.t_ff);
                file<<" "<<gas.Gamma_mer;
            }
        }


        if (mer) {
            if (i==0) {
                if (py) file<<" dMdt Gamma_mer Gamma_mer_th Gamma_mer_k";
                else file<<" dMdt Gamma_{mer} Gamma_{mer,th} Gamma_{mer,k}";
            }
            else file<<" "<<gas.dMdt<<" "<<gas.Gamma_mer<<" "<<gas.Gamma_mer_th<<" "<<gas.Gamma_mer_k;
        }
        file<<endl;
        i++;
    //printf("%3.2e\t %3.2e\n", Gamma_compr(gas.nH0,gas.T_K0,gas.f_Ma), Lambda_H2(gas.nH0,gas.T_K0,gas.y0));
    //cout<<gas.inDelay<<endl;
    }
    file.close();

    cout<<"\n\nff timescale nH=1.0 is : "<<t_freefall(1.)/Myr<<" Myr"<<endl;
    printf("R_J is %3.2e pc for n = %3.2e and T = %3.e K\n",gas.RJ/pc,gas.nH0,gas.T_K0);
    printf("evolve time is %3.2e Myr\n",gas.t_act/Myr);
    printf("final Ma in main is :  %f \n", gas.Ma);
    printf("v_tur in main is :  %f km/s\n", sqrt(gas.v_tur2/1.e10));
    printf("final de_tot by merger in main is :  %f erg\n", gas.de_tot);
    printf("final dT_tot by merger in main is :  %f K\n", gas.de_tot*(gamma_adb-1)*(mu*m_H)/k_B);
    printf("final dM_tot by merger in main is :  %5.3e Ms\n", gas.dM_tot/Ms);
    HALO halo (gas.Mh,gas.z);
    printf(" z from 35 to 10 %f Myr\n", (t_from_z(10) - t_from_z(35))/Myr );
    printf(" freefall timescale %f Myr\n", gas.t_ff0/Myr );
    printf("1.93/2.2*dTvir = %f K\n",1.89/2.2*(gas.MPs[gas.iMer].Tvir - gas.MPs[0].Tvir));
    return 0;
}