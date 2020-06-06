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

void evol(string treename, string fout, int MerMod, double Tbb, double J21, bool spec, bool Ma_on, int i_bsm){
    printf("################################################################################\n");
    printf("################################################################################\n");
    printf("################################################################################\n");
    int i =0;

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
    bool mer = false;
    double kform, rform, yequi, ycool;
    int itime=0;

    i = 0;
    //while (gas.t0<gas.t1 && gas.nH0<1.e15){
    while (gas.z>10 && gas.nH0<1.e11){
        if (i==0){
            if (py) file<<setw(16)<<"t"<<setw(16)<<"Dt"<<setw(16)<<"z"<<setw(16)<<"nH"<<setw(16)<<"T";
            else file <<" t Dt z nH T";
        }
        else file<<setw(16)<<gas.t_act/gas.t_ff0<<setw(16)<<gas.Dt/gas.t_ff0<<setw(16)<<gas.z<<setw(16)<<gas.nH0<<setw(16)<<gas.T_K0;
            
        //printf("nH=%3.2e, T_K=%3.2e\t k[15,+]=%3.2e, k[7,-]=%3.2e, y_H2=%3.2e\n",gas.nH0, gas.T_K0,gas.k[15],gas.k[7],gas.y0[2]);

        //cout<<z_ana(35,gas.t_act)<<gas.z<<endl; //not exactly the same...

        if (tscales) {            
            if (i==0) {
                if (py){
                    file<<setw(16)<<"t_ff"<<setw(16)<<"t_c"<<setw(16)<<"t_h"<<setw(16)<<"t_rcb"<<setw(16)<<"t_chem";
                    file<<setw(16)<<"t_ion"<<setw(16)<<"tc_H2"<<setw(16)<<"tc_H"<<setw(16)<<"i1"<<setw(16)<<"i2";
                }
                else file<<" t_{ff} t_{c} t_h t_{rcb} t_{chem} t_{ion} tc_{H2} tc_H 0 0";
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
                file<<setw(16)<<0<<setw(16)<<0;
                /* 
                file<<setw(16)<<gas.e0*gas.rho0/Lambda_H2(gas.nH0,gas.T_K0,gas.y0); //t_cH2
                file<<setw(16)<<gas.e0*gas.rho0/Lambda_H(gas.nH0,gas.T_K0,gas.y0[1],gas.y0[3],gas.k[1]); //t_cH
                 */
            }
        }

        if (fract){
            if (i==0){
                if (py){
                    file<<setw(16)<<"yH"<<setw(16)<<"yH2"<<setw(16)<<"ye"<<setw(16)<<"yH+"<<setw(16)<<"yH2+";
                    file<<setw(16)<<"yH-"<<setw(16)<<"yHe"<<setw(16)<<"yHe+"<<setw(16)<<"yHe++"<<setw(16)<<"i3"; //y_equi y_cool
                }
                else    file<<" y_H y_{H2} y_e y_{H+} y_{H2+} y_{H-} y_{He} y_{y_He+} y_{He++} 0";
            }
            else {
                for (int i=1; i<N_sp+1; i++) file<<setw(16)<<gas.y0[i];

                kform = gas.k[2]*gas.k[3]/(gas.k[3]+gas.k[7]/gas.nH0);
                rform = max(kform*gas.nH0*gas.y0[1]*gas.y0[3], gas.k[5]*pow(gas.nH0,2)*pow(gas.y0[1],3));
                            //H-, n小                              3b,  n>~10^9/cm^3
                yequi = min(rform/gas.k[6], rform/(gas.k[4]*gas.y0[1]*gas.nH0) );
                                // pd,  n<~100/cm^3      cd n大 
                //file<<setw(16)<<yequi;
                //file<<setw(16)<<gas.Dt*gas.rf[2]+gas.y0[2]; prediciton of yH2
                //ycool = 3.*k_B*gas.T_K0*gas.y0[2]/(2.*gas.t_ff*Lambda_H2(gas.nH0,gas.T_K0,gas.y0[2])*(mu*m_H));
                //file<<setw(16)<<ycool;
                file<<setw(16)<<0;

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
                itime += 1;
            }
        } */

        if (DM){
            HALO halo(gas.Mh,gas.z);
            double d0 = 200./3.*pow(halo.c,3)/(log(1+halo.c)-halo.c/(1+halo.c));
            if (i==0){
                if (py){
                    file<<setw(16)<<"nc_DM"<<setw(16)<<"z_Mh"<<setw(16)<<"t_Mh"<<setw(16)<<"Mh"<<setw(16)<<"Tvir";
                    file<<setw(16)<<"f_Ma"<<setw(16)<<"cs"<<setw(16)<<"v_tur"<<setw(16)<<"Vc"<<setw(16)<<"v_bsm";
                }
                else file <<" nc_{DM} z_{Mh} t_{Mh} Mh T_{vir} f_{Ma} cs v_{tur} Vc v_{bsm}";
            } //rho_c = d0*RHO_DM(gas.z)/m_H checked right. 
            else{
                file<<setw(16)<<gas.rhoc_DM/(mu*m_H)<<setw(16)<<gas.MPs[gas.iMer].z;
                file<<setw(16)<<gas.MPs[gas.iMer].t/gas.t_ff0<<setw(16)<<gas.Mh/Ms<<setw(16)<<gas.MPs[gas.iMer].Tvir;
                file<<setw(16)<<gas.f_Ma<<setw(16)<<gas.cs/km<<setw(16)<<sqrt(gas.v_tur2)/km<<setw(16)<<halo.Vc/km;
                file<<setw(16)<<gas.v_bsm/km;
            }
        }

        if (haloinfo){
            if (i==0) {
                if (py) file<<setw(16)<<"ievol"<<setw(16)<<"M_J"<<setw(16)<<"M_Jeff"<<setw(16)<<"Mg_intg"<<setw(16)<<"q";
                else file<<" ievol M_J M_{Jeff} Mg_{intg} q";
            }
            else file<<setw(16)<<gas.evol_stage<<setw(16)<<gas.MJ/Ms<<setw(16)<<gas.MJ_eff/Ms<<setw(16)<<gas.Mg_intg/Ms<<setw(16)<<gas.MPs[gas.iMer].mratio;
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
                if (py) file<<setw(16)<<"r_c"<<setw(16)<<"r_cH2"<<setw(16)<<"r_cH"<<setw(16)<<"r_hcompr"<<setw(16)<<"r_hmer";
                else file<<" r_c r_{c,H2} r_{c,H} r_{h,compr} r_{h,mer}";

            }
            else{
                file<<setw(16)<<gas.r_c;
                file<<setw(16)<<gas.r_cH2;
                file<<setw(16)<<gas.r_cH;
                file<<setw(16)<<Gamma_compr(gas.cs,gas.f_Ma,gas.t_ff);
                file<<setw(16)<<gas.Gamma_mer_th; //merger heating; kinetic energy input: *1/3
            }
        }

        if (mer) {
            if (i==0) {
                if (py) file<<setw(16)<<"dMdt"<<setw(16)<<"Gamma_mer"<<setw(16)<<"Gamma_mer_th"<<setw(16)<<"Gamma_mer_k";
                else file<<" dMdt Gamma_{mer} Gamma_{mer,th} Gamma_{mer,k}";
            }
            else file<<setw(16)<<gas.dMdt<<setw(16)<<gas.Gamma_mer<<setw(16)<<gas.Gamma_mer_th<<setw(16)<<gas.Gamma_mer_k;
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
    // printf("1.93/2.2*dTvir = %f K\n",1.89/2.2*(gas.MPs[gas.iMer].Tvir - gas.MPs[0].Tvir));
    printf("J21= %5.2f, z_col=%3.2f\n\n\n\n",J21,gas.z_col);
}

// resembling main_Jc1.cpp
double getT(double& zcol, bool write,int MerMod, double J, double Tb, string treename, bool spec, bool Ma_on, int i_bsm, double nH_tell){
    printf("\n**************************\nSTART in getT:\n");
    GAS gas(frac0,MerMod,J,Tb,treename,spec,Ma_on,i_bsm);

    int index=treename.find("fort.");
    string tree = treename.substr(index+5); // tree_id 输出

    int i=0;
    string fout = "tr"+tree+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+"J"+to_string(int(J))+".txt";
    fstream f1;
    f1.open(fout, ios::out | ios::trunc );
    f1<<setiosflags(ios::scientific)<<setprecision(5);

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
                cout<<"getT: writing file \" "<<fout<<"\"\n";
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
    f1.close();

    zcol = gas.z_col;
    printf("*******IN GET_T***********\n");
    cout<<"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0);
    printf("\tJ=%3.2f, z_col=%3.2f, nH_tell=%3.2e, T=%3.2e\n\n\n",J,zcol,nH_tell, gas.T_K0);
    return gas.T_K0;
}


void evol_Jc(string treename, string fout, double Tb, int MerMod, bool spec, bool Ma_on, int i_bsm){
    printf("############################\t JC_SOL\t ########################\n");
    double T_tell = 4000, nH_tell=1.e4;
    // boundary of bisection J21 
    double J0 = epE2, J1 = 1.e4; //包括所有的Tb所需range 没必要
    J0 = 500, J1 = 1.5e3;
    double T0, T1, T;
    double z0_col,z1_col,z_col;

    ofstream file;
    ifstream checkf_exist(fout.c_str());
    if (checkf_exist.good()) cout<<fout.c_str()<<" exist\n";
    else {
        file.open(fout, ios::out | ios::trunc);
        file<<setw(16)<<"tree"<<setw(16)<<"Tb"<<setw(16)<<"i_bsm"<<setw(16)<<"tur"<<setw(16)<<"z_col"<<setw(16)<<"Jc"<<endl;
        file.close();
    }
    file.open(fout, ios::out | ios::app);

    T0 = getT(z0_col,false,MerMod, J0, Tb, treename, spec, Ma_on, i_bsm,nH_tell); T1 = getT(z1_col,false,MerMod, J1, Tb, treename, spec, Ma_on, i_bsm,nH_tell);
    cout<<"***********log**************\n";
    if ( T0-T_tell>0 ) {
        while(T0-T_tell>0){
            printf("J0=%3.2f,T0=%3.2e,T_tell=%3.2e, wrong INITIAL LEFT boundary\n",J0,T0,T_tell);
            J1 = J0; T1 = T0; z1_col = z0_col;
            J0 *= 0.9;
            T0 = getT(z0_col,false,MerMod, J0, Tb, treename, spec, Ma_on, i_bsm,nH_tell);
        }
    }
    else if (T1-T_tell<0) {
        while (T1-T_tell<0){
            printf("J1=%3.2f,T1=%3.2e,T_tell=%3.2e, wrong INITIAL RIGHT boundary\n",J1,T1,T_tell);
            J0 = J1; T0 = T1; z0_col = z1_col;
            J1 *= 1.1;
            T1 = getT(z1_col,false,MerMod, J1, Tb, treename, spec, Ma_on, i_bsm,nH_tell); 
        }
    }

    while (J1-J0 > 0.01*J0){
        bool write = false;
        // if (J1-J0 > 0.05*J0) write = true;
        T = getT(z_col,write,MerMod, (J0+J1)/2., Tb, treename, spec, Ma_on, i_bsm,nH_tell);
        printf("#######\t########\t########\t#########");
        printf("J0=%3.2f, T0=%3.2f,J1=%3.2f,T1=%3.2f, Jmid=%4.3f\tTmid=%3.2e\n",J0,T0,J1,T1,(J0+J1)/2.,T);
        if (T-T_tell>=0) { J1 = (J0+J1)/2.; T1 = T; z1_col = z_col; }
        else { J0 = (J0+J1)/2.; T0 = T; z0_col = z_col; }
    }

    printf("J0=%3.2f\tT0=%3.2e\tJ1=%3.2f\tT1=%3.2e\n --> Jc_sol=%3.2f z_col=%3.2f\n",J0,T0,J1,T1, J1, z1_col);
    int index=treename.find("fort.");
    string tree = treename.substr(index+5); // tree_id 输出
    file<<setw(16)<<stoi(tree)<<setw(16)<<Tb<<setw(16)<<i_bsm<<setw(16)<< ((Ma_on)?1:0) <<setw(16)<<z1_col<<setw(16)<<J1<<endl;
    file.close();
}