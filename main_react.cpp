#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include<algorithm>

#include "reaction.h"
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
int n;
int main(){
    printf("################################################################################\n");
    printf("################################################################################\n");
    printf("################################################################################\n");
    char* treename = "../code_tree/fort.217"; //不行！！！会
    char* fout = "temp.txt";
    bool Ma_on = true;
    int MerMod = 0;
    double J21 = 60; //without mergers, WG11 shielding, Jcrit~50 for Tb=1.e4K, Jcrit~3000 for Tb=1.e5K
    double Tbb = 1.e4;
    double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-6, y_Hm = 1.0e-12, y_H2p = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;

//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
    double frac0[] = {0., y_H, y_H2, y_e, y_Hp, y_Hm};
    GAS gas(frac0,MerMod,J21,Tbb,treename,Ma_on);

    double t_ff0 = 1./C/sqrt(gas.nH0);
    double t1 = 1.9999*t_ff0; //to n = 1.e6/cm^3
    
    ofstream file;
    file.open(fout, ios::out | ios::trunc );
    bool py = true;
    bool DM = true;
    bool fract = true;
    bool react = false;
    bool tscales = true;
    bool haloinfo = true;
    bool heatingcooling = true;
    bool mer = true;
    double kform, rform, yequi, ycool;
    int itime=0;

    i = 0;
    int Ncount = 2000000;
    //while (gas.z>10 && gas.nH0<1.e11){
    for(int icount=0; icount<Ncount; icount++){
        printf("icount= %d\n",icount);
        if (i==0){
            if (py) file <<" t Dt z nH T";
            else file <<" t Dt z nH T";
        }
        else file<<" "<<gas.t_act<<" "<<gas.Dt<<" "<<gas.z<<" "<<gas.nH0<<" "<<gas.T_K0;
        //else file<<" "<<gas.t_act/t_ff0<<" "<<gas.Dt<<" "<<gas.z<<" "<<gas.nH0<<" "<<gas.T_K0;
        
        //cout<<z_ana(35,gas.t_act)<<gas.z<<endl; //not exactly the same...
        if (fract){
            if (i==0){
                if (py) file<<" yH yH2 ye yH+ yH- y_equi y_cool Ma v_tur MJ_eff";
                else    file<<" yH yH2 ye yH+ yH- y_{equi} y_{cool} Ma v_{tur} MJ_{eff}";
            }
            else {
                for (int i=1; i<N_sp+1; i++) file<<" "<<gas.y0[i];

                kform = gas.k[2]*gas.k[3]/(gas.k[3]+gas.k[7]/gas.nH0);
                rform = max(kform*gas.nH0*gas.y0[1]*gas.y0[3], gas.k[5]*pow(gas.nH0,2)*pow(gas.y0[1],3));
                            //H-, n小                              3b,  n>~10^9/cm^3
                yequi = min(rform/gas.k[6], rform/(gas.k[4]*gas.y0[1]*gas.nH0) );
                                // pd,  n<~100/cm^3      cd n大 
                file<<" "<<yequi;
                //file<<" "<<gas.Dt*gas.rf[2]+gas.y0[2]; prediciton of yH2
                ycool = 3.*k_B*gas.T_K0*gas.y0[2]/(2.*gas.t_ff*Lambda_H2(gas.nH0,gas.T_K0,gas.y0)*(mu*m_H));
                file<<" "<<ycool;

                file<<" "<<gas.Ma<<" "<<sqrt(gas.v_tur2);//( k_B*8000/(mu*m_H) ); // cs^2 = k_B*T_K/(mu*m_H)
                file<<" "<<gas.MJ_eff/Ms;

            }
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

        if (i!=0){
            gas.setMerger();
            gas.timescales(); 
            //gas.freefall(); 
            //cout<<"freefall done. "<<"n= "<<gas.nH0<<" /cc"<<endl;
            gas.react_sol(1); 
            //cout<<"react_sol done"<<endl;
            //gas.T_sol();
            //cout<<"T_sol done. "<<"T= "<<gas.T_K0<<" K"<<endl;
            gas.get_para();
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
    printf(" z from 35 to 10 %f Myr\n", (t_from_z(10) - t_from_z(35))/Myr );
    printf(" freefall timescale %f Myr\n", gas.t_ff0/Myr );
    return 0;
}