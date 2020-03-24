// include 1 merger tree
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include "reaction.h"
#include "class_gas.h"
#include "class_halo.h"
#include "my_linalg.h"
#include "Newton5.h"
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
    double T_K0=21; double nH0=4.5e-3;

    double z0 = 20;
    int MerMod = 0;
    double J21 = 0.; //without mergers, WG11 shielding, Jcrit~50 for Tb=1.e4K, Jcrit~3000 for Tb=1.e5K
    double Tbb = 1.e4;
    double t_ff0 = sqrt(1/nH0)/C;

    //double t1 = 1.999999*t_ff0; //to n = 1.e10/cm^3 (when n > 1.e9 /cm^3, 3body starts to act prominently
    double t1 = 1.9999*t_ff0; //to n = 1.e6/cm^3
    printf("t_ff = %e Myr\n",t_ff0/Myr);
    printf("t1 in main is :  %f t_ff\n", t1/t_ff0);

    double frac0[] = {0.,
                      1.-4.e-6-4.5e-3,
                      2.e-6,
                      4.5e-3,
                      4.5e-3+1.e-10,
                      1.e-10};
    GAS gas(z0,t1,T_K0,nH0,frac0,MerMod,J21,Tbb);
    ofstream f1;
    f1.open("data/atree.txt", ios::out | ios::trunc );
    for (i=1;i<gas.nMer;i++){
        if (i==1) f1<<"\tj\tt(Myr)\tdt(Myr)\tmhalo(10^8Ms)\tz\tTvir\n";
        f1<<gas.MPs[i].j<<" "<<gas.MPs[i].t/Myr<<" "<<gas.MPs[i].dt/Myr<<" "<<gas.MPs[i].mhalo/(1.e8*Ms)<<" "<<gas.MPs[i].z<<" "<<gas.MPs[i].Tvir<<endl;
    }
    f1.close();

    for (i=1; i<gas.nMer-1; i++){
        printf("time = %6.3e,\tdt = %6.3e\n", gas.MPs[i].t/Myr, (gas.MPs[i].t-gas.MPs[i-1].t)/Myr);
    }
    ofstream file;
    file.open("data/evolve.txt", ios::out | ios::trunc );

    bool fract = true;
    bool react = false;
    bool tscales = true;
    bool heatingcooling = true;
    bool dyn = false;
    double kform, rform, yequi, ycool;
    int itime=0;
    printf("R_J is %3.2e pc for n = %3.2e and T = %3.e K\n",gas.RJ/pc,gas.nH0,gas.T_K0);
    printf("M_J is %3.2e Ms for n = %3.2e and T = %3.e K\n",gas.MJ0/Ms,gas.nH0,gas.T_K0);

    i = 0;
    while (gas.t0<gas.t1 && gas.nH0<1.e12){
        if (i==0) file <<"\tt_{act} \tDt \tz \tnH \tT \trho_DM";    
        else file<<" "<<gas.t_act/t_ff0<<" "<<gas.Dt<<" "<<z_ana(z0,gas.t_act)<<" "<<gas.nH0<<" "<<gas.T_K0<<" "<<RHO_DM(gas.z0);

        gas.setMerger();
        gas.timescales(0); 
        gas.freefall();
        gas.react_sol(1); 
        gas.T_sol();
        gas.get_para();
        
        file<<endl;
        i++;
    }
    file.close();

    cout<<"ff timescale nH=1.0 is : "<<t_freefall(1.)/Myr<<" Myr"<<endl;
    cout<<"ff timescale nH = 0.1 is : "<<t_freefall(.1)/Myr<<" Myr"<<endl;
    cout<<"in Myr: Dt_m="<<gas.Dt_m/Myr<<"\tdt_m="<<gas.dt_m/Myr<<"\tt_delay="<<gas.t_delay/Myr<<endl;
    printf(" Merger: addMh = %e Ms in 20Myr\n",gas.addMh/Ms);
    printf("R_J is %3.2e pc for n = %3.2e and T = %3.e K\n",gas.RJ/pc,gas.nH0,gas.T_K0);
    printf("evolve time is %3.2e Myr\n",gas.t_act/Myr);

    return 0;
}