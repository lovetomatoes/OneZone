#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>

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


//#################################################################################################
//##################### BELOW COPY FROM EVOL. #####################################################
//#################################################################################################
double getT(int MerMod, double J, double Tb, char* treename, bool Ma_on, double nH_tell = 5.e5){
    double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-6, y_Hm = 1.0e-12, y_H2p = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;

    double frac0[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};

    GAS gas(frac0,MerMod,J,Tb,treename,Ma_on);
    while (gas.nH0<nH_tell){
        gas.setMerger();
        gas.timescales(); 
        gas.freefall(); 
        gas.react_sol(1); 
        gas.T_sol();
        gas.get_para();
    }
    return gas.T_K0;
}

void evol_Jc(char* treename, char* fout, double Tb, int MerMod, bool Ma_on){
    printf("################################################################################\n");
    printf("f_Ma is %d\n",(Ma_on)?1:0);
    double T_tell = 4000;
    // boundary of bisection J21
    double J0 = epE6, J1 = 1.e4;
    double T0, T1, T;
    char fout_malloc[100];
    sprintf(fout_malloc,fout); 
    ofstream file;
    if (FILE *f = fopen(fout_malloc, "r")) fclose(f);
    else {
        file.open(fout_malloc, ios::out | ios::trunc);
        file<<" Tb Jc"<<endl;
        file.close();
    }
    file.open(fout_malloc, ios::out | ios::app);


    T0 = getT(MerMod, J0, Tb, treename, Ma_on); T1 = getT(MerMod, J1, Tb, treename, Ma_on);
    cout<<"***********log**************\n";
    printf("T0 and T1: %3.2e\t%3.2e\n",T0-T_tell,T1-T_tell);
    if ( T0-T_tell>0 or T1-T_tell<0 ) cout<<"wrong INITIAL boundaries"<<endl;

    else{
        while (J1-J0 > 0.05*J0){
            T = getT(MerMod, (J0+J1)/2., Tb, treename, Ma_on);
            /* printf("#######\t########\t########\t#########");
            printf("J=%4.3f\tT=%3.2e\n",(J0+J1)/2.,T); */
            if (T-T_tell>=0) {J1 = (J0+J1)/2.; T1 = getT(MerMod, J1, Tb, treename, Ma_on);}
            else {J0 = (J0+J1)/2.; T0 = getT(MerMod, J0, Tb, treename, Ma_on);}
        }   
    }
    printf("J0=%4.3f\tT0=%3.2e\tJ1=%4.3f\tT1=%3.2e\n --> Jc_sol=%3.2e\n",J0,T0,J1,T1, (J0+J1)/2.);
    file<<" "<<Tb<<" "<<(J0+J1)/2.<<endl;
    file.close();
    
}

//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
int main(){
    printf("################################################################################\n");
    /* 
    char fout[100]; 
    sprintf(fout,"data/temp"); */
    int i =0;
    int n = 7;
    double Tbs[] = {8.e3, 1.e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};

    char* ftree = "../code_tree/fort.217"; //不行！！！会
    char* fout = "Jcs_cpp.txt";
    double Tb = 8.e3;
    bool Ma_on = true;
    int MerMod = 0;

    //evol_Jc(ftree, fout, Tb, MerMod, Ma_on);
    for (i=0;i<n;i++) evol_Jc(ftree, fout, Tbs[i], MerMod, Ma_on);

    return 0;
}