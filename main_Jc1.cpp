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
#include "Newton5.h"
#include "thermo.h"
#include "dyn.h"
#include "PARA.h"
using namespace std;
double min(double a, double b){
    return (a<=b)?a:b;
}
double max(double a, double b){
    return (a>=b)?a:b;
}
//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
int i =0;
int n;

double Tbs[] = {8.e3, 1.e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};
//  double Jcrits[] = {0.8, 20, 800, 1000, 1100, 1100, 1100}; //Sugimura 2014
double J0s[] = {3.12805, 52.002, 2054.69, 3054.69, 3390.62, 3484.38, 3515.62}; //对

//#################################################################################################
//##################### BELOW COPY FROM EVOL. #####################################################
//#################################################################################################
// resembling main_Jc1.cpp
// resembling main_Jc1.cpp
double getT(int MerMod, double J, double Tb, char* treename, bool Ma_on, double nH_tell = 5.e5){
    double frac0[] = {0.,
                        1.-4.e-6-4.5e-3,
                        2.e-6,
                        4.5e-3,
                        4.5e-3+1.e-10,
                        1.e-10};

    GAS gas(frac0,MerMod,J,Tb,treename,Ma_on);
    while (gas.nH0<nH_tell){
        gas.setMerger();
        gas.timescales(); 
        gas.freefall(); 
        gas.react_sol(0); 
        gas.T_sol();
        gas.get_para();
    }
    cout<<J<<"\t"<<gas.T_K0<<endl;
    return gas.T_K0;
}

void evol_Jc(char* treename, char* fout, double Tb, bool Ma_on){
    printf("################################################################################\n");
    int MerMod = 1; //printf("MODEL = %d\n",MerMod);
    double T_tell = 4000;
    //without mergers, WG11 shielding, Jcrit~50 for Tb=1.e4K, Jcrit~3000 for Tb=1.e5K
    double J0, J1;
    double T0, T1, T;
    char fout_malloc[100];
    // sprintf(fname, "../code_tree/fort.217"); 
    sprintf(fout_malloc,fout); 
    ofstream file;
    if (FILE *f = fopen(fout_malloc, "r")) fclose(f);
    else {
        file.open(fout_malloc, ios::out | ios::trunc);
        file<<" Tb Jc Tg"<<endl;
        file.close();
    }
    file.open(fout_malloc, ios::out | ios::app);

    //J0 = 0, J1 = 4000;
    //  WLI put, try for tree 217. check if J~[30,40]
    J0 = 30, J1 = 35;
    //  WLI put, try for tree 217. check if J~2.2*31.2=68
    //J0 = 50, J1 = 60;

    printf("T_tell = %3.2e\n",T_tell);
    T0 = getT(MerMod, J0, Tb, treename, Ma_on); T1 = getT(MerMod, J1, Tb, treename, Ma_on);
    cout<<"***********log**************\n";
    printf("T0 and T1: %3.2e\t%3.2e\n",T0-T_tell,T1-T_tell);
    if ( T0-T_tell>0 || T1-T_tell<0 ) cout<<"wrong INITIAL boundaries"<<endl;

    else{
        while (J1-J0 > 0.1*J0){
            T = getT(MerMod, (J0+J1)/2., Tb, treename, Ma_on);
            printf("#######\t########\t########\t#########");
            printf("J=%4.3f\tT=%3.2e\n",(J0+J1)/2.,T);
            if (T-T_tell>=0) {J1 = (J0+J1)/2.; T1 = getT(MerMod, J1, Tb, treename, Ma_on);}
            else {J0 = (J0+J1)/2.; T0 = getT(MerMod, J0, Tb, treename, Ma_on);}
        }   
    }
    printf("J0=%4.3f\tT0=%3.2e\tJ1=%4.3f\tT1=%3.2e\n",J0,T0,J1,T1);
    cout<<"Tb= "<<Tb<<"Jc_sol= "<<(J0+J1)/2.<<"Tg at nH_tell= "<<"is: "<<T<<endl;
    file<<" "<<Tb<<" "<<(J0+J1)/2.<<" "<<T<<endl;
    file.close();
    
}
//(Sugimura et al. 2014): 8000:0.8; 10000:20; 20000:800; 30000: 1000+; 50000: 1100; 100000:1100; 
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
int main(){
    printf("################################################################################\n");
    n = 7;
    /* char ftree[100];
    char fout[100]; 
    char* name1 = "../code_tree/fort.211";
    char* name2 = "fort.211";
    sprintf(ftree, name1);
    sprintf(fout,"data/temp"); */
    char* ftree = "../code_tree/fort.217"; //不行！！！会
    char* fout = "data/temp";
    double Tb = 1.e4;
    bool Ma_on = true;
    evol_Jc(ftree, fout, Tb, Ma_on);
    return 0;
}

//(Sugimura et al. 2014): 8000:0.8; 10000:20; 20000:800; 30000: 1000+; 50000: 1100; 100000:1100; 