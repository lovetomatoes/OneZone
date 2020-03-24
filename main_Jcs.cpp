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
#include <string> // for std::string s

//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
int MerMod = 1;

double T_tell = 4000;
double nH_tell = 5.e5;

//std::string s = std::to_string(MerMod);
double J0, J1;
double T0, T1, T;

double Tb, Tbs[] = {8.e3, 1.e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};
//  double Jcrits[] = {0.8, 20, 800, 1000, 1100, 1100, 1100}; //Sugimura 2014
double J0s[] = {3.12805, 52.002, 2054.69, 3054.69, 3390.62, 3484.38, 3515.62}; //对
/* // for arbitrary merger models 1 & 2
double J1s[] = {20.813, 128.418, 1246.09, 1167.97, 1175.78, 1175.78, 1175.78}; //right 
double J2s[] = {13.2446, 124.512, 1695.31, 1339.84, 2726.56, 2257.81, 2132.81}; */
/*     
// added merger heating from a tree
double J1s[] = {0.28, 5.1, 578, 968, 1093, 1156, 1156};
 */
int i=0;
int n;


double frac0[] = {0.,
                    1.-4.e-6-4.5e-3,
                    2.e-6,
                    4.5e-3,
                    4.5e-3+1.e-10,
                    1.e-10};
    
double getT(int MerMod,double J,double Tb){
    GAS gas(frac0,MerMod,J,Tb);
    while (gas.t0<gas.t1 && gas.nH0<nH_tell){
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

int main(){
    printf("################################################################################\n");

    printf("MODEL = %d\n",MerMod);
    n = 7;

    ofstream file;
    //file.open("data/T-Jcrit.txt", ios::out | ios::trunc );

    for(i=0;i<n;i++){
        cout<<"\n\ni= "<<i<<endl;
        Tb = Tbs[i];
        J0 = 0, J1 = 4000;
        T0 = getT(MerMod,J0,Tb); T1 = getT(MerMod,J1,Tb);
        cout<<"***********log**************\n";
        if ( T0-T_tell>0 | T1-T_tell<0 ) {cout<<"wrong INITIAL boundaries"<<endl; break;}
        else{
            T = getT(MerMod,(J0+J1)/2.,Tb);
            while (J1-J0 > 0.1*J0){
                T = getT(MerMod,(J0+J1)/2.,Tb);
                if (T-T_tell>=0) {J1 = (J0+J1)/2.; T1 = getT(MerMod,J1,Tb);}
                else {J0 = (J0+J1)/2.; T0 = getT(MerMod,J0,Tb);}
            }
            /* if (i==0) file<<"\tTb"<<"\tJ_{crit}"<<" \tJ_{"<<s<<"}"<<endl;
            if (i==0) file<<"\tTb"<<"\tJ_{0}"<<" \tJ_{1}"<<" \tJ_{2}"<<endl;
            file<<" "<<Tbs[i]<<" "<<J0s[i]<<" "<<J1s[i]<<" "<<J2s[i]<<endl;
 */
        }
        cout<<"****************************\n";
        cout<<" "<<Tbs[i]<<" "<<(J0+J1)/2.<<" "<<T<<endl;
        file<<" "<<Tbs[i]<<" "<<(J0+J1)/2.<<" "<<T<<endl;

    }
    cout<<" "<<"Tb"<<" "<<"(J0+J1)/2."<<" "<<"T"<<endl;
    //file<<" "<<"Tb"<<" "<<"(J0+J1)/2."<<" "<<"T"<<endl;
    //file.close();
    return 0;
}

//(Sugimura et al. 2014): 8000:0.8; 10000:20; 20000:800; 30000: 1000+; 50000: 1100; 100000:1100; 
