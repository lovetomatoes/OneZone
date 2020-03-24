#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include<algorithm>

#include "reaction.h"
#include "my_linalg.h"
#include "Newton.h"
#include "PARA.h"
using namespace std;

/* g++ -Wall -I/usr/local/include -c check_react.cpp 
g++  -L/usr/local/lib check_react.o Newton.o PARA.o my_linalg.o gsl_inverse.o reaction.o -lgsl -lgslcblas -lm -o check_react
./check_react
*/

//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************
int i =0, n;
int const N = N_sp+1;
int main(){
    for (i=0;i<3;i++)printf("################################################################################\n");
    ofstream myfile;
    myfile.open("../data/react_fine.txt", ios::out | ios::trunc );
    myfile<<" t yH yH2 ye yH+ yH2+ yH- yHe yHe+ yHe++\n";

    double J_LW = 60; //without mergers, WG11 shielding, Jcrit~50 for Tb=1.e4K, Jcrit~3000 for Tb=1.e5K
    double Tb = 1.e4;
    double tiny = 1.0e-20, yHe = 8.33333e-2, y_H2 = 1.0e-6, y_Hm = 1.0e-12, y_H2p = 1.0e-12;
    double y_Hp = 1.0e-4, y_H = 1.0 - 2.*y_H2 - 2.*y_H2p - y_Hm - y_Hp;
    double y_He = yHe - 2.*tiny, y_Hep = tiny, y_Hepp = tiny;
    double y_e = y_Hp + y_H2p - y_Hm + y_Hep + 2.*y_Hepp;
    double nH = 1.e8;
    double T_K = 8.e3;
    //* previous:     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *

    //* 1: H    2: H2    3: e    4: H+    5: H2+    6: H-    7: He    8: He+    9: He++ *

    double frac0[] = {0., y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp};

    //double t_ff0 = 1./C/sqrt(nH);
    double t_ff = 1.e15/sqrt(nH) * 1.e2;
    double dt = 1.e-3;
    double t0 = 0;
    
    int Ncount = 20000000, icount;
    int i=0, isp, jsp, iter0;
    double y0[N], y1[N], dy[N];

    // initially setting y0
    for (isp=0;isp<N;isp++) y0[isp] = frac0[isp];

    for(icount=0; icount<Ncount; icount++){
        //printf("icount= %d\n",icount); 
        //printf("time=%3.2e, ",t0);       
        if (t0>1.e19) break;
        if(dt<t_ff) dt *= 1.01;
        else dt = t_ff;
        
        for (isp=0; isp<N; isp++){
            dy[isp] = y0[isp]; y1[isp] = y0[isp];// checked right
            //printf("y0[%d]=%3.2e, y1=%3.2e, dy=%3.2e\n",isp,y0[isp],y1[isp],dy[isp]);
        }   
        iter0 = 0;
        while ( len_v(N, dy) > epE8*len_v(N, y1) ){
            //SOL_IMPLICIT(dy, y0, y1, dt, nH, T_K, J_LW, Tb);
            SOL_IMPLICIT(dy, y0, y1, dt, nH, T_K, 0, Tb);
        }

    myfile<<" "<<t0;
    for (isp=1; isp<N_sp+1; isp++) myfile<<" "<<y0[isp]; myfile<<endl;

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
    y0[1] = y_H;
    y0[2] = y_H2;
    y0[3] = y_e;
    y0[4] = y_Hp;
    y0[5] = y_H2p;
    y0[6] = y_Hm;
    y0[7] = y_He;
    y0[8] = y_Hep;
    y0[9] = y_Hepp;

    
    //for (isp=1;isp<N;isp++) y0[isp] = y1[isp];
    t0 += dt;
    }
    
    myfile.close();

    printf("n = %3.2e and T = %3.e K\n",nH,T_K);
    printf("t0 = %3.2e\n",t0);
    return 0;
}