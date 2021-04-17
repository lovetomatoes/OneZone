#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "evol.h"
#include "reaction.h"
#include "class_gas.h"
#include "LE_iso.h"
#include "class_halo.h"
#include "my_linalg.h"
#include "Newton.h"
#include "thermo.h"
#include "dyn.h"
#include "PARA.h"
#include "kpd.h"
using namespace std;
#include <time.h>

static int i,j;
static int i_bsm = 1, ntree=6, MerMod = 1,  itr; 
static double Tb = 2.e4;
static bool Ma_on = true, spec=false;
static clock_t t0, t1, t2;

int main(){
    printf("################################################################################\n");

    string ftree;
    double J21 = 0.;
    string fout, tree;
    string Jzname;
    int a[ntree] = {0,1,2,4479,14,23};
    fout = "./J_col"+to_string(i_bsm)+".txt";

    int Ni = 3, Nd = 15;
    int vi[Ni]; double vd[Nd];
    clock_t t0 = clock();

    for (i=0;i<ntree;i++){
        itr = a[i];
        ftree = "../treefiles/tree_"+to_string(itr);
        cout<<ftree<<endl;
        Jzname = "../Jpoints/tree_"+to_string(itr)+"_Jtrack"; // delta_z = 5

        // fout = "tr"+to_string(itr)+"_bsm"+to_string(i_bsm)+"Jdz5.txt";
        // evol(ftree, Jzname, fout, MerMod, Tb, J21,  spec, Ma_on, i_bsm);

        ofstream file;
        file<<setiosflags(ios::scientific)<<setprecision(3);
        ifstream checkf_exist(fout);
        if (checkf_exist.good()) {}
        else {
            file.open(fout, ios::out | ios::trunc);
            file<<setw(12)<<"tree"<<setw(12)<<"i_bsm"<<setw(12)<<"iso_col";
            file<<setw(12)<<"z_col"<<setw(12)<<"Mh_col"<<setw(12)<<"ng_col";
            file<<setw(12)<<"Tg_col"<<setw(12)<<"Mg_intg"<<setw(12)<<"MJ_col";
            file<<setw(12)<<"J_col"<<setw(12)<<"f_col"<<setw(12)<<"ng_loi";
            file<<setw(12)<<"Tg_loi"<<setw(12)<<"MJ_loi"<<setw(12)<<"f_loi";
            file<<setw(12)<<"Tg_1e4"<<setw(12)<<"f_1e4"<<setw(12)<<"Tg_max";
            file<<endl;
            file.close();
        }


        file.open(fout, ios::out | ios::app);

        evol_Jtrack(vi, vd, ftree, Jzname, Tb, MerMod, spec, Ma_on, i_bsm);

        for (int ii=0;ii<Ni;ii++) file<<setw(12)<<vi[ii];
        for (int id=0;id<Nd;id++) file<<setw(12)<<vd[id];
        file<<endl;
        file.close();
    }
    clock_t t1 = clock();
    printf("dt=%.2f\n", (double)(t1-t0)/CLOCKS_PER_SEC);

    return 0;
}
