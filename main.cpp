#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>

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
using namespace std;
#include <time.h>

static int i,j;
static int ntree = 10, MerMod = 1, i_bsm, itr; 
static double Tb = 2.e4;
static bool Ma_on = true, spec=false;
static clock_t t0, t1, t2;

static double evol_tiny = 1.0e-20, evol_yHe = 8.33333e-2, evol_y_H2 = 1.0e-5, evol_y_Hm = 1.0e-12, evol_y_H2p = 1.0e-12;
static double evol_y_Hp = 1.0e-4;
static double evol_y_H = 1.0 - 2.*evol_y_H2 - 2.*evol_y_H2p - evol_y_Hm - evol_y_Hp;
static double evol_y_He = evol_yHe - 2.*evol_tiny, evol_y_Hep = evol_tiny, evol_y_Hepp = evol_tiny;
static double evol_y_e = evol_y_Hp + evol_y_H2p - evol_y_Hm + evol_y_Hep + 2.*evol_y_Hepp;

static double frac0[] = {0., evol_y_H, evol_y_H2, evol_y_e, evol_y_Hp, evol_y_H2p, evol_y_Hm, evol_y_He, evol_y_Hep, evol_y_Hepp};

int main(){
    printf("################################################################################\n");

    string ftree = "../treefiles/tree_0"; //不行！！！会
    string fJtrack = "../treefiles/tree_0_Jtrack"; //不行！！！会
    double Tb = 2.e4;
    bool Ma_on = true, spec=false;
    int MerMod = 1;
    double J21 = 1000;
    int i_bsm;
    string fout, tree;
    int ntree = 10;
    ntree = 1; 
    string Jzname;

    t0 = clock();
    itr = 3460;
    tree = to_string(itr); // tree_id 输出
    ftree = "../treefiles/tree_"+to_string(itr); // cout<<ftree<<endl;
    Jzname = "../Jpoints/tree_"+to_string(itr)+"_Jtrack";

    // fout = "tree"+to_string(itr)+"Jtrack.txt";
    fout = "tree"+to_string(itr)+"J1e"+to_string(int(log10(J21)))+".txt";
    for (i_bsm=0; i_bsm<1; i_bsm++){
        // evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
        evol(ftree, Jzname, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
        // GAS gas(frac0,MerMod,J21,Tb,ftree, Jzname, spec,Ma_on,i_bsm);
    }


    t1 = clock();
    printf("time taken 1st part:%.2f \n", (double)(t1-t0)/CLOCKS_PER_SEC);

    // i = 2; tree = to_string(i);
    // ftree = "../treefiles/tree_"+to_string(i); 
    // cout<<ftree<<endl;

    // i_bsm = 2; J21 = 562.;
    // fout = "tr"+tree+"bsm"+to_string(i_bsm)+"J"+to_string(int(J21))+"_222.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

    t2 = clock();
    printf("time taken 2nd part:%.2f \n", (double)(t2-t1)/CLOCKS_PER_SEC);

// spec
    // fout = "Jcs_cpp.txt";
    // spec = false;
    // for (i=0;i<n;i++) evol_Jc(ftree, fout, Tbs[i], MerMod, spec, Ma_on);

    // fout = "Jcs_spec.txt";
    // spec = true;
    // evol_Jc(ftree, fout, Tbs[i], MerMod, true, Ma_on);

    return 0;
}
