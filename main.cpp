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
#include "kpd.h"
using namespace std;
#include <time.h>

static int i,j;
static int i_bsm = 1, ntree=10, MerMod = 1,  itr; 
static double Tb = 2.e4;
static bool Ma_on = true, spec=false;
static clock_t t0, t1, t2;

int main(){
    printf("################################################################################\n");

    string ftree;
    double J21 = 0.;
    string fout, tree;
    string Jzname;

    fout = "J_col"+to_string(i_bsm)+".txt";
    // int a[18] = {168,570,764,1186,2277,2669,2709,2797,3042,3195,4900,5062,7088,7242,7308,7586,8038,9503};

    clock_t t0 = clock();
    for (itr=0;itr<3;itr++){
        ftree = "../treefiles/tree_"+to_string(itr);
        cout<<ftree<<endl;

        // Jzname = "../Jpoints/tree_"+to_string(itr)+"_Jtrack"; // delta_z = 5
        // fout = "tr"+to_string(itr)+"_bsm"+to_string(i_bsm)+"Jdz5.txt";
        // evol(ftree, Jzname, fout, MerMod, Tb, 700.,  spec, Ma_on, i_bsm);

        Jzname = "../P_JLW/tree_"+to_string(itr)+"_Jtrack"; // delta_z = 1
        fout = "tr"+to_string(itr)+"_bsm"+to_string(i_bsm)+"Jdz1.txt";
        evol(ftree, Jzname, fout, MerMod, Tb, 700.,  spec, Ma_on, i_bsm);

        // evol_Jtrack(ftree, Jzname, fout, Tb, MerMod, spec, Ma_on, i_bsm);
    }
    clock_t t1 = clock();
    printf("dt=%.2f\n", (double)(t1-t0)/CLOCKS_PER_SEC);

    return 0;
}
