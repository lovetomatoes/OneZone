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
static int i_bsm = 1, MerMod = 1, itr; 
static int const ntree=3;
static double Tb = 2.e4;
static bool Ma_on = true, spec=false;
static clock_t t0, t1, t2;

int main(){
    printf("################################################################################\n");

    string ftree;
    double J21 = 0.;
    string fout, tree;
    string Jzname;
    int a[ntree] = {4479,14,23};

    clock_t t0 = clock();

    for (i=0;i<ntree;i++){
        itr = a[i];
        ftree = "../treefiles/tree_"+to_string(itr);
        fout = "./bsm"+to_string(i_bsm)+to_string(itr)+"evol.txt";
        Jzname = "../Jpoints/tree_"+to_string(itr)+"_Jtrack"; // delta_z = 5
        cout<<ftree<<endl;
        evol(ftree, Jzname, fout, MerMod, Tb, 0, spec, Ma_on, i_bsm);
    }
    clock_t t1 = clock();
    printf("dt=%.2f\n", (double)(t1-t0)/CLOCKS_PER_SEC);

    return 0;
}
