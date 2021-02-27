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
static int i_bsm = 0, ntree=100, MerMod = 1,  itr; 
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
    for (itr=0; itr<10; itr++){
        ftree = "../treefiles/tree_"+to_string(itr);
        Jzname = "../Jpoints/tree_"+to_string(itr)+"_Jtrack";
        cout<<ftree<<endl;
        cout<<Jzname<<endl;
        // evol_Jtrack(ftree, Jzname, fout, Tb, MerMod, spec, Ma_on, i_bsm);
    }

    return 0;
}
