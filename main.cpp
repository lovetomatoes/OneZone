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

int main(){
    printf("################################################################################\n");
    /* 
    char fout[100]; 
    sprintf(fout,"data/temp"); */
    int i =0;
    int n = 8;
    double Tbs[] = {8.e3, 1.e4, 1.5e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};

    char* ftree = "../code_tree/fort.217"; //不行！！！会
    char* fout = NULL;
    double Tb = 2.e4;
    bool Ma_on = true, spec=false;
    int MerMod = 1;
    double J21 = 5000;
    int i_bsm;

    ftree = "../tree_Hirano/fort.20";
    // i_bsm = 0;
    // fout = "20bsm0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 1;
    // fout = "20bsm1evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 2;
    // fout = "20bsm2evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 3;
    // fout = "20bsm3evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    Ma_on = false;
    i_bsm = 0;
    fout = "20bsm0tur0evolve_Vc.txt";
    evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

// bsm; turbulence
    // i_bsm = 0;
    // fout = "217bsm0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 1;
    // fout = "217bsm1evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 2;
    // fout = "217bsm2evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 3;
    // fout = "217bsm3evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // Ma_on = false;
    // i_bsm = 0;
    // fout = "217bsm0tur0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

// spec
    // fout = "Jcs_cpp.txt";
    // spec = false;
    // for (i=0;i<n;i++) evol_Jc(ftree, fout, Tbs[i], MerMod, spec, Ma_on);

    // fout = "Jcs_spec.txt";
    // spec = true;
    // evol_Jc(ftree, fout, Tbs[i], MerMod, true, Ma_on);

    return 0;
}