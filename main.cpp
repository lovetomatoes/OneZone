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
    char* fout = "Jcs_cpp.txt";
    double Tb = 2.e4;
    bool Ma_on = true;
    int MerMod = 0;
    double J21 = 50;

    fout = "evolve.txt";
    evol(ftree, fout, MerMod, Tb, J21, Ma_on);

    // fout = "Jcs_cpp.txt";
    // for (i=0;i<n;i++) evol_Jc(ftree, fout, Tbs[i], MerMod, Ma_on);

    return 0;
}