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

    string ftree = "../code_tree/fort.217"; //不行！！！会
    double Tb = 2.e4;
    bool Ma_on = true, spec=false;
    int MerMod = 1;
    double J21 = 0;
    int i_bsm;
    string fout;
    Ma_on = true;
    ftree = "../tree_Hirano/fort.20";

    i_bsm = 0; Ma_on = false;

    // double T0,z_col;
    // J21 = 1304;
    // fout = "J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // printf("--\n--\n--\n--\n--\n--\n");
    // J21 = 1312;
    // fout = "J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

    // fout = "Jc_bsm" + to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    // evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
//   //void evol_Jc(string treename, string fout, double Tb, int MerMod, bool spec, bool Ma_on, int i_bsm){

    Ma_on = false;
    J21 = 0;
    for (i_bsm=0; i_bsm<4; i_bsm++){
        string fout = "./evols/mer2J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
        evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    }
    J21 = 10;
    for (i_bsm=0; i_bsm<4; i_bsm++){
        string fout = "./evols/mer2J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
        evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    }

    Ma_on = true;
    J21 = 0;
    for (i_bsm=0; i_bsm<4; i_bsm++){
        string fout = "./evols/mer2J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
        evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    }
    J21 = 10;
    for (i_bsm=0; i_bsm<4; i_bsm++){
        string fout = "./evols/mer2J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
        evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    }

    // double T0,z_col,nH_tell=1.e4;
    // Ma_on = false;
    // J21 = 0;
    // for (i_bsm=0; i_bsm<4; i_bsm++){
    //     string fout = "./evols/mer10J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    //     evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // }
    // J21 = 10;
    // for (i_bsm=0; i_bsm<4; i_bsm++){
    //     string fout = "./evols/mer10J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    //     evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // }
    // Ma_on = true;
    // J21 = 0;
    // for (i_bsm=0; i_bsm<4; i_bsm++){
    //     string fout = "./evols/mer10J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    //     evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // }
    // J21 = 10;
    // for (i_bsm=0; i_bsm<4; i_bsm++){
    //     string fout = "./evols/mer10J"+to_string(int(J21))+"_bsm"+ to_string(i_bsm) + "tur"+ to_string((Ma_on)?1:0)+".txt";
    //     evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // }



    // i_bsm = 0;
    // fout = "20J1e1bsm0tur0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 1;
    // fout = "20J1e1bsm1tur0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 2;
    // fout = "20J1e1bsm2tur0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 3;
    // fout = "20J1e1bsm3tur0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

// bsm; turbulence
    // Ma_on = true;
    // i_bsm = 0;
    // fout = "20J1e1bsm0evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 1;
    // fout = "20J1e1bsm1evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 2;
    // fout = "20J1e1bsm2evolve_Vc.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);
    // i_bsm = 3;
    // fout = "20J1e1bsm3evolve_Vc.txt";
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