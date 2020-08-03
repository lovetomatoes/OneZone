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

int main(){
    printf("################################################################################\n");
    // cout<<"dt from z: "<<(t_from_z(20.) - t_from_z(20.2))/Myr<<endl;
    // cout<<"t_ff (n=300/cc): "<<t_freefall(300)/Myr<<endl;

    int i =0;
    int n = 8;
    double Tbs[] = {8.e3, 1.e4, 1.5e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};

    string ftree = "../code_tree/fort.217"; //不行！！！会
    double Tb = 2.e4;
    bool Ma_on = true, spec=false;
    int MerMod = 1;
    double J21 = 0;
    int i_bsm;
    string fout, tree;
    int ntree = 10;
    ntree = 3; 

    clock_t t0 = clock();
    // for (i=2;i<ntree;i++){
    //     tree = to_string(i); // tree_id 输出
    //     ftree = "../treefiles/tree_"+to_string(i); // cout<<ftree<<endl;
    //     fout = "Jcs250.txt";
    //     cout<<fout<<endl;
    //     for (i_bsm=2; i_bsm<3; i_bsm++){
    //         evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
    //     }
    // }

    clock_t t1 = clock();
    printf("time taken 1st part:%.2f\n", (double)(t1-t0)/CLOCKS_PER_SEC);

    i = 2; tree = to_string(i);
    ftree = "../treefiles/tree_"+to_string(i); cout<<ftree<<endl;

    i_bsm = 2; J21 = 562.;
    fout = "tr"+tree+"bsm"+to_string(i_bsm)+"J"+to_string(int(J21))+"_6.txt";
    evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

    clock_t t2 = clock();
    printf("time taken 2nd part:%.2f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    // i_bsm = 2; J21 = 753.;
    // fout = "tr"+tree+"bsm"+to_string(i_bsm)+"J"+to_string(int(J21))+"_1.txt";
    // // fout = "Hm0rcb.txt";
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

    // clock_t t3 = clock();
    // printf("time taken 3rd part:%.2f\n", (double)(t3-t2)/CLOCKS_PER_SEC);



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
