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
    ntree = 1; 

    for (i=0;i<ntree;i++){
        tree = to_string(i); // tree_id 输出
        ftree = "../treefiles/tree_"+to_string(i); // cout<<ftree<<endl;
        fout = "Jcs250.txt";
        cout<<fout<<endl;
        for (i_bsm=0; i_bsm<4; i_bsm++){
            evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
        }
    }

    // i = 0;
    // ftree = "../treefiles/tree_"+to_string(i); cout<<ftree;
    // // ftree = "../tree_Hirano/fort.120"; cout<<ftree;
    // fout = "tr"+tree+"Tb"+to_string(int(Tb))+"evol.txt"; cout<<fout;
    // evol(ftree, fout, MerMod, Tb, J21, spec, Ma_on, i_bsm);

    // Ma_on = false;
    // for (i_bsm=0; i_bsm<1; i_bsm++){
    //     evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
    // }
    // Ma_on = true;
    // for (i_bsm=0; i_bsm<4; i_bsm++){
    //     evol_Jc(ftree,fout,Tb,MerMod,spec,Ma_on,i_bsm);
    // }

    // double z_col, T, nH_tell=1.e4;
    // Ma_on = true;
    // J21 = 656;
    // i_bsm = 3;
    // T = getT(z_col, true, MerMod, J21, Tb, ftree, spec, Ma_on, i_bsm, nH_tell);


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
