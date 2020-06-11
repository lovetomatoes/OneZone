#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>

#include "class_halo.h"
#include "read_aTree.h"
#include "LE_adb.h"
#include "RK4.h"
#include "PARA.h"
#include "dyn.h"
#include <cmath>
using namespace std;
/* 
g++ read_aTree.cpp dyn.cpp PARA.cpp class_halo.o -o read.out && ./read.out
 */


static int i, num;
static double a;

void aTree(int& nMP, string treename, MainProgenitor* MPs){
    string filename = treename + "mer";
    string line, temp_str;
    ifstream readmer, readtree;
    readmer.open(filename.c_str());
    i = 0;

    if (readmer.good()){
        // cout<<filename<<"TREEmer exist\n";
        getline(readmer,temp_str);
        i = 0;
        while (getline(readmer,line)){
            stringstream ss(line);
            ss>>MPs[i].id_tree>>MPs[i].j>>MPs[i].mhalo>>MPs[i].z>>MPs[i].t>>MPs[i].c>>MPs[i].Tvir>>MPs[i].ng_adb>>MPs[i].mratio;
            MPs[i].t *= Myr; //time from universe age of 0
            MPs[i].mhalo *= Ms; // mhalo unit -> g
            // printf("IN READMER:\n MPs[%d].j=%d,mhalo=%3.2e,z=%3.2f,Tvir=%3.2e,c=%3.2e,n_adb=%3.2e\n",i,MPs[i].j,MPs[i].mhalo/(1.e11*Ms),MPs[i].z,MPs[i].Tvir,MPs[i].c,MPs[i].ng_adb );
            i++;
        }
        readmer.close();
        nMP = i;
        // printf("READMER:%d\n",nMP);
    }
// no mer file; calculate n_adb; write
    else{
        cout<<" no mer file\n";
        readtree.open(treename);
        int id_m;
        getline(readtree,temp_str);
        i = 0;
        while (getline(readtree,line)){
            // cout<<line;
            stringstream ss(line);
            ss>>MPs[i].j>>id_m>>MPs[i].mhalo>>MPs[i].z>>MPs[i].Tvir>>MPs[i].c>>MPs[i].id_tree; //current, for tree_Hirano/. file
            // ss>>MPs[i].j>>id_m>>MPs[i].mhalo>>MPs[i].z>>MPs[i].Tvir>>MPs[i].id_tree; //old, for code_tree/. files

            MPs[i].t = t_from_z(MPs[i].z); //time from universe age of 0
            MPs[i].mhalo *= (1.e11*Ms); // mhalo unit -> g
            // printf("IN READTREE:\n MPs[%d].j=%d,mhalo=%3.2e,z=%3.2f,Tvir=%3.2e,c=%3.2e,id_tree=%d \n",i,MPs[i].j,MPs[i].mhalo/(1.e11*Ms),MPs[i].z,MPs[i].Tvir,MPs[i].c,MPs[i].id_tree );
            i++;
        }
        readtree.close();
        nMP = i-1; //减1是因为analysis_trees.f 输出给fort.11的最后一行带\n
        // cout<<"tree: nMP= "<<nMP<<endl;

        // fort.211, fort.212得用mac来跑, haha上面报错segmentation fault

    // 调换位置 以num/2为轴 MPs[num]<->MPs[0]
        num = nMP - 1;
        double temp_dbl;
        // 原本到MPs[0]无效, 因为fort.11中第一行header: # j id_m mhalo_m/1.e11Ms z_m T_{vir,m} 赋给了MPs[0]
        for (i=0;i<num/2;i++){
            temp_dbl = MPs[num-i].z; MPs[num-i].z = MPs[i].z; MPs[i].z = temp_dbl;
            temp_dbl = MPs[num-i].t; MPs[num-i].t = MPs[i].t; MPs[i].t = temp_dbl;
            temp_dbl = MPs[num-i].mhalo; MPs[num-i].mhalo = MPs[i].mhalo; MPs[i].mhalo = temp_dbl;
            temp_dbl = MPs[num-i].c; MPs[num-i].c = MPs[i].c; MPs[i].c = temp_dbl;
            temp_dbl = MPs[num-i].Tvir; MPs[num-i].Tvir = MPs[i].Tvir; MPs[i].Tvir = temp_dbl;
            }
    // reset t starting point to the 1st progenitor; n_adb
        temp_dbl = MPs[0].t;
        for (i=0;i<nMP;i++){ 
            MPs[i].t -= temp_dbl; 
            HALO halo(MPs[i].mhalo,MPs[i].z);
            double ni=1; MPs[i].ng_adb = 1;
            Mg2N0_adb(MPs[i].ng_adb,ni,MPs[i].z,MPs[i].mhalo);
            MPs[i].dm = (i==nMP-1)? 0: MPs[i+1].mhalo - MPs[i].mhalo;
            MPs[i].mratio = (i==nMP-1)? 0:MPs[i].dm/MPs[i].mhalo;
        }

        ofstream f;
        f.open(filename, ios::out | ios::trunc);
        f<<setw(17)<<"id_tree"<<setw(17)<<"j"<<setw(17)<<"Mh_Ms";
        f<<setw(17)<<"z"<<setw(17)<<"t_Myr"<<setw(17)<<"c";
        f<<setw(17)<<"Tvir"<<setw(17)<<"n_adb"<<setw(17)<<"q"<<endl;
        f<<setiosflags(ios::scientific)<<setprecision(8);
        for (i=0;i<nMP;i++){
            f<<setw(17)<<MPs[i].id_tree<<setw(17)<<MPs[i].j<<setw(17)<<MPs[i].mhalo/Ms;
            f<<setw(17)<<MPs[i].z<<setw(17)<<MPs[i].t/Myr<<setw(17)<<MPs[i].c;
            f<<setw(17)<<MPs[i].Tvir<<setw(17)<<MPs[i].ng_adb<<setw(17)<<MPs[i].mratio<<endl;
        }
        f.close();
    }
// other attributes
    for (i=0;i<nMP;i++){
            HALO halo(MPs[i].mhalo,MPs[i].z);

            MPs[i].Tvir = halo.Tvir;
            MPs[i].c = halo.c;

            MPs[i].dt = (i==nMP-1)? 0: MPs[i+1].t - MPs[i].t;
            MPs[i].dm = (i==nMP-1)? 0: MPs[i+1].mhalo - MPs[i].mhalo;
            MPs[i].major = (MPs[i].mratio >= 0.25)?true:false;
            // printf("generation = %d, mratio = %3.2e, dm = %3.2e, dt = %3.2e, mhalo = %3.2e\n",MPs[i].j,MPs[i].mratio,MPs[i].dm/Ms,MPs[i].dt,MPs[i].mhalo/Ms);
        }
}


/*
rm ../tree_Hirano/fort.*mer
g++ read_aTree.cpp dyn.cpp PARA.cpp LE_adb.cpp RK4.o my_linalg.o class_halo.o  -o read.out && ./read.out
 */

/* int main(){
    MainProgenitor* MPs = NULL;
    MPs = new MainProgenitor [200];
    int nMP = 0;
    aTree(nMP, "../tree_Hirano/fort.20", MPs);
    int i=0;
    double n0_sol, ni=1;
    while(i<nMP){
        HALO halo(MPs[i].mhalo,MPs[i].z);
        // printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e\n",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
        // if (halo.Tvir<1.e4){
        //     printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e,\t",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
        //     printf("ng_adb:%3.2e\n",MPs[i].ng_adb);
        // }
        i++;      
    }

    delete [] MPs;
} */
