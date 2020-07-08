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
        cout<<filename<<"TREEmer exist\n";
        while (getline(readmer,line)){
            stringstream ss(line);
            ss>>MPs[i].id_tree>>MPs[i].j>>MPs[i].mhalo>>MPs[i].z>>MPs[i].t>>MPs[i].c>>MPs[i].Tvir>>MPs[i].ng_adb>>MPs[i].mratio;
            MPs[i].t *= Myr; //time from universe age of 0
            MPs[i].mhalo *= Ms; // mhalo unit -> g
            i++;
        }
        readmer.close();
        nMP = i-1;
        // printf("READMER:%d\n",nMP);
        // for (i=1;i<=nMP;i++){
        //     printf("IN READMER:\n i=%d,MPs[%d].j=%d,mhalo=%3.2e,z=%3.2f,Tvir=%3.2e,c=%3.2e,n_adb=%3.2e\n",
        //            i,MPs[i].j,MPs[i].mhalo/Ms,MPs[i].z,MPs[i].Tvir,MPs[i].c,MPs[i].ng_adb );
        // }
    }
// no mer file; calculate n_adb; write
    else{
        cout<<" no mer file\n";
        readtree.open(treename);
        int id_m;
        while (getline(readtree,line)){
            // cout<<line;
            stringstream ss(line);
            ss>>MPs[i].j>>id_m>>MPs[i].mhalo>>MPs[i].z>>MPs[i].Tvir>>MPs[i].c>>MPs[i].id_tree;
            MPs[i].t = t_from_z(MPs[i].z); // time from universe age of 0
            MPs[i].mhalo *= Ms; // mhalo unit -> g (for treefiles/. file)
            // printf("IN READTREE:\n MPs[%d].j=%d,mhalo=%3.2e,z=%3.2f,Tvir=%3.2e,c=%3.2e,id_tree=%d \n",
            //         i,MPs[i].j,MPs[i].mhalo/Ms,MPs[i].z,MPs[i].Tvir,MPs[i].c,MPs[i].id_tree );
            i++;
        }
        readtree.close();
        nMP = i-1;
        // cout<<"tree: nMP= "<<nMP<<endl;

    // 调换位置 以num/2为轴 MPs[1]<->MPs[nMP]
        num = nMP + 1; // cout<<"num= "<<num<<"num/2= "<<num/2<<endl;
        double temp_dbl;
        for (i=1;i<num/2;i++){
            temp_dbl = MPs[num-i].z; MPs[num-i].z = MPs[i].z; MPs[i].z = temp_dbl;
            temp_dbl = MPs[num-i].t; MPs[num-i].t = MPs[i].t; MPs[i].t = temp_dbl;
            temp_dbl = MPs[num-i].mhalo; MPs[num-i].mhalo = MPs[i].mhalo; MPs[i].mhalo = temp_dbl;
            temp_dbl = MPs[num-i].c; MPs[num-i].c = MPs[i].c; MPs[i].c = temp_dbl;
            temp_dbl = MPs[num-i].Tvir; MPs[num-i].Tvir = MPs[i].Tvir; MPs[i].Tvir = temp_dbl;
            }
    // reset t starting point to the 1st progenitor; n_adb
        temp_dbl = MPs[1].t;
        for (i=1;i<=nMP;i++){ 
            MPs[i].t -= temp_dbl; 
            HALO halo(MPs[i].mhalo,MPs[i].z);
            double ni=1; MPs[i].ng_adb = 1;
            // only calculate ng_adb for halo w/ Tvir <1.e5 K; wli
            if (MPs[i].Tvir<1.e5) Mg2N0_adb(MPs[i].ng_adb,ni,MPs[i].z,MPs[i].mhalo);
            MPs[i].dm = (i==nMP)? 0: MPs[i+1].mhalo - MPs[i].mhalo;
            MPs[i].mratio = (i==nMP)? 0:MPs[i].dm/MPs[i].mhalo;
        }

        ofstream f;
        f.open(filename, ios::out | ios::trunc);
        f<<setw(17)<<"id_tree"<<setw(17)<<"j"<<setw(17)<<"Mh_Ms";
        f<<setw(17)<<"z"<<setw(17)<<"t_Myr"<<setw(17)<<"c";
        f<<setw(17)<<"Tvir"<<setw(17)<<"n_adb"<<setw(17)<<"q"<<endl;
        f<<setiosflags(ios::scientific)<<setprecision(8);
        for (i=1;i<=nMP;i++){
            f<<setw(17)<<MPs[i].id_tree<<setw(17)<<MPs[i].j<<setw(17)<<MPs[i].mhalo/Ms;
            f<<setw(17)<<MPs[i].z<<setw(17)<<MPs[i].t/Myr<<setw(17)<<MPs[i].c;
            f<<setw(17)<<MPs[i].Tvir<<setw(17)<<MPs[i].ng_adb<<setw(17)<<MPs[i].mratio<<endl;
        }
        f.close();
    }

    printf("\t#\t#\t#\t tree_%d: nMP= %d\n",MPs[1].id_tree,nMP);
// reset some attributes
    for (i=1;i<=nMP;i++){
        HALO halo(MPs[i].mhalo,MPs[i].z);
    //consistent w/ class_gas
        MPs[i].Tvir = halo.Tvir;
        MPs[i].c = halo.c;

        MPs[i].dt = (i==nMP)? 0: MPs[i+1].t - MPs[i].t;
        MPs[i].dm = (i==nMP)? 0: MPs[i+1].mhalo - MPs[i].mhalo;
        MPs[i].major = (MPs[i].mratio >= 0.25)?true:false;
        // printf("generation = %d, mratio = %3.2e, dm = %3.2e, dt = %3.2e, mhalo = %3.2e\n",MPs[i].j,MPs[i].mratio,MPs[i].dm/Ms,MPs[i].dt,MPs[i].mhalo/Ms);
    }
}


/*
rm ../tree_Hirano/fort.*mer

rm ../treefiles/*mer
g++ read_aTree.cpp dyn.cpp PARA.cpp LE_adb.cpp RK4.o my_linalg.o class_halo.o  -o read.out && ./read.out
 */

/* int main(){
    MainProgenitor* MPs = NULL;
    MPs = new MainProgenitor [200];
    int nMP = 0;
    aTree(nMP, "../treefiles/tree_10", MPs);
    int i=1;
    double n0_sol, ni=1;
    cout<<"aa"<<3/2<<" "<<4/2<<endl;
    while(i<=nMP){
        // printf("%dth ");
        // HALO halo(MPs[i].mhalo,MPs[i].z);
        // printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e\n",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
        // if (halo.Tvir<1.e4){
        //     printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e,\t",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
        //     printf("ng_adb:%3.2e\n",MPs[i].ng_adb);
        // }
        i++;      
    }

    delete [] MPs;
} */