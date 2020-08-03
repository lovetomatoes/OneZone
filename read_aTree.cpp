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
        readtree.open(treename); cout<<treename<<"\n";
        int id_m;
        while (getline(readtree,line)){
            // cout<<line<<endl;
            stringstream ss(line);
            ss>>MPs[i].j>>MPs[i].mhalo>>MPs[i].z>>MPs[i].Tvir>>MPs[i].c>>MPs[i].id_tree;
            MPs[i].t = t_from_z(MPs[i].z); // time from universe age of 0
            MPs[i].mhalo *= Ms; // mhalo unit -> g (for treefiles/. file)
            // printf("IN READTREE:\n MPs[%d].j=%d,mhalo=%3.2e,z=%3.2f,Tvir=%3.2e,c=%3.2e,id_tree=%d \n",
            //         i,MPs[i].j,MPs[i].mhalo/Ms,MPs[i].z,MPs[i].Tvir,MPs[i].c,MPs[i].id_tree );
            i++;
        }
        readtree.close();
        nMP = i-1;
        cout<<"tree: nMP= "<<nMP<<endl;

    // 调换位置 以num/2为轴 MPs[1]<->MPs[nMP]
    // wli: 必须*0.5而非/2取整!!! 否则出错 (int, float)比大小先把int转化成float
        num = nMP + 1; // cout<<"num= "<<num<<"num/2= "<<num/2<<endl;
        double temp_dbl;
        for (i=1;i<0.5*num;i++){
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
            double ni=100.; MPs[i].ng_adb = 1; //初始尝试值不能太小 原本ni=1 牛顿迭代求解出错
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

    // printf("\t#\t#\t#\t tree_%d: nMP= %d\n",MPs[1].id_tree,nMP);
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
rm ../tree_Jc/tree*mer
g++ read_aTree.cpp dyn.cpp PARA.cpp LE_adb.cpp RK4.o my_linalg.o class_halo.o  -o read.out && ./read.out
 */

/* int main(){
    MainProgenitor* MPs = NULL;
    MPs = new MainProgenitor [300];
    int nMP = 0;
    aTree(nMP, "../tree_Jc/tree_0", MPs);
    // aTree(nMP, "../treefiles/tree_0", MPs);
    int i=1;
    double n0_sol, ni=1;
    cout<<"aa"<<3/2<<" "<<4/2<<endl;
    num = 5;
    for (i=1;i<0.5*num;i++) printf("i=%d, num=%d, int(num/2)=%d, num/2=%f\n",i, num, int(num/2),num/2);
    fstream f1;
    f1.open("ts.txt", ios::out | ios::trunc );
    f1<<setiosflags(ios::scientific)<<setprecision(5);
    f1<<setw(16)<<"z"<<setw(16)<<"Tvir"<<setw(16)<<"Mh";
    f1<<setw(16)<<"dt1000"<<setw(16)<<"t_ff"<<setw(16)<<"tg_ff";
    f1<<setw(16)<<"dt100"<<setw(16)<<"dt_col"<<setw(16)<<"t_dyn";
    f1<<setw(16)<<"ng_adb"<<setw(16)<<"nc_DM";
    f1<<endl;
    while(i<=nMP){
        HALO halo(MPs[i].mhalo,MPs[i].z);
        // printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e\n",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
        if (halo.Tvir<5.e4){
            // printf("%dth:z=%3.2f, Mh=%3.2e, Tv=%3.2e, Kv=%3.2e,\t",i,halo.z, halo.Mh/Ms, halo.Tvir, halo.Kvir);  
            // printf("nc_DM:%3.2e, ng_adb:%3.2e\n",halo.rho_c/(mu*m_H), MPs[i].ng_adb);
            // printf("t_ff= %3.2e, MPs.dt= %3.2e\n",t_freefall(halo.rho_c/(mu*m_H)+ MPs[i].ng_adb)/Myr, MPs[i].dt/Myr);
            f1<<setw(16)<<MPs[i].z<<setw(16)<<halo.Tvir<<setw(16)<<halo.Mh/Ms;
            f1<<setw(16)<<(t_from_z(MPs[i].z)-t_from_z(MPs[i].z+50./1000.))/Myr;
            f1<<setw(16)<<t_freefall(halo.rho_c/(mu*m_H)+ MPs[i].ng_adb)/Myr;
            f1<<setw(16)<<t_freefall(MPs[i].ng_adb)/Myr;
            f1<<setw(16)<<MPs[i].dt/Myr;
            f1<<setw(16)<< t_freefall(18.*pi*pi*RHO_crit(MPs[i].z)/(mu*m_H)) /Myr;
            f1<<setw(16)<<halo.t_dyn/Myr;
            f1<<setw(16)<<MPs[i].ng_adb;
            f1<<setw(16)<<halo.rho_c/(mu*m_H);
            f1<<endl;
        }
        i++;      
    }
    f1.close();
    delete [] MPs;
} */