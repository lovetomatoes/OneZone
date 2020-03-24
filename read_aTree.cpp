#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "class_halo.h"
#include "read_aTree.h"
#include "PARA.h"
#include "dyn.h"
#include <cmath>
using namespace std;
/* 
g++ read_aTree.cpp dyn.cpp PARA.cpp class_halo.o -o read && ./read
 */

void read_aTree (int& num, char* fname, MainProgenitor* MPs){

    //MainProgenitor MPs[200];

    ifstream infile;
    infile.open(fname);
    string line;
    int i = 0;
    int id_m = 0;
//    int id_tree;
    //printf("++++++++++++++++++++\n");
    while (getline(infile,line)){
        stringstream ss(line);
        ss>>MPs[i].j>>id_m>>MPs[i].mhalo>>MPs[i].z>>MPs[i].Tvir>>MPs[i].id_tree;
        MPs[i].t = t_from_z(MPs[i].z); //time from universe age of 0
        /* printf("MPs[%d].mhalo = %3.2e \n",i,MPs[i].mhalo);
        printf("MPs[%d].z = %3.2e \n",i,MPs[i].z);
        printf("MPs[%d].id_tree = %3.2e \n\n",i,MPs[i].id_tree); */
        i++;
    }
    num = i-2; //减1是因为i最后有一次++, 再减1是因为analysis_trees.f 输出给fort.11的最后一行带\n

    infile.close();
}


void aTree(int& nMer, MainProgenitor* MPs, char* filename){
    int num = 0, i = 0, j;
    /* char fname[100]; 
    sprintf(fname, "../code_tree/fort.217");  */
    // fort.211, fort.212得用mac来跑, haha上面报错segmentation fault
    //MainProgenitor* p = NULL;
    //read_aTree(num, fname, p);
    double temp_dbl;
    int temp_int;
    //read_aTree(num,fname,MPs);
    read_aTree(num,filename,MPs);
    // 调换位置，MPs[num]->MPs[1] 倒序赋值到MPs[0]->MPs[num-1]里去;
    // 原本到MPs[0]无效, 因为fort.11中第一行header: # j id_m mhalo_m/1.e11Ms z_m T_{vir,m} 赋给了MPs[0]
    for (i=0;i<(num+1)/2;i++){
        temp_dbl = MPs[num-i].z; MPs[num-i].z = MPs[i].z; MPs[i].z = temp_dbl;
        temp_dbl = MPs[num-i].t; MPs[num-i].t = MPs[i].t; MPs[i].t = temp_dbl;
        temp_dbl = MPs[num-i].mhalo; MPs[num-i].mhalo = MPs[i].mhalo; MPs[i].mhalo = temp_dbl;
        temp_dbl = MPs[num-i].Tvir; MPs[num-i].Tvir = MPs[i].Tvir; MPs[i].Tvir = temp_dbl;
        temp_int = MPs[num-i].j; MPs[num-i].j = MPs[i].j; MPs[i].j = temp_int;

        //MPs[i].t = t_from_z(MPs[i].z);
        /* in array MPs[num]: MPs[i].j generation, 
                           MPs[i].z redshift, 
                           MPs[i].t time,
                           MPs[i].mhalo halo mass grow from high redshift
                           MPs[i].Tvir virial temperature */
    }

    temp_dbl = MPs[0].t;
    for (i=0;i<num;i++){
        MPs[i].t -= temp_dbl; // reset t starting point to the 1st progenitor 
        MPs[i].mhalo *= (1.e11*Ms); // mhalo unit -> g
        //cout<<"j="<<MPs[i].j<<endl;
    }
    for (i=0;i<num;i++){
        HALO halo(MPs[i].mhalo,MPs[i].z);
        MPs[i].Tvir = halo.Tvir;

        MPs[i].dt = (i==num-1)? 0: MPs[i+1].t - MPs[i].t;
        MPs[i].dm = (i==num-1)? 0: MPs[i+1].mhalo - MPs[i].mhalo;
        MPs[i].mratio = (i==num-1)? 0:MPs[i].dm/MPs[i].mhalo;
        MPs[i].major = (MPs[i].mratio >= 0.25)?true:false;
        //if (MPs[i].major) cout<<"ratio = "<<MPs[i].mratio<<endl;
        //printf("generation = %d, mratio = %3.2e, dm = %3.2e, dt = %3.2e, mhalo = %3.2e\n",MPs[i].j,MPs[i].mratio,MPs[i].dm/Ms,MPs[i].dt,MPs[i].mhalo/Ms);
    }
    
    cout<<"\nnum = "<<num<<endl;
    nMer = num;
}


/* int main(){
    MainProgenitor* MPs = NULL;
    MPs = new MainProgenitor [200];
    int nMer = 0;
    aTree(nMer, MPs, "../code_tree/fort.211");
    delete [] MPs;
}
 */