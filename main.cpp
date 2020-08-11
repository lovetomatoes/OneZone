#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

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
#include "mpi.h"

#define ROOT    0
#define Ni      2
#define Nd      2
#define Ncore   2
#define Ntr     5
// Ncore = 4, use 4 cores, each claculate Ntr trees (bsm from 1 to 4)

static int i,j;
static int ntree = 10, MerMod = 1, i_bsm, itr; 
static double Tb = 2.e4;
static bool Ma_on = true, spec=false;

int main(int argc, char *argv[]){
    int size, rank;
    clock_t t0 = clock();

    string fout, tree, ftree;
    fout = "./Jcs_111.txt";
    fstream file;

    int data_i[Ni*Ntr], data_ri[Ncore*Ntr][Ni];
    double data_d[Nd*Ntr], data_rd[Ncore*Ntr][Nd];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank==ROOT){
        file.open(fout.c_str(), ios::out | ios::trunc );
        file<<setw(12)<<"tree"<<setw(12)<<"ibsm"<<setw(12)<<"Jc"<<setw(12)<<"zsol"<<endl;
    }
    
    i_bsm = 0;
    for (i=0;i<Ncore;i++){
        if (rank==i){
            for (itr=rank*Ntr; itr<(rank+1)*Ntr; itr++){
                ftree = "../treefiles/tree_"+to_string(itr);
                evol_Jc(&data_i[(itr-rank*Ntr)*Ni], &data_d[(itr-rank*Ntr)*Nd], itr, ftree, Tb, MerMod, spec, Ma_on, i_bsm);
                printf("rank=%d, data_i[0]=%d, data_d=%f\n",rank,data_i[0],data_d[0]);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(&data_i[0], Ni*Ntr, MPI_INT, &data_ri[0][0], Ni*Ntr, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Gather(&data_d[0], Nd*Ntr, MPI_DOUBLE, &data_rd[0][0], Nd*Ntr, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    if (rank==ROOT) {
        for (j=0;j<Ncore*Ntr;j++){
            file<<setw(12)<<data_ri[j][0]<<setw(12)<<data_ri[j][1];
            file<<setw(12)<<data_rd[j][0]<<setw(12)<<data_rd[j][1];
            file<<endl;
        }
        file.close();
    }
    MPI_Finalize();

    clock_t t1 = clock();
    if (rank==ROOT) printf("time taken 1st part:%.2f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
    return 0;
}
