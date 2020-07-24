#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Larson.h"
#include "RK4.h"
#include "PARA.h"
using namespace std;

/* g++ -c Larson.cpp RK4.cpp
g++ Larson.o class_halo.o dyn.o PARA.o RK4.o my_linalg.o -o larson
./larson
*/

//input para
//output: a profile of density distributions, Phi(x), x=r/a, Phi = log(rho/rho_g0)
// here x is time axis, x1: final time
// N: number of time steps 
extern double Rth = k_B/(mu*m_H); 
//
void colps(char* filename, double Tg, int const N){

    int const n = 2; // dimension of unknown vector
    int i,j;
    // set x1, largest xi 
    double x1 = 100;
    double* x = new double [N];
    
    double** y = new double* [N]; // y[0] = xi, y[1] = ln eta
    for (i=0;i<N;i++) y[i] = new double [n];

    double* dydx0 = new double [n];

    // dx
    double dx = x1/N;
    double err2=.01;
    // boundary conditions
    // // y[0] = xi, y[1] = ln eta
    //boundary x[0]=0, but 0 causes sigularity, using x[0]=dx/1000<<dx
    // 
    x[0] = dx/1.e8; y[0][0] = 0;
    double y1l = 0, y1r =  1;
    int N_y1 = 10000;
    double delta_y1 = (y1r-y1l)/N_y1, y10;

    int c=0; double* v=NULL; //no use, no other input parameters in *argv
    double min_err = 1000;
    for (j=0;j<N_y1;j++) {
        y[0][1] = y1l + j*delta_y1;
        for (i=1;i<N;i++){
            dx = x1/N;
            DyDx_L(x[i-1],y[i-1],dydx0,c,v);
            rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_L);
            check_conv(err2,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_L);
            x[i]=x[i-1]+dx;
            if (x[i]-y[i][0]>1+dx/2) break; // 后半段不要 算到x-xi=1处check
            if (x[i]-y[i][0] >=1-dx/2 && x[i]-y[i][0]<= 1+dx/2){
                //boundary condition: at x-xi=1, eta*x should =2
                if (pow(exp(y[i][1])*x[i]-2,2)< min_err){
                    printf("j=%d, eta=(%3.2e),\teta*x-2 =%3.2e,\t min_err=%3.2e\n",j,exp(y[0][1]), pow(exp(y[i][1])*x[i]-2,2), min_err);
                    y10 = y[0][1];
                    min_err = pow(exp(y[i][1])*x[i]-2,2);
                }      
            }
        }
    }
    printf("finala  eta =%3.2e,\t eta0 use: %3.2e\n",exp(y1l + j*delta_y1),exp(y10));
    fstream file;
//    file.open("data/profile.txt", ios::out | ios::trunc );
    file.open(filename, ios::out | ios::trunc );

    file<<"x xi lneta\n";
    y[0][1] = y10;
    for (i=1;i<N;i++){
        dx = x1/N;
        DyDx_L(x[i-1],y[i-1],dydx0,c,v);
        rk4(y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_L);
        check_conv(err2,y[i],x[i-1],dx,y[i-1],n,dydx0,c,v,DyDx_L);
        x[i]=x[i-1]+dx;
        ////output file<<"x xi lneta\n"
        file<<x[i]<<" "<<y[i][0]<<" "<<y[i][1]<<endl;
    }

    for (i=0;i<N;i++) delete [] y[i];
    delete [] x; delete [] y; delete [] dydx0;
    file.close();
}

void DyDx_L(double x, double* y, double* dydx, int argc, double* argv){
    //y[0] = xi; y[1] = log(eta)
    dydx[0] = (x-y[0])/x * (exp(y[1])*x*(x-y[0]) - 2) / (pow(x-y[0],2)-1);
    dydx[1] = (x-y[0])/x * (exp(y[1])*x -  2*(x-y[0])) / (pow(x-y[0],2)-1);
}


int main(){
// for N = 10000;  x1 = 10; 
    int const N = 100000;
    double Tg = 1.e4; 

    char* fname = "data/colps.txt";
    colps(fname,Tg,N);

    return 0;
}