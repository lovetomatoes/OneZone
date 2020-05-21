#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "my_linalg.h"
#include "RK4.h"
#include "PARA.h"
using namespace std;

// rk4 paras: x0, y(x0), dydx(x0), h = Delta_x, derivatives f(x,y,h)=dydx(x,y(x),h)
// then return f(x+h); may not be accurate enough, h*=0.5 and call it twice, compare y(x0+h) v.s. y(x0+h/2)

/* 
g++ -c RK4.cpp
g++ RK4.o my_linalg.o -o RK4
./RK4
 */

void DyDx_iso(double x, double* y, double* dydx, int argc, double* argv){
    double alpha = argv[0];
    double R = argv[1];
    dydx[0] = y[1]/(x*x);
    dydx[1] = x*x*exp(-y[0]); //only gas self-gravity
    dydx[1] = x*x/ (R * (alpha*x) * pow(1+alpha*x,2) ); //only dark matter density
    dydx[1] = x*x*exp(-y[0]) + x*x/ (R * (alpha*x) * pow(1+alpha*x,2) ); // DM + gas
    //printf("self grav=%3.2e\t, DM=%3.2e\n",x*x*exp(-y[0]),x*x/R/ (alpha*x) / pow(1+alpha*x,2));
}

void DyDx_adb_Kc(double x, double* y, double* dydx, int argc, double* argv){
    double n_adb = argv[0];
    double R = argv[1]; //printf("alpha=%3.2e, R=%3.2e\n", alpha, R);
    double alpha = argv[2];
    dydx[0] = -y[1]/(x*x);
    dydx[1] = x*x * pow(y[0], n_adb); //only gas self-gravity
    dydx[1] = x*x/ (R * (alpha*x) * pow(1+alpha*x,2) ); //only dark matter density
    dydx[1] = x*x * pow(y[0], n_adb) + x*x/ (R * (alpha*x) * pow(1+alpha*x,2) ); // DM + gas
}

void DyDx_adb_fit(double x, double* y, double* dydx, int argc, double* argv){
    double n_adb = argv[0];
    double R = argv[1]; //printf("alpha=%3.2e, R=%3.2e\n", alpha, R);
    double alpha = argv[2];
    double beta = argv[3];
    double gamma = 1. + 1./n_adb;
    // dydx[0] = y[1]*pow(y[0],2.-gamma)/ pow(x,2) / gamma; // K = Kcore const (+gas only checked; w/ K const y0=theta)
    // dydx[0] = ( y[1]*pow(y[0],2.-gamma)/pow(x,2) - y[0] )/ (x*gamma); // Kr; infinity; cannot check
    dydx[0] = ( y[1]*pow(y[0],2.-gamma)/pow(x,2) - y[0] )/ ((1+x)*gamma); // Kfit; finally into use

    dydx[1] = -beta * pow(x,2) * y[0]; //gas only 
    dydx[1] = -beta* pow(x,2) / (R * (alpha*x) * pow(1+alpha*x,2) ); // DM only 
    dydx[1] = -beta* pow(x,2)*( y[0] + 1./ (R * (alpha*x) * pow(1+alpha*x,2)) ) ; // gas+DM
}

void rk4(double* yout, double x, double h, double *y, int const n, double *dydx0,
            int argc, double* argv,
            void (*derivs)(double, double*, double*, int, double*)){
    int i;
    double xh, hh, h6;
    double *dym; double *dyt; double *yt;
    dym = new double [n]; dyt = new double [n]; yt = new double [n];
    hh = h*0.5;
    h6 = h/6.;
    xh = x + hh;
    //int c=1; double ** v=NULL;
    int c = argc;
    double* v = new double [c];
    for (i=0;i<c;i++) v[i]= argv[i];
    for (i=0;i<n;i++) yt[i] = y[i] + hh*dydx0[i];  // dydx0 = K1
    (*derivs)(xh,yt,dyt,c,v);                         // dyt = K2
    for (i=0;i<n;i++) yt[i] = y[i] + hh*dyt[i];
    (*derivs)(xh,yt,dym,c,v);                         // dym = K3
    for (i=0;i<n;i++){
        yt[i] = y[i] + h*dym[i];
        dym[i] += dyt[i];                         // dym += K2
    }
    (*derivs)(x+h,yt,dyt,c,v);                        // dyt = K4
    for (i=0;i<n;i++) yout[i] = y[i] + h6*(dydx0[i]+dyt[i]+2.*dym[i]);

    delete [] dym; delete [] dyt; delete [] yt; delete [] v;
}

void check_conv(double err_crit, double* yout, double x, double& h, double *y, int const n, double *dydx0,
           int argc, double* argv,
           void (*derivs)(double, double*, double*, int, double*) ){ // check_convergence...
        
    int i;
    double hh=h/2.;
    double x1;
    int c = argc;
    double* v = new double [c];
    for (i=0;i<c;i++) v[i]= argv[i];

    double* yout1; double* yout2; double* delta_y;

    yout1 = new double [n]; yout2 = new double [n]; delta_y = new double [n];

    int i_iter = 0, i_fold = 0;
    for (i=0;i<n;i++) delta_y[i] = yout[i];
    while (i_iter<=3) {

        (*derivs)(x,y,dydx0,c,v);
        rk4(yout1,x,hh,y,n,dydx0,c,v,(*derivs));
        (*derivs)(x+hh,yout1,dydx0,c,v);
        rk4(yout2,x+hh,hh,yout1,n,dydx0,c,v,(*derivs));
        for (i=0;i<n;i++) {
            delta_y[i] = yout2[i] - yout[i];
            yout[i] = yout2[i];
        }
        for (i=0;i<n;i++) if (!isfinite(yout[i])) {
            // throw "wrong";
            cout<<"INFINITY condition! ";
        }// 异常需处理...

        if ( len_v(n,delta_y)>err_crit*len_v(n,yout) ){
            // cout<<"REFINING GRID in check_conv!!!!!!!!!!!!!\n"; // report iteration times! 
            i_fold += 1;
            h *= 0.5;
            hh = h/2.;
            for (i=0;i<n;i++) yout[i] = yout1[i];
        }
        else break;
        // if (i_iter==3) cout<<"ITERATION TIMES > 3! REFINE GRID!!!!!!!!!!!!!!!\n"; 
        // wli: commented out. need to improve...
        i_iter ++;
    }
    delete [] yout1; delete [] yout2; delete [] delta_y; delete [] v;
}


void ydot(double x, double*y, double*dydx0, int argc, double** argv) {
    dydx0[0] = pow(x,3.) - y[0]/x;
}
void ydot_sin(double x, double*y, double*dydx0, int argc, double** argv) {
    dydx0[0] = sqrt(1-y[0]*y[0]);
    dydx0[1] = -y[0];
}

/* int main(){//改过函数参数设定
    void rk4 (double* yout, double x, double h, double* y, int const n, double* dydx0,   
              void (*derivs)(double, double*, double*, int argc, double** argv) );
    void ydot(double x, double*y, double*dydx0, int argc, double** argv);
    void ydot_sin(double x, double*y, double*dydx0, int argc, double** argv);
    void check_conv(double err_crit, double* yout, double x, double& h, double *y, int const n, double *dydx0,   
           void (*derivs)(double, double*, double*, int argc, double** argv) );

    // case 2: dydx for y[0]=sin(x), y[1]=cos(x)
    int const n=2;
    double y[n], dydx0[n], yout[n];
    int i;
    double x, h, hh;
    double x1;
    int c=1; double ** v=NULL;

    x = 0; h= 1;
    y[0] = 0; y[1] = 1;
    ydot_sin(x,y,dydx0,c,v);
    rk4(yout, x, h, y, n, dydx0, ydot_sin);
    for (i=0;i<n;i++) printf("%8.6f\n",yout[i]);
    check_conv(.01,yout,x,h,y,n,dydx0,ydot_sin);
    printf("new one: \n");
    for (i=0;i<n;i++) printf("%8.6f\n",yout[i]);
    printf("h=%3.2f\n",h);
}
 */
//prev position: line 57; y是一维矢量的情况
/* // case 1; dydx = x^3-y/x; y(1)=.4; h=.1
    int const n=1;
    double y[n], dydx0[n], yout[n];
    int i;
    double x, h;


    x = 1; h= .1;
    double fy[n], fdydx0[n], fyout[n];
    fy[0] = .4;
    ydot(x,fy,fdydx0);
    rk4(fyout, x, h, fy, n, fdydx0, ydot);
    for (i=0;i<n;i++ ) printf("%8.6f\n",fyout[i]); */