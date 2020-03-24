# include <iostream>
# include <stdio.h>
# include <cmath>

# include "my_linalg.h"
# include "gsl_inverse.h"
using namespace std;

/* to compile this file: 
g++ -c my_linalg.cpp && g++ my_linalg.o -o linalg && ./linalg
*/
void dot(double* vb, int N, double *mA, double* vx){
    double b[N];
    for (int i=0; i<N; i++){
        b[i] = 0.;
    }
    for (int i=0; i<N; i++){
        for (int j=0; j<N;j++){
            b[i] += mA[i*N+j]*vx[j];
        }
        *vb++ = b[i];
    }
}

double len_v(int N, double* v){
    double sum = 0.;
    for (int i=0; i<N; i++){
        sum += (*v) * (*v);
        v++;
    }
    return sqrt(sum);
}
/* 
int main(){
    void get_inverse();
    int const N = 6;

    double* vx = NULL; 
    double* vb = NULL;
    vb = new double[N];

    double* mA = NULL;
    // double aA [N*N] = {
    //    1,-1,0,0,0,-1,
    //    0,1,-1,-1,0,0,
    //    0,0,0,1,-1,1,
    //    6,6,4,0,0,0,
    //    0,0,-4,4,5,0,
    //    0,-6,0,4,0,6}; 
    double aA [N*N] = {
        1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1};
    double ax[N] = { 1, 2, 3, 4, 5, 6};
    mA = aA; vx = ax;
    dot(N,aA,vx,vb);
    cout<<"dot is: \n [\n";
    for (int i=0; i<N; i++){
        cout<<"\t"<<vb[i]<<endl;
    }
    cout<<"]"<<endl;
    delete [] vb;
    
    printf("#####################################\nlength_of_vector\n");
    cout<< len_v(6,vx)<<endl;
    
} */