#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/* g++ -c LaneEmden.cpp
g++ LaneEmden.o -o le
./le
 */

double n = 1.5; // adiabatic case
int i;

int main(){
    fstream file;
    file.open("Lane.txt", ios::out | ios::trunc );
    file<<"i ita theta\n";
    cout<<n<<endl;
    int N = 10000;
    double* ita = new double [N];
    double* theta = new double [N];
    double* DthetaDita = new double [N];
    // dx
    double dita = .001;
    // boundary conditions
    ita[0] = 0;
    theta[0] = 1;
    DthetaDita[0] = 0;
    // integration
    while (theta[i]>=0 && i<N-1){
         
        file<<i<<" "<<ita[i]<<" "<<theta[i]<<endl;
        ita[i+1] = ita[i] + dita; 
        DthetaDita[i+1] = ( pow(ita[i],2)*DthetaDita[i] - pow(theta[i],n)*pow(ita[i],2)*dita ) / pow(ita[i+1],2);
        theta[i+1] = theta[i] + DthetaDita[i+1]*dita;
        i+=1;
    }

    cout<<i<<endl;

    delete [] ita;
    delete [] theta;

    file.close();
    return 0;

}