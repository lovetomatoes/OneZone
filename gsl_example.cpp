#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
using namespace std;

int main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  //printf ("J0(%g) = %.18e\n", x, y);
  cout << "J0(5)= "<< y<<endl;
  return 0;
}