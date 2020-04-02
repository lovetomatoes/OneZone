// kpd for H2, H- and H2+
// g++ kpd.cpp PARA.cpp -o kpd.out && ./kpd.out
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "kpd.h"
#include "PARA.h"
using namespace std;


//int i,j; 
/* 与其他file里面的声明冲突了 
   报错 duplicate symbol _i in:
    /var/folders/q8/6kwqz_k539x_42fy_9xhrlhm0000gn/T/main-35c185.o
    kpd.o */
int const nnu = 45;
int const nT = 10;
int const nnu_spec = 1221;

int main(){
   double T_rad = 1.e5;
   double k_Hm=0, k_H2p=0;
   kpd_Hm_H2p(T_rad,k_Hm,k_H2p);
}

void kpd_Hm_H2p(double T_rad, double& k_Hm, double& k_H2p){
   int i,j;
   //double T_rad = 5.e4; 
   printf("IN KPD: T_rad=%5.2f\n", T_rad);
   double nurat0 = 1.005;
   double nu_min = 0.45*eV/h_p;
   double nurat = nurat0, nu_prev = nu_min;
   double nu_eV;
  
   double* nub = NULL; double* fluxb_cont = NULL; double* fluxb_line = NULL;
   nub = new double [nnu_spec]; fluxb_cont = new double [nnu_spec]; fluxb_line = new double [nnu_spec];
   string fspec_cont = "rspec_pop3_cont_fe00.txt";
   string fspec_line = "rspec_pop3_line_fe00.txt";
   read_spec(nub, fluxb_cont, fluxb_line, fspec_cont, fspec_line);

   double** sigmaa = new double* [nT];
   for (i=0;i<nT;i++) sigmaa[i] = new double[nnu];
   double* nua = new double [nnu];
   double Ta[nT] = {0,2520,3150,4200,5040,6300,8400,12600,16800,25200};
   read_sigma(nua,sigmaa);
   double sigma_Hm, sigma_H2p;
   double T_H2, sigma;
   
//     H2+ photodetachment
   T_H2 = 8.0e3;
   double E_ly = 12.4, E_end = 13.6;
   int k, NF = 1000, k_ly;
   double nu[NF], Planck[NF], Flux[NF];
   
   for (k=0;k<NF;k++){
      nu[k] = nu_prev*nurat;
      if (nu_prev<E_ly*eV/h_p and nu[k]>=E_ly*eV/h_p) {
         nu[k] = E_ly*eV/h_p;
         k_ly = k;
      }
      if (nu_prev<E_end*eV/h_p and nu[k]>=E_end*eV/h_p) {
         nu[k] = E_end*eV/h_p;
         //printf("E_end=%3.2e,  hnu=%3.2eV\n",E_end, h_p*nu[k]/eV);
         NF = k+1;
      }

      //printf(" hnu[%d]=%3.2e\n", k, h_p*nu[k]/eV);

      Planck[k] = pow(nu[k],3)/(exp(h_p*nu[k]/k_B/T_rad)-1.0);
      Flux[k] = pow(h_p*nu[k]/eV,-1.5);

      // interpolation from spectrum x:lambda in Angstron
      linear(nub,fluxb_cont,nnu_spec,c/nu[k]/Angstron,Flux[k]);
      //Flux[k] = Planck[k];
 
      nu_prev = nu[k];
   }

   
   printf("NF=%d\n",NF);// wli note: NF=684
   
   k_Hm=0.; k_H2p=0.;
   double dnu[NF-1];
   fstream f1;
   f1.open("../data/k_intg.txt", ios::out | ios::trunc );
   f1<<" k E Planck flux sigma_Hm sigma_H2p k_pd_Hm, k_pd_H2p\n";
   for (k=0;k<NF-1;k++){
      dnu[k] = nu[k+1]-nu[k];
      printf("dnu/nu*c=%3.2e\n",dnu[k]*c/nu[k]/km);
      if (nu[k] < E_end*eV/h_p){
         nu_eV = h_p*nu[k]/eV;
         Hm_CrossSec(sigma_Hm, nu_eV);
         H2p_bf_CrossSec(sigma_H2p, nu_eV, T_H2, Ta, nua, sigmaa);         
         k_Hm   += sigma_Hm*4*pi*Flux[k]/h_p/nu[k]*dnu[k];
         k_H2p += sigma_H2p*4*pi*Flux[k]/h_p/nu[k]*dnu[k];
         f1<<" "<<k<<" "<<h_p*nu[k]/eV<<" "<<Planck[k]<<" "<<Flux[k]<<" "<<sigma_Hm<<" "
           <<sigma_H2p<<" "<<k_Hm<<" "<<k_H2p<<endl;
      }

   }
   
   printf("norm for flux at 12.4eV: Planck:%3.2e Flux:%3.2e\n",Planck[k_ly],Flux[k_ly]);
   k_Hm *= (1.e-21/Flux[k_ly]);
   k_H2p *= (1.e-21/Flux[k_ly]);
   printf("k_Hm=%3.2e, k_H2p=%3.2e, kHm/kH2=%3.2e\n", k_Hm, k_H2p, k_Hm/1.39e-12);
   
   for(i=0; i<nT; i++) delete[] sigmaa[i];
   delete[] sigmaa;
   delete[] nua;
   delete[] nub; delete[] fluxb_cont; delete[] fluxb_line;
}

/* !     H- photodetachment
!     H-   +  ph.   ->   H    +    e
!     photon energy E in eV
!     ref. T.L.John (1988 A&A 193,189)
 */
void Hm_CrossSec(double& sigma, double E){
   double x,y_w,y_0,y_1,y,P,y_a,
          a1,s,a2,a3,
          rlmbd_0,rlmbd, C1, C2, C3, C4, C5, C6, f;
   
   rlmbd_0=1.6419;
   C1=152.519;
   C2=49.534;
   C3=-118.858;
   C4=92.536;
   C5=-34.194;
   C6=4.982;
   rlmbd=1.23985/E;
   if(rlmbd<=rlmbd_0) {
      x=sqrt(1./rlmbd-1./rlmbd_0);
      f=C1+x*(C2+x*(C3+x*(C4+x*(C5+x*C6))));
      sigma=1.e-18*pow(rlmbd,3)*pow(x,3)*f;
   }
   else sigma = 0.0;
   if(rlmbd<=0.09) sigma = 0.0;
}

/* C     cross section of H2+ b-f
C     H2+   +  ph.   ->   H   +   H+
C     frequency rnu in eV
C     ref. P.C.Stancil (1994 ApJ 430,360)
 */
void H2p_bf_CrossSec(double& sigma, double nu_eV, double T_K, double* Ta, double* nua, double** sigmaa){
//         H2p_bf_CrossSec(sigma_H2p, nu[k], T_H2, Ta, nua, sigmaa);
  /*  cout<<"IN H2P CROSS_SEC:  ";
   for(j=0; j<nnu; j++) cout<<nua[j]<<"\t";
   cout<<endl<<endl; */
   int i,j;
   double s1, s2, lnu, lnu1, lnu2, ls, ls1, ls2; 
   double sa[nT];
   T_K=min(T_K,Ta[nT-1]); //T_K = 8000K
   if(nu_eV<nua[0]){
      for (i=0;i<nT;i++) sa[i]=sigmaa[i][0];
      linear(Ta,sa,nT,T_K,s1);
      for (i=0;i<nT;i++) sa[i]=sigmaa[i][1];
      linear(Ta,sa,nT,T_K,s2);
      lnu=log(nu_eV);
      lnu1=log(nua[0]); lnu2=log(nua[1]);
      ls1=log(s1); ls2=log(s2);
      ls = ((lnu-lnu1)/(lnu2-lnu1))*ls2 + ((lnu2-lnu)/(lnu2-lnu1))*ls1; // linear interpolation
      sigma = exp(ls);
   }
   else if (nu_eV>nua[nnu-1]) sigma=0.0;
   else bilinear(Ta,nua,sigmaa,nT,nnu,T_K,nu_eV,sigma);
   
}

void linear(double* xa, double* ya, int m, double x, double& y){
   int ms;
   int i,j;
   double y1, y2, t;
   if( x<xa[0] or x>xa[m-1]) printf("xa_min=%3.2e, x=%3.2e, xa_max=%3.2e", xa[0],x,xa[m-1]);
   for (i=0;i<m;i++){
      if (x-xa[i]<=0.) {
         ms = i;
         break;
      }
   }
   if (ms==0) ms=1; //x=xa[0]
   y1 = ya[ms-1];
   y2 = ya[ms];
   t=(x-xa[ms-1])/(xa[ms]-xa[ms-1]);
   y=(1.-t)*y1+t*y2;
}

void bilinear(double* x1a, double* x2a, double** ya, int m, int n, double x1, double x2, double& y){
   double t, u, y1, y2, y3, y4;
   int ms, ns;
   int i,j;

   if((x1<x1a[0]) or (x1>x1a[m-1]) or (x2<x2a[0]) or (x2>x2a[n-1])) cout<<"IN BILINEAR, WRONG BOUNDARY"<<endl;
   else{ 
      for (i=0;i<m;i++){
         if (x1-x1a[i]<=0.0){
            ms=i;
            break;
         }
      }
      for (i=0;i<n;i++){
         if (x2-x2a[i]<=0.0){
            ns=i;
            break;
         }
      }
      if(ms==0) ms=1;
      if(ns==0) ns=1;
      y1=ya[ms-1][ns-1]; 
      y2=ya[ms][ns-1]; 
      y3=ya[ms][ns];
      y4=ya[ms-1][ns];

      t=(x1-x1a[ms-1])/(x1a[ms]-x1a[ms-1]);
      u=(x2-x2a[ns-1])/(x2a[ns]-x2a[ns-1]);
      y=(1.-t)*(1.-u)*y1 + t*(1.-u)*y2 + t*u*y3 + (1.-t)*u*y4;
   }
}

void read_sigma(double* nua, double** sigmaa){
   string line, nu_str, sigma_str, name;
   int i,j;
   int const nfile = nT;
   int T[nfile] = {0,2520,3150,4200,5040,6300,8400,12600,16800,25200};
   for (i=0;i<nfile;i++){
      name = "./H2p/sigma_"+to_string(T[i])+"K"; //cout<<name<<endl;
      ifstream inFile(name); if (!inFile) cout<<"read error\n";
      j = 0;
      while (getline(inFile, line)){
         istringstream isolated(line);
         isolated>>nu_str>>sigma_str;
         if (i==0) nua[j] = stod(nu_str);
         sigmaa[i][j] = stod(sigma_str);
         // cout<<j<<"\t"<<sigmaa[i][j]<<"\n";
         j++;
      }
      if (j!=nnu) printf("line count error\n");
      inFile.close();
   }

   /*  ifstream fa("./H2p/a"); if (!fa) cout<<"H2p/a read error\n";
   for(i=0;i<148;i++){
      getline(fa,line);
      istringstream isolated(line);
      string Taa_str, aa_str;
      isolated>>Taa_str>>aa_str;
      lnT[i] = log(stod(Taa_str));
      lna[i] = log(stod(aa_str));
   }
   fa.close(); */
}

void read_spec(double* nub, double* fluxb_cont, double* fluxb_line, string fspec_cont, string fspec_line){
   string line, nu_str, flux1_str, flux2_str, flux3_str, flux4_str, name;
   int i,j;
   //in file is actually wavelength, here just call it nu 
   name = "./spec_send/"+fspec_cont; //cout<<name<<endl;
   ifstream inFile(name); if (!inFile) cout<<"read error\n";
   j = 0;
   while (getline(inFile, line)){
      istringstream isolated(line);
      isolated>>nu_str>>flux1_str>>flux2_str>>flux3_str>>flux4_str;
      nub[j] = stod(nu_str);
      fluxb_cont[j] = stod(flux1_str);
      j++;
      //if (j==1 or j==nnu_spec) cout<<nu_str<<" "<<flux4_str<<" "<<flux2_str<<" "<<flux3_str<<" "<<flux4_str<<endl;
   }
   
   if (j!=nnu_spec) printf("line count error\n j=%d\n",j);
   inFile.close();

   name = "./spec_send/"+fspec_line; //cout<<name<<endl;
   inFile.open(name); if (!inFile) cout<<"read error\n";
   j = 0;
   while (getline(inFile, line)){
      istringstream isolated(line);
      isolated>>nu_str>>flux1_str>>flux2_str>>flux3_str>>flux4_str;
      //nub[j] = stod(nu_str);
      fluxb_line[j] = stod(flux1_str);
      j++;
      //if (j==1 or j==nnu_spec) cout<<nu_str<<" "<<flux4_str<<" "<<flux2_str<<" "<<flux3_str<<" "<<flux4_str<<endl;
   }
   if (j!=nnu_spec) printf("line count error\n j=%d\n",j);
   inFile.close();
}