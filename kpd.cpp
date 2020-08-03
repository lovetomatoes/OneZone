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
int const nnu_H2p = 45;
int const nT = 10;
int const nnu_spec = 1221;
int const nnu_Hm = 14000;

/* int main(){
   double T_rad = 1.e5;
   double k_Hm=0, k_H2p=0;
   double T0 = 10, Trat = 1.001, lgT; 
   double kra17=0, kra79;
   kpd_Hm_H2p(T_rad,k_Hm,k_H2p,true);
   double Ta[40], ka[40];
   read_k(40,Ta,ka);
   fstream f;
   f.open("Hm/kra.txt", ios::out | ios::trunc);
   f<<"kra17 kra79\n";
   while (T0<5.e4){
      kra(kra17, T0, 40, Ta, ka);
      lgT = log10(T0);
      if (T0<=6.e3) kra79 = pow(10., -17.845 + 0.762*lgT + 0.1523*pow(lgT,2) - 0.03274*pow(lgT,3) );
      else kra79 = pow(10., -16.4199 + 0.1998*pow(lgT,2) - 5.447e-3*pow(lgT,4) + 4.0415e-5*pow(lgT,6) );
      f<<T0<<" "<<kra17<<" "<<kra79<<endl;
      T0 *= Trat;
   }
   
   f.close();
}
 */
void kpd_Hm_H2p(double T_rad, double& k_Hm, double& k_H2p, bool spec){
   int i,j;
   // //double T_rad = 5.e4; 

   // if (spec) printf("IN KPD: using real spectrum\t");
   // else printf("IN KPD: T_rad=%5.2f\t", T_rad);
   double nurat0 = 1.0005;
   double nu_min = 0.45*eV/h_p;
   double nurat = nurat0, nu_prev = nu_min;
   double nu_eV;
  
   double* nub = NULL; double* fluxb_cont = NULL; double* fluxb_line = NULL;
   nub = new double [nnu_spec]; fluxb_cont = new double [nnu_spec]; fluxb_line = new double [nnu_spec];
   string fspec_cont = "rspec_pop3_cont_fe00.txt";
   string fspec_line = "rspec_pop3_line_fe00.txt";
   if (spec) read_spec(nub, fluxb_cont, fluxb_line, fspec_cont, fspec_line);

   double* nuc = NULL; double* sigmac = NULL;
   nuc = new double [nnu_Hm]; sigmac = new double [nnu_Hm];

   double** sigmaa = new double* [nT];
   for (i=0;i<nT;i++) sigmaa[i] = new double[nnu_H2p];
   double* nua = new double [nnu_H2p];
   double Ta[nT] = {0,2520,3150,4200,5040,6300,8400,12600,16800,25200};
   read_sigma(nua,sigmaa, nuc,sigmac);
   double sigma_Hm, sigma_H2p;
   double T_H2, sigma;
   
//     H2+ photodetachment
   T_H2 = 8.0e3;
   double E_ly = 12.4, E_lylim = 13.6, E_end = 15;
   int k, NF = 10000, k_ly;
   double nu[NF], Planck[NF], Flux[NF], Flux_c[NF], Flux_l[NF];
   
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
      Flux[k] = Planck[k]; Flux_c[k] = Planck[k];

      // interpolation from spectrum x:lambda in Angstron
      if (spec){
         linear(nub,fluxb_cont,nnu_spec,c/nu[k]/Angstron,Flux_c[k]);
         linear(nub,fluxb_line,nnu_spec,c/nu[k]/Angstron,Flux_l[k]);
         Flux[k] = Flux_c[k] + Flux_l[k];
         //Flux[k] = Flux_c[k];
      }

      nu_prev = nu[k];
   }

   
   //printf("NF=%d\n",NF);// wli note: NF=6820
   
   k_Hm=0.; k_H2p=0.;
   double dnu[NF-1];
   fstream f1;
   f1.open("../data/k_intg.txt", ios::out | ios::trunc );
   f1<<" k E Planck flux sigma_Hm sigma_H2p k_pd_Hm, k_pd_H2p\n";
   for (k=0;k<NF-1;k++){
      dnu[k] = nu[k+1]-nu[k];
      //printf("dnu/nu*c=%3.2e\n",dnu[k]*c/nu[k]/km); // = 150km/s
      if (nu[k] < E_lylim*eV/h_p){
         nu_eV = h_p*nu[k]/eV;
         Hm_CrossSec(sigma_Hm, nu_eV, nuc, sigmac);
         H2p_bf_CrossSec(sigma_H2p, nu_eV, T_H2, Ta, nua, sigmaa);         
         k_Hm   += sigma_Hm*4*pi*Flux[k]/h_p/nu[k]*dnu[k];
         k_H2p += sigma_H2p*4*pi*Flux[k]/h_p/nu[k]*dnu[k];
         f1<<" "<<k<<" "<<h_p*nu[k]/eV<<" "<<Planck[k]<<" "<<Flux[k]<<" "<<sigma_Hm<<" "
           <<sigma_H2p<<" "<<k_Hm<<" "<<k_H2p<<endl;
      }

   }
   
   //printf("norm for flux at 12.4eV: Planck:%3.2e Flux:%3.2e\n",Planck[k_ly],Flux[k_ly]);
   k_Hm *= (1.e-21/Flux_c[k_ly]);
   k_H2p *= (1.e-21/Flux_c[k_ly]);
   //printf("k_Hm=%3.2e, k_H2p=%3.2e, kHm/kH2=%3.2e\n", k_Hm, k_H2p, k_Hm/1.39e-12);
   // printf("kHm/kH2=%3.2e\n", k_Hm/1.39e-12);

   for(i=0; i<nT; i++) delete[] sigmaa[i];
   delete[] sigmaa;
   delete[] nua;
   delete[] nuc; delete[] sigmac;
   delete[] nub; delete[] fluxb_cont; delete[] fluxb_line;
}

/* !     H- photodetachment
!     H-   +  ph.   ->   H    +    e
!     photon energy E in eV
!     McLaughlin et al 2017
 */

void Hm_CrossSec(double& sigma, double nu_eV, double* nuc, double* sigmac){
   if (nu_eV<nuc[0]) sigma = 0;
   else if (nu_eV>13.6) sigma = 0;
   else linear(nuc,sigmac,nnu_Hm,nu_eV,sigma);
   sigma = sigma*Mb;
   
   /* //John (1988 A&A 193,189)
   //printf("%3.2e %3.2e ",nu_eV,sigma);
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
   rlmbd=1.23985/nu_eV;
   if(rlmbd<=rlmbd_0) {
      x=sqrt(1./rlmbd-1./rlmbd_0);
      f=C1+x*(C2+x*(C3+x*(C4+x*(C5+x*C6))));
      sigma=1.e-18*pow(rlmbd,3)*pow(x,3)*f;
   }
   else sigma = 0.0;
   if(rlmbd<=0.09) sigma = 0.0;
   //printf("%3.2e\n",sigma); */
}

/* // ref. T.L.John (1988 A&A 193,189)
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
}*/

/* C     cross section of H2+ b-f
C     H2+   +  ph.   ->   H   +   H+
C     frequency rnu in eV
C     ref. P.C.Stancil (1994 ApJ 430,360)
 */
void H2p_bf_CrossSec(double& sigma, double nu_eV, double T_K, double* Ta, double* nua, double** sigmaa){
//         H2p_bf_CrossSec(sigma_H2p, nu[k], T_H2, Ta, nua, sigmaa);
  /*  cout<<"IN H2P CROSS_SEC:  ";
   for(j=0; j<nnu_H2p; j++) cout<<nua[j]<<"\t";
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
   else if (nu_eV>nua[nnu_H2p-1]) sigma=0.0;
   else bilinear(Ta,nua,sigmaa,nT,nnu_H2p,T_K,nu_eV,sigma);
   
}

void linear(double* xa, double* ya, int m, double x, double& y){
   int ms;
   int i,j;
   double y1, y2, t;
   //if( x<xa[0] or x>xa[m-1]) printf("xa_min=%3.2e, x=%3.2e, xa_max=%3.2e\n", xa[0],x,xa[m-1]);
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

void read_sigma(double* nua, double** sigmaa, double* nuc, double* sigmac){
   string line, nu_str, sigma_str, name;
   int i,j;
// H2p
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
      if (j!=nnu_H2p) printf("line count error\n");
      inFile.close();
   }

// Hm
   name = "./Hm/JPBaa6c1fSuppdata1.dat";
   ifstream inFile(name); if (!inFile) cout<<"read error\n";
   for (j=0;j<4;j++) getline(inFile, line);
   for (j=0;j<nnu_Hm;j++){
      getline(inFile, line);
      istringstream isolated(line);
      isolated>>nu_str>>sigma_str;
      nuc[j] = stod(nu_str);
      sigmac[j] = stod(sigma_str);    
   }
   inFile.close();
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

///**************************************************************************************
///                         H- from radiative attatchment  H + e- --> H- + ph  

void read_k(int n_ra, double* Ta, double* ka){
   string line, T_str, k_str, name;
   int i,j;
// Hm
   name = "./Hm/JPBaa6c1fSuppdata2.dat";
   ifstream inFile(name); if (!inFile) cout<<"read error\n";
   for (j=0;j<3;j++) getline(inFile, line);
   for (j=0;j<n_ra;j++){
      getline(inFile, line);
      istringstream isolated(line);
      isolated>>T_str>>k_str;
      //cout<<j<<" "<<T_str<<" "<<k_str<<endl;
      Ta[j] = stod(T_str);
      ka[j] = stod(k_str);
   }
   inFile.close();
}

void kra(double& y, double T, int n_ra, double* Ta, double* ka){
   linear(Ta, ka, n_ra, T, y);
}
