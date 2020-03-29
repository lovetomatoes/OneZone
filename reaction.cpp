// calculate the reaction coefficients  reaction rates
// 6 major reactions, 5 species
// update to 9 species; 42+2 reactions (last 2: ph dissociation of H2  Hm)

#include <iostream>
#include <stdio.h>
#include <cmath>
#include<algorithm>

#include "reaction.h"
#include "PARA.h"
using namespace std;

void react_coef(double *k, double nH, double y_H, double y_H2, double T_K, double J_LW, double Tb){
    /* double T_eV = 8.61735e-5*T_K;
    double lnT_eV = log(T_eV);
    double T300 = T_K/300.;

//*********************************************************************************************************
// GLOVER  ABEL 2008
//  1)   H+    +   e     ->   H     +   ph.
//     Glover  Abel 2008 reation (13), Table A1
//     Ferland et al. 1992; case B
    k[1] = 2.753e-14*pow(315614./T_K, 1.500)*pow((1.+pow(115188./T_K, 0.407)),-2.242);

//  2)   H     +   e     ->   H-    +   ph.
//     Glover  Abel 2008, reation (1), Table A1
    if (T_K<=6000) k[2] = pow(10, -17.845 + 0.762*lgT + 0.1523*pow(lgT,2) - 0.03274*pow(lgT,3) );
    else k[2] = pow(10, -16.4199 + 0.1998*pow(lgT, 2) - 5.447e-3*pow(lgT,4) + 4.0415e-5*pow(lgT,6) );

//  3)   H-    +   H     ->   H2    +   e              
//     Glover  Abel 2008 reation (2), Table A1; Eq(5)
    k[3] = 1.3e-9;

//  4)   H2  +   H    <-> 3 H  
//     Glover  Abel 2008, reation (9) in Table A1; 
//     He therein not included
    double ncrH, ncrH2, ncr;
    double kLTE, kv0;
    ncrH = pow(10, 3-0.416*log10(T_K/1.e4)-0.327*pow(log10(T_K/1.e4),2));
    ncrH2 = pow(10, 4.845-1.3*log10(T_K/1.e4)+1.62*pow(log10(T_K/1.e4), 2));
    ncr = 1./ (y_H/ncrH+ y_H2/ncrH2); // Eq(14)
    kv0  = 6.67e-12*sqrt(T_K)*exp(-1-63593./T_K);
    kLTE = 3.52e-9*exp(-43900./T_K);
    k[4] = pow(10, nH/ncr/(1+nH/ncr)*log10(kLTE) + 1./(1+nH/ncr)*log10(kv0));
//    k[4] = 0;

//  5)   3 H       ->     H2     +      H
//    Glover  Abel 2008 reaction (30), Table A1
//    in text follow Palla, Salpeter  Stahler (1983)
    k[5] = 5.5e-29/T_K;
    //k[5] = 0.;
    // k(15)=k(7)*zH2/pow(zH,2)*1.493e-20/pow(T_K,1.5)*exp(chi_H2/T_K);

//  6)    H2    +     J21        ->         2 H  
// Wolcott-Green  Haiman 2011
    // mu = 1.2 in PARA.h
    double R_core = sqrt(pi*k_B*T_K/(G*m_H*nH*mu*m_H));
    double fsh;
//WG 11:
    double x, b5;
    x = nH*y_H2*R_core/5.e14;
    b5 = sqrt(k_B*T_K/m_H)/1.e5;
    fsh = 0.965/pow(1+x/b5,1.1) + 0.035/sqrt(1+x)*exp(-8.5e-4*sqrt(1+x)); 
//DB 96: 
    //fsh = min(1, pow(nH*y_H2*R_core/1.e14,-0.75)); 
//   tried fsh = 1, fsh严重影响pd给出的y_H2,equi
    double kh2pd = 1.39e-12 * J_LW;
    k[6] = kh2pd * fsh;
    //k[6] = 0;

// 7)    H-    +     ph       ->   H    +     e
    double ratioT = 0;
    if (Tb == 8.e3) ratioT = 8.7e5;
    if (Tb == 1.e4) ratioT = 4.6e4; // T=10^4K, ratioT = 4.6e4; 
    if (Tb == 2.e4) ratioT = 2.1e2;
    if (Tb == 3.e4) ratioT = 4.6e1;
    if (Tb == 5.e4) ratioT = 1.7e1;
    if (Tb == 1.e5) ratioT = 10; //T=10^5K, ratioT = 10.
    if (Tb == 2.e5) ratioT = 8.1;

    k[7] = ratioT*kh2pd;
    
    //k[7] = 0;

// 8)    H    +     e-       ->   H+    +     2e-
    // Abel et al. (1997)
    k[8] = exp(-32.71396786+(13.536556
           +(-5.73932875+(1.56315498+(-0.2877056
           +(3.48255977e-2+(-2.63197617e-3
           +(1.11954395e-4-2.03914985e-6*lnT_eV)*lnT_eV)
           *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV);
    //k[8] = 0.;

 */

////////////////////////////////////////////////
//     Chemical reaction rate coefficients    //
//     mostly from GLover  Abel (2008)        //
////////////////////////////////////////////////
    double gamma_h1, gamma_l1, gamma_h2, gamma_l2, nc_1, nc_2, 
           pp, gamma, xk_L, xk_H, n_cr, a, n_cr_H, n_cr_H2, 
           zH, zH2, chi_H2, z_n, xk_rr_A, xk_rr_B, xk_di;
    double T_eV = 8.61735e-5*T_K;
    double lnT_eV = log(T_eV);
    double lgT = log10(T_K);
    double lgT4 = lgT - 4.;
    double T300 = T_K/300.;
    int n_max=5, n;

    chi_H2  =  5.1965e4; //binding energy of H2

    n_cr_H  = pow( 10., 3.0 - 0.416*lgT4 - 0.327*pow(lgT4,2) );
    n_cr_H2 = pow( 10., 4.845 - 1.3*lgT4 + 1.62*pow(lgT4,2) );
    n_cr = 1.0/(y_H/n_cr_H + 2.0*y_H2/n_cr_H2); // Eq(14) GA08
    a=1./(1.+nH/n_cr);

    zH=0.0;
    for (n=1; n<=n_max; n++){
        z_n = 2.*pow(n,2)*exp(-1.57798e5*(1.-1./pow(n,2))/T_K);
        zH = zH+z_n;
    }

    zH2 = pow( 10., 2.20859 -1.8089*lgT + 0.451858* pow(lgT,2.0) ); //
    //printf("IN REACTION.CPP, n_cr=%3.2e, zH=%3.2e, zH2=%3.2e\n",n_cr,zH,zH2);
//================= Hydrogen (1-20) =================//
//
//  (1)   H     +   e     ->   H+    + 2 e   
//     Abel et al. (1997): GA08 TableA1-12
    k[1]= exp(-32.71396786 + (13.536556
        +(-5.73932875 + (1.56315498+(-0.2877056
        +(3.48255977e-2 + (-2.63197617e-3
        +(1.11954395e-4 - 2.03914985e-6*lnT_eV)*lnT_eV)
        *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV);

//  (2)   H+    +   e     ->   H     +   ph.
//     Glover  Jappsen (2007): GA08 TableA1-13
    //     (Case A)
    //       k[2]=1.269d-13*((315614d0/T_K)**(1.503d0))
    //           *((1.0d0+(604625d0/T_K)**(0.47d0))**(-1.923d0))
    //     (Case B)
    k[2] = 2.753e-14* pow( 315614./T_K , 1.5) * pow( 1.0 + pow(115188./T_K, 0.407) ,-2.242);

//  (3)   H     +   e     ->   H-    +   ph. 
//     Wishart (1979): GA08 TableA1-1
    if (T_K<=6.e3) k[3] = pow(10., -17.845 + 0.762*lgT + 0.1523*pow(lgT,2) - 0.03274*pow(lgT,3) );
    else k[3] = pow(10., -16.4199 + 0.1998*pow(lgT,2) - 5.447e-3*pow(lgT,4) + 4.0415e-5*pow(lgT,6) );
    //k[3] *= 0.5;
    k[3] = 1.4e-18*pow(T_K,0.928)*exp(-T_K/16200.); //Galli Palla 1998

//  (4)   H-    +   H     ->   H2    +   e
//     Kreckel et al. (2010): 
    //     The most recent experiment
    k[4] = 1.35e-9 * ( pow(T_K,9.8493e-2)
           + 0.32852*pow(T_K,0.5561) 
           + 2.771e-7*pow(T_K,2.1826))
           /(1.0 + 6.1910e-3*pow(T_K,1.0461)
           + 8.9712e-11*pow(T_K,3.0424)  
           + 3.2576e-14*pow(T_K,3.7741));
      
// 5)   H     +   H+    ->   H2+   +   ph.
//     Coppola et al. (2011): // wli: did not find... but plotted checked
    if(T_K<=30.0) k[5] = 2.1e-20/pow( T_K/30.0, 0.15); 
    else k[5] = pow(10., (-18.20 - 3.194*lgT + 1.786* pow(lgT,2) - 0.2072*pow(lgT,3) ) );

//  (6)   H2+   +   H     ->   H2    +   H+         
//     Galli  Palla (1998): GA08 TableA1-4
    k[6]=6.4e-10;
      
//  (7)   H2    +   H     -> 3 H    // cd1     
//     GA08 TableA1-9
    xk_L  = 6.67e-12*sqrt(T_K)/exp(1 + 63593./T_K);
    xk_H = 3.52e-9*exp(-43900./T_K);
//     connect between v=0 and LTE  
    if(a==1.) k[7] = xk_L;
    else if (a==0.) k[7] = xk_H;
    else k[7] = pow(xk_H, 1.-a) * pow(xk_L,a);
// xk(7) below same as in Kohei's react_coef.f: got same result w/ fort.10; not used
/* //     Martin et al. (1996)
    gamma_h1=-1.784239e2-6.842243e1*log10(T_K)
        +4.320243e1*pow(lgT,2)-4.633167*pow(lgT,3)
        +6.970086e1*log10(1.0+4.087038e4/T_K);
    gamma_h2=-2.370570e4/T_K;
    gamma_l1=1.288953e2-5.391334e1*lgT
        +5.315517*pow(lgT,2)
        -1.973427e1*log10(1.0+1.678095e4/T_K);
    gamma_l2=-2.578611e4/T_K;
    nc_1=1.482123e1-4.890915*lgT
        +4.749030e-1*pow(lgT,2)-1.338283e2/T_K;
    nc_2=-1.164408+nc_1;
    pp=8.227443e-1+5.864073e-1*exp(-T_K/1850)
        -2.056313*exp(-T_K/440);
    nc_1=pow(1.e1,nc_1);
    nc_2=pow(1.0e1,nc_2);
    gamma=gamma_h1
        -(gamma_h1-gamma_l1)/(1.0+pow(nH/nc_1,pp))
        +gamma_h2
        -(gamma_h2-gamma_l2)/(1.0+pow(nH/nc_2,pp));
    //k[7]=pow(1.0e1,gamma);
 */

//  (8)   H2    +   H+    ->   H2+   +   H        
//     Savin (2004) low density: GA08 TableA1-7 // wli in table should be loge(T_K) not log10(T_K)
    xk_L = (-3.3232183e-7 + (3.3735382e-7 + (-1.4491368e-7 
         + (3.4172805e-8  + (-4.7813720e-9 + (3.9731542e-10 
         + (-1.8171411e-11 + 3.5311932e-13*log(T_K))
         *log(T_K))*log(T_K))*log(T_K))*log(T_K))*log(T_K))*log(T_K) )
         *exp(-21237.15/T_K);
//     Coppola et al.(2011) high density Table3-4
    xk_H=exp(-33.081+6.3173e-5*T_K-2.3478e4/T_K
        -1.8691e-9*pow(T_K,2) )*1.e6;
// wli if T_K <30.0 xk_L wrong...
    if (xk_L<0.) xk_L = -xk_L; 
//     connect between v=0 and LTE
    if(a==1.) k[8] = xk_L;
    else if (a==0.) k[8] = xk_H;
    else k[8] = pow(xk_H, 1.-a) * pow(xk_L,a);

//  (9)   H2    +   e     -> 2 H     +   e
//     Trevisan  Tennyson (2002): GA08 TableA1-8
    xk_L = 4.49e-9*pow(T_K,0.11) *exp(-101858./T_K);
    xk_H = 1.91e-9*pow(T_K,0.136)*exp(-53407.1/T_K);
//     connect between v=0 and LTE
    if(a==1.) k[9] = xk_L;
    else if (a==0.) k[9] = xk_H;
    else k[9] = pow(xk_H, 1.-a) * pow(xk_L,a);

//  (10)   H-    +   e     ->   H     + 2 e              
//     Janev et al. (1987): GA08 TableA1-7
    k[10]=exp(-18.01849334+(2.3608522
        +(-0.2827443+(1.62331664e-2+(-3.36501203e-2
        +(1.17832978e-2+(-1.65619470e-3+(1.0682752e-4
        -2.63128581e-6*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)
        *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV);

//  (11)   H-    +   H+    -> 2 H
//     Croft, Dickinson  Gadea (1999): GA08 TableA1-5
    k[11] = 2.4e-6/sqrt(T_K)*(1.0+T_K/2.e4);

//  (12)   H-    +   H+    ->   H2+   +   e          
//     Poulaert et al. (1978): GA08 TableA1-16
    if(T_K<=8.e3) k[12] = 6.9e-9/pow(T_K,0.35);
    else k[12] = 9.6e-7/pow(T_K,0.9);

//  (13)   H2+   +   e     -> 2 H                
//     Savin et al. (2004): GA08 TableA1-6
    if(T_K<=617.0) k[13]=1.e-8;
    else k[13]=1.32e-6*pow(T_K,-0.76);

//  (14)   H2+   +   H-    ->   H2    +   H         
//     Dalgarno  Lepp (1987): GA08 TableA1-21
    k[14] = 1.40e-7/sqrt(T_K/300.);

//  (15) 3 H               ->   H2    +   H                
//     inverse rate of dissociation //Omukai 01: Table2 not used;
    // k[15] = k[7]*zH2/pow(zH,2)*1.493e-20 /pow(T_K,1.5)*exp(chi_H2/T_K); abandoned since k[7] large error at low T
//     GA08 TableA1-30 text following Abel 2002, lowest rates among literature; also default value in Glover & Savin 2009 section 3.1.7
    if (T_K<=300.) k[15] = 1.14e-31*pow(T_K,-0.38);
    else k[15] = 3.9e-30/T_K;
//  (16) 2 H2           -> 2 H    +   H2
//     Martin+1998; Shapiro&Kang 1987 GA08 TableA1-10 
    xk_L=5.996e-30*pow(T_K,4.1881) / pow(1.0+6.761e-6*T_K, 5.6881) *exp(-54657.4/T_K);
    xk_H=1.3e-9*exp(-53300.0/T_K);
    //     connect between v=0 and LTE
    if(a==1.) k[16] = xk_L;
    else if(a==0.) k[16] = xk_H;
    else k[16] = pow(xk_H, 1.-a) * pow(xk_L,a);

//  (17) 2 H    +   H2  -> 2 H2
//     inverse rate of dissociation //Omukai 01: Table2; not used afraid error at low T like collision w/ H;
    //k[17] = k[16]*zH2/pow(zH,2)*1.493e-20/pow(T_K,1.5)*exp(chi_H2/T_K);
// Palla Salpeter Stahler 1983 
    k[17] = 6.9e-30/T_K;
//     GA08 TableA1-31 following Palla, Salpeter& Stahler 1983 actually quoting Jacobs+1967
//    k[17] = k[15]/8.;
//  (18) H-   +   H  -> 2 H   +  e 
//     Janev et al. (1987): GA08 TableA1-15
    if(T_eV<=0.1) k[18] = 2.5634e-9 * pow(T_eV,1.78186);
    else
        k[18] = exp( -2.0372609e1 
              + (1.13944933+(-1.4210135e-1 + (8.4644554e-3
              +(-1.4327641e-3+(2.0122503e-4+(8.6639632e-5
              +(-2.5850097e-5 + (2.4555012e-6 - 8.0683825e-8*lnT_eV)
              *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV );

//  (19) H-   +   H2+  -> 3 H
//     Dalgarno  Lepp (1987): GA08 TableA1-22
    k[19] = 1.4e-7/sqrt(T_K/3.e2);

//  (20) H2   +   e   ->  H-  +  H
//     Schulz  Asundi (1967): GA08 TableA1-23
    k[20] = 2.7e-8/pow(T_K,1.27)*exp(-43000./T_K);

//  (21) H2   +   ph  -> 2 H
// Wolcott-Green  Haiman 2011
    double R_core = sqrt(pi*k_B*T_K/(G*m_H*nH*mu*m_H));
    double fsh;
//WG 11 shielding factor:
    double x, b5;
    x = nH*y_H2*R_core/5.e14;
    b5 = sqrt(k_B*T_K/m_H)/1.e5;
    fsh = 0.965/pow(1+x/b5,1.1) + 0.035/sqrt(1+x)*exp(-8.5e-4*sqrt(1+x)); 
//DB 96 shielding factor: 
    //fsh = min(1, pow(nH*y_H2*R_core/1.e14,-0.75)); 
    //fsh = 1; // fsh严重影响pd给出的y_H2,equi
    double kh2pd = 1.39e-12 * J_LW;
    k[21] = kh2pd * fsh;

//  (22) H-   +   ph  ->  H   +  e
    double ratioT = 0;
    if (Tb == 8.e3) ratioT = 8.7e5;
    if (Tb == 1.e4) ratioT = 4.6e4; // T=10^4K, ratioT = 4.6e4; 
    if (Tb == 2.e4) ratioT = 2.1e2;
    if (Tb == 3.e4) ratioT = 4.6e1;
    if (Tb == 5.e4) ratioT = 1.7e1;
    if (Tb == 1.e5) ratioT = 10; //T=10^5K, ratioT = 10.
    if (Tb == 2.e5) ratioT = 8.1;
    //printf("from RATIO = %3.2e\n",ratioT*kh2pd);

//  (23) H2+  +   ph  ->  H   +  H+

//  (24-30) blank:
    for (n=24;n<=30;n++) k[n]=0.0;

    
//================= Helium (31-) =================// 

//  (31) He   +   e   ->   He+   +   2 e
//     GA08 TableA1-17
    k[31] = exp(-4.409864886e1 + (2.391596563e1 + (-1.07532302e1 
          + (3.05803875 + (-5.68511890e-1 + (6.79539123e-2 
          + (-5.00905610e-3 + (2.06723616e-4 - 3.64916141e-6*lnT_eV)
          *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV );

// 32) He+  +   e   ->   He   +   ph.
//     GA08 TableA1-19 //wli: but not knowing why 0.68:0.32
    xk_rr_A = 1.e-11/sqrt(T_K)
        *( 12.72 - 1.615*lgT
        -  0.3162*pow(lgT,2)
        +  0.0493*pow(lgT,3));
    xk_rr_B = 1.e-11/sqrt(T_K)
        *( 11.19 - 1.676*lgT
        -  0.2852*pow(lgT,2)
        +  4.433e-2*pow(lgT,3));
    xk_di = 1.9e-3/pow(T_K,1.5)*exp(-473421./T_K)
        *(1.0 + 0.3*exp(-94684/T_K));

    k[32] = 0.68*xk_rr_A + 0.32*xk_rr_B + xk_di;

//  (33) He+  +   e   ->   He++  +   2 e
//     GA08 TableA1-18
    k[33] = exp(-6.87104099e1
        + (4.393347633e1 + (-1.848066990e1 + (4.701626490 + (-7.6924663e-1 
        + (8.113042e-2 + (-5.32402063e-3 + (1.97570531e-4 - 3.16558106e-6*lnT_eV)
        *lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV)*lnT_eV  );

// 34) He++  +   e   ->   He+  +   ph.  
//     GA08 TableA1-20 //wli: not knowing why Case B
    k[34] = 5.506e-14*pow(1262456.0/T_K,1.5)
        /pow(1.0+pow(460752.0/T_K,0.407),2.242);

//  (35) H2   +   He   ->   He   +  H  +  H
//     GA08 TableA1-11
    xk_L = pow(10., -27.029+3.801*lgT-29487.0/T_K);
    xk_H = pow(10., -2.729-1.75*lgT-23474.0/T_K);
    if(a==1.) k[35]=xk_L;
    else if(a==0.) k[35]=xk_H;
    else k[35]= pow(xk_H, 1.-a) * pow(xk_L,a);


//  (36) H2   +   He+   ->   He   +  H  +  H+
//     GA08 TableA1-24
    k[36] = 3.7e-14*exp(35.0/T_K);

//  (37) H2   +   He+   ->   H2+   +   He
//     GA08 TableA1-25
    k[37] = 7.2e-15;

//  (38) He+  +   H   ->   He   +   H+
//     GA08 TableA1-26
    k[38] = 1.2e-15*pow(T_K/3.0e2,0.25);

//  (39) He   +   H+  ->   He+  +   H
//     GA08 TableA1-27
    if(T_K<=1.0e4) k[39] = 1.26e-9/pow(T_K,0.75)*exp(-127500./T_K);
    else k[39] = 4.e-37*pow(T_K,4.74);
    
//  (40) He+  +   H-  ->   He   +   H
//     GA08 TableA1-28
    k[40] = 2.32e-7/pow(T_K/300.,0.52) * exp(T_K/22400.);

//  (41) He   +   H-  ->   He   +   H   +   e
//     GA08 TableA1-29
    k[41] = 4.1e-17*pow(T_K,2)*exp(-19870.0/T_K);

//  (42) H  +  H  +  He  ->  H2  +  He
//     GA08 TableA1-32
    k[42] = 6.9e-32/pow(T_K,0.4);

    for (int i=0;i<N_react1;i++) {
        if (not isfinite(k[i])) k[i]=0.;
        if (k[i]<0.) k[i]=0.;
    }
}

void react_rat(double *r_f_tot, double *y, double *k, double nH, double T_K){
//****************************************************************
//*     dy(i)/dt=r_f_tot(i)                                      * 
//*     This subroutine returns reaction rate for each spieces,  *
//*     r_f_tot(i).                                              *
//****************************************************************
//      N_sp = number of spiecies
//      N_react = number of reactions

    double r_f[N_react1][N_sp1];
    double rate;    
    int isp, ire;

    for (int isp=1; isp<N_sp1; isp++){
        for (int ire=1; ire<N_react1; ire++){
            r_f[ire][isp] = 0.;
        }        
    }

/****************************************************************
*     dy(i)/dt=r_f_tot(i)                                      * 
*     This subroutine returns reaction rate for each spieces,  *
*     r_f_tot(i).                                              *
****************************************************************/
/****************************************************
*     SPECIES                                      *
*     1 : H      2 : H2     3 : e      4 : H+      *
*     5 : H2+    6 : H-                            *
*     7 : He     8 : He+    9 : He++               *
****************************************************/

//  1)   H     +   e     ->   H+    +  2 e
//       1         3          4        2 * 3
    rate  = k[1] * y[1] * y[3] * nH;
    r_f[1][1] = -rate;
    r_f[1][3] = rate;
    r_f[1][4] = rate;

//  2)   H+    +   e     ->   H     +   ph.
//       4         3          1
    rate  = k[2] * y[3] * y[4] * nH;
    r_f[2][1] = rate;
    r_f[2][3] = -rate;
    r_f[2][4] = -rate;

//  3)   H     +   e     ->   H-    +   ph.
//       1         3          6
    rate  = k[3] * y[1] * y[3] * nH;
    r_f[3][1] = -rate;
    r_f[3][3] = -rate;
    r_f[3][6] = rate;

//  4)   H-    +   H     ->   H2    +   e
//       6         1          2         3
    rate  = k[4] * y[1] * y[6] * nH;
    r_f[4][1] = -rate;
    r_f[4][2] = rate;
    r_f[4][3] = rate;
    r_f[4][6] = -rate;

//  5)   H     +   H+    ->   H2+   +   ph.
//       1         4          5
    rate  = k[5] * y[1] * y[4] * nH;
    r_f[5][1] = -rate;
    r_f[5][4] = -rate;
    r_f[5][5] = rate;

//  6)   H2+   +   H     ->   H2    +   H+
//       5         1          2         4
    rate  = k[6] * y[1] * y[5] * nH;
    r_f[6][1] = -rate;
    r_f[6][2] = rate;
    r_f[6][4] = rate;
    r_f[6][5] = -rate;

//  7)   H2    +   H     -> 3 H
//       2         1        3 * 1
    rate  = k[7] * y[1] * y[2] * nH;
    r_f[7][1] = 2.*rate;
    r_f[7][2] = -rate;

//  8)   H2    +   H+    ->   H2+   +   H
//       2         4          5         1
    rate  = k[8] * y[2] * y[4] * nH;
    r_f[8][1] = rate;
    r_f[8][2] = -rate;
    r_f[8][4] = -rate;
    r_f[8][5] = rate;

//  9)   H2    +   e     -> 2 H     +   e
//       2         3        2 * 1       3
    rate  = k[9] * y[2] * y[3] * nH;
    r_f[9][1] = 2.*rate;
    r_f[9][2] = -rate;

//  10)  H-    +   e     ->   H     +  2 e
//       6         3          1        2 * 3
    rate  = k[10] * y[3] * y[6] * nH;
    r_f[10][1] = rate;
    r_f[10][3] = rate;
    r_f[10][6] = -rate;

//  11)  H-    +   H+    -> 2 H
//       6         4        2 * 1
    rate  = k[11] * y[4] * y[6] * nH;
    r_f[11][1] = 2.*rate;
    r_f[11][4] = -rate;
    r_f[11][6] = -rate;

//  12)  H-    +   H+    ->   H2+   +   e
//       6         4          5         3
    rate  = k[12] * y[4] * y[6] * nH;
    r_f[12][3] = rate;
    r_f[12][4] = -rate;
    r_f[12][5] = rate;
    r_f[12][6] =-rate;

//  13)  H2+   +   e     -> 2 H
//       5         3        2 * 1
    rate  = k[13] * y[3] * y[5] * nH;
    r_f[13][1] = 2.*rate;
    r_f[13][3] = -rate;
    r_f[13][5] = -rate;

//  14)  H2+   +   H-    ->   H2    +   H
//       5         6          2         1
    rate  = k[14] * y[5] * y[6] * nH;
    r_f[14][1] = rate;
    r_f[14][2] = rate;
    r_f[14][5] = -rate;
    r_f[14][6] = -rate;

//  15) 3 H              ->   H2    +   H
//      3 * 1                 2         1
    rate  = k[15] * y[1] * y[1] * y[1] * nH * nH;
    r_f[15][1] = -2.*rate;
    r_f[15][2] = rate;

//  16) 2 H2             -> 2 H     +   H2
//      2 * 2               2 * 1       2
    rate  = k[16] * y[2] * y[2] * nH;
    r_f[16][1] = 2.*rate;
    r_f[16][2] = -rate;

//  17) 2 H    +   H2    -> 2 H2
//      2 * 1      2        2 * 2
    rate  = k[17] * y[1] * y[1] * y[2] * nH * nH;
    r_f[17][1] = -2.*rate;
    r_f[17][2] = rate;

//  18)  H-    +   H     -> 2 H     +   e
//       6         1        2 * 1       3
    rate  = k[18] * y[1] * y[6] * nH;
    r_f[18][1] = rate;
    r_f[18][3] = rate;
    r_f[18][6] = -rate;

//  19)  H-    +   H2+   -> 3 H
//       6         5        3 * 1
    rate  = k[19] * y[5] * y[6] * nH;
    r_f[19][1] = 3.*rate;
    r_f[19][5] = -rate;
    r_f[19][6] = -rate;

//  20)  H2    +   e     ->   H-    +   H
//       2         3          6         1
    rate  = k[20] * y[2] * y[3] * nH;
    r_f[20][1] = rate;
    r_f[20][2] = -rate;
    r_f[20][3] = -rate;
    r_f[20][6] = rate;

//  21)  H2   +   ph   -> 2 H
//       2                  1
    rate  = k[21] * y[2];
    r_f[21][1] = 2.*rate;
    r_f[21][2] = -rate;

//  22)  H-   +   ph  ->  H   +  e
//       6                1      3
    rate  = k[22] * y[6]; // wli  rate = 0.;
    r_f[22][1] = rate;
    r_f[22][3] = rate;
    r_f[22][6] = -rate;

//  23)  H2+  +   ph  ->  H   +  H+
//       5                1      4
    rate  = k[23] * y[5];
    r_f[23][1] = rate;
    r_f[23][4] = rate;
    r_f[23][5] = -rate;

//  31)  He    +   e     ->   He+   +  2 e
//       7         3          8        2 * 3
    rate  = k[31] * y[3] * y[7] * nH;
    r_f[31][3] = rate;
    r_f[31][7] = -rate;
    r_f[31][8] = rate;

//  32)  He+   +   e     ->   He    +   ph.
//       8         3          7
    rate  = k[32] * y[3] * y[8] * nH;
    r_f[32][3] = -rate;
    r_f[32][7] = rate;
    r_f[32][8] = -rate;

//  33)  He+   +   e     ->   He++  +  2 e
//       8         3          9        2 * 3
    rate  = k[33] * y[3] * y[8] * nH;
    r_f[33][3] = rate;
    r_f[33][8] = -rate;
    r_f[33][9] = rate;

//  34)  He++  +   e     ->   He+   +   ph.
//       9         3          8
    rate  = k[34] * y[3] * y[9] * nH;
    r_f[34][3] = -rate;
    r_f[34][8] = rate;
    r_f[34][9] = -rate;

//  35)  H2    +   He    ->   He    +  2 H
//       2         7          7        2 * 1
    rate  = k[35] * y[2] * y[7] * nH;
    r_f[35][1] = 2.*rate;
    r_f[35][2] = -rate;

//  36)  H2    +   He+   ->   He    +   H   +   H+
//       2         8          7         1       4
    rate  = k[36] * y[2] * y[8] * nH;
    r_f[36][1] = rate;
    r_f[36][2] = -rate;
    r_f[36][4] = rate;
    r_f[36][7] = rate;
    r_f[36][8] = -rate;

//  37)  H2    +   He+   ->   H2+   +   He
//       2         8          5         7
    rate  = k[37] * y[2] * y[8] * nH;
    r_f[37][2] = -rate;
    r_f[37][5] = rate;
    r_f[37][7] = rate;
    r_f[37][8] = -rate;

//  38)  He+   +   H     ->   He    +   H+
//       8         1          7         4
    rate  = k[38] * y[1] * y[8] * nH;
    r_f[38][1] = -rate;
    r_f[38][4] = rate;
    r_f[38][7] = rate;
    r_f[38][8] = -rate;

//  39)  He    +   H+    ->   He+   +   H
//       7         4          8         1
    rate  = k[39] * y[4] * y[7] * nH;
    r_f[39][1] = rate;
    r_f[39][4] = -rate;
    r_f[39][7] = -rate;
    r_f[39][8] = rate;

//  40)  He+   +   H-    ->   He    +   H
//       8         6          7         1
    rate  = k[40] * y[6] * y[8] * nH;
    r_f[40][1] = rate;
    r_f[40][6] = -rate;
    r_f[40][7] = rate;
    r_f[40][8] = -rate;

//  41)  He    +   H-    ->   He    +   H   +   e
//       7         6          7         1       3
    rate  = k[41] * y[6] * y[7] * nH;
    r_f[41][1] = rate;
    r_f[41][3] = rate;
    r_f[41][6] = -rate;

//  42) 2 H    +   He    ->   H2    +   He
//      2 * 1      7          2         7
    rate  = k[42] * y[1] * y[1] * y[7] * nH * nH;
    r_f[42][1] =-2.*rate;
    r_f[42][2] = rate;

//************************************************************
    for (int isp=1; isp<N_sp1; isp++){
        r_f_tot[isp] = 0.;
        for (int ire=1; ire<N_react1; ire++){
            r_f_tot[isp] += r_f[ire][isp];
        }
    }
}



//****************************************************************
//*     SPECIES                                                  *
//*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
//****************************************************************

    //y_Hm=k(2)*y[3]/k[3] 没明白这里H-为什么不算在species 里面

//********* primordial gas reactions *********
/* 
//   1)   H+    +   e     ->   H     +   ph.
//         (4         3          1)
    rate = k[1] * y[4] * y[3] * nH;
    r_f[1][3] = -rate;
    r_f[1][4] = -rate;
    r_f[1][1] = rate;
//   2)   H     +   e     ->   H-    +   ph. 
//         (1         3          x)
    rate = k[2] * y[1] * y[3] * nH;
    r_f[2][1] = -rate;
    r_f[2][3] = -rate;
    r_f[2][5] = rate;
//   3)   H-    +   H     ->   H2    +   e       
//         (x         1          2         3)
    rate = k[3] * y[1] * y[5] * nH;
    r_f[3][1] = -rate;
    r_f[3][5] = -rate;
    r_f[3][2] = rate;
    r_f[3][3] = rate;
    //cout<<" H- : \t"<<rate<<endl;

//   4)   H2    +   H     -> 3 H    
//         (2         1        3*1)
    rate = k[4] * y[2] * y[1] * nH;
    r_f[4][2] = -rate;
    r_f[4][1] = 2. * rate;
//   5) 3 H               ->   H2    +   H  
//       (3*1                    2         1)
    rate = k[5] * pow(y[1],3) * pow(nH,2);
    r_f[5][1] = -2. * rate;
    r_f[5][2] = rate;
    //cout<<" 3body: \t"<<rate<<endl;

// 6)    H2    +     J_LW        ->         2 H  
//  J_LW in units of J21
    rate = k[6] * y[2];
    r_f[6][1] = 2. * rate;
    r_f[6][2] = -rate;
// 7)    H-    +     ph       ->   H    +     e
    rate = k[7] * y[5];
    r_f[7][1] = rate;
    r_f[7][3] = rate;
    r_f[7][5] = -rate;
// 8)    H    +     e-       ->   H+    +     2e-
    rate = k[8] * y[1] * y[3] * nH;
    r_f[8][1] = -rate;
    r_f[8][3] = rate;
    r_f[8][4] = rate;

//************************************************************
    for (int isp=0; isp<N_sp+1; isp++){
        r_f_tot[isp] = 0.;
        for (int ire=1; ire<N_react+1; ire++){
            r_f_tot[isp] += r_f[ire][isp];
        }
    }     
 */
