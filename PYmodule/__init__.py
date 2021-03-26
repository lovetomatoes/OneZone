import numpy as np
from ctypes import * # c 类型库
import struct
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import os

#os.system("g++ evol.cpp -L/usr/local/lib class_gas.o LE_iso.o read_aTree.o class_halo.o dyn.o thermo.o reaction.o Newton5.o my_linalg.o gsl_inverse.o RK4.o -lgsl -lgslcblas -lm -o cc.so -shared -fPIC")
#libc = CDLL('cc.so') # 装入动态链接库 ## 居然必须放在这里

global G, h0, H0, Omega_m0, Omega_L0, m_H, mu, Ms, pi, km, pc, Myr, alpha_T
G, c, k_B, m_H = 6.67408e-8, 2.9979245e10, 1.38064852e-16, 1.66053904e-24
pi = 3.141593
mu = 1.2
Ms = 2.e33
pc = 3.e18
Mpc = 1.e6*pc
km = 1.e5
Myr = 1.e6*(365*24*3600)
Omega_m0 = 0.311 
Omega_L0 = 1 - Omega_m0
h0 = .677
H0 = h0*100*km/Mpc

alpha_T = 2.324e4

def Omega_mz(z):
    return Omega_m0*(1+z)**3 /(Omega_m0*(1+z)**3 + Omega_L0)

def Hz(z):
    return H0*np.sqrt( Omega_m0*(1+z)**3 + Omega_L0 ) 
 
def RHO_crit(z):
    return 3*pow(H0,2)/(8*pi*G)*(1+z)**3*Omega_m0/Omega_mz(z) 

class HALO:
    def __init__(self,M,z0):
        self.Mh = M
        self.z = z0
        self.c = 18*pow(self.Mh/(1.e11*Ms), -0.13)/(1+self.z) #concentration parameter c from Dekel & Birnboim 2006 Eq(22)
        c, z = self.c, self.z
        self.d = Omega_mz(z) - 1 
        d = self.d
        self.Delta_crit = 18.0*pi*pi + 82*d - 39*d*d  # Delta_crit ~ 200, overdensity
        Delta_crit = self.Delta_crit

        self.delta0 = self.Delta_crit/3.*pow(c,3)/(-c/(1+c) + np.log(1+c)) # characteristic overdensity parameter 
        delta0 = self.delta0
        self.rho_crit = RHO_crit(z)  # mean density of DM at z
        self.rho_c = self.rho_crit * delta0 

        self.Rvir = pow( self.Mh/(4./3*pi*Delta_crit*self.rho_crit),1./3. ) 
        self.Rs = self.Rvir/self.c 
        self.Vc = np.sqrt(G*self.Mh/self.Rvir) 

        self.t_dyn = self.Rvir/self.Vc 
        self.Tvir = G*self.Mh*(mu*m_H)/(2.*k_B*self.Rvir)
        self.gc = 2*c/(np.log(1+c) - c/(1+c)) 
        self.alpha = self.Tvir/self.Mh**(2./3)

    def Rho_r(self, r):
        rho_crit, delta0, Rvir = self.rho_crit, self.delta0, self.Rvir
        c, x = self.c, r/Rvir
        return rho_crit*delta0/( c*x * (1+c*x)**2 )

    # x = r/Rvir  c = Rvir/Rs
    def F_NFW(self,x):
        c = self.c
        return -c*x/(1+c*x) + np.log(1+c*x)
        
    def M_enc(self,r):
        rho_crit, delta0, Rs, Rvir = self.rho_crit, self.delta0, self.Rs, self.Rvir
        M_r = 4*pi*rho_crit*delta0*pow(Rs,3)*self.F_NFW(r/Rvir)
        return M_r


    def Phi(self, r):
        # lim r -> 0
        #return -4*pi*G*rho_crit*delta0*Rs*Rs 
        rho_crit, delta0, Rs = self.rho_crit,  self.delta0, self.Rs
        return -4*pi*G*rho_crit*delta0*(Rs**3)/r*np.log(1+r/Rs) 




class iso_gas:
    def __init__(self,z,T):
        self.z = z; self.halo_T = T
        self.Mh = 1.e8*Ms*pow(T/alpha_T*11/(1+z),1.5)
        halo = HALO(self.Mh,self.z)
        self.rs = halo.Rs; self.rvir = halo.Rvir; self.rhoc = halo.rho_c
        print("z, Mh, rs, rvir, rhoc")
        print(z,self.Mh/Ms, self.rs, self.rvir, self.rhoc)

        self.Tg = 1.e4
        beta = (4*np.pi*G*mu*m_H *self.rhoc)/(k_B*self.Tg)
        self.R_EQ = 9/4*beta*self.rs**2

    def a(self,R):
        return (k_B*self.Tg/ (4*np.pi*G*mu*m_H*R*self.rhoc) )**.5
    def Req(self,red=.1):#0.1 nice...
        return self.R_EQ
    def Mg(self,R,power_a):
        self.rho_g0 = self.rhoc*R
        Rmax = self.rvir
    #------------------------------------------------------------------
    # Rcore 由 R v.s. Req 决定; 之后都-2 profile
        if R>self.Req(): # gas dominate, Rcore = a, Rout = rvir
            Rcore = self.a(R)
        else: # DM dominate, only a core within r1.
            r1 = min(self.rvir, self.a(self.Req()))#外围最大积分到rvir
            Rcore = r1
        if power_a==3:
            return ( 4*np.pi/3 *self.rho_g0 *Rcore**3 + 4*np.pi* self.rho_g0* Rcore**3*np.log(Rmax/Rcore) )/Ms
        else:
            return ( 4*np.pi/3 *self.rho_g0 *Rcore**3 + 4*np.pi* self.rho_g0 /(3-power_a)* Rcore**power_a* (Rmax**(3-power_a)-Rcore**(3-power_a)) )/Ms
    #--------------------------------------------------------------------
    # 不连续 Rcore 由 R v.s. Req 决定;
    # 如果是 DM dominate, 只用一个core的质量
        # if R>self.Req(): # gas dominate, Rcore = a, Rout = rvir
        #     Rcore = self.a(R)
        #     if power_a==3:
        #         return ( 4*np.pi/3 *self.rho_g0 *Rcore**3 + 4*np.pi* self.rho_g0* Rcore**3*np.log(Rmax/Rcore) )/Ms
        #     else:
        #         return ( 4*np.pi/3 *self.rho_g0 *Rcore**3 + 4*np.pi* self.rho_g0 /(3-power_a)* Rcore**power_a* (Rmax**(3-power_a)-Rcore**(3-power_a)) )/Ms
    
        # else: # DM dominate, only a core within r1.
        #     r1 = min(self.rvir, self.a(self.Req()))#外围最大积分到rvir
        #     Rcore = r1
        #     return ( 4*np.pi/3 *self.rho_g0 *Rcore**3 )/Ms
