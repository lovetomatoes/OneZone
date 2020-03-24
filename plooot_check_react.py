from astropy.io import ascii
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.table import Table
global G, m_H
G = 6.67408e-8; m_H = 1.66053904e-24
nHI = 1.; T_K = 100
#****************************************************************
#*     SPECIES                                                  *
#*     1 : H      2 : H2     3 : e      4 : H+      5 : H-      *
#****************************************************************
y0 = np.array([[1.-2*1.e-6-1.e-4-1.e-10], [1.e-6], [1.e-4], [1.e-4+1.e-10], [1.e-10]])
t_ff = np.sqrt(3/(4*np.pi*G*m_H*nHI));

Tc = ascii.read("./asol.txt"); Tp = ascii.read("./pyasol.txt")
Tc = Table(Tc,names=["t",'HI','H2','e','H+','H-']);
Tp = Table(Tp,names=["t",'HI','H2','e','H+','H-'])

""" plt.plot(Tc['t']/t_ff,Tc['HI']/Tc['HI'][0],label = "cHI");# plt.plot(Tp['t'],1.0001*Tp['H2'],label='p')
plt.plot(Tc['t']/t_ff,Tc['H2']/Tc['H2'][0],label = "cH2");
plt.plot(Tc['t']/t_ff,Tc['e']/Tc['e'][0],label = "ce");
plt.plot(Tc['t']/t_ff,Tc['H+']/Tc['H+'][0],label = "cH+");
plt.plot(Tc['t']/t_ff,Tc['H-']/Tc['H-'][0],label = "cH-"); 

plt.plot(Tp['t']/t_ff,Tp['HI']/Tp['HI'][0],label = "pHI");# plt.plot(Tp['t'],1.0001*Tp['H2'],label='p')
plt.plot(Tp['t']/t_ff,Tp['H2']/Tp['H2'][0],label = "pH2");
plt.plot(Tp['t']/t_ff,Tp['e']/Tp['e'][0],label = "pe");
plt.plot(Tp['t']/t_ff,Tp['H+']/Tp['H+'][0],label = "pH+");
plt.plot(Tp['t']/t_ff,Tp['H-']/Tp['H-'][0],label = "pH-"); """

plt.plot(Tc['t']/t_ff,Tc['HI'],label = "cHI");# plt.plot(Tp['t'],1.0001*Tp['H2'],label='p')
plt.plot(Tc['t']/t_ff,Tc['H2'],label = "cH2");
plt.plot(Tc['t']/t_ff,Tc['e'],label = "ce");
plt.plot(Tc['t']/t_ff,Tc['H+'],label = "cH+");
plt.plot(Tc['t']/t_ff,Tc['H-'],label = "cH-"); 

plt.plot(Tp['t']/t_ff,Tp['HI'],label = "pHI");# plt.plot(Tp['t'],1.0001*Tp['H2'],label='p')
plt.plot(Tp['t']/t_ff,Tp['H2'],label = "pH2");
plt.plot(Tp['t']/t_ff,Tp['e'],label = "pe");
plt.plot(Tp['t']/t_ff,Tp['H+'],label = "pH+");
plt.plot(Tp['t']/t_ff,Tp['H-'],label = "pH-");
plt.xscale("log"); plt.yscale("log")
plt.legend(loc="best")
plt.show()
