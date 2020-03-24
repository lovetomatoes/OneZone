# pwd --> ~/Desktop/Kohei-project/C_pp_chemical
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table

global G, c, k_B, m_H, gamma

G, c, k_B, m_H = 6.67408e-8, 2.9979245e10, 1.38064852e-16, 1.66053904e-24
gamma = 5./3
Ms = 2.e33

#T1 = ascii.read('data/evolve_Mer_delaydilu.txt', guess=False,delimiter=' ')
#T1 = ascii.read('data/evolve216.txt', guess=False,delimiter=' ')
#T2 = ascii.read('data/evolve217.txt', guess=False,delimiter=' ')
#T10 = ascii.read('data/evolve216nomer.txt', guess=False,delimiter=' ')
#T20 = ascii.read('data/evolve217nomer.txt', guess=False,delimiter=' ')
T1 = ascii.read('data/evolve50_mer.txt', guess=False,delimiter=' ')
T2 = ascii.read('data/evolve50_nomer.txt', guess=False,delimiter=' ')
T1 = ascii.read('data/evolve.txt', guess=False,delimiter=' ')

print(T1.info)

t = T1['t']
nH = T1['nH']; T_K = T1['T']
Tvir = T1['Tvir']
#inDelay = T1['inDelay']
inMer = T1['inMer']
yH2 = T1['yH2']
ye = T1['ye']
Mh = T1['Mh']
z = T1['z']
gamma_adb = 5/3

T1['Taccu'] = np.zeros(len(T1))
T1['Taccu'][0] = T1['T'][0]
imer = 0
for i in range(1,len(T1)):
    if T1['iMer'][i] > imer:
        T1['Taccu'][i] = T1['Taccu'][i-1] + T1['dMh'][i-1]/T1['Mh'][i-1]*(gamma_adb-1)*T1['Tvir'][i-1]
        T1['Taccu'][i] = T1['Taccu'][i-1] + T1['dMh'][i-1]/T1['Mh'][i]*(gamma_adb-1)*T1['Tvir'][i-1]
        #T1['Taccu'][i] = T1['Taccu'][i-1]* (T1['Mh'][i]/T1['Mh'][i-1] )**(2/3)
        imer += 1
    else:
        T1['Taccu'][i] = T1['Taccu'][i-1]
plt.figure(figsize=(12,12),dpi=200)
#plt.plot(T1['z'], T1['T'], label=r'$T_{g}$');
plt.plot(T1['z'], T1['T'], label=r'$T_{g}$');
plt.plot(T1['z'], T1['Tvir'], label=r'$T_{vir}$');
#plt.plot(T1['z'], T1['Tvir_Mh'], label=r'$T_{vir,Mh}$');

#plt.plot(T1['z'], T1['Taccu'], linewidth=2, label=r'$T_{accu}$');
plt.xlim(35,10.25)
#plt.xlim(20,15);plt.ylim(2000,5000)
#plt.xlim(35,32); plt.ylim(250,500)
plt.legend(loc='upper left',fontsize=30)
plt.xlabel('z',fontsize=30)
plt.ylabel('T (K)',fontsize=30);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.savefig('figs/z_T_Tvir.png')
plt.show()

""" 
# T1: evolve.txt, 
T1['Taccu'] = np.zeros(len(T1))
T1['Taccu'][0] = T1['T'][0]
imer = 0
for i in range(1,len(T1)):
    if T1['iMer'][i] > imer:
        T1['Taccu'][i] = T1['Taccu'][i-1] + T1['dMh'][i-1]/T1['Mh'][i-1]*(gamma_adb-1)*T1['Tvir'][i-1]
        imer += 1
    else:
        T1['Taccu'][i] = T1['Taccu'][i-1]
plt.figure(figsize=(10,8),dpi=200)
plt.plot(T1['z'], T1['T'], linewidth=2, label=r'$T_{g}$');
plt.plot(T1['z'], T1['Tvir'], linewidth=2, label=r'$T_{vir}$');
plt.plot(T1['z'], T1['Taccu'], linewidth=2, label=r'$T_{accu}$');

plt.xlim(35,30)
plt.legend(loc='lower right',fontsize=30)
plt.xlabel('z',fontsize=30)
plt.ylabel('T (K)',fontsize=30);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.savefig('figs/z_T_Tvir.png')
plt.show()
 """

""" 
plt.figure(figsize=(10,8),dpi=200)
plt.plot(T1['nH'], T1['T'], linewidth=2, label=r'$T_{g}$');
plt.plot(T1['nH'], T1['Tvir'], linewidth=2, label=r'$T_{vir}$');
#plt.ylim(bottom=100)
plt.legend(loc='lower right',fontsize=30)
plt.xlabel(r'$nH (cm^{-3})$',fontsize=30)
plt.ylabel('T (K)',fontsize=30);
plt.xscale('log')
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.savefig('figs/T_Tvir.png')
plt.show()

T3 = np.loadtxt('data/Jaaaa.txt',delimiter=' ')
print(T3[:,0],T3[:,1],T3[:,2])
plt.figure(figsize=(8,6),dpi=200)
plt.loglog(T3[:,0],T3[:,1], linewidth=2, label=r'$T_{no\,merger}$',c='orange');
plt.loglog(T3[:,0],T3[:,2], linewidth=2, label=r'$T_{mer1}$',c='purple');
plt.loglog(T3[:,0],T3[:,3], linewidth=2, label=r'$T_{mer2}$',c='g');
plt.loglog(T3[:,0],T3[:,4], linewidth=2, label=r'$T_{mer3}$',c='cyan');
#plt.scatter(T3[:,0],T3[:,2], s=5, label=r'$T_{merger}$');
plt.legend(loc='best',fontsize=20)
plt.xlabel(r'$T_{rad} (K)$',fontsize=20);
plt.ylabel(r'$J_{21,\, crit}$',fontsize=20)
plt.xscale('log');plt.yscale('log')
plt.xticks(fontsize=15);plt.yticks(fontsize=15)
plt.savefig('figs/Jc.png')
plt.show()
 """

""" plt.figure(figsize=(8,6),dpi=200)
plt.loglog(T1['nH'], T1['T'], linewidth=2, label=r'$T_{merger}$');
plt.loglog(T2['nH'], T2['T'], linewidth=2, label=r'$T_{no \, merger}$');
plt.text(1.e1,2.e2,r'$J_{LW}=50, T_{rad}=10^4K$',fontsize=15)
plt.ylim(bottom=100)
plt.legend(loc='lower right',fontsize=20)
plt.xlabel(r'$nH (cm^{-3})$',fontsize=20)
plt.ylabel('T (K)',fontsize=20);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.savefig('figs/n-T_2mer.png')
plt.show()

plt.figure(figsize=(8,6),dpi=200)
plt.loglog(T1['t'], T1['T'], linewidth=2, label=r'$T_{merger}$');
plt.loglog(T2['t'], T2['T'], linewidth=2, label=r'$T_{no \, merger}$');
#plt.loglog(T10['t'], T10['T'], label=r'$T_{smooth,nomer}$');
#plt.loglog(T20['t'], T20['T'], label=r'$T_{violent,nomer}$');
plt.scatter(t,1000*inDelay, s=.5, c='g', label = 'in Delay')
#plt.scatter(T2['t'],800*T2['inDelay'], s=.05, label = 'in Delay')
plt.ylim(bottom=100)
plt.text(1.e-3,2.e2,r'$J_{LW}=50, T_{rad}=10^4K$',fontsize=20)
plt.legend(loc='best',fontsize=20)
plt.xlabel(r'$t \, t_{ff,0})$',fontsize=20)
plt.ylabel('T (K)',fontsize=20);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.savefig('figs/t-T_2mer.png')
plt.show()

plt.figure(figsize=(8,6),dpi=200)
plt.plot(t,T_K, label=r'$T_{gas}$');
plt.plot(t,Tvir, label=r'$T_{vir}$');
plt.scatter(t,1000*inDelay, s=.05, c='g', label = 'in Delay')

plt.plot(T2['t'],T2['T'], label=r'$T_{gas}$');
plt.plot(T2['t'],T2['Tvir'], label=r'$T_{vir}$');
plt.scatter(T2['t'],800*T2['inDelay'], s=.05, label = 'in Delay')

plt.xlabel(r'$t \, t_{ff,0})$')
plt.ylabel('T (K)');
plt.xscale('log');plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/n-T.png')
plt.show()

plt.figure(figsize=(8,6),dpi=200)

Ms = 2.e33;
plt.plot(t,T2['Mh']/Ms, linewidth=.5, label=r'$M_h$')
plt.plot(t,T2['Mh_major']/Ms, linewidth=.5, label=r'$M_{h,major}$')
plt.scatter(t,1.e5*T2['inDelay'], s=.5, c='g', label=r'in Delay');
plt.xlabel(r'$t \, (t_{ff,0})$')
plt.ylabel(r'$M_h (M_\odot)$');
plt.xscale('log');plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/Mh.png')
plt.show()

plt.figure(figsize=(8,6),dpi=200)
plt.scatter(t,ye, s=.5, c='g', label = r'$y_e$')
plt.xlabel(r'$t \, (t_{ff,0})$')
plt.ylabel(r'$y_e$');
plt.xscale('log');plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/t-ye.png')
plt.show()

plt.figure(figsize=(8,6),dpi=200)
plt.scatter(nH,ye, s=.5, c='g', label = r'$y_e$')
plt.xlabel(r'$nH \, (cm^{-3})$')
plt.ylabel(r'$y_e$');
plt.xscale('log');plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/t-ye.png')
plt.show()
 """



""" 
plt.scatter(nH,1000*T1['inDelay'], s=.5, c='g', label = 'in Delay')
plt.xlabel(r'$nH \, (t_{ff,0})$')
plt.ylabel('T (K)');
plt.xscale('log');plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/t-T.png')
plt.show()


plt.figure(figsize=(8,6),dpi=200)
#plt.plot(T1['t'], T1['T'], label=r'$T_{gas}$');
#plt.plot(T1['t'],T1['Tvir'], label=r'$T_{vir}$');
plt.plot(z, nH, label=r'$T_{gas}$');

#plt.scatter(T1['t'],1000*T1['inDelay'], s=.5, c='g', label = 'in Delay')
#plt.xlabel(r'$t \, (t_{ff,0})$')
#plt.ylabel('T (K)');
#plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
#plt.savefig('figs/t-T.png')
plt.show()

fig, ax1 = plt.subplots(figsize=(8,6),dpi=200)
ax2 = ax1.twiny()
ax1.plot(t,T_K,'g-')
ax2.plot(nH,T_K,'b-')
 
ax1.set_xlabel("time",color='g')
ax1.set_ylabel("ye")
 
ax2.set_xlabel("nH",color='b')
plt.xscale('log');plt.yscale('log')
plt.show()
 """



'''plt.figure(dpi=200)
plt.plot(nHIs,t_heatings, label = 'Heating')
plt.plot(nHIs,t_coolings, label = 'Cooling')
plt.plot(nHIs,t_ffs, label = 'ff')

plt.xlabel('nHI (cm^-3)'); plt.ylabel('timescale(s)')

plt.xscale('log'); plt.yscale('log')
plt.legend(loc='best')
plt.savefig('../figs/totN.png')
plt.show()'''