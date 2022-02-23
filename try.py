from PYmodule import *

halo_d0 = 0
halo_rvir = 0
halo_rs = 0
halo_c = 0
halo_Dc = 0
Tgas = 1e4


T=ascii.read('profile_Nsol.txt', guess=False,delimiter=' ')

print(T.info)
r_pc = T['r_pc'][:-1]
dPdr_over_rho = -1./(mu*m_H*T['ng'][1:])*k_B*Tgas*(T['ng'][1:]-T['ng'][:-1])/(T['r_pc'][1:]-T['r_pc'][:-1])/pc
g_r = T['g_r'][:-1]
for i in range(len(T)-1):
    print(g_r[i]-dPdr_over_rho[i])


plt.figure(figsize=(10,8),dpi=200)

plt.plot(r_pc,((g_r-dPdr_over_rho)/g_r)**2)
plt.scatter(r_pc,((g_r-dPdr_over_rho)/g_r)**2)
plt.xlabel('pc',fontsize=30)
plt.ylabel(r'$cm s^{-1}$',fontsize=30);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.legend(loc='best',fontsize=15)
plt.xscale('log'); plt.yscale('log')
plt.savefig('force.png')
