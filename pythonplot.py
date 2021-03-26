from PYmodule import *

def Tv(Mh,z):
    return 2.324e4* (Mh/1.e8)**(2./3.) *  (1+z)/11.

Ts = []
Ts_iso = [] # n_tell = 1.e4, T_tell = 4000K
Ts_H2 = []
Ts_isofail = []
Ts_isoOK = []

T_tell = 4000

J1000_min = 1.e7; Jcol_min = 1.e7
J1000_max = 1.e-7; Jcol_max = 1.e-7

for i_bsm in range(4):
    T=ascii.read('Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ')
    T['Tv_col'] = Tv(T['Mh_col'],T['z_col'])
    T['Tv_1000'] = Tv(T['Mh_1000'],T['z_1000'])
    Ts_iso.append(T[T['Tg_col']>T_tell])
    Ts_isofail.append(T[np.logical_and(T['Tg_col']>T_tell, T['iso_col']==0)])
    Ts_isoOK.append(T[np.logical_and(T['Tg_col']>T_tell, T['iso_col']==1)])
    Ts_H2.append(T[np.logical_or(T['iso_col']==0, T['Tg_col']<T_tell)])
    Ts.append(T)
    if np.max(T['J_1000'])>J1000_max:
        J1000_max = np.max(T['J_1000'])
    if np.min(T['J_1000'])<J1000_min:
        J1000_min = np.min(T['J_1000'])
    if np.max(T['J_col'])>Jcol_max:
        Jcol_max = np.max(T['J_col'])
    if np.min(T['J_col'])<Jcol_min:
        Jcol_min = np.min(T['J_col'])
    print("max Mh_col:",np.max(Ts_iso[i_bsm]['Mh_col']))
    print("min z_col:",np.min(Ts_iso[i_bsm]['z_col']))
    print("max z_col:",np.max(Ts_iso[i_bsm]['z_col']))

T = Ts_iso[1]
id = np.argmax(T['J_1000'])
print(len(T['tree'][T['J_1000']>1.e4]))
# print('tree id: ',T['tree'][id],'tree Tvir: ', T['Tv_1000'][id], 'J_1000: ',T['J_1000'][id])
# print('Tg_1000: ', T['Tg_1000'][id])
# print('z_1000: ', T['z_1000'][id],'Mh_1000: ', T['Mh_1000'][id])

# plt.figure(figsize=(10,8),dpi=200)
# for i in range(len(Ts)):
#     T_iso = Ts_iso[i]
#     T_H2 = Ts_H2[i]
#     plt.scatter(T_iso['Tv_1000'],T_iso['J_1000'],facecolors='none', edgecolors='C'+str(2*i),label='iso v_bsm='+str(i))
#     plt.scatter(T_H2['Tv_1000'],T_H2['J_1000'],marker='s',facecolors='none', edgecolors='C'+str(2*i+1),label='H2')
# plt.xscale('log')
# plt.yscale('log')
# plt.ylim(J1000_min/2, J1000_max*2)
# plt.xlabel('Tv_1000',fontsize=30)
# plt.ylabel('J_1000',fontsize=30);
# plt.xticks(fontsize=20);plt.yticks(fontsize=20)
# plt.legend(loc='best',fontsize=15)
# plt.xscale('log')
# plt.yscale('log')
# plt.savefig('Tv_J_1000.png')

plt.figure(figsize=(10,8),dpi=200)
for i in range(100):
    T = ascii.read('../treefiles/tree_'+str(i)+'mer')
    plt.plot(T['z'], T['Mh_Ms'],linewidth=1,c='black',alpha=.5)

for i in range(len(Ts)):
    T_iso = Ts_isoOK[i]
    T_H2 = Ts_H2[i]
    plt.scatter(T_iso['z_1000'],T_iso['Mh_1000'], s=20,facecolors='none', edgecolors='C'+str(2*i),label='iso v_bsm='+str(i))
    # plt.scatter(T_H2['z_1000'],T_H2['Mh_1000'],marker='s',facecolors='none', edgecolors='C'+str(2*i+1),label='H2')
    T_iso = Ts_isofail[i]
    T_H2 = Ts_H2[i]
    plt.scatter(T_iso['z_1000'],T_iso['Mh_1000'], marker='x', s=20,c='C'+str(2*i),label='fail v_bsm='+str(i))

plt.yscale('log')
plt.xlim(5,50)
plt.xlabel('z',fontsize=30)
plt.ylabel('Mh_1000',fontsize=30);
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.legend(loc='best',fontsize=15)
plt.yscale('log')
plt.savefig('z_Mh_1000_isoOK+tree.png')

# plt.figure(figsize=(10,8),dpi=200)
# for i in range(len(Ts)):
#     T_iso = Ts_iso[i]
#     T_H2 = Ts_H2[i]
#     plt.scatter(T_iso['z_1000'],T_iso['Tv_1000'],facecolors='none', edgecolors='C'+str(2*i),label='iso v_bsm='+str(i))
#     plt.scatter(T_H2['z_1000'],T_H2['Tv_1000'],marker='s',facecolors='none', edgecolors='C'+str(2*i+1),label='H2')
# plt.yscale('log')
# plt.xlim(5,50)
# plt.xlabel('z',fontsize=30)
# plt.ylabel('Tv_1000',fontsize=30);
# plt.xticks(fontsize=20);plt.yticks(fontsize=20)
# plt.legend(loc='best',fontsize=15)
# plt.yscale('log')
# plt.savefig('z_Tv_1000.png')

# plt.figure(figsize=(10,8),dpi=200)
# for i in range(len(Ts)):
#     T_iso = Ts_iso[i]
#     T_H2 = Ts_H2[i]
#     plt.scatter(T_iso['z_1000'],T_iso['J_1000'],facecolors='none', edgecolors='C'+str(2*i),label='iso v_bsm='+str(i))
#     plt.scatter(T_H2['z_1000'],T_H2['J_1000'],marker='s',facecolors='none', edgecolors='C'+str(2*i+1),label='H2')
# plt.yscale('log')
# plt.xlim(5,50); plt.ylim(J1000_min/2, J1000_max*2)
# plt.xlabel('z',fontsize=30)
# plt.ylabel('J_1000',fontsize=30);
# plt.xticks(fontsize=20);plt.yticks(fontsize=20)
# plt.legend(loc='best',fontsize=15)
# plt.yscale('log')
# plt.savefig('z_J_1000.png')