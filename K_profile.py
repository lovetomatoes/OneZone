from PYmodule import *

Rvir = 1.
Rcore = 0.1*Rvir
Kism = 1.
Kvir = 1.
# r in physical units
def K(r):
    return max(Kism, .1*Kvir) + Kvir*(r/Rvir)
x = np.logspace(-2,1.5,num=20)
y = K(x)

# plt.figure(figsize=(10,8),dpi=200)
plt.loglog(x,y)

# plt.xlabel('r',fontsize=30)
# plt.ylabel('K',fontsize=30);
# plt.xticks(fontsize=20);plt.yticks(fontsize=20)
# plt.legend(loc='best',fontsize=15)

# plt.savefig('z_J_1000.png')
plt.show()