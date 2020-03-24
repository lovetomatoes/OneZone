from PYmodule import *
# 先 make evol 再运行evol.py

# tree files e.g.: "../code_tree/fort.217"
finprefix = "../code_tree/fort."
trs = [211,212,216,217]
Nfile = 4

Mermode = 1
Tb = 1.e4


Ma_on = bool(1)
J21 = 60

for i in range(1):
    
    treename = finprefix+ str(trs[i]);print(treename)
    foutprefix = "../data/evolveTb4J"
    foutprefix = foutprefix + str(J21) + "_"
    foutprefix = foutprefix+"fMaon" if Ma_on else foutprefix+"fMaoff"
    foutprefix = '../data/temp'
    libc.evol(c_char_p(bytes(treename,encoding='utf-8')), c_char_p(bytes(foutprefix+str(trs[i])+".txt",encoding='utf-8')), c_int(Mermode), c_double(J21), Ma_on)
    print("J21= ",J21, foutprefix+'.txt',treename)

# Ma_on = bool(1)
# J21 = 60
# for i in range(Nfile):
#     if i==1 :#or i==1:# or  i==3:
#         treename = finprefix+ str(trs[i]);print(treename)
#         foutprefix = "../data/evolveTb4J"
#         foutprefix = foutprefix + str(J21) + "_"
#         foutprefix = foutprefix+"fMaon" if Ma_on else foutprefix+"fMaoff"
#         #libc.evol(c_char_p(bytes(treename,encoding='utf-8')), c_char_p(bytes(foutprefix+str(trs[i])+".txt",encoding='utf-8')), c_int(Mermode), c_double(J21), Ma_on)
#         fout = "../data/Jcs.txt"
#         libc.evol_Jc(c_char_p(bytes(treename,encoding='utf-8')), c_char_p(bytes(fout,encoding='utf-8')), c_double(Tb),Ma_on)
#         print("tree= ",treename)


# Tbs[] = {8.e3, 1.e4, 2.e4, 3.e4, 5.e4, 1.e5, 2.e5};
# Jcrits[] = {0.8, 20, 800, 1000, 1100, 1100, 1100}; #Sugimura 2014
# J0s[] = {3.12805, 52.002, 2054.69, 3054.69, 3390.62, 3484.38, 3515.62}; #right...
