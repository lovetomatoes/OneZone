double getT(int MerMod, double J, double Tb, char* treename, double nH_tell, bool spec, bool Ma_on);
extern "C" void evol(char* treename, char* fout, int MerMod, double Tbb, double J21, bool spec, bool Ma_on);
extern "C" void evol_Jc(char* treename, char* fout, double Tb, int MerMod, bool spec, bool Ma_on);