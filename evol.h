double getT(int MerMod, double J, double Tb, char* treename, double nH_tell, bool Ma_on);
extern "C" void evol(char* treename, char* fout, int mode, double Jlw, bool Ma_on);
extern "C" void evol_Jc(char* treename, char* fout, double Tb, bool Ma_on);