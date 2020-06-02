double getT(double& zcol, int MerMod, double J, double Tb, std::string treename, bool spec, bool Ma_on, int i_bsm, double nH_tell);
extern "C" void evol(std::string treename, std::string fout, int MerMod, double Tbb, double J21, bool spec, bool Ma_on, int i_bsm);
extern "C" void evol_Jc(std::string treename, std::string fout, double Tb, int MerMod, bool spec, bool Ma_on, int i_bsm);