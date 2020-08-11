double getT(int argc, double* argv, int MerMod, double J, double Tb, std::string treename, bool spec, bool Ma_on, int i_bsm, double nH_tell);
extern "C" void evol(std::string treename, std::string fout, int MerMod, double Tbb, double J21, bool spec, bool Ma_on, int i_bsm);
void evol_Jc(int* vi, double* vd, int itr, std::string treename, double Tb, int MerMod, bool spec, bool Ma_on, int i_bsm);