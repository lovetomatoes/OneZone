extern "C" void profile(std::string filename, double Tg, double R, double z, double Mh);
void BOUNDARY(double& N_VIR, double& MG_VIR, double Tg, double R, double z, double Mh);
void Nvir2N0(double& n_sol, double& nvir_max, double ni, double Tg, double z, double Mh);