extern "C" double Mg_max(double Tg, double z, double Mh);
double R_EQ(double Tg, double rhoc, double rs);
double A_TR(double Tg, double rhoc, double R);
extern "C" void profile(std::string filename, double Tg, double R, double z, double Mh);
void BOUNDARY(double& N_VIR, double& MG_VIR, double Tg, double R, double z, double Mh);
void Nvir2N0(double& n_sol, double& nvir_max, double ni, double Tg, double z, double Mh);