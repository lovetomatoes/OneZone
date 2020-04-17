extern "C" double Mg_max(double Tg, double z, double Mh);
double R_EQ(double Tg, double rhoc, double rs);
double A_TR(double Tg, double rhoc, double R);
extern "C" double Mg(char* filename, double Tg, double R, double z, double Mh);
extern "C" double N_VIR(double Tg, double R, double z, double Mh);
void N_VIR_max(double& Rm, double& n_max, double Tg, double z, double Mh);
void Nvir2N0(double& n_sol, double& nvir_max, double ni, double Tg, double z, double Mh);