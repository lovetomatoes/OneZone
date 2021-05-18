double Lambda_H2(double nH, double T_K, double* y);
double Lambda_H(double nH, double T_K, double y_H,double y_e, double k_ion);
double Lambda_Hep(double nH, double T_K, double y_Hep, double y_e, double y_He, double k_Heion);
double Lambda_metal(double nH, double T_K, double y_H, double y_H2, double y_e, double Z, double f_C);
double Gamma_chem(double nH, double T_K, double* y, double* k);
double Gamma_compr(double cs2_eff, double t_ff);
//double Gamma_dyn(double z, double Mh_in_g);
//double Gamma_merger(double z, double addMh);