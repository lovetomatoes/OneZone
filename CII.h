double CIIcool(double nH, double T_K, double Nc_CII, double y_m, double y_a, double y_e,
               double esc_10, double tau_c);
void CIIpop(double &f_0, double& f_1, double &esc_10);
void twolevel(double &f_0, double &f_1, double *esc,double *error);
double Q_bg(double T_nu);
double beta_esc(double tau_L, double tau_C);