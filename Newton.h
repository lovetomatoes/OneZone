void FX(double* F, double* x, double* x0, double dt, double* k, double nH);
void J(double *j, double *x, double dt, double *k, double nH);
void J_INV(double *j_inv, double *x, double dt, double *k, double nH);
void SOL_IMPLICIT(double* dy, double* x0, double* y1, double dt, double nH, double T_K, double J_LW, double Tb);