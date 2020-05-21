void DyDx_iso(double x, double* y, double* dydx, int argc, double* argv);
void DyDx_adb_Kc(double x, double* y, double* dydx, int argc, double* argv);
void DyDx_adb_fit(double x, double* y, double* dydx, int argc, double* argv);
void rk4(double* yout, double x, double h, double *y, int const n, double *dydx0,
            int argc, double* argv,
            void (*derivs)(double, double*, double*, int, double*) );

void check_conv(double err_crit, double* yout, double x, double& h, double *y, int const n, double *dydx0,
            int argc, double* argv,
            void (*derivs)(double, double*, double*, int argc, double* argv) );