
void read_sigma(double* nua, double** sigmaa);
void linear(double* xa, double* ya, int m, double x, double& y);
void bilinear(double* x1a, double* x2a, double** ya, int m, int n, double x1, double x2, double& y);
void Hm_CrossSec(double& sigma, double E);
void H2p_bf_CrossSec(double& sigma, double nu, double T_K, double* Ta, double* nua, double** sigmaa);
void kpd_Hm_H2p(double T_rad, double& k_Hm, double& k_H2p);