//const double fn=4;

class HALO
{
    public:
        //******************************************
        // member function definition
        HALO(double M0, double z0);
        ~HALO();
        void update();
        double F_NFW(double x);
        double Rho_r(double r);
        double M_enc(double r);
        void Update();
        double Phi(double r);
        double MerHeating(double dMhdt, double rin);
        double V_tur(double n_gascore);
        double V_inf2(double rout, double rin);
    //private:
        int i,j,k,id;
        double z, t;
        double Mh;
        double c, delta0;
        double Delta_crit;
        double rho_crit, rho_c, rho_vir, Tvir, Rvir, Vc, Rs;
        double n200;
        double gc;
        double t_dyn;
        double alpha;

};

double Mh_Tz(double Tvir, double z);
