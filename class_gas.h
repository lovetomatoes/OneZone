#include "read_aTree.h"

using namespace std; //std::ofstream

class GAS
{
    public:
        //******************************************
        void add_Nt(int);
        void get_para();
        void a_react_sol(bool);
        void react_sol(bool);
        void ana_e();
        void ana_H2();

        void setMerger();
        void timescales();
        void freefall();
        void T_sol();
        static void read_Jz(string, double*, double*, int&);

        // member function definition -- constructor
        GAS(double *frac0, int MergerModel, double J21, double Tbb, string treefile, string Jzfile, bool spec, bool Ma_turn, int bsm);
        ~GAS();
    //private:
        int N, Nt;
        int i_LE; int MerMod;
        int nMer, iMP;
        double z0, z;
        double t0, t1, Dt, t_act, t_ff0, t_delay;
        double dMdt, dEdt;

        double T_K0, nH0, rho0, P0, e0, S0;
        double rhoc_DM;
        double J_LW, Tb;
        double kappa_Hm, kappa_H2p;
        double *y0, *y1, *ys,  *k, *rf;
        double yequi, ypd, ycd, ycool, ycool_crit;
        double Jc_pd, Jc_cd, Jc_pred, Jc_pred_max;
        double delta_H2_compr, delta_H2_compr_min;
        double fMa_H2crit, gMa_H2crit, n_H2crit, z_H2crit;
        double a, b, c, d, e;
        double r_h, r_c, r_cH, r_cH2, r_cMetal;
        double Z;
        double t_ff, t_h, t_c, t_rcb, t_chem, t_ion;
        int evol_stage, i_bsm;
        bool intoequi;
        double n_iso, Mh_prev, t_prev, dt_iso;
        // evol_stage: label density evolve: MerMod==0 直接freefall; MerMod==1, evol_stage=1 等熵,2 Eli saturated core,3 等温,4 freefall
        double Gamma_mer, Gamma_mer_th, Gamma_mer_k;
        double cs, MJ, RJ, MJ_eff, Mgas, M_BE, Mg_intg;
        double M_major, M_minor, MJ0, Ma, f_Ma, g_Ma, v_tur2, reduction, ncore, nb200;
        double de_tot, dM_tot;
        double v_tur2_eff, cs2_eff, v_bsm; // v_tur2_eff = Ptur/rho; cs2_eff = Ptot/rho
        bool inMer, inDelay, not_adb, Ma_on;
        double Mh;
        MainProgenitor* MPs; //structure defined in read_aTree.h
        //HALO halo (1,1); must have initial parameters...
        ofstream file_ingas;
        int n_za;
        double z_col, ng_col, Tg_col, Mh_col, J_col, f_col, MJ_col;
        double z_1e4, Tg_1e4, Mh_1e4, J_1e4, f_1e4, MJ_1e4;
        double Tg_max;
        double z_loi, ng_loi, Tg_loi, Mh_loi, J_loi, f_loi, MJ_loi;

    private :
        int static const n_ra = 40;
        double* Ta, *ka;
        double* Ja, *za;

};

const double fraction = 0.75; // fraction of merger heating turned into thermal energy