extern "C" double Mg_max(double Tg, double z, double Mh);
double R_EQ(double Tg, double rhoc, double rs);
double A_TR(double Tg, double rhoc, double R);
extern "C" double Mg(char* filename, double Tg, double R, double z, double Mh);
extern "C" double Mg2ng(double Mg, double ni, double Tg, double z, double Mh);