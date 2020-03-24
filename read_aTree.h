struct MainProgenitor
{
    int j;
    double mhalo;
    double dm;
    double z;
    double t;
    double dt;
    double Tvir;
    int id_tree;
    double mratio;
    bool major;
};

void read_aTree (int& num, char* fname);
void aTree(int& nMer, MainProgenitor* MPs, char* filename);