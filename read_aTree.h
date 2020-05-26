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
    double ng_adb;
};

void read_aTree (int& num, char* fname, MainProgenitor* MPs);
void aTree(int& nMer, char* filename, MainProgenitor* MPs);