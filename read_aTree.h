struct MainProgenitor
{
    int id_tree;
    int j;
    double mhalo;
    double z;
    double t;
    double dm;
    double dt;
    double c;
    double Tvir;
    double mratio;
    bool major;
    double ng_adb;
};

void aTree(int& nMer, std::string filename, MainProgenitor* MPs);