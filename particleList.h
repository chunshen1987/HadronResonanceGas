#ifndef particleList_H
#define particleList_H

#include "particle.h"
#include<string>
#include<fstream>
#include<vector>
using namespace std;

class particleList
{
    private:
       string particleListFilename;
       vector<particle*> partList;
       double edSystem, sdSystem, pressureSys;

    public:
       particleList(string particleTableName);
       ~particleList();

       void readParticlelistTable(string tableName);
       void calculateSystemEOS(double mu_B=0.0, double mu_S=0.0);
       void calculate_particle_mu(double mu_B, double mu_S);
       void calculate_particle_yield(double Temperature, double mu_B=0.0, double mu_S=0.0);
       void calculateSystemenergyDensity(double Temperature, double mu_B=0.0, double mu_S=0.0);
       void calculateSystemPressure(double Temperature, double mu_B=0.0, double mu_S=0.0);
       void calculateSystementropyDensity(double Temperature, double mu_B=0.0, double mu_S=0.0);
       int get_particle_idx(int particle_monval);
       int getParticlelistSize() {return(partList.size());};
};


#endif
