#ifndef particleList_H
#define particleList_H

#include "parameters.h"
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

    public:
       particleList(string particleTableName);
       ~particleList();

       void readParticlelistTable();
       void calculate_particle_mu(double mu_B, double mu_S);
       void calculate_particle_yield(double Temperature);
       int get_particle_idx(int particle_monval);
};


#endif
