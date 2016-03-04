#ifndef particleList_H
#define particleList_H

#include "particle.h"
#include "Chemical_potential.h"
#include<string>
#include<fstream>
#include<vector>
using namespace std;

class particleList
{
    private:
       string particleListFilename;
       vector<particle*> partList;
       double edSystem, sdSystem, pressureSys, net_baryon_density;

    public:
       particleList(string particleTableName);
       ~particleList();

       int getParticlelistSize() {return(partList.size());};

       void readParticlelistTable(string tableName);
       void calculate_particle_decay_probability(Chemical_potential* mu_tb);
       void calculate_particle_chemical_potential2(
                       double Temperature, Chemical_potential* mu_tb);
       void calculateSystemEOS(double mu_B=0.0, double mu_S=0.0);
       void calculate_particle_mu(double mu_B, double mu_S);
       void calculate_particle_yield(double Temperature, double mu_B=0.0, 
                                     double mu_S=0.0);
       void calculateSystemenergyDensity(double Temperature);
       void calculateSystemPressure(double Temperature);
       void calculateSystementropyDensity(double Temperature);
       void calculateSystemnetbaryonDensity(double Temperature);
       int get_particle_idx(int particle_monval);
       void calculate_particle_chemical_potential(
                       double Temperature, Chemical_potential* mu_tb);
       void calculateSystemEOS_and_output_in_2D();
       void output_particle_chemical_potentials(Chemical_potential* mu_tb);
};


#endif
