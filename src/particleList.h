// Copyright 2016 Chun Shen
#ifndef SRC_PARTICLELIST_H_
#define SRC_PARTICLELIST_H_

#include <string>
#include <fstream>
#include <vector>

#include "./particle.h"
#include "./Chemical_potential.h"

using namespace std;

class particleList {
 private:
    string particleListFilename;
    vector<particle*> partList;
    vector<particle*> stable_particle_list;
    double edSystem, sdSystem, pressureSys, net_baryon_density;

 public:
    particleList(string particleTableName);
    ~particleList();
    int getParticlelistSize() {return(partList.size());}
    void readParticlelistTable(string tableName);
    void sort_particle_list_in_particle_mass();
    void calculate_particle_decay_probability();
    void calculate_particle_decays();
    void calculate_particle_chemical_potential2(
                    double Temperature, Chemical_potential* mu_tb);
    void calculateSystemEOS(double mu_B = 0.0, double mu_S = 0.0);
    void calculate_particle_mu(double mu_B, double mu_S);
    void calculate_particle_yield(double Temperature, double mu_B = 0.0,
                                  double mu_S = 0.0);
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


#endif  // SRC_PARTICLELIST_H_
