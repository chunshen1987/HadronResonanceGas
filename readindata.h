#ifndef READINDATA_H
#define READINDATA_H

#include<fstream>
#include<vector>
using namespace std;

typedef struct
{
   int daughter_monval;
   string daughter_name;
   vector<int> parent_monval;
   vector<double> parent_mass;
   vector<double> decay_contribution;
}decay_contribution_table;

typedef struct
{
  int monval;			// Montecarlo number according PDG
  string name;
  double mass;
  double width;
  int gspin;			// spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin;			// isospin degeneracy
  int charge;
  int decays;			// amount of decays listed for this resonance
  int stable;			// defines whether this particle is considered as stable
  double yield;
  double decay_contribution;
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  double mu;
  int sing;       //Bose-Einstein or Dirac-Fermi statistics
}particle_info;

int read_resonance(particle_info* particle);
void read_decdat_mu(int N_stable, double* particle_mu);
void calculate_particle_mu(int Nparticle, particle_info* particle, double* particle_mu);
void calculate_particle_yield(int Nparticle, particle_info* particle, double Temperature);
int get_particle_idx(particle_info* particle, int Nparticle, int particle_monval);

#endif
