#ifndef particle_H
#define particle_H

#include<string>
using namespace std;

class particle
{
    private:
       int monval;			// Monte-Carlo number according PDG
       string name;                 // particle name
       double mass;                 // particle mass (GeV)
       double width;                // decay width
       int gspin;			      // spin degeneracy
       int baryon;                  // baryon number
       int strange;                 // strangeness
       int charm;                   // charmness
       int bottom;                  // bottomness
       int gisospin;			// isospin degeneracy
       int charge;                  // charge
       int decays;			// amount of decays listed for this resonance
       int stable;			// defines whether this particle is considered as stable
       double mu;                   // chemical potential
       int sign;                    // Bosons or Fermions

       double yield;                // particle yield at given T and mu
     
       int NdecayChannel;           // number of decay channels
       double* decays_branchratio;  // branching ratio of each decay channel
       int* decays_Npart;           // number of daughter particles of each decay channel
       int** decays_part;           // identity of daughter particles

       int channelIdx;
    public:
       particle(int monval_in, string name_in, double mass_in, double width_in, int gspin_in, int baryon_in, int strange_in, int charm_in, int bottom_in, int gisospin_in, int charge_in, int NdecayChannel_in);
       ~particle();
       
       void addResonancedecays(double branchratio, int Npart, int* decays_part);
       int getAntiparticleMonval();
       int getMonval() {return(monval);};
       string getName() {return(name);};
       double getMass() {return(mass);};
       int getBaryon() {return(baryon);};
       int getCharge() {return(charge);};
       int getSpinfactor() {return(gspin);};
       double getMu() {return(mu);};
       int getSign() {return(sign);};

};

#endif
