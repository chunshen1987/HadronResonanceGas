// Copyright 2016 Chun Shen
#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#include<string>

class particle {
 private:
    double hbarC;           //! GeV*fm
    int monval;             //! Monte-Carlo number according PDG
    std::string name;       //! particle name
    double mass;            //! particle mass (GeV)
    double width;           //! decay width
    int gspin;              //! spin degeneracy
    int baryon;             //! baryon number
    int strange;            //! strangeness
    int charm;              //! charmness
    int bottom;             //! bottomness
    int gisospin;           //! isospin degeneracy
    int charge;             //! charge
    int decays;             //! amount of decays listed for this resonance
    int stable;             //! defines whether this particle is stable
    double mu;              //! chemical potential
    int sign;               //! Bosons or Fermions

    double yield;                //! particle thermal yield at given T and mu
    double stable_yield;         //! particle yield after decays

    double ed, sd, pressure;     //! thermodynamic quantities at given T and mu

    int NdecayChannel;           //! number of decay channels
    double* decays_branchratio;  //! branching ratio of each decay channel
    // number of daughter particles of each decay channel
    int* decays_Npart;
    int** decays_part;           //! identity of daughter particles

    // array to record particle decay probability
    // for the particle decays into stable ones
    double* decay_probability;

    int channelIdx;
    int trunOrder;               // truncated order in the summation

 public:
    particle(int monval_in, std::string name_in, double mass_in,
             double width_in, int gspin_in, int baryon_in, int strange_in,
             int charm_in, int bottom_in, int gisospin_in, int charge_in,
             int NdecayChannel_in);
    ~particle();

    void addResonancedecays(double branchratio, int Npart, int* decays_part);
    void create_decay_probability_array(int dimension) {
        decay_probability = new double[dimension];
        for (int i = 0; i < dimension; i++)
            decay_probability[i] = 0.0;
    }

    void calculateChemicalpotential(double mu_B, double mu_S, double mu_Q);
    void calculateParticleYield(double Temperature, double mu_B,
                                double mu_S, double mu_Q);
    double calculateEnergydensity(double Temperature);
    double calculatePressure(double Temperature);
    double calculateEntropydensity(double Temperature);

    int getAntiparticleMonval();

    int getMonval() {return(monval);}
    std::string getName() {return(name);}
    double getMass() {return(mass);}
    int getBaryon() const {return(baryon);}
    int getStrangeness() const {return(strange);}
    int getElectricCharge() const {return(charge);}
    int getSpinfactor() {return(gspin);}
    double getMu() {return(mu);}
    int getSign() {return(sign);}
    double getParticleYield() {return(yield);}
    double getEnergyDensity() {return(ed);}
    void set_particle_stable_yield(double yield_in) {stable_yield = yield_in;}
    double get_particle_stable_yield() {return(stable_yield);}
    int getNdecays() {return(decays);}
    int getNdecayChannel() {return(NdecayChannel);}
    int getdecaysNpart(int i) {return(decays_Npart[i]);}
    int getdecays_part(int i, int j) {return(decays_part[i][j]);}
    double getdecays_branchratio(int i) {return(decays_branchratio[i]);}
    double getdecay_probability(int i) {return(decay_probability[i]);}
    int getStable() {return(stable);}
    void setStable(int s) {stable = s;}
    void setMu(double chem) {mu = chem;}
    void set_decay_probability(int i, double val) {
        decay_probability[i] = val;
    }
    double calculate_dsdmu(double Temperature);
    double calculate_deoverTdmu(double Temperature);
    double calculate_dPoverTdmu(double Temperature);
    double calculate_dndmu(double Temperature);

    void getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients);
};

#endif  // SRC_PARTICLE_H_
