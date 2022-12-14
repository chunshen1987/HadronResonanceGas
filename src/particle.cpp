// Copyright 2016 Chun Shen
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

#include "./particle.h"

using namespace std;

particle::particle(int monval_in, string name_in, double mass_in,
                   double width_in, int gspin_in, int baryon_in,
                   int strange_in, int charm_in, int bottom_in,
                   int gisospin_in, int charge_in, int NdecayChannel_in) {
    hbarC = 0.19733;
    trunOrder = 5;
    monval = monval_in;
    name = name_in;
    mass = mass_in;
    width = width_in;
    gspin = gspin_in;
    baryon = baryon_in;
    strange = strange_in;
    charm = charm_in;
    bottom = bottom_in;
    gisospin = gisospin_in;
    charge = charge_in;
    NdecayChannel = NdecayChannel_in;
    stable = 0;
    stable_yield = 0.0;

    // determine Bose/Fermi statistic for particle
    if (baryon == 0)
       sign = -1;
    else
       sign = 1;

    decays_branchratio = new double[NdecayChannel];
    decays_Npart = new int[NdecayChannel];
    decays_part = new int* [NdecayChannel];

    channelIdx = 0;
}


particle::~particle() {
    delete [] decays_branchratio;
    delete [] decays_Npart;
    for (int i = 0; i < NdecayChannel; i++)
       delete [] decays_part[i];
    delete [] decays_part;
    delete [] decay_probability;
}


void particle::addResonancedecays(double branchratio, int Npart,
                                  int* decayChannelparts) {
// add resonance decay channel for particle
    if (channelIdx > NdecayChannel-1) {
        cout << "Warning: channelidx exceed number of decay channels! "
             << "Please check" << endl;
        exit(1);
    }
    decays_branchratio[channelIdx] = branchratio;
    decays_Npart[channelIdx] = Npart;
    decays_part[channelIdx] = new int[Npart];
    for (int i = 0; i < Npart; i++)
       decays_part[channelIdx][i] = decayChannelparts[i];
    channelIdx++;

    if (NdecayChannel == 1 && decays_Npart[0] == 1
       && decays_part[0][0] == monval
       && fabs(branchratio - 1.0) < 1e-15) {
       // particle is stable
        stable = 1;
    }
}


int particle::getAntiparticleMonval() {
// return the corresponding anti-particle' Monte-Carlo value
    if (baryon == 0 && charge == 0 && strange == 0)
        return(monval);
    else
        return(-monval);
}


void particle::calculateChemicalpotential(double mu_B, double mu_S,
                                          double mu_Q) {
    mu = mu_B*baryon + mu_S*strange + mu_Q*charge;
}


void particle::calculateParticleYield(double Temperature, double mu_B,
                                      double mu_S, double mu_Q) {
    // this function compute particle thermal yield at given T, mu
    //int sf_expint_truncate_order = 1;

    double N_eq = 0.0;                  // equilibrium contribution
    //double deltaN_bulk_term1 = 0.0;     // contribution from bulk delta f
    //double deltaN_bulk_term2 = 0.0;     // contribution from bulk delta f
    //double deltaN_qmu_term1 = 0.0;      // contribution from baryon diffusion
    //double deltaN_qmu_term2 = 0.0;      // contribution from baryon diffusion

    double beta = 1./Temperature;
    calculateChemicalpotential(mu_B, mu_S, mu_Q);
    double lambda = exp(beta*mu);  // fugacity factor

    double prefactor = gspin/(2*M_PI*M_PI)/hbarC/hbarC/hbarC;
    double mbeta = mass*beta;

    for (int n = 1; n < trunOrder; n++) {
        double arg = n*mbeta;  // argument inside bessel functions
        double theta = pow(-sign, n-1);
        double fugacity = pow(lambda, n);
        double K_2 = gsl_sf_bessel_Kn(2, arg);
        //double K_1 = gsl_sf_bessel_K1(arg);
        // cout << arg << "  " << K_1 << "  " << K_2;

        // equilibrium contribution
        N_eq += theta/n*fugacity*K_2;

        // bulk viscous contribution
        //deltaN_bulk_term1 += theta*fugacity*(mass*beta*K_1 + 3./n*K_2);
        //deltaN_bulk_term2 += theta*fugacity*K_1;

        // baryon diffusion contribution
        //deltaN_qmu_term1 += theta/n*fugacity*K_2;
        //double I_1_1 = exp(-arg)/arg*(2./(arg*arg) + 2./arg - 1./2.);
        //double I_1_2 = 3./8.*gsl_sf_expint_E2(arg);
        //double I_1_n = I_1_1 + I_1_2;

        //double double_factorial = 1.;  // record (2k-5)!!
        //double factorial = 2.;         // record k! start with 2!
        //double factor_2_to_k_power = 4.;      // record 2^k start with 2^2
        //for (int k = 3; k < sf_expint_truncate_order; k++) {
        //    double_factorial *= (2*k - 5);
        //    factorial *= k;
        //    factor_2_to_k_power *= 2;
        //    double I_1_k = (3.*double_factorial/factor_2_to_k_power/factorial
        //                    *gsl_sf_expint_En(2*k-2, arg));
        //    I_1_n += I_1_k;
        //}
        //I_1_n = -(mbeta*mbeta*mbeta)*I_1_n;
        //deltaN_qmu_term2 += n*theta*fugacity*I_1_n;
    }

    // equilibrium contribution
    double prefactor_Neq = mass*mass*Temperature;
    N_eq = prefactor*prefactor_Neq*N_eq;

    // bulk viscous contribution
    //deltaN_bulk_term1 = mass*mass/beta*deltaN_bulk_term1;
    //deltaN_bulk_term2 = mass*mass*mass/3.*deltaN_bulk_term2;

    // baryon diffusion contribution
    //deltaN_qmu_term1 = mass*mass/(beta*beta)*deltaN_qmu_term1;
    //deltaN_qmu_term2 = 1./(3.*beta*beta*beta)*deltaN_qmu_term2;

    yield = N_eq;
}


double particle::calculateEnergydensity(double Temperature) {
    double results;
    results = 0.0;
    double prefactor = gspin/(2*M_PI*M_PI)*pow(mass, 4);
    for (int j = 0; j <trunOrder; j++) {
        double arg = (j+1)*mass/Temperature;
        double lambda = exp(mu/Temperature);
        results += (pow((-1.0)*sign, j)*pow(lambda, j+1)
                    *(3.*gsl_sf_bessel_Kn(2, arg)/(arg*arg)
                      + gsl_sf_bessel_Kn(1, arg)/arg));
    }
    results = results*prefactor;
    ed = results/pow(hbarC, 3);    // unit: GeV/fm^3
    return(ed);
}


double particle::calculatePressure(double Temperature) {
    double results;
    results = 0.0;
    double prefactor = gspin/(2*M_PI*M_PI)*pow(mass, 2)*pow(Temperature, 2);
    for (int j = 0; j < trunOrder; j++) {
        double arg = (j+1)*mass/Temperature;
        double lambda = exp(mu/Temperature);
        results += (pow((-1.0)*sign, j)/pow(j+1., 2)
                    *pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg));
    }
    results = results*prefactor;
    pressure = results/pow(hbarC, 3);    // unit : GeV/fm^3
    return(pressure);
}


//! calculate the entropy density using the first law of thermodynamics
//! at give T and mu
double particle::calculateEntropydensity(double Temperature) {
    sd = (ed + pressure - mu*yield)/Temperature;    // unit : 1/fm^3
    return(sd);
}


//! calculate the first order derivative dn/dmu [1/GeV]
double particle::calculate_dndmu(double Temperature) {
    double results;
    results = 0.0;
    double prefactor = gspin/(2*M_PI*M_PI)*mass;
    for (int j = 0; j < trunOrder; j++) {
        double arg = (j+1)*mass/Temperature;
        double lambda = exp(mu/Temperature);
        results +=
            pow((-1.0)*sign, j)*pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg);
    }
    results = results*prefactor;
    return(results);
}


//! calculate the first order derivative d(P/T)/dmu [1/(GeV fm^3)]
double particle::calculate_dPoverTdmu(double Temperature) {
    double results;
    results = 0.0;
    double prefactor = gspin/(2*M_PI*M_PI)*mass*mass;
    for (int j = 0; j < trunOrder; j++) {
        double arg = (j+1)*mass/Temperature;
        double lambda = exp(mu/Temperature);
       results += (pow((-1.0)*sign, j)*pow(lambda, j+1)/(j+1)
                   *gsl_sf_bessel_Kn(2, arg));
    }
    results = results*prefactor/pow(hbarC, 3);
    return(results);
}


//! calculate the first order derivative d(e/T)/dmu [1/(GeV fm^3)]
double particle::calculate_deoverTdmu(double Temperature) {
    double results;
    results = 0.0;
    double prefactor = gspin/(2*M_PI*M_PI)*pow(mass, 4);
    for (int j = 0; j < trunOrder; j++) {
        double arg = (j+1)*mass/Temperature;
        double lambda = exp(mu/Temperature);
        results += (pow((-1.0)*sign, j)*pow(lambda, j+1)*(j+1)
                    *(3*gsl_sf_bessel_Kn(2, arg)/(arg*arg)
                      + gsl_sf_bessel_Kn(1, arg)/arg));
    }
    results = results*prefactor/pow(hbarC, 3);
    return(results);
}


//! calculate the first order derivative ds/dmu [1/(GeV fm^3)]
double particle::calculate_dsdmu(double Temperature) {
    double dPoverTdmu = calculate_dPoverTdmu(Temperature);
    double deoverTdmu = calculate_deoverTdmu(Temperature);
    double dndmu = calculate_dndmu(Temperature);
    double dsdmu = dPoverTdmu - yield - mu*dndmu + deoverTdmu;
    return(dsdmu);
}
