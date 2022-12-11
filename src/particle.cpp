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
    return;
}

int particle::getAntiparticleMonval() {
// return the corresponding anti-particle' Monte-Carlo value
    if (baryon == 0 && charge == 0 && strange == 0)
        return(monval);
    else
        return(-monval);
}

void particle::calculateChemicalpotential(double mu_B, double mu_S) {
    mu = mu_B*baryon + mu_S*strange;
    return;
}

void particle::calculateParticleYield(double Temperature, double mu_B,
                                      double mu_S) {
    // this function compute particle thermal yield at given T, mu
    //int sf_expint_truncate_order = 1;

    double N_eq = 0.0;                  // equilibrium contribution
    //double deltaN_bulk_term1 = 0.0;     // contribution from bulk delta f
    //double deltaN_bulk_term2 = 0.0;     // contribution from bulk delta f
    //double deltaN_qmu_term1 = 0.0;      // contribution from baryon diffusion
    //double deltaN_qmu_term2 = 0.0;      // contribution from baryon diffusion

    double beta = 1./Temperature;
    calculateChemicalpotential(mu_B, mu_S);
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
    return;
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

double particle::calculateEntropydensity(double Temperature) {
// calculate the entropy density using the first law of thermodynamics
// at give T and mu
    sd = (ed + pressure - mu*yield)/Temperature;    // unit : 1/fm^3
    return(sd);
}

double particle::calculate_dndmu(double Temperature) {
// calculate the first order derivative dn/dmu [1/GeV]
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

double particle::calculate_dPoverTdmu(double Temperature) {
// calculate the first order derivative d(P/T)/dmu [1/(GeV fm^3)]
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

double particle::calculate_deoverTdmu(double Temperature) {
// calculate the first order derivative d(e/T)/dmu [1/(GeV fm^3)]
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

double particle::calculate_dsdmu(double Temperature) {
// calculate the first order derivative ds/dmu [1/(GeV fm^3)]
    double dPoverTdmu = calculate_dPoverTdmu(Temperature);
    double deoverTdmu = calculate_deoverTdmu(Temperature);
    double dndmu = calculate_dndmu(Temperature);
    double dsdmu = dPoverTdmu - yield - mu*dndmu + deoverTdmu;
    return(dsdmu);
}

// void particle::getbulkvisCoefficients(
//                 double Tdec, double* bulkvisCoefficients) {
// // get transport coefficient for bulk delta f
//     double Tdec_fm = Tdec/hbarC;  // [1/fm]
//     double Tdec_fm_power[11];     // cache the polynomial power of Tdec_fm
//     Tdec_fm_power[1] = Tdec_fm;
//     for (int ipower = 2; ipower < 11; ipower++)
//         Tdec_fm_power[ipower] = Tdec_fm_power[ipower-1]*Tdec_fm;
//     if (bulk_deltaf_kind == 1) {    // relaxation type
//         // parameterization from JF
//         // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
//         // Both fits are reliable between T=100 -- 180 MeV,
//         // do not trust it beyond
//         bulkvisCoefficients[0] = (642096.624265727
//                                   - 8163329.49562861*Tdec_fm_power[1]
//                                   + 47162768.4292073*Tdec_fm_power[2]
//                                   - 162590040.002683*Tdec_fm_power[3]
//                                   + 369637951.096896*Tdec_fm_power[4]
//                                   - 578181331.809836*Tdec_fm_power[5]
//                                   + 629434830.225675*Tdec_fm_power[6]
//                                   - 470493661.096657*Tdec_fm_power[7]
//                                   + 230936465.421*Tdec_fm_power[8]
//                                   - 67175218.4629078*Tdec_fm_power[9]
//                                   + 8789472.32652964*Tdec_fm_power[10]);
// 
//         bulkvisCoefficients[1] = (1.18171174036192
//                                   - 17.6740645873717*Tdec_fm_power[1]
//                                   + 136.298469057177*Tdec_fm_power[2]
//                                   - 635.999435106846*Tdec_fm_power[3]
//                                   + 1918.77100633321*Tdec_fm_power[4]
//                                   - 3836.32258307711*Tdec_fm_power[5]
//                                   + 5136.35746882372*Tdec_fm_power[6]
//                                   - 4566.22991441914*Tdec_fm_power[7]
//                                   + 2593.45375240886*Tdec_fm_power[8]
//                                   - 853.908199724349*Tdec_fm_power[9]
//                                   + 124.260460450113*Tdec_fm_power[10]);
//     } else if (bulk_deltaf_kind == 2) {
//         // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
//         // Both fits are reliable between T=100 -- 180 MeV
//         // do not trust it beyond
//         bulkvisCoefficients[0] = (21091365.1182649
//                                   - 290482229.281782*Tdec_fm_power[1]
//                                   + 1800423055.01882*Tdec_fm_power[2]
//                                   - 6608608560.99887*Tdec_fm_power[3]
//                                   + 15900800422.7138*Tdec_fm_power[4]
//                                   - 26194517161.8205*Tdec_fm_power[5]
//                                   + 29912485360.2916*Tdec_fm_power[6]
//                                   - 23375101221.2855*Tdec_fm_power[7]
//                                   + 11960898238.0134*Tdec_fm_power[8]
//                                   - 3618358144.18576*Tdec_fm_power[9]
//                                   + 491369134.205902*Tdec_fm_power[10]);
// 
//         bulkvisCoefficients[1] = (4007863.29316896
//                                   - 55199395.3534188*Tdec_fm_power[1]
//                                   + 342115196.396492*Tdec_fm_power[2]
//                                   - 1255681487.77798*Tdec_fm_power[3]
//                                   + 3021026280.08401*Tdec_fm_power[4]
//                                   - 4976331606.85766*Tdec_fm_power[5]
//                                   + 5682163732.74188*Tdec_fm_power[6]
//                                   - 4439937810.57449*Tdec_fm_power[7]
//                                   + 2271692965.05568*Tdec_fm_power[8]
//                                   - 687164038.128814*Tdec_fm_power[9]
//                                   + 93308348.3137008*Tdec_fm_power[10]);
//     } else if (bulk_deltaf_kind == 3) {
//         bulkvisCoefficients[0] = (160421664.93603
//                                   - 2212807124.97991*Tdec_fm_power[1]
//                                   + 13707913981.1425*Tdec_fm_power[2]
//                                   - 50204536518.1767*Tdec_fm_power[3]
//                                   + 120354649094.362*Tdec_fm_power[4]
//                                   - 197298426823.223*Tdec_fm_power[5]
//                                   + 223953760788.288*Tdec_fm_power[6]
//                                   - 173790947240.829*Tdec_fm_power[7]
//                                   + 88231322888.0423*Tdec_fm_power[8]
//                                   - 26461154892.6963*Tdec_fm_power[9]
//                                   + 3559805050.19592*Tdec_fm_power[10]);
//         bulkvisCoefficients[1] = (33369186.2536556
//                                   - 460293490.420478*Tdec_fm_power[1]
//                                   + 2851449676.09981*Tdec_fm_power[2]
//                                   - 10443297927.601*Tdec_fm_power[3]
//                                   + 25035517099.7809*Tdec_fm_power[4]
//                                   - 41040777943.4963*Tdec_fm_power[5]
//                                   + 46585225878.8723*Tdec_fm_power[6]
//                                   - 36150531001.3718*Tdec_fm_power[7]
//                                   + 18353035766.9323*Tdec_fm_power[8]
//                                   - 5504165325.05431*Tdec_fm_power[9]
//                                   + 740468257.784873*Tdec_fm_power[10]);
//     } else if (bulk_deltaf_kind == 4) {
//         bulkvisCoefficients[0] = (1167272041.90731
//                                   - 16378866444.6842*Tdec_fm_power[1]
//                                   + 103037615761.617*Tdec_fm_power[2]
//                                   - 382670727905.111*Tdec_fm_power[3]
//                                   + 929111866739.436*Tdec_fm_power[4]
//                                   - 1540948583116.54*Tdec_fm_power[5]
//                                   + 1767975890298.1*Tdec_fm_power[6]
//                                   - 1385606389545*Tdec_fm_power[7]
//                                   + 709922576963.213*Tdec_fm_power[8]
//                                   - 214726945096.326*Tdec_fm_power[9]
//                                   + 29116298091.9219*Tdec_fm_power[10]);
//         bulkvisCoefficients[1] = (5103633637.7213
//                                   - 71612903872.8163*Tdec_fm_power[1]
//                                   + 450509014334.964*Tdec_fm_power[2]
//                                   - 1673143669281.46*Tdec_fm_power[3]
//                                   + 4062340452589.89*Tdec_fm_power[4]
//                                   - 6737468792456.4*Tdec_fm_power[5]
//                                   + 7730102407679.65*Tdec_fm_power[6]
//                                   - 6058276038129.83*Tdec_fm_power[7]
//                                   + 3103990764357.81*Tdec_fm_power[8]
//                                   - 938850005883.612*Tdec_fm_power[9]
//                                   + 127305171097.249*Tdec_fm_power[10]);
//     }
//     return;
// }
