#include<iostream>
#include<cmath>

#include<gsl/gsl_sf_bessel.h>

#include "particle.h"
using namespace std;

particle::particle(int monval_in, string name_in, double mass_in, double width_in, int gspin_in, int baryon_in, int strange_in, int charm_in, int bottom_in, int gisospin_in, int charge_in, int NdecayChannel_in)
{
    hbarC = 0.19733;
    trunOrder = 10;
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

    //determine Bose/Fermi statistic for particle
    if(baryon == 0)
       sign = -1;
    else
       sign = 1;

    decays_branchratio = new double [NdecayChannel];
    decays_Npart = new int [NdecayChannel];
    decays_part = new int* [NdecayChannel];
    
    channelIdx = 0;
}

particle::~particle()
{  
    delete [] decays_branchratio;
    delete [] decays_Npart;
    for(int i = 0; i < NdecayChannel; i++)
       delete [] decays_part[i];
    delete [] decays_part;
}

void particle::addResonancedecays(double branchratio, int Npart, int* decayChannelparts)
//add resonance decay channel for particle
{
    if(channelIdx > NdecayChannel-1)
    {
        cout << "Warning: channelidx exceed number of decay channels! Please check" << endl;
        exit(1);
    }
    decays_branchratio[channelIdx] = branchratio;
    decays_Npart[channelIdx] = Npart;
    decays_part[channelIdx] = new int [Npart];
    for(int i = 0; i < Npart; i++)
       decays_part[channelIdx][i] = decayChannelparts[i];
    channelIdx++;
    if(NdecayChannel == 1 && decays_Npart[0] == 1) stable = 1;

    return;
}

int particle::getAntiparticleMonval()
//return the corresponding anti-particle' Monte-Carlo value
{
   if(baryon == 0 && charge == 0 && strange == 0)
      return(monval);
   else
      return(-monval);
}

void particle::calculateChemicalpotential(double mu_B, double mu_S)
{
   mu = mu_B*baryon + mu_S*strange;
   return;
}

void particle::calculateParticleYield(double Temperature)
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*mass*mass;
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)/(j+1)*pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg);
   }
   results = results*prefactor;
   yield = results;
   return;
}

double particle::calculateEnergydensity(double Temperature)
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*pow(mass,4);
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)*pow(lambda, j+1)*(3.*gsl_sf_bessel_Kn(2, arg)/(arg*arg) + gsl_sf_bessel_Kn(1, arg)/arg);
   }
   results = results*prefactor;
   ed = results/pow(hbarC,3);    // unit: GeV/fm^3
   return(ed);
}

double particle::calculatePressure(double Temperature)
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*pow(mass,2)*pow(Temperature, 2);
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)/pow(j+1.,2)*pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg);
   }
   results = results*prefactor;
   pressure = results/pow(hbarC, 3);    // unit : GeV/fm^3
   return(pressure);
}

double particle::calculateEntropydensity(double Temperature)
// calculate the entropy density using the first law of thermodynamics at give T and mu
{
   sd = (ed + pressure - mu*yield)/Temperature;    // unit : 1/fm^3
   return(sd);
}

double particle::calculate_dndmu(double Temperature)
// calculate the first order derivative dn/dmu [1/GeV]
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*mass;
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)*pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg);
   }
   results = results*prefactor;
   return(results);
}

double particle::calculate_dPoverTdmu(double Temperature)
// calculate the first order derivative d(P/T)/dmu [1/(GeV fm^3)]
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*mass*mass;
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)*pow(lambda, j+1)/(j+1)*gsl_sf_bessel_Kn(2, arg);
   }
   results = results*prefactor/pow(hbarC, 3);
   return(results);
}

double particle::calculate_deoverTdmu(double Temperature)
// calculate the first order derivative d(e/T)/dmu [1/(GeV fm^3)]
{
   double results;
   results = 0.0;
   double prefactor = gspin/(2*M_PI*M_PI)*pow(mass, 4);
   for(int j=0; j<trunOrder; j++)
   {
      double arg = (j+1)*mass/Temperature;
      double lambda = exp(mu/Temperature);
      results += pow((-1.0)*sign, j)*pow(lambda, j+1)*(j+1)*(3*gsl_sf_bessel_Kn(2, arg)/(arg*arg) + gsl_sf_bessel_Kn(1, arg)/arg);
   }
   results = results*prefactor/pow(hbarC, 3);
   return(results);
}

double particle::calculate_dsdmu(double Temperature)
// calculate the first order derivative ds/dmu [1/(GeV fm^3)]
{
   double dPoverTdmu = calculate_dPoverTdmu(Temperature);
   double deoverTdmu = calculate_deoverTdmu(Temperature);
   double dndmu = calculate_dndmu(Temperature);
   double dsdmu = dPoverTdmu - yield - mu*dndmu + deoverTdmu;
   return(dsdmu);
}
