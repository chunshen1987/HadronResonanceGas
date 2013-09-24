#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include<gsl/gsl_sf_bessel.h>

#include "particleList.h"
using namespace std;


particleList::particleList(string particleTableName)
{
   particleListFilename = particleTableName; //filename of pdg data file
}

particleList::~particleList()
{

}

void particleList::readParticlelistTable()
//read in particle information from pdg data file
{
   cout << "Reading in particle resonance decay table...";
   ifstream resofile(particleListFilename.c_str());
   int monval;
   string name;
   double mass, width;
   int gspin, gisospin;
   int baryon, strange, charm, bottom, charge;
   int decays, decayNpart;
   double decayBranchratio;
   int decayPart[5] = {0,0,0,0,0};

   int dummy_int;
   while(1)
   {
      resofile >> monval;
      if(resofile.eof()) break;
      resofile >> name;
      resofile >> mass;
      resofile >> width;
      resofile >> gspin;
      resofile >> baryon;
      resofile >> strange;
      resofile >> charm;
      resofile >> bottom;
      resofile >> gisospin;
      resofile >> charge;
      resofile >> decays;
      partList.push_back(new particle(monval, name, mass, width, gspin, baryon, strange, charm, bottom, gisospin, charge, decays));
      if(baryon == 1)
      {
         ostringstream antiname;
         antiname << "Anti-" << name;
         partList.push_back(new particle(-monval, antiname.str(), mass, width, gspin, -baryon, -strange, -charm, -bottom, gisospin, -charge, decays));
      }
      for (int j = 0; j < decays; j++)
      {
         resofile >> dummy_int;
         resofile >> decayNpart;
         resofile >> decayBranchratio;
         resofile >> decayPart[0];
         resofile >> decayPart[1];
         resofile >> decayPart[2];
         resofile >> decayPart[3];
         resofile >> decayPart[4];
         decayNpart = abs(decayNpart);
         int* tempptr = new int [decayNpart];
         for(int ipart = 0; ipart < decayNpart; ipart++)
            tempptr[ipart] = decayPart[ipart];

         if(baryon == 0)
            partList.back()->addResonancedecays(decayBranchratio, decayNpart, tempptr);
         else
         {
            partList.at(partList.size()-2)->addResonancedecays(decayBranchratio, decayNpart, tempptr);
            for(int ipart = 0; ipart < decayNpart; ipart++)
            {
               int particleId = get_particle_idx(tempptr[ipart]);
               tempptr[ipart] = partList[particleId]->getAntiparticleMonval();
            }
            partList.back()->addResonancedecays(decayBranchratio, decayNpart, tempptr);
         }
         delete [] tempptr;
      }
   }
   resofile.close();
   partList.erase(partList.begin());  //delete gamma
   cout << "done! Antiparticles are added!" << endl;
   cout << "There are totally " << partList.size() << " particles." << endl;
   return;
}

int particleList::get_particle_idx(int particle_monval)
// return the idx in particleList for given particle Monte-Carlo number
{
   for(int i = 0; i < partList.size(); i++)
      if(partList[i]->getMonval() == particle_monval)
         return(i);
   cout << "Warning: can not fine particle " << particle_monval << endl;
   exit(1);
}

void particleList::calculate_particle_mu(double mu_B, double mu_S)
// calculate particle chemical potentials
// need to add support for partial chemical equilibrium
{
   for(int i = 0; i < partList.size(); i++)
      partList[i]->calculateChemicalpotential(mu_B, mu_S);
   return;
}

void particleList::calculate_particle_yield(double Temperature, double mu_B, double mu_S)
//calculate particle yield
{
   calculate_particle_mu(mu_B, mu_S);
   for(int i = 0; i < partList.size(); i++)  
   {
      partList[i]->calculateParticleYield(Temperature);
      cout << partList[i]->getName() << " : " << partList[i]->getParticleYield() << endl;
   }
   return;
}

void particleList::calculateSystemenergyDensity(double Temperature, double mu_B, double mu_S)
//calculate the energy density of the system at given T and mu
{
   double result = 0.0e0;
   calculate_particle_mu(mu_B, mu_S);
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculateEnergydensity(Temperature);
   edSystem = result;
   return;
}

void particleList::calculateSystemPressure(double Temperature, double mu_B, double mu_S)
//calculate the pressure of the system at given T and mu
{
   double result = 0.0e0;
   calculate_particle_mu(mu_B, mu_S);
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculatePressure(Temperature);
   pressureSys = result;
   return;
}

void particleList::calculateSystementropyDensity(double Temperature, double mu_B, double mu_S)
//calculate the entropy density of the system at given T and mu
{
   double result = 0.0e0;
   calculate_particle_mu(mu_B, mu_S);
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculateEntropydensity(Temperature);
   sdSystem = result;
   cout << sdSystem << endl;
   return;
}
