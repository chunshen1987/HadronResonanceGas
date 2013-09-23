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

/*
void calculate_particle_yield(int Nparticle, particle_info* particle, double Temperature)
{
   double results;
   int order = 10;
   for(int i=1; i<Nparticle; i++)
   {
      results = 0.0;
      double prefactor = particle[i].gspin/(2*M_PI*M_PI)*particle[i].mass*particle[i].mass;
      for(int j=0; j<order; j++)
      {
         double arg = (j+1)*particle[i].mass/Temperature;
         double lambda = exp(particle[i].mu/Temperature);
         results += pow((-1.0)*particle[i].sing, j)/(j+1)*pow(lambda, j+1)*gsl_sf_bessel_Kn(2, arg);
      }
      results = results*prefactor;
      particle[i].yield = results;
   }
   return;
}

void calculate_particle_mu(int Nparticle, particle_info* particle, double* particle_mu)
{
   for(int i=0; i<Nparticle; i++) particle[i].mu = 0.0;
   return;
}

void read_decdat_mu(int N_stable, double* particle_mu)
{
  cout<<" -- Read chemical potential for stable particles...";
  ostringstream decdat_mu_stream;
  double dummy;
  decdat_mu_stream << "results/decdat_mu.dat";
  ifstream decdat_mu(decdat_mu_stream.str().c_str());

  decdat_mu >> dummy;  //not used in the code plz ignore it
  for(int i=0; i<N_stable; i++)
     decdat_mu >> particle_mu[i];
  cout<<"done" << endl;
  return;
}

int get_particle_idx(particle_info* particle, int Nparticle, int particle_monval)
{
   int idx = 0;
   int i;
   for(i=0; i<Nparticle; i++)
      if(particle[i].monval == particle_monval)
      {
         idx = i;
         return(idx);
      }
   if(i == Nparticle)
   {
      cout << "get_particle_idx: Error : can not find particle index in the particle list." << endl;
      exit(0);
   }
   return(0);
}*/
