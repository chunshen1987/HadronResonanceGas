#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include<gsl/gsl_sf_bessel.h>

#include "readindata.h"
using namespace std;

int read_resonance(particle_info* particle)
{
   int Nparticle=0; 
   cout << "Reading in particle resonance decay table...";
   ifstream resofile("EOS/pdg.dat");
   int local_i = 0;
   int dummy_int;
   while (!resofile.eof())
   {
      resofile >> particle[local_i].monval;
      resofile >> particle[local_i].name;
      resofile >> particle[local_i].mass;
      resofile >> particle[local_i].width;
      resofile >> particle[local_i].gspin;	      //spin degeneracy
      resofile >> particle[local_i].baryon;
      resofile >> particle[local_i].strange;
      resofile >> particle[local_i].charm;
      resofile >> particle[local_i].bottom;
      resofile >> particle[local_i].gisospin;     //isospin degeneracy
      resofile >> particle[local_i].charge;
      resofile >> particle[local_i].decays;
      for (int j = 0; j < particle[local_i].decays; j++)
      {
         resofile >> dummy_int;
         resofile >> particle[local_i].decays_Npart[j];
         resofile >> particle[local_i].decays_branchratio[j];
         resofile >> particle[local_i].decays_part[j][0];
         resofile >> particle[local_i].decays_part[j][1];
         resofile >> particle[local_i].decays_part[j][2];
         resofile >> particle[local_i].decays_part[j][3];
         resofile >> particle[local_i].decays_part[j][4];
      }
      
      //decide whether particle is stable under strong interactions
      if(particle[local_i].decays_Npart[0] == 1) 	      
         particle[local_i].stable = 1;
      else
         particle[local_i].stable = 0;

      //add anti-particle entry
      if(particle[local_i].baryon == 1)
      {
         local_i++;
         particle[local_i].monval = -particle[local_i-1].monval;
         ostringstream antiname;
         antiname << "Anti-" << particle[local_i-1].name;
         particle[local_i].name = antiname.str();
         particle[local_i].mass = particle[local_i-1].mass;
         particle[local_i].width = particle[local_i-1].width;
         particle[local_i].gspin = particle[local_i-1].gspin;
         particle[local_i].baryon = -particle[local_i-1].baryon;
         particle[local_i].strange = -particle[local_i-1].strange;
         particle[local_i].charm = -particle[local_i-1].charm;
         particle[local_i].bottom = -particle[local_i-1].bottom;
         particle[local_i].gisospin = particle[local_i-1].gisospin;
         particle[local_i].charge = -particle[local_i-1].charge;
         particle[local_i].decays = particle[local_i-1].decays;
         particle[local_i].stable = particle[local_i-1].stable;
         for (int j = 0; j < particle[local_i].decays; j++)
         {
            particle[local_i].decays_Npart[j]=particle[local_i-1].decays_Npart[j];
            particle[local_i].decays_branchratio[j]=particle[local_i-1].decays_branchratio[j];
            for (int k=0; k< Maxdecaypart; k++)
            {
               int idx = 0;  
               for(int ii=0; ii < local_i; ii++) // find the index for decay particle
               {
                  if(particle[local_i-1].decays_part[j][k] == particle[ii].monval)
                  {
                     idx = ii;
                     break;
                  }
               }
               if(idx == local_i-1 && particle[local_i-1].stable == 0)  // check
               {
                  cout << "Error: can not find decay particle index for anti-baryon!" << endl;
                  cout << "particle monval : " << particle[local_i-1].decays_part[j][k] << endl;
                  exit(1);
               }
               if(particle[idx].baryon == 0 && particle[idx].charge == 0 && particle[idx].strange == 0)
                  particle[local_i].decays_part[j][k]= particle[local_i-1].decays_part[j][k];
               else
                  particle[local_i].decays_part[j][k]= -particle[local_i-1].decays_part[j][k];
            }
         }
       }
       local_i++;	// Add one to the counting variable "i" for the meson/baryon
   }
   resofile.close();
   Nparticle=local_i-1; //take account the final fake one
   for(int i=0; i < Nparticle; i++)
   {
      if(particle[i].baryon==0)
         particle[i].sing=-1;
      else
         particle[i].sing=1;
   }
   cout << "done! Antiparticles are added!" << endl;
   return(Nparticle);
}

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
}
