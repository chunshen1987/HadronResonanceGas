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
   readParticlelistTable(particleListFilename);
}

particleList::~particleList()
{

}

void particleList::readParticlelistTable(string tableName)
//read in particle information from pdg data file
{
   cout << "Reading in particle resonance decay table...";
   ifstream resofile(tableName.c_str());
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

void particleList::output_particle_chemical_potentials(Chemical_potential* mu_tb)
{
   ofstream particle_mu_table("chemical_potentials.dat");
   //output names
   particle_mu_table << scientific << setw(18) << setprecision(8)
                     << 0.0 << "    ";
   for(int i = 0; i < partList.size(); i++)
   {
      particle_mu_table << scientific << setw(18) << setprecision(8)
                        << partList[i]->getName() << "    ";
   }
   particle_mu_table << endl;
   //output mass
   particle_mu_table << scientific << setw(18) << setprecision(8)
                     << 0.0 << "    ";
   for(int i = 0; i < partList.size(); i++)
   {
      particle_mu_table << scientific << setw(18) << setprecision(8)
                        << partList[i]->getMass() << "    ";
   }
   particle_mu_table << endl;
   //output baryon number
   particle_mu_table << scientific << setw(18) << setprecision(8)
                     << 0.0 << "    ";
   for(int i = 0; i < partList.size(); i++)
   {
      particle_mu_table << scientific << setw(18) << setprecision(8)
                        << partList[i]->getBaryon() << "    ";
   }
   particle_mu_table << endl;

   double Ti = 0.1;
   double Tf = 0.2;
   double dT = 0.001;
   int nT = (Tf - Ti)/dT + 1;
   for(int i = 0 ; i < nT; i++)
   {
      double T_local = Ti + i*dT;
      calculate_particle_chemical_potential(T_local, mu_tb);
      particle_mu_table << scientific << setw(18) << setprecision(8)
                        << T_local << "    " ;
      for (int j = 0; j < partList.size(); j++)
          particle_mu_table << scientific << setw(18) << setprecision(8)
                            << partList[j]->getMu() << "    ";
      particle_mu_table << endl;
   }
   particle_mu_table.close();
}

void particleList::calculate_particle_chemical_potential(double Temperature, Chemical_potential* mu_tb)
{
   int N_mu = mu_tb->get_Nstable();
   double *mu_stable = new double [N_mu];
   mu_tb->output_stable_mu(Temperature, mu_stable);
   
   int Nstable_particle;
   int Idummy;
   char cdummy[256];
   ifstream particletable("EOS/EOS_particletable.dat");
   particletable >> Nstable_particle;
   if(N_mu != Nstable_particle)
   {
      cout << "chemical potential table is not compatible with EOS_particletable.dat" << endl;
      exit(1);
   }
   double *stable_particle_monval = new double [Nstable_particle];
   for(int i=0; i<Nstable_particle; i++)
   {
       particletable >> Idummy >> stable_particle_monval[i];
       particletable.getline(cdummy, 256);
   }
   particletable.close();
      
   for(int i=0; i<Nstable_particle; i++)
      for(int j=0; j<partList.size(); j++)
         if(partList[j]->getMonval() == stable_particle_monval[i])
         {
            partList[j]->setStable(1);
            partList[j]->setMu(mu_stable[i]);
            break;
         }

   for(int i=0; i < partList.size() ; i++)
   {
      double mu_temp = 0.0;
      if(partList[i]->getStable() == 0)
      {
         for(int j=0; j < partList[i]->getNdecayChannel(); j++)
         {
            for(int k=0; k < abs(partList[i]->getdecaysNpart(j)); k++)
            {
               for(int l=0; l < partList.size(); l++)
               {
                  if(partList[i]->getdecays_part(j,k) == partList[l]->getMonval())
                  {
                     mu_temp += partList[i]->getdecays_branchratio(j)*partList[l]->getMu();
                     break;
                  }
                  if(l == (partList.size() - 1) && partList[i]->getdecays_part(j,k) != 22)
                     cout<<"warning: can not find particle " <<  partList[i]->getdecays_part(j,k) << endl;
               }
            }
         }
         partList[i]->setMu(mu_temp); 
      }
   }

   //for (int i = 0; i < partList.size(); i++)
   //   cout << partList[i]->getMonval() << "   " << partList[i]->getMu() << endl;

   delete [] stable_particle_monval;
   delete [] mu_stable;
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
      partList[i]->calculateParticleYield(Temperature);
   return;
}

void particleList::calculateSystemenergyDensity(double Temperature, double mu_B, double mu_S)
//calculate the energy density of the system at given T and mu
{
   double result = 0.0e0;
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculateEnergydensity(Temperature);
   edSystem = result;
   return;
}

void particleList::calculateSystemPressure(double Temperature, double mu_B, double mu_S)
//calculate the pressure of the system at given T and mu
{
   double result = 0.0e0;
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculatePressure(Temperature);
   pressureSys = result;
   return;
}

void particleList::calculateSystementropyDensity(double Temperature, double mu_B, double mu_S)
//calculate the entropy density of the system at given T and mu
{
   double result = 0.0e0;
   for(int i = 0; i < partList.size(); i++)
      result += partList[i]->calculateEntropydensity(Temperature);
   sdSystem = result;
   return;
}

void particleList::calculateSystemEOS(double mu_B, double mu_S)
//calculate the EOS of given system, e,p,s as functions of T at given mu_B and mu_S
{ 
   cout << "calculate the EOS of the system with muB = " << mu_B << " GeV and mu_S = " << mu_S << " GeV .... " << endl;
   int nT = 200;
   double T_i = 0.01;        // unit: (GeV)
   double T_f = 0.2;          // unit: (GeV)
   double dT = (T_f - T_i)/(nT - 1);
   double* temp_ptr = new double [nT];
   double* ed_ptr = new double [nT];
   double* sd_ptr = new double [nT];
   double* pressure_ptr = new double [nT];
   double* cs2_ptr = new double [nT];
   for(int i = 0; i < nT; i++)
   {
      temp_ptr[i] = T_i + i*dT;
      calculate_particle_yield(temp_ptr[i], mu_B, mu_S);
      calculateSystemenergyDensity(temp_ptr[i]);
      calculateSystemPressure(temp_ptr[i]);
      calculateSystementropyDensity(temp_ptr[i]);
      ed_ptr[i] = edSystem;
      sd_ptr[i] = sdSystem;
      pressure_ptr[i] = pressureSys;
   }
   //calculate speed of sound cs^2 = dP/de
   for(int i = 0; i < nT - 1; i++)
      cs2_ptr[i] = (pressure_ptr[i+1] - pressure_ptr[i])/(ed_ptr[i+1] - ed_ptr[i]+1e-30);
   cs2_ptr[nT-1] = cs2_ptr[nT-2];

   //output EOS table
   ostringstream EOSfilename;
   EOSfilename << "./EOS_muB_" << mu_B << "_muS_" << mu_S << ".dat";
   ofstream output(EOSfilename.str().c_str());
   for(int i = 0; i < nT; i++)
      output << scientific << setw(20) << setprecision(8)
             << temp_ptr[i] << "   "  << ed_ptr[i] 
             << "   " << sd_ptr[i] << "   " << pressure_ptr[i] 
             << "   " << cs2_ptr[i] << endl;
   output.close();
   delete [] temp_ptr;
   delete [] ed_ptr;
   delete [] sd_ptr;
   delete [] pressure_ptr;
   delete [] cs2_ptr;
   return;
}
