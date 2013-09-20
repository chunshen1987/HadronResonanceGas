//===============================================================================
//  calculate the HBT radii from VISH2+1
//
//
//  Programmer: Chun Shen
//       Email: shen.201@asc.ohio-state.edu
//
//        Date: 11/22/11
//  The HBT part of the code is reorganized in the class style for easy 
//  maintainance and furture extension.  --- 04/15/2012
//
//===============================================================================


#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>


#include<gsl/gsl_sf_bessel.h>
#include "parameters.h"
#include "readindata.h"
#include "Stopwatch.h"

using namespace std;

int main()
{
   Stopwatch sw;
   sw.tic();
   
   double Temperature;
   cin >> Temperature;
   
   //read particle resonance decay table
   particle_info *particle = new particle_info [Maxparticle];
   int Nparticle=read_resonance(particle);
   cout <<"read in total " << Nparticle << " particles!" << endl;

   // read in stable particle table
   int Nstable_particle;
   int Idummy;
   char cdummy[256];
   ifstream particletable("EOS/EOS_particletable.dat");
   particletable >> Nstable_particle;
   double *stable_particle_monval = new double [Nstable_particle];
   double *particle_mu = new double [Nstable_particle];
   for(int i=0; i<Nstable_particle; i++)
   {
       particletable >> Idummy >> stable_particle_monval[i];
       particletable.getline(cdummy, 256);
   }
   particletable.close();
   cout << "read in data finished!" << endl;

   read_decdat_mu(Nstable_particle, particle_mu);
   calculate_particle_mu(Nparticle, particle, particle_mu);
   
   calculate_particle_yield(Nparticle, particle, Temperature);

   decay_contribution_table* decayTb = new decay_contribution_table [Nstable_particle];
   
   //calculate decay contributions for each speices stable particle
   for(int m=0; m<Nstable_particle; m++) 
   {
      decayTb[m].daughter_monval = stable_particle_monval[m];
      int ii = get_particle_idx(particle, Nparticle, stable_particle_monval[m]);
      decayTb[m].daughter_name = particle[ii].name;

      for(int i=0; i<Nparticle; i++)
         particle[i].decay_contribution = 0.0e0;

      for(int i=0; i<Nparticle; i++) // loop over entire particle list
      {
         if(particle[i].decays_Npart[0]!=1) //determine whether particle[i] will decay or not; particle[i].decays_Npart[0] = 1 for stable particles
         {
            double direct_decay = 0.0;
            double indirect_decay = 0.0;
            bool decay_flag = false;
            for(int j=0; j < particle[i].decays; j++) // loop over each decay channel for particle[i]
            {
               bool direct_flag = false;
               bool indirect_flag = false;
               int direct_Npart = 0;
               double indirect_temp = 0.0;
               for(int k=0; k < abs(particle[i].decays_Npart[j]); k++) // loop over each daughter particle for given decay channel with index j 
               {
                  // direct decay contribution
                  if(particle[i].decays_part[j][k] == stable_particle_monval[m])
                  {
                     direct_Npart++;
                     direct_flag = true;
                  }
                  else // take account for secondary decay contributions
                  {
                     int idx = get_particle_idx(particle, Nparticle, particle[i].decays_part[j][k]);
                     if(particle[idx].decays_Npart[0]!=1)
                     {
                        if(fabs(particle[idx].decay_contribution)>1e-12)
                        {
                           indirect_temp += particle[idx].decay_contribution;
                           indirect_flag = true;
                        }
                     }
                  }
               }
               if(direct_flag)
               {
                  direct_decay += particle[i].decays_branchratio[j]*direct_Npart*particle[i].yield;
                  decay_flag = true;
               }
               if(indirect_flag)
               {
                  indirect_decay += particle[i].decays_branchratio[j]*indirect_temp*particle[i].yield;
                  decay_flag = true;
               }
            }
            if(decay_flag)
            {
                double total_decay = direct_decay + indirect_decay;
                particle[i].decay_contribution = total_decay/particle[i].yield;

                decayTb[m].parent_monval.push_back(particle[i].monval);
                decayTb[m].parent_mass.push_back(particle[i].mass);
                decayTb[m].decay_contribution.push_back(total_decay);
            }
            else
                particle[i].decay_contribution = 0.0;
         }
      }
   }
   
   for(int l=0; l<Nstable_particle; l++)
   {
      ostringstream output_stream;
      output_stream << "reso_" << decayTb[l].daughter_monval << "_dNdy.dat";
      ofstream output(output_stream.str().c_str());
      cout << decayTb[l].daughter_name << endl;
      for(int i=0; i<decayTb[l].parent_monval.size(); i++)
      {
         string name;
         for(int ii = 0; ii<Nparticle; ii++)
            if(particle[ii].monval == decayTb[l].parent_monval[i])
            {
               name = particle[ii].name;
               break;
            }
         output << decayTb[l].parent_monval[i] << "   " << decayTb[l].parent_mass[i] << "   " << name << "   " << decayTb[l].decay_contribution[i] << endl;
      }
   }

   sw.toc();
   cout << "Program totally finished in " << sw.takeTime() << " sec." << endl;
   return 0;
}
