//===============================================================================
//  calculate the thermodynamic quantities of hadron resonance gas
//
//
//  Programmer: Chun Shen
//       Email: shen.201@asc.ohio-state.edu
//
//  Date: 04/22/13
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
#include "particle.h"

using namespace std;

int main()
{
   particle testparticle(211, "pion", 0.14, 0.0, 1, 0, 0, 0, 0, 1, 1, 1, -1, 2);
   int* decayPart = new int [2];
   decayPart[0] = 111; decayPart[1] = 21;
   testparticle.addResonancedecays(0.5, 2, decayPart);
   delete [] decayPart;
   int* decayPart1 = new int [3];
   decayPart1[0] = 111; decayPart1[1] = 21; decayPart1[2] = 32;
   testparticle.addResonancedecays(0.5, 3, decayPart1);
   delete [] decayPart1;

   exit(0);

   Stopwatch sw;
   sw.tic();

   double Temperature = 0.12;
   
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

   calculate_particle_mu(Nparticle, particle, particle_mu);
   
   calculate_particle_yield(Nparticle, particle, Temperature);

   for(int i=0; i<Nparticle; i++)
   {
      cout << particle[i].yield << endl;
   }

   sw.toc();
   cout << "Program totally finished in " << sw.takeTime() << " sec." << endl;
   return 0;
}
