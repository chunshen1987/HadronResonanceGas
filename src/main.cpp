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
#include "Stopwatch.h"
#include "particle.h"
#include "particleList.h"
#include "Chemical_potential.h"

using namespace std;

int main()
{
   Stopwatch sw;
   sw.tic();

   particleList hadronList("EOS/pdg.dat");
   Chemical_potential chem_table;
   chem_table.readin_chempotential_table(
                  "chemical_potential_tb/s95p/s95p-v1/s95p-v1-CE_chemvsT.dat");
   //chem_table.readin_stable_particle_list("EOS/EOS_particletable.dat");
   //hadronList.calculate_particle_decay_probability(&chem_table);
   //hadronList.output_particle_chemical_potentials(&chem_table);

   hadronList.calculateSystemEOS_and_output_in_2D();
   //hadronList.calculateSystemEOS(0.02);
   //hadronList.calculateSystemEOS(0.4);
   
   sw.toc();
   cout << "Program totally finished in " << sw.takeTime() << " sec." << endl;
   return 0;
}
