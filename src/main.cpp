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
   chem_table.readin_chempotential_table("chemical_potential_tb/s95p/s95p-PCE165-v0/s95p-v0-PCE165_chemvsT.dat");
   hadronList.output_particle_chemical_potentials(&chem_table);
   cout << "done!" << endl;

   hadronList.calculateSystemEOS();
   hadronList.calculateSystemEOS(0.02);
   hadronList.calculateSystemEOS(0.4);
   
   sw.toc();
   cout << "Program totally finished in " << sw.takeTime() << " sec." << endl;
   return 0;
}
