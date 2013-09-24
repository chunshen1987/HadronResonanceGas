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

using namespace std;

int main()
{
   Stopwatch sw;
   sw.tic();

   particleList hadronList("EOS/pdg.dat");
   hadronList.calculateSystemEOS();
   
   sw.toc();
   cout << "Program totally finished in " << sw.takeTime() << " sec." << endl;
   return 0;
}
