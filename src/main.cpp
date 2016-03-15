//==============================================================================
//  calculate the thermodynamic quantities of hadron resonance gas
//
//  Copyright 2016 Chun Shen
//  Email: shen.201@asc.ohio-state.edu
//==============================================================================

#include<sys/time.h>
#include<gsl/gsl_sf_bessel.h>

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>

#include "./Stopwatch.h"
#include "./particle.h"
#include "./particleList.h"
#include "./Chemical_potential.h"

int main() {
    Stopwatch sw;
    sw.tic();

    particleList hadronList("EOS/pdg.dat");
    hadronList.calculate_particle_yield(0.125, 0.3005, 0.0);

    sw.toc();
    std::cout << "Program totally finished in " << sw.takeTime() << " sec."
              << std::endl;
    return 0;
}
