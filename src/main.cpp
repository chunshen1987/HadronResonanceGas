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

    //particleList hadronList("EOS/pdg.dat");
    particleList hadronList("EOS/pdg-urqmd_v3.3+.dat");

    //hadronList.calculateSystemEOS(0., 0., 0.);
    //hadronList.calculateSystemEOS2D(0., 0.);
    //hadronList.calculateSystemEOS2DNS(0., 0.);
    //hadronList.calculateSystemEOS2DNSNQ(0., 0.4);
    //hadronList.calculateSystemEOS3DNS(0.);

    // read in external tables
    string filename = "Latin_HPC_points4D.txt";
    std::ifstream filein(filename.c_str());
    if (!filein.good()) {
        std::cout << "Error: can not find the file: " << filename << std::endl;
        exit(1);
    }
    std::vector<std::vector<double>> points;
    double T, muB, muQ, muS;
    filein >> T >> muB >> muQ >> muS;
    while(!filein.eof()) {
        std::vector<double> line;
        line.push_back(T);
        line.push_back(muB);
        line.push_back(muQ);
        line.push_back(muS);
        points.push_back(line);
        filein >> T >> muB >> muQ >> muS;
    }
    std::cout << "size = " << points.size() << std::endl;
    hadronList.calculateSystemEOSfromTable(points);

    sw.toc();
    std::cout << "Program totally finished in " << sw.takeTime() << " sec."
              << std::endl;
    return 0;
}
