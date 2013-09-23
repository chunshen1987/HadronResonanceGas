#include<iostream>

#include "particle.h"
using namespace std;

particle::particle(int monval_in, string name_in, double mass_in, double width_in, int gspin_in, int baryon_in, int strange_in, int charm_in, int bottom_in, int gisospin_in, int charge_in, int NdecayChannel_in)
{
    monval = monval_in;
    name = name_in;
    mass = mass_in;
    width = width_in;
    gspin = gspin_in;
    baryon = baryon_in;
    strange = strange_in;
    charm = charm_in;
    bottom = bottom_in;
    gisospin = gisospin_in;
    charge = charge_in;
    NdecayChannel = NdecayChannel_in;
    stable = 0;

    //determine Bose/Fermi statistic for particle
    if(baryon == 0)
       sign = -1;
    else
       sign = 1;

    decays_branchratio = new double [NdecayChannel];
    decays_Npart = new int [NdecayChannel];
    decays_part = new int* [NdecayChannel];
    
    channelIdx = 0;
}

particle::~particle()
{  
    delete [] decays_branchratio;
    delete [] decays_Npart;
    for(int i = 0; i < NdecayChannel; i++)
       delete [] decays_part[i];
    delete [] decays_part;
}

void particle::addResonancedecays(double branchratio, int Npart, int* decayChannelparts)
//add resonance decay channel for particle
{
    if(channelIdx > NdecayChannel-1)
    {
        cout << "Warning: channelidx exceed number of decay channels! Please check" << endl;
        exit(1);
    }
    decays_branchratio[channelIdx] = branchratio;
    decays_Npart[channelIdx] = Npart;
    decays_part[channelIdx] = new int [Npart];
    for(int i = 0; i < Npart; i++)
       decays_part[channelIdx][i] = decayChannelparts[i];
    channelIdx++;
    if(NdecayChannel == 1 && decays_Npart[0] == 1) stable = 1;

    return;
}

int particle::getAntiparticleMonval()
//return the corresponding anti-particle' Monte-Carlo value
{
   if(baryon == 0 && charge == 0 && strange == 0)
      return(monval);
   else
      return(-monval);
}
