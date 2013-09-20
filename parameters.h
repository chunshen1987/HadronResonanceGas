#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
using namespace std;

const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

const string path = "results";

const double tol = 1e-15;  //tolarence
const int flagneg = 1;     //neglect all points that are negative

#endif
