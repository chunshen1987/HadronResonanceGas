#include <iostream>
#include "Chemical_potential.h"
#include "Arsenal.h"

using namespace std;

Chemical_potential::Chemical_potential()
{

}

Chemical_potential::~Chemical_potential()
{
    delete [] T;
    for(int i = 0; i < Nstable; i++)
        delete [] mu_table[i];
    delete [] mu_table;
    delete EOS_Mu_Table_ptr;
}

void Chemical_potential::readin_chempotential_table(string filename)
{
    EOS_Mu_Table_ptr = new Table2D(filename);
    Tb_length = EOS_Mu_Table_ptr->getTbsizeX();
    Nstable = EOS_Mu_Table_ptr->getTbsizeY() - 1;
    T = new double [Tb_length];
    mu_table = new double* [Nstable];
    for(int i = 0; i < Nstable; i++)
        mu_table[i] = new double [Tb_length];
    for(int i = 0; i < Tb_length; i++)
    {
        T[i] = EOS_Mu_Table_ptr->getTbdata(i, 0);
        for(int j = 0; j < Nstable; j++)
            mu_table[j][i] = EOS_Mu_Table_ptr->getTbdata(i, j+1);
    }
}


void Chemical_potential::output_stable_mu(double Temperature, double* mu)
{
    for(int i = 0; i < Nstable; i++)
        interpolation1D_linear(T, mu_table[i], &Temperature, &mu[i], Tb_length);
}
