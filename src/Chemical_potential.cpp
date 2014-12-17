#include <iostream>
#include <fstream>
#include "Chemical_potential.h"
#include "Arsenal.h"

using namespace std;

Chemical_potential::Chemical_potential()
{
    read_in_flag = 0;
}

Chemical_potential::~Chemical_potential()
{
    delete [] T;
    for(int i = 0; i < Nstable; i++)
        delete [] mu_table[i];
    delete [] mu_table;
    if(read_in_flag == 99)
        delete EOS_Mu_Table_ptr;
}

void Chemical_potential::readin_stable_particle_list(string filename)
{
    read_in_flag = 1;
    ifstream in_file(filename.c_str());
    in_file >> Nstable;
    stable_particle_monval_list = new int [Nstable];
    int dummy;
    string temp;
    for(int i = 0; i < Nstable; i++)
    {
        in_file >> dummy >> stable_particle_monval_list[i];
        getline(in_file, temp);
    }
    in_file.close();

    T_chem = 0.165;
    double T_min = 0.05;
    nT = 200;
    double dT = (T_chem - T_min)/(nT - 1);
    T = new double [nT];
    for(int i = 0; i < nT; i++)
        T[i] = T_chem - i*dT;

    mu_table = new double* [Nstable];
    for(int i = 0; i < Nstable; i++)
    {
        mu_table[i] = new double [nT];
        for(int j = 0; j < nT; j++)
            mu_table[i][j] = 0.0;
    }
}

void Chemical_potential::readin_chempotential_table(string filename)
{
    read_in_flag = 99;
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

void Chemical_potential::get_stable_mu_table(int iT, double* mu)
{
    for(int i = 0; i < Nstable; i++)
        mu[i] = mu_table[i][iT];
    return;
}

void Chemical_potential::output_stable_mu(double Temperature, double* mu)
{
    for(int i = 0; i < Nstable; i++)
        interpolation1D_linear(T, mu_table[i], &Temperature, &mu[i], Tb_length);
}
