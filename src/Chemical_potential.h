#ifndef CHEMICAL_POTENTIAL_H
#define CHEMICAL_POTENTIAL_H

#include "Table2D.h"

class Chemical_potential
{
   private: 
      int Tb_length, Nstable;
      Table2D* EOS_Mu_Table_ptr;
      double* T;
      double T_chem;
      int nT;
      double** mu_table;
      int read_in_flag;
      int *stable_particle_monval_list;
      
   public:
      Chemical_potential();
      ~Chemical_potential();

      void readin_chempotential_table(string filename);
      void readin_stable_particle_list(string filename);

      int get_Tblength() {return(Tb_length);};
      int get_Nstable() {return(Nstable);};
      int get_stable_particle_monval(int i) {return(stable_particle_monval_list[i]);};
      double get_T(int i) {return(T[i]);};

      void get_stable_mu_table(int iT, double* mu);
      void output_stable_mu(double Temperature, double* mu); 

};


#endif
