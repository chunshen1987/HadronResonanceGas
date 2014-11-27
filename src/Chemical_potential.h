#ifndef CHEMICAL_POTENTIAL_H
#define CHEMICAL_POTENTIAL_H

#include "Table2D.h"

class Chemical_potential
{
   private: 
      int Tb_length, Nstable;
      Table2D* EOS_Mu_Table_ptr;
      double* T;
      double** mu_table;
      
   public:
      Chemical_potential();
      ~Chemical_potential();

      void readin_chempotential_table(string filename);
      void Set_chemical_potential();

      int get_Tblength() {return(Tb_length);};
      int get_Nstable() {return(Nstable);};
      double get_T(int i) {return(T[i]);};

      void output_stable_mu(double Temperature, double* mu); 

};


#endif
