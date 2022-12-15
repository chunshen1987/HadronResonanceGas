// Copyright 2016 Chun Shen

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "particleList.h"

using namespace std;

particleList::particleList(string particleTableName) {
    particleListFilename = particleTableName;  // filename of pdg data file
    readParticlelistTable(particleListFilename);
    sort_particle_list_in_particle_mass();
    calculate_particle_decay_probability();
}

particleList::~particleList() {
    partList.clear();
    stable_particle_list.clear();
}

void particleList::readParticlelistTable(string tableName) {
// read in particle information from pdg data file
    cout << "Reading in particle resonance decay table...";
    ifstream resofile(tableName.c_str());
    int monval;
    string name;
    double mass, width;
    int gspin, gisospin;
    int baryon, strange, charm, bottom, charge;
    int decays, decayNpart;
    double decayBranchratio;
    int decayPart[5] = {0, 0, 0, 0, 0};

    int dummy_int;
    while (1) {
        resofile >> monval;
        if (resofile.eof())
            break;
        resofile >> name;
        resofile >> mass;
        resofile >> width;
        resofile >> gspin;
        resofile >> baryon;
        resofile >> strange;
        resofile >> charm;
        resofile >> bottom;
        resofile >> gisospin;
        resofile >> charge;
        resofile >> decays;
        partList.push_back(
            new particle(monval, name, mass, width, gspin, baryon, strange,
                         charm, bottom, gisospin, charge, decays));
        if (baryon == 1) {
            ostringstream antiname;
            antiname << "Anti-" << name;
            partList.push_back(
                new particle(-monval, antiname.str(), mass, width, gspin,
                             -baryon, -strange, -charm, -bottom, gisospin,
                             -charge, decays));
        }
        for (int j = 0; j < decays; j++) {
            resofile >> dummy_int;
            resofile >> decayNpart;
            resofile >> decayBranchratio;
            resofile >> decayPart[0];
            resofile >> decayPart[1];
            resofile >> decayPart[2];
            resofile >> decayPart[3];
            resofile >> decayPart[4];
            decayNpart = abs(decayNpart);

            int* tempptr = new int[decayNpart];
            for (int ipart = 0; ipart < decayNpart; ipart++)
                tempptr[ipart] = decayPart[ipart];

            if (baryon == 0) {
                partList.back()->addResonancedecays(decayBranchratio,
                                                    decayNpart, tempptr);
            } else {
                partList.at(partList.size()-2)->addResonancedecays(
                                        decayBranchratio, decayNpart, tempptr);
                for (int ipart = 0; ipart < decayNpart; ipart++) {
                    int particleId = get_particle_idx(tempptr[ipart]);
                    tempptr[ipart] =
                            partList[particleId]->getAntiparticleMonval();
                }
                partList.back()->addResonancedecays(decayBranchratio,
                                                    decayNpart, tempptr);
            }
            delete [] tempptr;
        }
    }
    resofile.close();
    partList.erase(partList.begin());  // delete gamma
    cout << "done! Antiparticles are added!" << endl;
    cout << "There are totally " << partList.size() << " particles." << endl;
    return;
}

void particleList::sort_particle_list_in_particle_mass() {
    for (int i = 1; i < partList.size(); i++) {
        int k = i;
        int j = i - 1;
        while (partList[k]->getMass() < partList[j]->getMass() && j >= 0) {
            particle* temp = partList[j];
            partList[j] = partList[k];
            partList[k] = temp;
            k--;
            j--;
        }
    }
    // for (int i = 0; i < partList.size(); i++) {
    //     cout << partList[i]->getMass() << endl;
    // }
    // exit(0);
}

void particleList::calculate_particle_decay_probability() {
    for (int i = 0; i < partList.size(); i++) {
        if (partList[i]->getStable() == 1) {
            stable_particle_list.push_back(partList[i]);
        }
    }
    int Nstable = stable_particle_list.size();
    // cout << "there are " << Nstable << " stable particles." << endl;
    // for (int i = 0; i < Nstable; i++) {
    //     cout << "stable particle name: " << stable_particle_list[i]->getName()
    //          << endl;
    // }

    for (int i = 0; i < partList.size(); i++) {
        partList[i]->create_decay_probability_array(Nstable);
        for (int j = 0; j < Nstable; j++) {
            if (partList[i]->getMonval()
                == stable_particle_list[j]->getMonval()) {
                partList[i]->setStable(1);
                partList[i]->set_decay_probability(j, 1.0);
                break;
            }
        }
    }

    for (int i = 0; i < partList.size(); i++) {
        if (partList[i]->getStable() == 0) {  // unstable resonances
            double* temp = new double[Nstable];
            for (int j = 0; j < Nstable; j++)
                temp[j] = 0.0;
            for (int j = 0; j < partList[i]->getNdecayChannel(); j++) {
                for (int k = 0; k < abs(partList[i]->getdecaysNpart(j)); k++) {
                    if (partList[i]->getdecays_part(j, k) == 22)
                        continue;
                    for (int l = 0; l < partList.size(); l++) {
                        if (partList[i]->getdecays_part(j, k)
                                   == partList[l]->getMonval()) {
                            for (int m = 0; m < Nstable; m++)
                                temp[m] += (
                                    partList[i]->getdecays_branchratio(j)
                                    *partList[l]->getdecay_probability(m));
                            break;
                        }
                        if (l > i) {
                                cout << "decay: " << partList[l]->getMonval()
                                     << "   " << partList[i]->getMonval()
                                     << "   "
                                     << partList[i]->getdecays_part(j, k)
                                     << endl;
                        }
                    }
                }
            }
            for (int j = 0; j < Nstable; j++)
                partList[i]->set_decay_probability(j, temp[j]);
            delete [] temp;
        }
    }
}

void particleList::calculate_particle_decays() {
    for (int j = 0; j < stable_particle_list.size(); j++) {
        double stable_yield = 0.;
        for (int k = 0; k < partList.size(); k++) {
            stable_yield += (partList[k]->getdecay_probability(j)
                             *partList[k]->getParticleYield());
        }
        stable_particle_list[j]->set_particle_stable_yield(stable_yield);
    }
    //for (int i = 0; i < stable_particle_list.size(); i++) {
    //    cout << "stable particle name: " << stable_particle_list[i]->getName()
    //         << ", yield: "
    //         << stable_particle_list[i]->get_particle_stable_yield()
    //         << ", thermal yield:" 
    //         << stable_particle_list[i]-> getParticleYield()
    //         << endl;
    //}
}

void particleList::output_particle_chemical_potentials(
                                                Chemical_potential* mu_tb) {
    ofstream particle_mu_table("chemical_potentials.dat");
    // output names
    particle_mu_table << scientific << setw(18) << setprecision(8)
                      << 0.0 << "    ";
    for (int i = 0; i < partList.size(); i++) {
        particle_mu_table << scientific << setw(18) << setprecision(8)
                         << partList[i]->getName() << "    ";
    }
    particle_mu_table << endl;
    // output mass
    particle_mu_table << scientific << setw(18) << setprecision(8)
                      << 0.0 << "    ";
    for (int i = 0; i < partList.size(); i++) {
        particle_mu_table << scientific << setw(18) << setprecision(8)
                          << partList[i]->getMass() << "    ";
    }
    particle_mu_table << endl;
    // output baryon number
    particle_mu_table << scientific << setw(18) << setprecision(8)
                      << 0.0 << "    ";
    for (int i = 0; i < partList.size(); i++) {
        particle_mu_table << scientific << setw(18) << setprecision(8)
                          << partList[i]->getBaryon() << "    ";
    }
    particle_mu_table << endl;

    double Ti = 0.1;
    double Tf = 0.2;
    double dT = 0.001;
    int nT = (Tf - Ti)/dT + 1;
    for (int i = 0 ; i < nT; i++) {
        double T_local = Ti + i*dT;
        calculate_particle_chemical_potential(T_local, mu_tb);
        // calculate_particle_chemical_potential2(T_local, mu_tb);
        particle_mu_table << scientific << setw(18) << setprecision(8)
                          << T_local << "    ";
        for (int j = 0; j < partList.size(); j++)
            particle_mu_table << scientific << setw(18) << setprecision(8)
                              << partList[j]->getMu() << "    ";
        particle_mu_table << endl;
    }
    particle_mu_table.close();
}

void particleList::calculate_particle_chemical_potential2(
                            double Temperature, Chemical_potential* mu_tb) {
    int Nstable = mu_tb->get_Nstable();
    double *mu_stable = new double[Nstable];
    mu_tb->output_stable_mu(Temperature, mu_stable);
    for (int i = 0; i < partList.size() ; i++) {
        double mu_temp = 0.0;
        for (int j = 0; j < Nstable; j++)
            mu_temp += partList[i]->getdecay_probability(j)*mu_stable[j];
        partList[i]->setMu(mu_temp);
    }
    delete [] mu_stable;
}

void particleList::calculate_particle_chemical_potential(
                        double Temperature, Chemical_potential* mu_tb) {
    int N_mu = mu_tb->get_Nstable();
    double *mu_stable = new double[N_mu];
    mu_tb->output_stable_mu(Temperature, mu_stable);

    int Nstable_particle;
    int Idummy;
    char cdummy[256];
    ifstream particletable("EOS/EOS_particletable.dat");
    particletable >> Nstable_particle;
    if (N_mu != Nstable_particle) {
       cout << "chemical potential table is not compatible with "
            << "EOS_particletable.dat" << endl;
       exit(1);
    }
    double *stable_particle_monval = new double[Nstable_particle];
    for (int i = 0; i < Nstable_particle; i++) {
        particletable >> Idummy >> stable_particle_monval[i];
        particletable.getline(cdummy, 256);
    }
    particletable.close();

    for (int i = 0; i < Nstable_particle; i++)
        for (int j = 0; j < partList.size(); j++)
            if (partList[j]->getMonval() == stable_particle_monval[i]) {
                partList[j]->setStable(1);
                partList[j]->setMu(mu_stable[i]);
                break;
            }
    for (int i = 0; i < partList.size() ; i++) {
        double mu_temp = 0.0;
        if (partList[i]->getStable() == 0) {
            for (int j=0; j < partList[i]->getNdecayChannel(); j++) {
                for (int k=0; k < abs(partList[i]->getdecaysNpart(j)); k++) {
                    for (int l=0; l < partList.size(); l++) {
                        if (partList[i]->getdecays_part(j, k)
                            == partList[l]->getMonval()) {
                            mu_temp += (
                                    partList[i]->getdecays_branchratio(j)
                                    *partList[l]->getMu());
                            break;
                        }
                    if (l == (partList.size() - 1)
                                   && partList[i]->getdecays_part(j, k) != 22)
                        cout << "warning: can not find particle "
                             << partList[i]->getdecays_part(j, k) << endl;
                    }
                }
            }
            partList[i]->setMu(mu_temp);
        }
    }

    delete [] stable_particle_monval;
    delete [] mu_stable;
}

int particleList::get_particle_idx(int particle_monval) {
// return the idx in particleList for given particle Monte-Carlo number
    for (int i = 0; i < partList.size(); i++)
        if (partList[i]->getMonval() == particle_monval)
            return(i);
    cout << "Warning: can not find particle " << particle_monval << endl;
    exit(1);
}


//! calculate particle chemical potentials
//! need to add support for partial chemical equilibrium
void particleList::calculate_particle_mu(double mu_B, double mu_S,
                                         double mu_Q) {
    for (int i = 0; i < partList.size(); i++)
        partList[i]->calculateChemicalpotential(mu_B, mu_S, mu_Q);
}


//! calculate particle yield
void particleList::calculate_particle_yield(double Temperature, double mu_B,
                                            double mu_S, double mu_Q) {
    //calculate_particle_mu(mu_B, mu_S, mu_Q);
    for (int i = 0; i < partList.size(); i++)
        partList[i]->calculateParticleYield(Temperature, mu_B, mu_S, mu_Q);
    //calculate_particle_decays();
}


//! calculate particle yield when only mu changed
void particleList::calculate_particle_yieldFugacity(
        double Temperature, double mu_B, double mu_S, double mu_Q) {
    for (int i = 0; i < partList.size(); i++) {
        partList[i]->calculateParticleYieldFugacity(Temperature, mu_B,
                                                    mu_S, mu_Q);
    }
}


double particleList::calculateSystemNetbaryonDensity() {
    double result = 0.0;
    for (int i = 0; i < partList.size(); i++) {
        result += partList[i]->getBaryon()*partList[i]->getParticleYield();
    }
    return(result);
}


double particleList::calculateSystemNetStrangenessDensity() {
    double result = 0.0;
    for (int i = 0; i < partList.size(); i++) {
        result += (partList[i]->getStrangeness()
                   *partList[i]->getParticleYield());
    }
    return(result);
}


double particleList::calculateSystemNetElectricChargeDensity() {
    double result = 0.0;
    for (int i = 0; i < partList.size(); i++) {
        result += (partList[i]->getElectricCharge()
                   *partList[i]->getParticleYield());
    }
    return(result);
}


//! calculate the energy density of the system at given T and mu
double particleList::calculateSystemenergyDensity(double Temperature) {
    double result = 0.0e0;
    for (int i = 0; i < partList.size(); i++)
        result += partList[i]->calculateEnergydensity(Temperature);
    return(result);
}


//! calculate the pressure of the system at given T and mu
double particleList::calculateSystemPressure(double Temperature) {
    double result = 0.0e0;
    for (int i = 0; i < partList.size(); i++)
        result += partList[i]->calculatePressure(Temperature);
    return(result);
}


//! calculate the entropy density of the system at given T and mu
double particleList::calculateSystementropyDensity(double Temperature) {
    double result = 0.0e0;
    for (int i = 0; i < partList.size(); i++)
        result += partList[i]->calculateEntropydensity(Temperature);
    return(result);
}


double particleList::findMuSFromNS(double T, double muB, double muQ,
                                   double nS) {
    // the limit was set by Lambda mass = 1.116 GeV
    // muS + muB < m_Lambda and muB_max = 0.8 GeV
    double muS_min = -0.15;
    double muS_max = 0.35;
    calculate_particle_yieldFugacity(T, muB, muS_min, muQ);
    double nS_min = calculateSystemNetStrangenessDensity();
    calculate_particle_yieldFugacity(T, muB, muS_max, muQ);
    double nS_max = calculateSystemNetStrangenessDensity();
    if ((nS - nS_min)*(nS - nS_max) > 0) {
        cout << "No solution for nS = " << nS << " 1/fm^3, with T = "
             << T << " GeV, muB = " << muB << " GeV, muQ = " << muQ
             << " GeV. nS_min(muS=-0.15GeV) = " << nS_min << " 1/fm^3, "
             << "nS_max(muS=0.40GeV) = " << nS_max << " 1/fm^3." << endl;
        exit(1);
    }
    double muS_mid = 0;
    double nS_mid = 0;
    double nB_mid = 0;
    int iter = 0;
    do {
        muS_mid = (muS_min + muS_max)/2.;
        calculate_particle_yieldFugacity(T, muB, muS_mid, muQ);
        nS_mid = calculateSystemNetStrangenessDensity();
        nB_mid = std::abs(calculateSystemNetbaryonDensity());
        if ((nS - nS_mid)*(nS - nS_max) > 0) {
            muS_max = muS_mid;
            nS_max = nS_mid;
        } else {
            muS_min = muS_mid;
            nS_min = nS_mid;
        }
        iter++;
        if (iter > 1000) {
            cout << "check : " << muS_mid << ", nS_mid = " << nS_mid
                 << "muS_min = " << muS_min << ", muS_max = "
                 << muS_max << ", nB = " << nB_mid << endl;
        }
    } while (std::abs(muS_max - muS_min) > 1e-10
             && (std::abs(nS_mid - nS) > 1e-8
                 || std::abs((nS_mid - nS)/(nB_mid + 1e-16)) > 1e-4));
    return(muS_mid);
}


void particleList::findMuQMuQFromNSNQ(double T, double muB, double nS,
                                      double nQovernBRatio,
                                      double &muS, double &muQ) {
    if (std::abs(muB) < 1e-4) {
        muS = 0.;
        muQ = 0.;
        return;
    }
    double muQ_min = -0.1;
    double muQ_max =  0.1;
    double muS_min = findMuSFromNS(T, muB, muQ_min, nS);
    double muS_max = findMuSFromNS(T, muB, muQ_max, nS);
    calculate_particle_yieldFugacity(T, muB, muS_min, muQ_min);
    double nQ_min = calculateSystemNetElectricChargeDensity();
    double nB_min = std::abs(calculateSystemNetbaryonDensity());
    calculate_particle_yieldFugacity(T, muB, muS_max, muQ_max);
    double nQ_max = calculateSystemNetElectricChargeDensity();
    double nB_max = std::abs(calculateSystemNetbaryonDensity());

    double r_min = nQ_min/(nB_min + 1e-16);
    double r_max = nQ_max/(nB_max + 1e-16);
    if ((nQovernBRatio - r_min)*(nQovernBRatio - r_max) > 0) {
        cout << "No solution for nS = " << nS << " 1/fm^3, nQ/nB = "
             << nQovernBRatio << ", with T = "
             << T << " GeV, muB = " << muB << " GeV. "
             << "nQ/nB(muQ=-0.1GeV) = " << r_min
             << ", nQ = " << nQ_min << ", nB = " << nB_min
             << ", nQ/nB(muQ=0.1GeV) = " << r_max
             << ", nQ = " << nQ_max << ", nB = " << nB_max
             << endl;
        exit(1);
    }
    double muQ_mid = 0;
    double muS_mid = 0;
    double r_mid = 0;
    int iter = 0;
    do {
        muQ_mid = (muQ_min + muQ_max)/2.;
        muS_mid = findMuSFromNS(T, muB, muQ_mid, nS);
        calculate_particle_yieldFugacity(T, muB, muS_mid, muQ_mid);
        double nQ_mid = calculateSystemNetElectricChargeDensity();
        double nB_mid = std::abs(calculateSystemNetbaryonDensity());
        r_mid = nQ_mid/(nB_mid + 1e-16);
        if ((nQovernBRatio - r_mid)*(nQovernBRatio - r_max) > 0) {
            muQ_max = muQ_mid;
            r_max = r_mid;
        } else {
            muQ_min = muQ_mid;
            r_min = r_mid;
        }
        iter++;
        if (iter > 1000) {
            cout << "ratio = " << r_mid << ", muQ = " << muQ_mid
                 << ", nQ = " << nQ_mid << ", nB = " << nB_mid << endl;
        }
    } while (std::abs(r_mid - nQovernBRatio) > 1e-4);
    muS = muS_mid;
    muQ = muQ_mid;
}


void particleList::calculateSystemEOS(double mu_B, double mu_S, double mu_Q) {
// calculate the EOS of given system, e,p,s as functions of T
// at given mu_B and mu_S and mu_Q
    cout << "calculate the EOS of the system with muB = " << mu_B
         << " GeV and mu_S = " << mu_S << " GeV .... " << endl;

    int nT = 181;
    double T_i = 0.01;         // unit: (GeV)
    double T_f = 0.19;          // unit: (GeV)
    double dT = (T_f - T_i)/(nT - 1);
    std::vector<double> temp_ptr(nT, 0.);
    std::vector<double> ed_ptr(nT, 0.);
    std::vector<double> sd_ptr(nT, 0.);
    std::vector<double> pressure_ptr(nT, 0.);
    std::vector<double> net_baryon_ptr(nT, 0.);
    std::vector<double> cs2_ptr(nT, 0.);
    for (int i = 0; i < nT; i++) {
        temp_ptr[i] = T_i + i*dT;
        calculate_particle_yield(temp_ptr[i], mu_B, mu_S, mu_Q);
        net_baryon_ptr[i] = calculateSystemNetbaryonDensity();
        ed_ptr[i] = calculateSystemenergyDensity(temp_ptr[i]);
        pressure_ptr[i] = calculateSystemPressure(temp_ptr[i]);
        sd_ptr[i] = calculateSystementropyDensity(temp_ptr[i]);
    }
    // calculate speed of sound cs^2 = dP/de
    for (int i = 0; i < nT - 1; i++)
        cs2_ptr[i] = ((pressure_ptr[i+1] - pressure_ptr[i])
                      /(ed_ptr[i+1] - ed_ptr[i] + 1e-30));
    cs2_ptr[nT-1] = cs2_ptr[nT-2];

    // output EOS table
    ostringstream EOSfilename;
    EOSfilename << "./EOS_muB_" << mu_B << "_muS_" << mu_S << ".dat";
    ofstream output(EOSfilename.str().c_str());
    output << "# T [GeV]  e [GeV/fm^3]  nB [1/fm^3]  s [1/fm^3]  P [GeV/fm^3]"
           << "  cs^2" << endl;
    for (int i = 0; i < nT; i++)
        output << scientific << setw(20) << setprecision(8)
               << temp_ptr[i] << "   "
               << ed_ptr[i] << "   " << net_baryon_ptr[i] << "  "
               << sd_ptr[i] << "   " << pressure_ptr[i] << "   "
               << cs2_ptr[i] << endl;
    output.close();
}


void particleList::calculateSystemEOS2D(double mu_S, double mu_Q) {
    // calculate the EOS of given system, e, p, s as functions of T and mu_B
    // at given mu_S and mu_Q
    cout << "calculate the EOS of the system with mu_S = " << mu_S
         << " GeV and mu_Q = " << mu_Q << " GeV .... " << endl;

    int nT = 181;
    double T_i = 0.01;         // unit: (GeV)
    double T_f = 0.19;          // unit: (GeV)
    double dT = (T_f - T_i)/(nT - 1);

    int nmuB = 801;
    double muB_i = 0.;
    double muB_f = 0.8;
    double dmuB = (muB_f - muB_i)/(nmuB - 1);

    double HBARC = 0.19733;
    double unitFac = HBARC*HBARC*HBARC;

    std::vector<double> temp_ptr(nT, 0.);
    std::vector<double> muB_ptr(nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        temp_ptr[i] = T_i + i*dT;
    }
    for (int i = 0; i < nmuB; i++) {
        muB_ptr[i] = muB_i + i*dmuB;
    }

    std::vector<double> ed_ptr(nT*nmuB, 0.);
    std::vector<double> sd_ptr(nT*nmuB, 0.);
    std::vector<double> pressure_ptr(nT*nmuB, 0.);
    std::vector<double> net_baryon_ptr(nT*nmuB, 0.);
    std::vector<double> net_strangeness_ptr(nT*nmuB, 0.);
    std::vector<double> net_electric_ptr(nT*nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        calculate_particle_yield(temp_ptr[i], 0., 0., 0.);
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            calculate_particle_yieldFugacity(temp_ptr[i], muB_ptr[j],
                                             mu_S, mu_Q);
            net_baryon_ptr[idx] = calculateSystemNetbaryonDensity();
            ed_ptr[idx] = calculateSystemenergyDensity(temp_ptr[i]);
            pressure_ptr[idx] = calculateSystemPressure(temp_ptr[i]);
            sd_ptr[idx] = calculateSystementropyDensity(temp_ptr[i]);
            net_strangeness_ptr[idx] = calculateSystemNetStrangenessDensity();
            net_electric_ptr[idx] = calculateSystemNetElectricChargeDensity();
        }
    }

    // output EOS table
    ostringstream EOSfilename;
    EOSfilename << "./EOS2D_muS_" << mu_S << "_muQ_" << mu_Q << ".dat";
    ofstream output(EOSfilename.str().c_str());
    output << "# T [GeV]  muB [GeV]  e/T^4  nB/T^3  s/T^3  P/T^3  "
           << "nS/T^3  nQ/T^3" << endl;
    for (int i = 0; i < nT; i++) {
        double T3 = pow(temp_ptr[i], 3);
        double T4 = pow(temp_ptr[i], 4);
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            output << scientific << setw(20) << setprecision(8)
                   << temp_ptr[i] << "  " << muB_ptr[j] << "  "
                   << ed_ptr[idx]/T4*unitFac << "  "
                   << net_baryon_ptr[idx]/T3*unitFac << "  "
                   << sd_ptr[idx]/T3*unitFac << "  "
                   << pressure_ptr[idx]/T4*unitFac << "  "
                   << net_strangeness_ptr[idx]/T3*unitFac << "  "
                   << net_electric_ptr[idx]/T3*unitFac << endl;
        }
    }
    output.close();
}


void particleList::calculateSystemEOS2DNS(double n_S, double mu_Q) {
    // calculate the EOS of given system, e, p, s as functions of T and mu_B
    // at given n_S and mu_Q
    cout << "calculate the EOS of the system with n_S = " << n_S
         << " GeV and mu_Q = " << mu_Q << " GeV .... " << endl;

    int nT = 181;
    double T_i = 0.01;         // unit: (GeV)
    double T_f = 0.19;          // unit: (GeV)
    double dT = (T_f - T_i)/(nT - 1);

    int nmuB = 801;
    double muB_i = 0.;
    double muB_f = 0.8;
    double dmuB = (muB_f - muB_i)/(nmuB - 1);

    double HBARC = 0.19733;
    double unitFac = HBARC*HBARC*HBARC;

    std::vector<double> temp_ptr(nT, 0.);
    std::vector<double> muB_ptr(nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        temp_ptr[i] = T_i + i*dT;
    }
    for (int i = 0; i < nmuB; i++) {
        muB_ptr[i] = muB_i + i*dmuB;
    }

    std::vector<double> ed_ptr(nT*nmuB, 0.);
    std::vector<double> sd_ptr(nT*nmuB, 0.);
    std::vector<double> pressure_ptr(nT*nmuB, 0.);
    std::vector<double> net_baryon_ptr(nT*nmuB, 0.);
    std::vector<double> muS_ptr(nT*nmuB, 0.);
    std::vector<double> net_electric_ptr(nT*nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        calculate_particle_yield(temp_ptr[i], 0., 0., 0.);
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            double mu_S = findMuSFromNS(temp_ptr[i], muB_ptr[j], mu_Q, n_S);
            calculate_particle_yieldFugacity(temp_ptr[i], muB_ptr[j],
                                             mu_S, mu_Q);
            net_baryon_ptr[idx] = calculateSystemNetbaryonDensity();
            ed_ptr[idx] = calculateSystemenergyDensity(temp_ptr[i]);
            pressure_ptr[idx] = calculateSystemPressure(temp_ptr[i]);
            sd_ptr[idx] = calculateSystementropyDensity(temp_ptr[i]);
            muS_ptr[idx] = mu_S;
            net_electric_ptr[idx] = calculateSystemNetElectricChargeDensity();
            //double nS = calculateSystemNetStrangenessDensity();
            //cout << "check: T = " << temp_ptr[i] << " GeV, muB = "
            //     << muB_ptr[j] << " GeV, muS = " << muS_ptr[idx] << " GeV, "
            //     << "nS = " << nS << endl;
        }
    }

    // output EOS table
    ostringstream EOSfilename;
    EOSfilename << "./EOS2D_nS_" << n_S << "_muQ_" << mu_Q << ".dat";
    ofstream output(EOSfilename.str().c_str());
    output << "# T [GeV]  muB [GeV]  e/T^4  nB/T^3  s/T^3  P/T^3  "
           << "muS [GeV]  nQ/T^3" << endl;
    for (int i = 0; i < nT; i++) {
        double T3 = pow(temp_ptr[i], 3);
        double T4 = pow(temp_ptr[i], 4);
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            output << scientific << setw(20) << setprecision(8)
                   << temp_ptr[i] << "  " << muB_ptr[j] << "  "
                   << ed_ptr[idx]/T4*unitFac << "  "
                   << net_baryon_ptr[idx]/T3*unitFac << "  "
                   << sd_ptr[idx]/T3*unitFac << "  "
                   << pressure_ptr[idx]/T4*unitFac << "  "
                   << muS_ptr[idx] << "  "
                   << net_electric_ptr[idx]/T3*unitFac << endl;
        }
    }
    output.close();
}


void particleList::calculateSystemEOS2DNSNQ(double n_S, double nQovernB) {
    // calculate the EOS of given system, e, p, s as functions of T and mu_B
    // at given n_S and nQ/nB ratio
    cout << "calculate the EOS of the system with n_S = " << n_S
         << " GeV and nQ/nB = " << nQovernB << " .... " << endl;

    int nT = 181;
    double T_i = 0.01;         // unit: (GeV)
    double T_f = 0.19;          // unit: (GeV)
    double dT = (T_f - T_i)/(nT - 1);

    int nmuB = 801;
    double muB_i = 0.;
    double muB_f = 0.8;
    double dmuB = (muB_f - muB_i)/(nmuB - 1);

    double HBARC = 0.19733;
    double unitFac = HBARC*HBARC*HBARC;

    std::vector<double> temp_ptr(nT, 0.);
    std::vector<double> muB_ptr(nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        temp_ptr[i] = T_i + i*dT;
    }
    for (int i = 0; i < nmuB; i++) {
        muB_ptr[i] = muB_i + i*dmuB;
    }

    std::vector<double> ed_ptr(nT*nmuB, 0.);
    std::vector<double> sd_ptr(nT*nmuB, 0.);
    std::vector<double> pressure_ptr(nT*nmuB, 0.);
    std::vector<double> net_baryon_ptr(nT*nmuB, 0.);
    std::vector<double> muS_ptr(nT*nmuB, 0.);
    std::vector<double> muQ_ptr(nT*nmuB, 0.);
    for (int i = 0; i < nT; i++) {
        calculate_particle_yield(temp_ptr[i], 0, 0, 0);
        cout << "computing for T = " << temp_ptr[i] << " GeV ..." << endl;
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            double mu_S = 0.;
            double mu_Q = 0.;
            findMuQMuQFromNSNQ(temp_ptr[i], muB_ptr[j], n_S, nQovernB,
                               mu_S, mu_Q);
            calculate_particle_yieldFugacity(temp_ptr[i], muB_ptr[j],
                                             mu_S, mu_Q);
            net_baryon_ptr[idx] = calculateSystemNetbaryonDensity();
            ed_ptr[idx] = calculateSystemenergyDensity(temp_ptr[i]);
            pressure_ptr[idx] = calculateSystemPressure(temp_ptr[i]);
            sd_ptr[idx] = calculateSystementropyDensity(temp_ptr[i]);
            muS_ptr[idx] = mu_S;
            muQ_ptr[idx] = mu_Q;
            //double nS = calculateSystemNetStrangenessDensity();
            //double nQ = calculateSystemNetElectricChargeDensity();
            //cout << "check: T = " << temp_ptr[i] << " GeV, muB = "
            //     << muB_ptr[j] << " GeV, muS = " << muS_ptr[idx] << " GeV, "
            //     << "muQ = " << muQ_ptr[idx] << " GeV. nB = "
            //     << net_baryon_ptr[idx] << ", nS = " << nS << ", nQ = "
            //     << nQ << ", nQ/nB = " << nQ/(net_baryon_ptr[idx] + 1e-16)
            //     << endl;
        }
    }

    // output EOS table
    ostringstream EOSfilename;
    EOSfilename << "./EOS2D_nS_" << n_S << "_nQovernB_" << nQovernB << ".dat";
    ofstream output(EOSfilename.str().c_str());
    output << "# T [GeV]  muB [GeV]  e/T^4  nB/T^3  s/T^3  P/T^3  "
           << "muS [GeV]  muQ [GeV]" << endl;
    for (int i = 0; i < nT; i++) {
        double T3 = pow(temp_ptr[i], 3);
        double T4 = pow(temp_ptr[i], 4);
        for (int j = 0; j < nmuB; j++) {
            int idx = i*nmuB + j;
            output << scientific << setw(20) << setprecision(8)
                   << temp_ptr[i] << "  " << muB_ptr[j] << "  "
                   << ed_ptr[idx]/T4*unitFac << "  "
                   << net_baryon_ptr[idx]/T3*unitFac << "  "
                   << sd_ptr[idx]/T3*unitFac << "  "
                   << pressure_ptr[idx]/T4*unitFac << "  "
                   << muS_ptr[idx] << "  " << muQ_ptr[idx] << endl;
        }
    }
    output.close();
}


void particleList::printParticleContributions() {
    for (unsigned int i = 0; i < partList.size(); i++) {
        cout << "particle " << partList[i]->getMonval() << ": "
             << partList[i]->getParticleYield() << "  "
             << partList[i]->getEnergyDensity() << endl;
    }
}
