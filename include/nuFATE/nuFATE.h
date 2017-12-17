#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <memory>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

#ifndef NUFATE_H
#define NUFATE_H

namespace nufate{

struct Result {
  double* eval;
  double* evec;
  double* ci;
  std::vector<double> energy_nodes_;
  std::vector<double> phi_0_;
};

class nuFATE {
  private:
    const double GF = 1.16e-5;
    const double hbarc = 1.97e-14;
    const double GW = 2.085;
    const double MW = 80.385e0;
    const double mmu = 0.106e0;
    const double me = 511.e-6;
    const double pi = 3.14159265358979323846;
    const double MZ = 91.18;
    const double s2t = 0.23;
    const double gL =  s2t-0.5;
    const double REarth = 6371.;
    const double gR = s2t;
  private:
    int newflavor_;
    double newgamma_;
    std::string newh5_filename_;
  public:
    nuFATE(int, double, std::string);
    Result get_eigs();
    double get_t_earth(double theta); //calculates column density for a given angle in radians
    int getFlavor() const;
    double getGamma() const;
    std::string getFilename() const;
    double getNumNodes() const;
  protected:
    void set_glashow_total();
    void set_glashow_partial();
    void set_RHS_matrices(std::shared_ptr<double> RMatrix_, std::shared_ptr<double> dxs_array);
    static double rho_earth(double theta, void * p); //returns density for a given angle in radians along x distance
    double readDoubleAttribute(hid_t, std::string) const;
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::vector<double> logspace(double min,double max,unsigned int samples) const;
  private:
    std::string grpdiff_;
    std::string grptot_;
    hid_t file_id_;
    hid_t root_id_;
    unsigned int NumNodes_;
    double Emax_;
    double Emin_;
    int dxsdim_[2];
    std::vector<double> energy_nodes_;
    std::vector<double> sigma_array_;
    std::vector<double> DeltaE_;
    std::vector<double> phi_0_;
    std::vector<double> glashow_total_;
    std::shared_ptr<double> dxs_array_;
    std::shared_ptr<double> tau_array_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> A_;
    std::shared_ptr<double> RHregen_;
    std::shared_ptr<double> glashow_partial_;
    std::shared_ptr<double> Enuin_;
    std::shared_ptr<double> Enu_;
    std::shared_ptr<double> selectron_;
    std::shared_ptr<double> den_;
    std::shared_ptr<double> t1_;
    std::shared_ptr<double> t2_;
    std::shared_ptr<double> t3_;

};

}

#endif
