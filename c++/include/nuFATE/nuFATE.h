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

#ifndef NUFATE_H
#define NUFATE_H

namespace nufate {

struct Result {
  double* eval;
  double* evec;
  double* ci;
  std::vector<double> energy_nodes_;
  std::vector<double> phi_0_;
};

class nuFACE {
  private:
    double const GF = 1.16e-5;
    double const hbarc = 1.97e-14;
    double const GW = 2.085;
    double const MW = 80.385e0;
    double const mmu = 0.106e0;
    double const me = 511.e-6;
    double const pi = 3.14159265358979323846;
    double const MZ = 91.18;
    double const s2t = 0.23;
    double const gL =  s2t-0.5;
    double const gR = s2t;
  private:
    int newflavor_;
    double newgamma_;
    std::string newh5_filename_;
  public:
    nuFACE(int, double, std::string);
    Result get_eigs();
    int getFlavor() const;
    double getGamma() const;
    std::string getFilename() const;
  protected:
    void set_glashow_total(unsigned int NumNodes_, std::vector<double> energy_nodes_);
    void set_glashow_partial(unsigned int NumNodes_, std::vector<double> energy_nodes_);
    void set_RHS_matrices(std::shared_ptr<double> RMatrix_, unsigned int NumNodes, std::vector<double> energy_nodes, std::vector<double> sigma_array_, std::shared_ptr<double> dxs_array_);
    double readDoubleAttribute(hid_t, std::string) const;
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::vector<double> logspace(double min,double max,unsigned int samples) const;
  private:
    std::string grpdiff;
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


#endif
