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


#ifndef NUFATE_H
#define NUFATE_H

namespace nufate {

struct Result {
  std::shared_ptr<double> eval;
  std::shared_ptr<double> evec;
  std::shared_ptr<double> ci;
  std::shared_ptr<double> energy_nodes;
  std::shared_ptr<double> phi_0_;
};

class nuFACE {
  private:
    const double GD = 1.;
  private:
    int newflavor_;
    double newgamma_;
    std::string newh5_filename_;
  public:
    nuFACE(int flavor, double gamma, std::string h5_file);
    Result get_eigs();
    int getFlavor() const;
    double getGamma() const;
    std::string getFilename() const;
  protected:
    void set_glashow_total(unsigned int NumNodes, std::shared_ptr<double> energy_nodes);
    std::shared_ptr<double> get_RHS_matrices(unsigned int NumNodes, std::shared_ptr<double> energy_nodes, std::shared_ptr<double> sigma_array_, std::shared_ptr<double> dxs_array_);
    double readDoubleAttribute(hid_t, std::string) const;
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::shared_ptr<double> logspace(double min,double max,unsigned int samples) const;
  private:
    unsigned int size_;
    int dxsdim_[2];
    std::shared_ptr<double> phi_0_;
    std::shared_ptr<double> energy_nodes_;
    std::shared_ptr<double> logpoints_;
    std::shared_ptr<double> sigma_array_;
    std::shared_ptr<double> dxs_array_;
    std::shared_ptr<double> tau_array_;
    std::shared_ptr<double> DeltaE_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> RHregen_;
    std::shared_ptr<double> glashow_total_;
    unsigned int NumNodes_;
};

} // close namespace

#endif
