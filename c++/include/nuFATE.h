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

struct _result{
        std::shared_ptr<double> eval;
        std::shared_ptr<double> evec;
        std::shared_ptr<double> ci;
        std::shared_ptr<double> energy_nodes;
        std::shared_ptr<double> phi_0_;
};    

typedef struct _result result;


class nuFACE {

public:
    
    result result1(double* eval, double* evec, double* ci, double* energy_nodes,double*phi_0_);
    double readDoubleAttribute(hid_t, std::string);
    unsigned int readUIntAttribute(hid_t, std::string);
    std::shared_ptr<double> logspace(double min,double max,unsigned int samples);
    std::shared_ptr<double> get_glashow_total(unsigned int NumNodes, std::shared_ptr<double> energy_nodes);
    std::shared_ptr<double> get_RHS_matrices(unsigned int NumNodes, std::shared_ptr<double> energy_nodes, std::shared_ptr<double> sigma_array_, std::shared_ptr<double> dxs_array_);
    result1 get_eigs(int flavor, double gamma, std::string h5_filename);
    int getFlavor() const;
    double getGamma() const;
    std::string getFilename() const;

private:
    unsigned int size;
    int newflavor;
    int dxsdim[2];
    std::shared_ptr<double> phi_0_;
    double newgamma;
    std::string newh5_filename;
    std::shared_ptr<double> energy_nodes;
    std::shared_ptr<double> logpoints_;
    std::shared_ptr<double> sigma_array_;
    std::shared_ptr<double> dxs_array_; 
    std::shared_ptr<double> tau_array_;
    std::shared_ptr<double> DeltaE_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> RHregen_;
    std::shared_ptr<double> glashow_total_;
    unsigned int NumNodes;

};

#endif
