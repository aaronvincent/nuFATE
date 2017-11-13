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

#ifndef NUFACE_SECS_H
#define NUFACE_SECS_H

struct _result{
        double* eval;
        double* evec;
        double* ci;
        double* energy_nodes;
        double* phi_0_;
};    

typedef struct _result result;

class nuFACE_secs {

public:
    
    result result1(double* eval, double* evec, double* ci, double* energy_nodes,double*phi_0_);
    double readDoubleAttribute(hid_t, std::string);
    unsigned int readUIntAttribute(hid_t, std::string);
    double* logspace(double min,double max,unsigned int samples);
    double* get_glashow_total(unsigned int NumNodes, double* energy_nodes);
    double* get_RHS_matrices(unsigned int NumNodes, double* energy_nodes, double* sigma_array_, double* sig3_array_, double* dxs_array_, double* sec_array_, double* regen_array_);
    result get_eigs(int flavor, double gamma, std::string h5_filename);

private:

    unsigned int NumNodes;
    double* energy_nodes;
    double* sigma_array_;
    double* dxs_array_;
    double* phi_0_;
    double* sig3_array_;
    double* sec_array_;
    double* regen_array_;

    double* RHSMatrix_;
    double* RHSMatrix1_;
    double* RHSMatrix2_;
    double* RHSMatrix3_;
    double* RHSMatrix4_;

    double* glashow_total_;
    double* glashow_piece_;

    int newflavor;
    double newgamma;
    std::string newh5_filename;

    int dxsdim[2];

};

#endif
