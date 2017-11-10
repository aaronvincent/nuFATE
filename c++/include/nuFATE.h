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


#ifndef NUFATE_H
#define NUFATE_H

class nuFACE {

public:

    double readDoubleAttribute(hid_t, std::string);

    unsigned int readUIntAttribute(hid_t, std::string);

    double* logspace(double min,double max,unsigned int samples);

    double* get_glashow_total(unsigned int NumNodes, double* energy_nodes);

    double* get_RHS_matrices(unsigned int NumNodes, double* energy_nodes, double* sigma_array_, double* dxs_array_);

    get_eigs(int, double, std::string);


    int getFlavor() const;

    double getGamma() const;

    std::string getFilename() const;

private:
    unsigned int size;
    int newflavor;
    int dxsdim[2];
    double* phi_0_;
    double newgamma;
    std::string newh5_filename;
    double* energy_nodes;
    double* sigma_array_;
    double* dxs_array_; 
    double* tau_array_;
    double* DeltaE_;
    double* RHSMatrix_;
    double* RHregen_;
    double* glashow_total_;
    unsigned int NumNodes;

};

#endif
