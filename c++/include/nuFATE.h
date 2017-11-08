#include <string>
#include <math.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"


#ifndef NUFATE.H
#define NUFATE.H

class nuFACE {

public:


    double readDoubleAttribute(hid_t, std::string);

    double* logspace(double min,double max,unsigned int samples);

    double* get_RHS_matrices(double dim, double* nodes, std::shared_ptr<double> sigma, double* dxs);

    get_eigs(int flavor, double gamma, std::string h5_filename);

    int getFlavor() const;

    double getGamma() const;

    std::string getFilename() const;

private:

    //dimensions of differential cross section array
    double dxsdim[2];
    //spectral index
    double newgamma;
    //filename containing cross sections
    std::string newh5_filename;
    double* energy_nodes;
    std::shared_ptr<double> sigma_array_;
    //differential cross section 2d array
    double* dxs_array_;
    //tau regeneration
    double* tau_array_;
    
    double* DeltaE_;
    
    double* RHSMatrix_;
    
    double* RHregen_;
    
    double* glashow_total_;
    
    int newflavor;
}
