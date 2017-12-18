#include "nuFATE.h"
#include <iostream>

namespace nufate{

nuFATE::nuFATE(int flavor, double gamma, std::string h5_filename) : newflavor_(flavor), newgamma_(gamma), newh5_filename_(h5_filename) {
    //open h5file containing cross sections (xsh5)
    file_id_ = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
    grptot_ = "/total_cross_sections";
    grpdiff_ = "/differential_cross_sections";
    hid_t group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);
    //assign some important variables
    Emax_ = readDoubleAttribute(group_id, "max_energy");
    Emin_ = readDoubleAttribute(group_id, "min_energy");
    NumNodes_ = readUIntAttribute(group_id, "number_energy_nodes");
    // allocate memory
    AllocateMemoryForMembers(NumNodes_);
    //set energy nodes and deltaE
    energy_nodes_ = logspace(Emin_, Emax_, NumNodes_);
    // calculate and set energy bin widths
    SetEnergyBinWidths();
    // load cross sections from file
    LoadCrossSectionFromHDF5();
    // set the initial flux
    SetInitialFlux();
}

nuFATE::nuFATE(int flavor, double gamma, std::vector<double> energy_nodes, std::vector<double> sigma_array, std::vector<std::vector<double>> dsigma_dE):
  newflavor_(flavor), newgamma_(gamma), energy_nodes_(energy_nodes), sigma_array_(sigma_array)
{
  NumNodes_ = energy_nodes_.size();
  Emax_ = energy_nodes_.back();
  Emin_ = energy_nodes_.front();
  if(sigma_array.size() != NumNodes_)
    throw std::runtime_error("nuFATE::nuFATE Total cross section array does not match energy nodes size.");
  if(dsigma_dE.size() != NumNodes_ or dsigma_dE.front().size() != NumNodes_)
    throw std::runtime_error("nuFATE::nuFATE Differential cross section array does not match energy nodes size.");
  AllocateMemoryForMembers(NumNodes_);
  SetEnergyBinWidths();
  SetInitialFlux();
  SetCrossSectionsFromInput(dsigma_dE);
}

void nuFATE::SetCrossSectionsFromInput(std::vector<std::vector<double>> dsigma_dE){
    for(unsigned int i = 0; i<NumNodes_; i++){
        for(unsigned int j=0; j<NumNodes_; j++)
        *(dxs_array_.get()+i*NumNodes_+j) = dsigma_dE[i][j];
    }
    total_cross_section_set_ = true;
    differential_cross_section_set_ = true;
}

void nuFATE::AllocateMemoryForMembers(unsigned int NumNodes){
  // if the energy nodes are different or if memory has not been allocated. Allocate it.
  if(NumNodes_ != NumNodes or (not memory_allocated_)){
    NumNodes_ = NumNodes;
    //allocate memory that will be used in functions below
    glashow_total_ = std::vector<double>(NumNodes_);
    glashow_partial_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    A_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    Enuin_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    Enu_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    den_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    selectron_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    dxs_array_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  }
  memory_allocated_ = true;
}

void nuFATE::SetEnergyBinWidths(){
  if(energy_nodes_.size() == 0)
    throw std::runtime_error("nuFATE::SetEnergyBinWidths Energy nodes need to be set before the widths can be set.");
  DeltaE_ = std::vector<double>(NumNodes_);
  for(unsigned int i = 0; i < NumNodes_-1;i++){
      DeltaE_[i] = log(energy_nodes_[i+1]) - log(energy_nodes_[i]);
  }
}

//reads an attribute of type double from h5 object
double nuFATE::readDoubleAttribute(hid_t object, std::string name) const{
    double target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}

//reads an attribute of type unsigned int from h5 object
unsigned int nuFATE::readUIntAttribute(hid_t object, std::string name) const{
    unsigned int target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}

//given a minimum, maximum, and number of points, this function returns equally spaced array in logspace
std::vector<double> nuFATE::logspace(double Emin,double Emax,unsigned int div) const {
    if(div==0)
        throw std::length_error("number of samples requested from logspace must be nonzero");
    std::vector<double> logpoints(div);
    double Emin_log,Emax_log;
    Emin_log = log10(Emin);
    Emax_log = log10(Emax);
    double step_log = (Emax_log - Emin_log)/double(div-1);
    logpoints[0]=Emin;
    double EE = Emin_log+step_log;
    for(unsigned int i=1; i<div-1; i++, EE+=step_log)
        logpoints[i] = std::pow(10,EE);
    logpoints[div-1]=Emax;
    return logpoints;
}

//sets the contribution from glashow
void nuFATE::set_glashow_total(){
    for(unsigned int i=0; i<NumNodes_; i++){
        glashow_total_[i] = 2.*me*energy_nodes_[i];
        double x = glashow_total_[i];
        glashow_total_[i] = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return;
}

void nuFATE::set_glashow_partial(){
    for (unsigned int i =0; i<NumNodes_;i++){
        for(unsigned int j = 0; j<NumNodes_;j++){
            *(Enuin_.get()+i*NumNodes_+j) = energy_nodes_[i];
            *(Enu_.get()+i*NumNodes_+j) = energy_nodes_[j];
            *(Enu_.get()+i*NumNodes_+j) = 1 - *(Enu_.get()+i*NumNodes_+j)/ *(Enuin_.get()+i*NumNodes_+j);
            *(selectron_.get()+i*NumNodes_+j) = 2.*me* *(Enuin_.get()+i*NumNodes_+j);
            *(den_.get()+i*NumNodes_+j) = std::pow(1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2),2) + std::pow(GW,2)/std::pow(MW,2);
            *(t1_.get()+i*NumNodes_+j) = std::pow(gR,2)/std::pow((1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)),2);
            *(t2_.get()+i*NumNodes_+j) = gL/(1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)) + (1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2))/ *(den_.get()+i*NumNodes_+j);
            *(t3_.get()+i*NumNodes_+j) = GW/MW/ *(den_.get()+i*NumNodes_+j);
            if (*(Enu_.get()+i*NumNodes_+j) > 0.){
                *(glashow_partial_.get()+i*NumNodes_+j) = (std::pow(GF,2)* *(selectron_.get()+i*NumNodes_+j)/pi*(*(t1_.get()+i*NumNodes_+j)+ (std::pow(*(t2_.get()+i*NumNodes_+j) , 2) + std::pow(*(t3_.get()+i*NumNodes_+j),2)) * std::pow((1.-*(Enu_.get()+i*NumNodes_+j)),2))*std::pow(hbarc,2))/ *(Enuin_.get()+i*NumNodes_+j);
            } else {
                *(glashow_partial_.get()+i*NumNodes_+j) = 0.;
            }
        }
    }
    return;
}

void nuFATE::set_RHS_matrices(std::shared_ptr<double> RMatrix, std::shared_ptr<double> dxsarray) {
    for(unsigned int i = 0; i < NumNodes_; i++)
    {
        for(unsigned int j= i+1; j < NumNodes_; j++){
            double e1 = 1./ energy_nodes_[j];
            double e2 = energy_nodes_[i] * energy_nodes_[i];
            *(RMatrix.get()+i*NumNodes_+j) = DeltaE_[j - 1] * *(dxsarray.get()+j * dxsdim_[1]+i) * e1 * e2;
        }
    }
    RHS_set_ = true;
    return;
}

void nuFATE::LoadCrossSectionFromHDF5(){
    hid_t group_id;
    group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);

    if (newflavor_ == -1) {
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_.data());
    }  else if (newflavor_ == -2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "numubarxs", sigma_array_.data());
    }  else if (newflavor_ == -3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nutaubarxs", sigma_array_.data());
    }  else if (newflavor_ == 1){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nuexs", sigma_array_.data());
    }  else if (newflavor_ == 2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
        H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "numuxs", sigma_array_.data());
    }  else if (newflavor_ == 3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nutauxs", sigma_array_.data());
    }

    hsize_t dxarraysize[2];
    group_id = H5Gopen(root_id_, grpdiff_.c_str(), H5P_DEFAULT);
   if (newflavor_ > 0){
        H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        if((unsigned int)(dxsdim_[0]) != NumNodes_ or (unsigned int)(dxsdim_[1]) != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Secondaries arrays do not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "dxsnu", dxs_array_.get());
    } else {
        H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        if((unsigned int)(dxsdim_[0]) != NumNodes_ or (unsigned int)(dxsdim_[1]) != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Secondaries arrays do not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "dxsnubar", dxs_array_.get());
    }

   total_cross_section_set_ = true;
   differential_cross_section_set_ = true;
}

void nuFATE::SetInitialFlux(){
    phi_0_ = std::vector<double>(NumNodes_);
    for (unsigned int i = 0; i < NumNodes_; i++){
        phi_0_[i] = std::pow(energy_nodes_[i],(2.-newgamma_));
    }
    initial_flux_set_ = true;
}

void nuFATE::AddSecondaryTerms(){
    hid_t group_id;
    if (newflavor_ == -3 and add_tau_regeneration_){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", tau_array_.get());
        RHregen_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
        set_RHS_matrices(RHregen_, tau_array_);
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++)
            *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
        }
    } else if (newflavor_ == 3 and add_tau_regeneration_){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", tau_array_.get());
        RHregen_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
        set_RHS_matrices(RHregen_, tau_array_);
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++)
            *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
        }
    } else if (newflavor_ == -1 and add_glashow_term_){
        set_glashow_total();
        for (unsigned int i = 0; i < NumNodes_; i++){
            sigma_array_[i] = sigma_array_[i] + glashow_total_[i]/2.;
            *(RHSMatrix_.get() +i) = *(RHSMatrix_.get() +i) + *(glashow_partial_.get() + i)/2.;
        }
    }
}

Result nuFATE::getEigensystem(){
    if(not initial_flux_set_)
      throw std::runtime_error("nuFATE::getEigensystem initial flux not set.");
    if(not total_cross_section_set_)
      throw std::runtime_error("nuFATE::getEigensystem total cross section not set.");
    if(not differential_cross_section_set_)
      throw std::runtime_error("nuFATE::getEigensystem differential cross section not set.");

    set_RHS_matrices(RHSMatrix_, dxs_array_);

    if(add_secondary_term_)
      AddSecondaryTerms();

    for (unsigned int i = 0; i < NumNodes_; i++){
        *(RHSMatrix_.get()+i*NumNodes_+i) = *(RHSMatrix_.get()+i*NumNodes_+i) + sigma_array_[i];
    }

    gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_.get(), NumNodes_, NumNodes_);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (NumNodes_);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (NumNodes_, NumNodes_);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (NumNodes_);

    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

    int s;
    gsl_vector *ci = gsl_vector_alloc(NumNodes_);
    gsl_permutation * p = gsl_permutation_alloc(NumNodes_);
    gsl_vector_view b = gsl_vector_view_array (&phi_0_.front(), NumNodes_);
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, ci);

    //free unneeded memory
    gsl_permutation_free (p);
    gsl_eigen_nonsymmv_free (w);

    struct Result r1;
    r1.eval = eval->data;
    r1.evec = evec->data;
    r1.ci = ci->data;
    r1.energy_nodes_ = energy_nodes_;
    r1.phi_0_ = phi_0_;

    return r1;
}

struct rho_earth_params{double theta;};

double nuFATE::rho_earth(double x, void * p){
    double RE = 6371.;
    struct rho_earth_params * params = (struct rho_earth_params *)p;
    double theta = (params->theta);
    double xmax = 2.*abs(RE*cos(theta));
    double R = std::pow(RE,2) + std::pow((xmax-x),2) + 2.*RE*(xmax-x)*cos(theta);
    double r = std::pow(R,0.5);
    double p1;
    double p2;
    double p3;

    if (r<1221.){
        p1 = -0.0002177;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
    } else if (r<3480.){
        p1 = -0.0002409;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
    } else if (r<5721.){
        p1 = -3.764e-05;
        p2 = -0.1876;
        p3 = 6664.;
    } else if (r<5961.){
        p1 = 0.;
        p2 = -1.269;
        p3 = 1.131e+04;
    } else if (r<6347.){
        p1 = 0.;
        p2 = -0.725;
        p3 = 7887.;
    } else if (r<6356.){
        p1 = 0.;
        p2 = 0.;
        p3 = 2900.;
    } else if (r<6368.){
        p1 = 0.;
        p2 = 0.;
        p3 = 2600.;
    } else {
        p1 = 0.;
        p2 = 0.;
        p3 = 1020.;
    }

    double rho = p1*std::pow(r,2)+p2*r+p3;
    return rho*1.0e-3;
}

double nuFATE::getEarthColumnDensity(double theta){
    double t;
    if (theta < pi/2.){
       t = 0;
    } else {
      double kmtocm = 1.0e5;
      double result, error;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      gsl_function F;
      struct rho_earth_params params = {theta};
      F.function = &nuFATE::rho_earth;
      F.params = &params;
      double xmax = 2.*abs(REarth*cos(theta));

      gsl_integration_qags(&F, 0, xmax, 1.0e-18, 1.0e-3, 1000, w, &result, &error);
      t = result*kmtocm;
    }

   return t;
}


int nuFATE::getFlavor() const {
    return newflavor_;
}

double nuFATE::getGamma() const {
    return newgamma_;
}

std::string nuFATE::getFilename() const {
    return newh5_filename_;
}

double nuFATE::getNumNodes() const {
    return NumNodes_;
}

} //close namespace
