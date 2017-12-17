#include "nuFATE_secs.h"

namespace nufate{

nuFATE_secs::nuFATE_secs(int flavor, double gamma, std::string h5_filename) : newflavor_(flavor), newgamma_(gamma), newh5_filename_(h5_filename) {
    //open h5file containing cross sections (xsh5)
    hid_t group_id;
    file_id_ = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
    grptot_ = "/total_cross_sections";
    grpdiff_ = "/differential_cross_sections";
    group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);
    //assign some important variables
    Emax_ = readDoubleAttribute(group_id, "max_energy");
    Emin_ = readDoubleAttribute(group_id, "min_energy");   
    NumNodes_ = readUIntAttribute(group_id, "number_energy_nodes");
    //allocate memory that will be used in functions below
    energy_nodes_.resize(NumNodes_);
    DeltaE_ = std::vector<double>(NumNodes_);
    glashow_total_ = std::vector<double>(NumNodes_);
    glashow_partial_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free); 
    RHSMatrix1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix4_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);  
    Enuin_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    Enu_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    den_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    selectron_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    energy_nodes_ = logspace(Emin_, Emax_, NumNodes_);
    for(unsigned int i = 0; i < NumNodes_-1;i++){
        DeltaE_[i] = log10(energy_nodes_[i+1]) - log10(energy_nodes_[i]);
    }


}

double nuFATE_secs::readDoubleAttribute(hid_t object, std::string name) const{
        double target;
        hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
        if(status<0)
            throw std::runtime_error("Failed to read attribute '"+name+"'");
        H5Aclose(attribute_id);
        return target;
}

unsigned int nuFATE_secs::readUIntAttribute(hid_t object, std::string name) const{
  unsigned int target;
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return target;
}

std::vector<double> nuFATE_secs::logspace(double Emin,double Emax,unsigned int div) const {
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

void nuFATE_secs::set_glashow_total(){
    for(unsigned int i=0; i<NumNodes_; i++){
        glashow_total_[i] = 2.*me*energy_nodes_[i];
        double x = glashow_total_[i];
        glashow_total_[i] = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return;
}

void nuFATE_secs::set_glashow_partial(){
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


void nuFATE_secs::set_RHS_matrices(){
    for(unsigned int i=0; i<NumNodes_; i++){
        for(unsigned int j=0; j<NumNodes_; j++){
                *(RHSMatrix1_.get()+i*NumNodes_+j) = 0.;
                *(RHSMatrix2_.get()+i*NumNodes_+j) = 0.;
                *(RHSMatrix3_.get()+i*NumNodes_+j) = 0.;
                *(RHSMatrix4_.get()+i*NumNodes_+j) = 0.;
            }
    }
    for (unsigned int i=0; i<NumNodes_; i++){
        for(unsigned int j=i+1; j<NumNodes_; j++){
            *(RHSMatrix1_.get()+i*NumNodes_+j) = DeltaE_[j-1] * *(dxs_array_.get()+j*NumNodes_+i) * (1./energy_nodes_[j]);
            *(RHSMatrix2_.get()+i*NumNodes_+j) = DeltaE_[j-1] * *(sec_array_.get()+j*NumNodes_+i) * (1./energy_nodes_[j]) * std::pow(energy_nodes_[i],2);
            *(RHSMatrix4_.get()+i*NumNodes_+j) = DeltaE_[j-1] * *(dxs_array_.get()+j*NumNodes_+i) * *(regen_array_.get()+j*NumNodes_+i) * (1./energy_nodes_[j]) * std::pow(energy_nodes_[i],2);
        }
    }
    for (unsigned int i = 0; i < NumNodes_; i++){
        *(RHSMatrix1_.get()+i*NumNodes_+i) = *(RHSMatrix1_.get()+i*NumNodes_+i) - sigma_array_[i];    
        *(RHSMatrix4_.get()+i*NumNodes_+i) = *(RHSMatrix4_.get()+i*NumNodes_+i) - sig3_array_[i];
    }
    rsize_ = 2*NumNodes_;
    for (unsigned int i = 0; i<rsize_; i++){
        if(i<NumNodes_){
            for (unsigned int j = 0; j<rsize_;j++){
                if (j<NumNodes_){
                    *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix1_.get()+i*NumNodes_+j);
                } else {
                        *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix3_.get()+i*NumNodes_+j);
                }
            }
        } else {
            for (unsigned int j = 0; j<rsize_;j++){
                if (j<NumNodes_){
                    *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix2_.get()+i*NumNodes_+j);
                } else {
                        *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix4_.get()+i*NumNodes_+j);
                 }
            }
        }
    }

    return;
}

Result nuFATE_secs::get_eigs() {

    hid_t group_id;
    group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);

    if (newflavor_ == -1) {
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]); 
        H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_.data());
    }  else if (newflavor_ == -2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]); 
        H5LTread_dataset_double(group_id, "numubarxs", sigma_array_.data());
    }  else if (newflavor_ == 1){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nuexs", sigma_array_.data());
    }  else if (newflavor_ == 2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "numuxs", sigma_array_.data());
    }

    if (newflavor_ > 0){        
        hsize_t sarraysize[1];
        size_t dim1, dim2;
        H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
        sig3_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nutauxs", sig3_array_.data());
        
        hsize_t dxarraysize[2];
        group_id = H5Gopen(root_id_, grpdiff_.c_str(), H5P_DEFAULT);        
        H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);    
        dim1 = dxarraysize[0];
        dim2 = dxarraysize[1];
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        dxs_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "dxsnu", dxs_array_.get());
        
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
        dim1 = tauarraysize[0];
        dim2 = tauarraysize[1];
        regen_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tfull", regen_array_.get());
        hsize_t secarraysize[2];
        H5LTget_dataset_info(group_id,"secfull", secarraysize,NULL,NULL);
        sec_array_ = std::shared_ptr<double>((double *)malloc(secarraysize[0]*secarraysize[1]*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "secfull", sec_array_.get());
    } else {
        hsize_t sarraysize[1];
        size_t dim1, dim2;
        H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
        sig3_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nutaubarxs", sig3_array_.data());

        hsize_t dxarraysize[2];
        group_id = H5Gopen(root_id_, grpdiff_.c_str(), H5P_DEFAULT);        
        H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
        dim1 = dxarraysize[0];
        dim2 = dxarraysize[1];
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        dxs_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "dxsnubar", dxs_array_.get());

        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
        dim1 = tauarraysize[0];
        dim2 = tauarraysize[1];
        regen_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", regen_array_.get());

        hsize_t secarraysize[2];
        H5LTget_dataset_info(group_id,"secbarfull", secarraysize,NULL,NULL);
        sec_array_ = std::shared_ptr<double>((double *)malloc(secarraysize[0]*secarraysize[1]*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "secbarfull", sec_array_.get());
    }

    set_RHS_matrices();

    if (newflavor_ == -1){
        set_glashow_total();
        set_glashow_partial();

        for (unsigned int i = 0; i<NumNodes_;i++){
            *(glashow_partial_.get()+i*NumNodes_+i) = (*(glashow_partial_.get()+i*NumNodes_+i) + glashow_total_[i])/2.;
        }
        for(unsigned int i=0; i<NumNodes_;i++){
            for(unsigned int j =0; j<NumNodes_;j++){
                *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(glashow_partial_.get()+i*NumNodes_+j);
            }
        }  
    }
    phi_0_ = std::vector<double>(2*NumNodes_);
    for (unsigned int i = 0; i < NumNodes_; i++){
        phi_0_[i] = std::pow(energy_nodes_[i],(2.-newgamma_));
        phi_0_[i+NumNodes_] = std::pow(energy_nodes_[i],(2.-newgamma_));
    }

    gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_.get(), rsize_, rsize_);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (rsize_);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (rsize_, rsize_);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (rsize_);

    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

    int s;
    gsl_vector *ci = gsl_vector_alloc(rsize_);
    gsl_permutation *p = gsl_permutation_alloc(rsize_);
    gsl_vector_view b = gsl_vector_view_array (&phi_0_.front(), rsize_);
    
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

struct rho_earth_params {double theta;};

double nuFATE_secs::rho_earth(double x, void * p){
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

double nuFATE_secs::get_t_earth(double theta){
    double t;
    if (theta < pi/2.){
       t = 0;
    } else {
      double kmtocm = 1.0e5;
      double result, error;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      gsl_function F;
      struct rho_earth_params params = {theta};
      F.function = &nuFATE_secs::rho_earth;
      F.params = &params;
      double xmax = 2.*abs(REarth*cos(theta));

      gsl_integration_qags(&F, 0, xmax, 1.0e-18, 1.0e-3, 1000, w, &result, &error);
      t = result*kmtocm;
    }

   return t;
}




} //close namespace
