#include "nuFATE.h"

namespace nufate {

nuFACE::nuFACE(int flavor, double gamma, std::string h5_file):
  newflavor_(flavor), newgamma_(gamma), newh5_filename_(h5file), GF(1.16e-5) {

};

double nuFACE::readDoubleAttribute(hid_t object, std::string name){
        double target;
        hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
        if(status<0)
            throw std::runtime_error("Failed to read attribute '"+name+"'");
        H5Aclose(attribute_id);
        return target;
}

unsigned int nuFACE::readUIntAttribute(hid_t object, std::string name) const{
  unsigned int target;
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return target;
}

std::vector<double> nuFACE::logspace(double Emin,double Emax,unsigned int div) const {
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
        logpoints[i] = exp(EE);
    logpoints[div-1]=Emax;

    return logpoints;
}

void nuFACE::set_glashow_total(unsigned int NumNodes, std::vector<double> energy_nodes){
    double GF = 1.16e-5;
    double hbarc = 1.97e-14;
    double GW = 2.085;
    double MW = 80.385e0;
    double mmu = 0.106e0;
    double me = 511.e-6;
    double pi = 3.14159265358979323846;

    glashow_total_ = std::maked_shared<double>(malloc(NumNodes*sizeof(double)),free);

    for(int i=0; i<NumNodes; i++){
        double x = *(glashow_total_+i);
        *(glashow_total_ +i) = 2.*me*(*(energy_nodes+i));
        *(glashow_total_ +i) = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return glashow_total_;
}

double* nuFACE::get_glashow_partial(unsigned int NumNodes, double* energy_nodes){

    Enuin_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    for (int i =0; i<NumNodes;i++){
        for(int j = 0; j<NumNodes;j++){
            *(Enuin_+i*NumNodes+j) = *(energy_nodes+i);
        }
    }

    Enu_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    for (int i =0; i<NumNodes;i++){
        for(int j = 0; j<NumNodes;j++){
            *(Enu_+i*NumNodes+j) = *(energy_nodes+j);
            *(Enu_+i*NumNodes+j) = 1 - *(Enu_+i*NumNodes+j)/ *(Enuin_+i*NumNodes+j);
        }
    }

    double GF = 1.16e-5;
    double hbarc=1.97e-14;
    double GW = 2.085;
    double MW = 80.385;
    double MZ = 91.18;
    double me=511.e-6;
    double s2t = 0.23;
    double gL =  s2t-0.5;
    double gR = s2t;
    double pi = 3.14159265358979323846;


    selectron_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    den_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    t1_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    t2_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    t3_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    dsig_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    for (int i =0; i<NumNodes;i++){
        for(int j = 0; j<NumNodes;j++){
            *(selectron_+i*NumNodes+j) = 2.*me* *(Enuin_+i*NumNodes+j);
            *(den_+i*NumNodes+j) = std::pow(1. - *(selectron_+i*NumNodes+j)/std::pow(MW,2),2) + std::pow(GW,2)/std::pow(MW,2);
            *(t1_+i*NumNodes+j) = std:pow(gR,2)/std::pow((1.+ *(Enu_+i*NumNodes+j)* *(selectron_+i*NumNodes+j)/std:pow(MZ,2)),2);
            *(t2_+i*NumNodes+j) = gL/(1.+ *(Enu_+i*NumNodes+j)* *(selectron_+i*NumNodes+j)/std:pow(MZ,2)) + (1. - *(selectron_+i*NumNodes+j)/std::pow(MW,2))/ *(den_+i*NumNodes+j);
            *(t3_+i*NumNodes+j) = GW/MW/ *(den_+i*NumNodes+j);
            if (*(Enu_+i*NumNodes+j) > 0.){
                *(dsig_+i*NumNodes+j) = (std::pow(GF,2)* *(selectron_+i*NumNodes+j)/pi*(*(t1_+i*NumNodes+j)+ (std:pow(*(t2_+i*NumNodes+j),2)+std:pow(*(t3_+i*NumNodes+j),2))*std::pow((1-*(Enu_+i*NumNodes+j),2)))*std::pow(hbarc,2))/ *(Enuin_+i*NumNodes+j);
            } else {
                *(dsig_+i*NumNodes+j) = 0.;
            }
        }
    }
    return dsig_;
}

double* nuFACE::get_RHS_matrices(unsigned int NumNodes, double* energy_nodes, double* sigma_array_, double* dxs_array_){
        size = NumNodes;
        DeltaE_ = (double *)malloc(NumNodes*sizeof(double));
        for(int i = 0; i < NumNodes-1;i++){
            *(DeltaE_ + i) = log10(*(energy_nodes+i+1)) - log10(*(energy_nodes+i));
        }

        RHSMatrix_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
        
        //double RHSMatrix[size][size] = {};
        
        for(int i = 0; i < NumNodes; i++) 
        {
            for(int j= i+1; j < NumNodes; j++){
                double e1 = 1./ *(energy_nodes+j);
                double e2 = *(energy_nodes+i) * *(energy_nodes+i);
                *(RHSMatrix_+i*NumNodes+j) = *(DeltaE_ + j - 1) * *(dxs_array_+j * dxsdim[1]+i) * e1 * e2;
            }
        }

        //double* RHSMatrix_ = &RHSMatrix[0][0];
        return RHSMatrix_;
}       


nuFACE::get_eigs(int flavor, double gamma, std::string h5_filename) {

        newflavor = flavor;
        newgamma = gamma;
        newh5_filename = h5_filename;

        //open h5file containing cross sections (xsh5)

        hid_t file_id,group_id,root_id;

        file_id = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        root_id = H5Gopen(file_id, "/", H5P_DEFAULT);

        std::string grptot = "/total_cross_sections";
        std::string grpdiff = "/differential_cross_sections";
        group_id = H5Gopen(root_id, grptot.c_str(), H5P_DEFAULT);

        double Emin = readDoubleAttribute(group_id, "max_energy");
        double Emax = readDoubleAttribute(group_id, "min_energy");   
        unsigned int NumNodes = readUIntAttribute(group_id, "number_energy_nodes");
        energy_nodes = logspace(Emin, Emax, NumNodes);


        if (flavor == -1) {
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double)); 
            H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_);
        }  else if (flavor == -2){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double)); 
            H5LTread_dataset_double(group_id, "numubarxs", sigma_array_);
        }  else if (flavor == -3){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "nutaubarxs", sigma_array_);            
        }  else if (flavor == 1){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "nuexs", sigma_array_);
        }  else if (flavor == 2){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "numuxs", sigma_array_);
        }  else if (flavor == 3){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
            sigma_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "nutauxs", sigma_array_);
        }
            

        hsize_t dxarraysize[2];
        group_id = H5Gopen(root_id, grpdiff.c_str(), H5P_DEFAULT);
        
        if (flavor > 0){
            H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);    
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "dxsnu", dxs_array_);
        } else {
            H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "dxsnu", dxs_array_);
        }

        RHSMatrix_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, dxs_array_);

        if (flavor = -3){
            std::string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = tauarraysize[1];
            tau_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "tbarfull", tau_array_);
            RHregen_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, tau_array_);
            for (int i = 0; i<NumNodes; i++){
                for(int j=0; j<NumNodes;j++)
                *(RHSMatrix_+i*NumNodes+j) = *(RHSMatrix_+i*NumNodes+j) + *(RHregen_+i*NumNodes+j);
            }
        } else if(flavor = 3){
            std::string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = tauarraysize[1];
            tau_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "tbarfull", tau_array_);
            RHregen_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, tau_array_);
            for (int i = 0; i<NumNodes; i++){
                for(int j=0; j<NumNodes;j++)
                *(RHSMatrix_+i*NumNodes+j) = *(RHSMatrix_+i*NumNodes+j) + *(RHregen_+i*NumNodes+j);
            }
        } else if(flavor = -1){
            double* glashow_total_ = get_glashow_total(NumNodes,energy_nodes);
            for (int i = 0; i < NumNodes; i++){
                *(sigma_array_+i) = *(sigma_array_+i) + *(glashow_total_ + i)/2.; 
                *(RHSMatrix_ +i) = *(RHSMatrix_ +i) + *(glashow_partial_ + i)/2.;

            }
        }
        phi_0_ = (double *)malloc(NumNodes*sizeof(double));
        for (int i = 0; i < NumNodes; i++){
            *(phi_0_ + i) = std::pow(*(energy_nodes +i),(2-gamma));
        }
        for (int i = 0; i < NumNodes; i++){
            *(RHSMatrix_+i*NumNodes+i) = *(RHSMatrix_+i*NumNodes+i) + *(sigma_array_+i);    
        }

        //compute eigenvalues and eigenvectors

        gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_, NumNodes, NumNodes);
        gsl_vector *eval = gsl_vector_alloc (NumNodes);
        gsl_matrix *evec = gsl_matrix_alloc (NumNodes, NumNodes);
        gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (NumNodes);

        gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

        int s;
        gsl_vector *ci = gsl_vector_alloc(NumNodes);
        gsl_permutation *p = gsl_permutation_alloc(NumNodes);
        
        gsl_linalg_LU_decomp (&m.matrix, p, &s);
        gsl_linalg_LU_solve (&m.matrix, p, phi_0_, ci);

        //free unneeded memory
        gsl_permutation_free (p);
        gsl_eigen_nonsymmv_free (w);

        return eval, evec, ci, energy_nodes, phi_0_;



}

int nuFACE::getFlavor() const {
    return newflavor;
}

double nuFACE::getGamma() const {
    return newgamma;
}

std::string nuFACE::getFilename() const {
    return newh5_filename;
}

} // close namespace
