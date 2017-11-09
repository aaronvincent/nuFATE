#include "nuFATE.h"



double nuFACE::readDoubleAttribute(hid_t object, std::string name){
        double target;
        hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
        if(status<0)
            throw std::runtime_error("Failed to read attribute '"+name+"'");
        H5Aclose(attribute_id);
        return target;
}

double* nuFACE::logspace(double Emin,double Emax,unsigned int div){
        if(div==0)
            throw std::length_error("number of samples requested from logspace must be nonzero");
        double logpoints[div];
        double Emin_log,Emax_log;
        Emin_log = log(Emin);
        Emax_log = log(Emax);
        
        double step_log = (Emax_log - Emin_log)/double(div-1);
        
        logpoints[0]=Emin;
        double EE = Emin_log+step_log;
        for(unsigned int i=1; i<div-1; i++, EE+=step_log)
            logpoints[i] = exp(EE);
        logpoints[div-1]=Emax;
        double* logpoints_ = &logpoints[0];
        return logpoints_;
}

double* nuFACE::get_glashow_total(double NumNodes, double* energy_nodes){
    double GF = 1.16e-5
    double hbarc = 1.97e-14
    double GW = 2.085
    double MW = 80.385e0
    double mmu = 0.106e0
    double me = 511.e-6
    double pi = 3.14159265358979323846
    glashow_total_ = (double *)malloc(NumNodes*sizeof(double));

    for(int i=0; i<NumNodes; i++){
        double x = *(glashow_total_+i);
        *(glashow_total_ +i) = 2.*me*(*(energy_nodes+i));
        *(glashow_total_ +i) = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return glashow_total_;
}


double* nuFACE::get_RHS_matrices(double NumNodes, double* energy_nodes, std::shared_ptr<double> sigma_array_, double* dxs_array_){
        
        DeltaE_ = (double *)malloc(NumNodes*sizeof(double));
        for(int i = 0; i < NumNodes-1;i++){
            *(DeltaE_ + i) = log10(*(energy_nodes+i+1)) - log10(*(energy_nodes+i));
        }

        double RHSMatrix[NumNodes][NumNodes] = {};
        
        for(int i = 0; i < NumNodes; i++) 
        {
            for(int j= i+1; j < NumNodes )
            {
                double e1 = 1./ *(energy_nodes+j);
                double e2 = *(energy_nodes+i) * *(energy_nodes+i);
                RHSMatrix[i][j] = *(DeltaE_ + j - 1) * *(dxs_array_+j * dxsdim[1]+i) * e1 * e2;
            }
        }

        double* RHSMatrix_ = &RHSMatrix[0][0];
        return RHSMatrix_;
}       

nuFACE::get_eigs(int flavor, double gamma, string h5_filename) {

        newflavor = flavor;
        newgamma = gamma;
        newh5_filename = h5_filename;

        //open h5file containing cross sections

        hid_t file_id,group_id,root_id;

        file_id = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        root_id = H5Gopen(file_id, "/", H5P_DEFAULT);

        std::string grptot = "/total_cross_sections";
        std::string grpdiff = "/differential_cross_sections";
        group_id = H5Gopen(root_id, grptot.c_str(), H5P_DEFAULT);
        
        //Get energy information
        double Emin = readDoubleAttribute(group_id, 'max_energy')
        double Emax = readDoubleAttribute(group_id, 'min_energy')
        double NumNodes = readDoubleAttribute(group_id, 'number_energy_nodes')
        energy_nodes = logspace(Emin, Emax, Numnodes)

        //Get sigma_array 
        
        if (flavor == -1) {
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_);
        }  else if (flavor == -2){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "numubarxs", sigma_array_);
        }  else if (flavor == -3){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "nutaubarxs", sigma_array_);            
        }  else if (flavor == 1){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "nuexs", sigma_array_);
        }  else if (flavor == 2){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "numuxs", sigma_array_);
        }  else if (flavor == 3){
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
            sigma_array_ = std::make_shared<double>(sarraysize[0]); 
            H5LTread_dataset_double(group_id, "nutauxs", sigma_array_);
        }
            
        //Get differential cross sections
        
        hsize_t dxarraysize[2];
        group_id = H5Gopen(root_id, grpdiff.c_str(), H5P_DEFAULT);
        
        if (flavor > 0){
            H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);    
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset(group_id, "dxsnu", dxs_array_);
        } else {
            H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset(group_id, "dxsnu", dxs_array_);
        }

        //Find RHS matrix
        
        RHSMatrix_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, dxs_array_);

        //Account for tau regeneration/Glashow resonance

        if (flavor == -3){
            std:string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = taurraysize[1];
            tau_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset(group_id, "tbarfull", tau_array_);
            RHregen_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, tau_array_);
            for (int i = 0; i<NumNodes; i++){
                for(int j=0; j<NumNodes;j++)
                *(RHSMatrix_+i*NumNodes+j) = *(RHSMatrix_+i*NumNodes+j) + *(RHregen_+i*NumNodes+j);
            }
        } else if(flavor == 3){
            std:string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = taurraysize[1];
            tau_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset(group_id, "tbarfull", tau_array_);
            RHregen_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, tau_array_);
            for (int i = 0; i<NumNodes; i++){
                for(int j=0; j<NumNodes;j++)
                *(RHSMatrix_+i*NumNodes+j) = *(RHSMatrix_+i*NumNodes+j) + *(RHregen_+i*NumNodes+j);
            }
        } else if(flavor == -1){
            double* glashow_total_ = get_glashow_total(energy_nodes);
            for (int i = 0; i < sarraysize[0]; i++){
                *(sigma_array_+i) = *(sigma_array_+i) + *(glashow_total_ + i)/2.; 
                *(RHSMatrix_ +i) = *(RHSMatrix_ +i) + *(glashow_partial_ + i)/2.;

            }
        }
        phi_0_ = (double *)malloc(NumNodes*sizeof(double));
        for (int i = 0; i < NumNodes; i++){
            *(phi_0_ + i) = std::pow(*(energy_nodes +i),(2-gamma));
        }
        for (int i = 0; i < sarraysize[0]; i++){
            *(RHSMatrix_+i*sarraysize[0]+i) = *(RHSMatrix_+i*sarraysize[0]+i) + *(sigma_array_+i);    
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

string nuFACE::getFilename() const {
    return newh5_filename;
}
