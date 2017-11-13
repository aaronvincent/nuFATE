#include "nuFACE_secs.h"

double nuFACE_secs::readDoubleAttribute(hid_t object, std::string name){
        double target;
        hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
        if(status<0)
            throw std::runtime_error("Failed to read attribute '"+name+"'");
        H5Aclose(attribute_id);
        return target;
}

unsigned int nuFACE_secs::readUIntAttribute(hid_t object, std::string name){
  unsigned int target;
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return target;
}


double* logspace(double Emin,double Emax,unsigned int div){
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
        double* logpoints_ = logpoints;
        return logpoints_;
}

double* nuFACE::get_glashow_total(unsigned int NumNodes, double* energy_nodes){
    double GF = 1.16e-5;
    double hbarc = 1.97e-14;
    double GW = 2.085;
    double MW = 80.385e0;
    double mmu = 0.106e0;
    double me = 511.e-6;
    double pi = 3.14159265358979323846;
    glashow_total_ = (double *)malloc(NumNodes*sizeof(double));

    for(int i=0; i<NumNodes; i++){
        double x = *(glashow_total_+i);
        *(glashow_total_ +i) = 2.*me*(*(energy_nodes+i));
        *(glashow_total_ +i) = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return glashow_total_;
}

double* nuFACE_secs::get_RHS_matrices(unsigned int NumNodes, double* energy_nodes, double* sigma_array_, double* sig3_array_, double* dxs_array_, double* sec_array_, double* regen_array_){

    DeltaE_ = (double *)malloc(NumNodes*sizeof(double));
    for(int i = 0; i < NumNodes-1;i++){
        *(DeltaE_ + i) = log10(*(energy_nodes+i+1)) - log10(*(energy_nodes+i));
    }

    RHSMatrix1_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    RHSMatrix2_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    RHSMatrix3_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    RHSMatrix4_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
    RHSMatrix_ =  (double *)malloc(4*NumNodes*NumNodes*sizeof(double));
        for(int i=0; i<NumNodes, i++){
            for(int j=0; j<NumNodes, j++){
                *(RHSMatrix1_+i*NumNodes+j) = 0.;
                *(RHSMatrix2_+i*NumNodes+j) = 0.;
                *(RHSMatrix3_+i*NumNodes+j) = 0.;
                *(RHSMatrix4_+i*NumNodes+j) = 0.;
            }
        }
    //matrix 1: nue or numu NC:
    for (int i=0; i<NumNodes, i++){
        for(int j=i+1; j<NumNodes, j++){
            *(RHSMatrix1_+i*NumNodes+j) = *(DeltaE_+j-1) * *(dxs_array_+j*NumNodes+i) * (1./ *(energy_nodes+j));
            *(RHSMatrix2_+i*NumNodes+j) = *(DeltaE_+j-1) * *(sec_array_+j*NumNodes+i) * (1./ *(energy_nodes+j)) * (*(energy_nodes+i)* *(energy_nodes+i));
            *(RHSMatrix4_+i*NumNodes+j) = *(DeltaE_+j-1) * *(dxs_array_+j*NumNodes+i) * *(regen_array_+j*NumNodes+i) * (1./ *(energy_nodes+j)) * (*(energy_nodes+i)* *(energy_nodes+i));
        }
    }
    for (int i = 0; i < NumNodes; i++){
        *(RHSMatrix1_+i*NumNodes+i) = *(RHSMatrix1_+i*NumNodes+i) - *(sigma_array_+i);    
        *(RHSMatrix4_+i*NumNodes+i) = *(RHSMatrix4_+i*NumNodes+i) - *(sig3_array_+i);
    }

    //matrix2:
    /*for (int i=0; i<NumNodes, i++){
        for(int j=i+1; j<NumNodes, j++){
            *(RHSMatrix2_+i*NumNodes+j) = *(DeltaE_+j-1) * *(sec_array_+j*NumNodes+i) * (1./ *(energy_nodes+j)) * (*(energy_nodes+i)* *(energy_nodes+i));
        }
    }
    //matrix4 (matrix3 is zero)
    for (int i=0; i<NumNodes, i++){
        for(int j=i+1; j<NumNodes, j++){
            *(RHSMatrix4_+i*NumNodes+j) = *(DeltaE_+j-1) * *(dxs_array_+j*NumNodes+i) * *(regen_array_+j*NumNodes+i) * (1./ *(energy_nodes+j)) * (*(energy_nodes+i)* *(energy_nodes+i));*/
    int rsize = 2*NumNodes;
    
    for (int i = 0; i<rsize; i++){
        if(i<NumNodes){
            for (int j = 0; j<rsize;j++){
                if (j<NumNodes){
                    *(RHSMatrix_+i*rsize+j) = *(RHSMatrix1_+i*NumNodes+j);
                    else {
                        *(RHSMatrix_+i*rsize+j) = *(RHSMatrix3_+i*NumNodes+j);
                    }
                }     
            }
        } else {
            for (int j = 0; j<rsize;j++){
                if (j<NumNodes){
                    *(RHSMatrix_+i*rsize+j) = *(RHSMatrix2_+i*NumNodes+j);
                } else {
                        *(RHSMatrix_+i*rsize+j) = *(RHSMatrix4_+i*NumNodes+j);
                 }
            }     
        }
    }

    return RHSMatrix_;
}


result nuFACE_secs:get_eigs(int flavor, double gamma, std::string h5_filename) {

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
        }

        if (flavor > 0){
            
            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
            sig3_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "nutauxs", sig3_array_);
            
            hsize_t dxarraysize[2];
            group_id = H5Gopen(root_id, grpdiff.c_str(), H5P_DEFAULT);        
            H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);    
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "dxsnu", dxs_array_);

            std::string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = tauarraysize[1];
            regen_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "tfull", regen_array_);

            hsize_t secarraysize[2];
            H5LTget_dataset_info(group_id,"secfull", secarraysize,NULL,NULL);
            sec_array_ = (double *)malloc(secarraysize[0]*secarraysize[1]*sizeof(double));
            H5LTread_dataset_double(group_id, "secfull", sec_array_);

        } else {

            hsize_t sarraysize[1];
            H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
            sig3_array_ = (double *)malloc(sarraysize[0]*sizeof(double));
            H5LTread_dataset_double(group_id, "nutaubarxs", sig3_array_);

            hsize_t dxarraysize[2];
            group_id = H5Gopen(root_id, grpdiff.c_str(), H5P_DEFAULT);        
            H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
            size_t dim1 = dxarraysize[0];
            size_t dim2 = dxarraysize[1];
            dxsdim[0] = dxarraysize[0];
            dxsdim[1] = dxarraysize[1];
            dxs_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "dxsnubar", dxs_array_);

            std::string grptau = "/tau_decay_spectrum";
            group_id = H5Gopen(root_id, grptau.c_str(), H5P_DEFAULT);
            hsize_t tauarraysize[2];
            H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
            size_t dim1 = tauarraysize[0];
            size_t dim2 = tauarraysize[1];
            regen_array_ = (double *)malloc(dim1*dim2*sizeof(double));
            H5LTread_dataset_double(group_id, "tbarfull", regen_array_);

            hsize_t secarraysize[2];
            H5LTget_dataset_info(group_id,"secbarfull", secarraysize,NULL,NULL);
            sec_array_ = (double *)malloc(secarraysize[0]*secarraysize[1]*sizeof(double));
            H5LTread_dataset_double(group_id, "secbarfull", sec_array_);

        }

        RHSMatrix_ = get_RHS_matrices(NumNodes, energy_nodes, sigma_array_, sig3_array_, dxs_array_, sec_array_, regen_array_);

        if (flavor == -1){
            total_glashow_ = get_glashow_total(NumNodes, energy_nodes);
            glashow_piece_ = (double *)malloc(NumNodes*NumNodes*sizeof(double));
            glashow_piece_ =  get_glashow_partial(NumNodes, energy_nodes);

            for (int i = 0; i<NumNodes;i++){
                *(glashow_piece_+i*NumNodes+i) = (*(glashow_piece_+i*NumNodes+i) + *(total_glashow_+i))/2.;
                *(phi_0_ + i) = std::pow(*(energy_nodes +i),(2-gamma)); 
            }
            for(int i=0; i<NumNodes;i++){
                for(int j =0; j<NumNodes;j++){
                    *(RHSMatrix_+i*NumNodes+j) = *(RHSMatrix_+i*NumNodes+j) + *(glashow_piece_+i*NumNodes+j);
                }
            }  
        }

        phi_0_ = (double *)malloc(rsize*sizeof(double));
        for (int i = 0; i < NumNodes; i++){
            *(phi_0_ + i) = std::pow(*(energy_nodes +i),(2-gamma));
            *(phi_0_ + i + NumNodes) = std::pow(*(energy_nodes +i),(2-gamma));
        }

        gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_, rsize, rsize);
        gsl_vector_complex *eval = gsl_vector_complex_alloc (rsize);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (rsize, rsize);
        gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (rsize);

        gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

        int s;
        gsl_vector *ci = gsl_vector_alloc(rsize);
        gsl_permutation *p = gsl_permutation_alloc(rsize);
        
        gsl_linalg_LU_decomp (&m.matrix, p, &s);
        gsl_linalg_LU_solve (&m.matrix, p, phi_0_, ci);

        //free unneeded memory
        gsl_permutation_free (p);
        gsl_eigen_nonsymmv_free (w);

        struct result r1;
        r1.eval = eval;
        r1.evec = evec;
        r1.ci = ci;
        r1.energy_nodes = energy_nodes;
        r1.phi_0_ = phi_0_;

        return r1;

}
