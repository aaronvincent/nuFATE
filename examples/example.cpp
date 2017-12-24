#include <iostream>
#include "nuFATE.h"


int main(){

/*
 * flavor: Select a neutrino flavor (1, 2, 3, for electron, muon, and tau, respectively. the number is negative for anti-neutrinos)
 *  gamma: Spectral index of the incoming flux
 *  file: path to file containing the cross sections.
*/
    int flavor = -2;
    double gamma = 2.2;
    std::string file = "/Users/wipacuser/Code/nuFATE/nuFATE/nuFATE/c++/examples/NuFATECrossSections.h5";

    //Initialize an instance of the nuFATE class with these three parameters.

    nufate::nuFATE object(flavor, gamma, file);

    //Result is a struct that stores the solution to the cascade equation.

    nufate::Result result;

    //get_eigs() solves the cascade equation
    result = object.getEigensystem();
    int NumNodes;
    NumNodes = object.getNumNodes();

    std::vector<double> eval = result.eval;
    std::shared_ptr<double> evec = result.evec;
    std::vector<double> ci = result.ci;
    std::vector<double> energy_nodes = result.energy_nodes_;
//    for(int i =0; i<NumNodes; i++){
//      for(int j=0; j<NumNodes;j++){
        //int i1 = 2*i;
        //int j1 = 2*j;
//        std::cout << *(evec+i1*NumNodes+j1) << std::endl;
//      }
//    }
    double Na = 6.0221415e23;
    double zenith = 2.2689280276;
    double t;

    std::vector<double> abs;
    std::vector<double> phi_sol;
//Calculate earth column density for a given zenith
    t = object.getEarthColumnDensity(zenith) * Na;
//    std::cout << t << std::endl;
    for(int i=0; i<NumNodes; i++){
      //steps of 2i because gsl saves evals and evecs as a complex and a real piece (complex piece is 0 here)
      //double w = *(eval+2*i);
      //abs.push_back(exp(-t*w));
      double sum = 0.;
      for (int j=0; j<NumNodes;j++){
        abs.push_back(ci[j] * exp(-t*eval[j]));
        sum+= abs[j] *  *(evec.get()+i*NumNodes+j) * (std::pow(energy_nodes[i],-2) / std::pow(energy_nodes[i],-gamma));
      }
      phi_sol.push_back(sum);
    }

    std::cout << "Solution = " << std::endl;
    for(int i =0; i<NumNodes; i++){
      std::cout << phi_sol[i] << std::endl;
    }

   return 0;
}
