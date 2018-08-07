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
    bool include_secondaries = false;
    std::string file = "/Users/wipacuser/Code/nuFATE/nuFATE/nuFATE/c++/examples/NuFATECrossSections.h5";
    //Initialize an instance of the nuFATE class with these parameters.
    nufate::nuFATE object(flavor, gamma, file, include_secondaries);
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
    std::vector<double> phi_0 = result.phi_0_;
    //Calculate earth column density for a given zenith
    double Na = 6.0221415e23;
    double zenith = 2.2689280276;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;

    std::vector<double> abs;
    std::vector<double> phi_sol;
    //Get Attenuation!
    if(not include_secondaries){
      //Without Secondaries
      for(int i=0; i<NumNodes; i++){
        double sum = 0.;
        for (int j=0; j<NumNodes;j++){
          abs.push_back(ci[j] * exp(-t*eval[j]));
          sum+= abs[j] *  *(evec.get()+i*NumNodes+j) * (std::pow(energy_nodes[i],-2) / std::pow(energy_nodes[i],-gamma));
        }
        phi_sol.push_back(sum);
      }
      //Print Solution
      std::cout << "Solution = " << std::endl;
      for(int i =0; i<NumNodes; i++){
        std::cout << phi_sol[i] << std::endl;
      }

    } else{
      //With Secondaries
        int rsize = 2*NumNodes;
        for(int i=0; i<rsize; i++){
          double sum = 0.;
          abs.clear();
          for (int j=0; j<rsize;j++){
            abs.push_back(ci[j] * exp(-t*eval[j]));
            sum+= (abs[j] * *(result.evec.get()+i*rsize+j)) / phi_0[i];
          }
        phi_sol.push_back(sum);
        }

        //Print Solution
        std::cout << "Solution Including secondaries= " << std::endl;
        for(int i =0; i<NumNodes; i++){
          std::cout << phi_sol[i] << std::endl;
        }
    }

    return 0;
}
