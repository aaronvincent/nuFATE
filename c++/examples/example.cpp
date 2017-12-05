#include <iostream>
#include "nuFATE.h"

int main(){
    int flavor = -1;
    double gamma = 2.2;
    std::string file = "../../data/NuFATECrossSections.h5";
    nufate::nuFACE object(flavor, gamma, file);

    nufate::Result result;

    result = object.get_eigs();
    double NumNodes;
    NumNodes = object.getNumNodes();
    for(int i =0; i<NumNodes; i++){
        double* eval = result.eval;
        std::cout << *(eval+i) << std::endl;
    }

    return 0;
}
