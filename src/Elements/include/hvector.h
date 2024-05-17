//
// Created by Yanjia on 2024/5/9.
//

#ifndef INC_05_OVERALLPROCESS_HVECTOR_H
#define INC_05_OVERALLPROCESS_HVECTOR_H


#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

using namespace arma;
using namespace std;

struct HVector {

    mat H3i;
    mat H2i;
    mat H1i;

};

HVector hVector(imat Eid,mat Nid,double nn);
arma::mat nodeNormal(imat Eid,mat Nid,double nn);



#endif //INC_05_OVERALLPROCESS_HVECTOR_H
