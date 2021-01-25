//
//  invisible_e.hpp
//  
//
//  Created by Dipyaman Pramanik on 24/01/21.
//

#ifndef invisible_e_hpp
#define invisible_e_hpp

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include "numerical.hpp"
#include "probability.hpp"

class invisible_electron_prob
{
private:
    double Propagation(double);
    Probability prob;
public:
    int Init_prob(std::string,int);
    double Tan_th12,Th13,Dm21,Tau,L;
    double Calculate_probability(double);
};

#endif /* invisible_e_hpp */
