//
//  convers_prob.hpp
//  
//
//  Created by Dipyaman Pramanik on 26/01/21.
//

#ifndef convers_prob_hpp
#define convers_prob_hpp

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include "numerical.hpp"
#include "probability.hpp"

class converse
{
private:
    Probability prob;
public:
    int Init_prob(std::string, int);
    double Tan_th12,Th13,Dm21,L,conv;
    double Calculate_probability(double);
};

#endif /* convers_prob_hpp */
