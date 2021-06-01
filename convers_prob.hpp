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

#define SOL_YES 1
#define SOL_NO 0

class converse
{
private:
    Probability prob;
    std::string flux_file;
    neutrino_data flux;
    Cubic_interpolator flux_interpolator;
public:
    int oscillation;
    int Init_prob(std::string, int);
    double Tan_th12,Th13,Dm21,L,conv;
    int init_interpolate_flux();
    double Calculate_probability(double);
    double Calculate_flux(double E);
};

#endif /* convers_prob_hpp */
