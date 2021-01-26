//
//  invisible_e.cpp
//  
//
//  Created by Dipyaman Pramanik on 24/01/21.
//

#include "invisible_e.hpp"

int invisible_electron_prob::Init_prob(std::string which_type,int num_of_params)
{
    prob.osc_params.push_back(atan(Tan_th12));
    prob.osc_params.push_back(Th13);
    prob.osc_params.push_back(0);
    prob.osc_params.push_back(0);
    prob.osc_params.push_back(Dm21);
    prob.osc_params.push_back(0);
    
    prob.Init_probability_engine();
    
    return 0;
}

double invisible_electron_prob::Calculate_probability(double E)
{
    prob.Calculate_probability(E,0,0);
    
    double Prob = square(cos(atan(Tan_th12)))*Propagation(E)*prob.pday[4];
    
    return Prob;
}

double invisible_electron_prob::Propagation(double E)
{
    double Gamma_i = 1/(E*Tau);
    double res = (1.0-exp(-Gamma_i*L*4.96e-6));
    
    return res;
}
