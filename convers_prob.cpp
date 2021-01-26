//
//  convers_prob.cpp
//  
//
//  Created by Dipyaman Pramanik on 26/01/21.
//

#include "convers_prob.hpp"

int converse::Init_prob(std::string which_type,int num_of_params)
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

double converse::Calculate_probability(double E)
{
    prob.Calculate_probability(E,0,0);
    
    double Prob = square(cos(atan(Tan_th12)))*prob.pday[4]*conv;
    
    return Prob;
    
}
