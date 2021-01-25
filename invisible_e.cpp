//
//  invisible_e.cpp
//  
//
//  Created by Dipyaman Pramanik on 24/01/21.
//

#include "invisible_e.hpp"

int invisible_electron_prob::Init_prob(std::string which_type,int num_of_params)
{
    prob.osc_params[0] = atan(Tan_th12);
    prob.osc_params[1] = Th13;
    prob.osc_params[2] = 0;
    prob.osc_params[3] = 0;
    prob.osc_params[4] = Dm21;
    prob.osc_params[5] = 0;
    
    prob.interpolate_data();
    
    return 0;
}

double invisible_electron_prob::Calculate_probability(double E)
{
    prob.Calculate_probability(E,0,0);
    
    double prob = square(cos(atan(Tan_th12)))*Propagation(E)*double(prob.E_p_day[4]);
    
    
    return prob;
}

double invisible_electron_prob::Propagation(double E)
{
    double exponent = L/(Tau*E);
    double res = exp(-exponent);
    
    return res;
}
