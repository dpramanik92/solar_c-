//
//  interactions.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#include "interactions.hpp"

weighted_differential::weighted_differential(std::string type_of_coupling)
{
    which_type = type_of_coupling;
    
    
}

double weighted_differential::weighted_rate(double delta,double E_alpha,double E_beta)
{
    double w;
    
    if(delta!=0)
    {
        double factor2 = 1.0/E_alpha;
        double factor3 = (1+sqr(delta)-E_beta/E_alpha-sqr(delta)*E_alpha/E_beta);
        
        if(which_type=="Scalar")
        {
            double factor1 = 1.0/((1-sqr(delta))*sqr(1+delta));
            
            w = factor1*factor2*factor3;
        }
        else if(which_type=="Pseudo")
        {
            double factor1 = 1.0/((1-sqr(delta))*sqr(1-delta));
            
            w = factor1*factor2*factor3;
            
        }
        else if(which_type=="mixed")
        {
            double factor1 = 1.0/((1-sqr(delta*delta)));
            
            w = factor1*factor2*factor3;
        }
        
        
    }
    
    
    
    return w;
}


double sqr(double x)
{
    return x*x;

}
