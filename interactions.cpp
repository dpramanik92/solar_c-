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
        double factor3 = (1.0+sqr(delta)-E_beta/E_alpha-sqr(delta)*E_alpha/E_beta);
        
        
        if(which_type=="Scalar")
        {
            double factor1 = 1.0/((1.0-sqr(delta))*sqr(1.0+delta));
            
            w = factor1*factor2*factor3;
            
            //~ std::cout<<E_alpha/E_beta<<"\t"<<factor3<<"\n";
        }
        else if(which_type=="Pseudo")
        {
            double factor1 = 1.0/((1-sqr(delta))*sqr(1-delta));
            
            w = factor1*factor2*factor3;
            
        }
        else if(which_type=="Mixed")
        {
            double factor1 = 1.0/((1-sqr(delta*delta)));
            
            w = factor1*factor2*factor3;
        }
        else
        {
            std::cerr<<"ERROR!! Invalid option\n";
            
        }
        
    }
    
    
    
    return w;
}



Branching_ratio::Branching_ratio(std::string type_of_coupling)
{
    which_type = type_of_coupling;
    
    
}

double Branching_ratio::calc_branching(double delta)
{
    double Br;
    
    if(delta!=0)
    {
        
        if(which_type=="Scalar")
        {
            double term1 = (1.0+sqr(delta))/(2.0*sqr(1+delta));
            double term2 = 2.0*sqr(delta)*log(delta)/(sqr(1+delta)*(1-sqr(delta)));
            
            Br = term1 + term2;
        }
        else if(which_type=="Pseudo")
        {
            double term1 = (1.0+sqr(delta))/(2.0*sqr(1-delta));
            double term2 = 2.0*sqr(delta)*log(delta)/(sqr(1-delta)*(1-sqr(delta)));
            
            Br = term1 + term2;
            
        }
        else if(which_type=="Mixed")
        {
            double term1 = 0.5;
            double term2 = 2.0*sqr(delta)*log(delta)/(1-sqr(delta*delta));
            
            Br = term1 + term2;

        }
        
        
    }
    
    
    
    return Br;
}


double sqr(double x)
{
    return x*x;

}
