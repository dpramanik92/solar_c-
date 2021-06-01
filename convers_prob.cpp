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
    init_interpolate_flux();

    prob.Init_probability_engine();
    
    return 0;
    
    
}

int converse::init_interpolate_flux()
{
    flux.read_flux("flux/b8spec-2006.dat");

    vec flux_x,flux_y;

    for(int i=0;i<flux.Flux.data[0].size();i++)
    {
        flux_x.push_back(flux.Flux.data[0][i]);
        flux_y.push_back(flux.Flux.data[1][i]);
    }

    flux_interpolator.set_cubic_spline(flux_x,flux_y,false);

    return 0;


}



double converse::Calculate_flux(double E)
{
    double f = flux_interpolator.interpolate(E)*Calculate_probability(E);
        

    return f;
}


double converse::Calculate_probability(double E)
{
    double Prob;

    if(oscillation==SOL_YES)
    {   
        prob.Calculate_probability(E,0,0);
    
        Prob = square(cos(atan(Tan_th12)))*prob.pday[4]*conv;
    

    }
    else if(oscillation==SOL_NO)
    {
        Prob = conv;
    }
    else
    {
        std::cerr<<"ERROR!!INVALID OPTION!"<<std::endl;
        exit(-1);
    }

    return Prob;
    
}
