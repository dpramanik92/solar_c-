//
//  dec_probability.hpp
//  Solar_neutrino
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#ifndef dec_probability_hpp
#define dec_probability_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include "read_files.hpp"
#include "numerical.hpp"
#include "interactions.hpp"



typedef std::vector<double> vec;

#define SOL_NO 0
#define SOL_YES 1


class dec_prob
{
private:
    std::string which_type;
    vec *ff;
    double *rhosol;
    int init_interpolate_flux();
    Cubic_interpolator flux_interpolator;
    neutrino_data flux;
    double E_max;
    double integrate();
    double integrand(double);
    int prob_inside_sun(double);
    int prod_data_n_row;
    double Propagation(double E);
    double survival(double E);
    int calculate_total_flux();
    double x_max;
    double life;
    int particle;
    int _channel;


public:
    double Tan_Th12,Th13,Dm21,Delta,Tau1,Tau2;
    vec E_vec;
    int n_of_params;
    double *osc_params;
    std::string file_path;
    std::string outfile;
    double Energy;
    double P_ij;
    double pday[8],p_after_decay[8],pnight[8][9],pmed[8],pe1[8],pe2[8];
    vec Prob;
    double L;
    int Init_prob(std::string type_name,int);
    ~dec_prob();
    int interpolate_data();
    int Calculate_decayed_flux(vec);
    double Calculate_decayed_flux(double,int,int,int);
    int Decayed_Flux_curve(double,double,double,bool);
    int free_data();
    int regeneration;
    int Print_probability();
    int Wrap_oscparams();
    double integrate_flux();
    
    
};




#endif /* visiible_anti_hpp */
