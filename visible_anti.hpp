//
//  visiible_anti.hpp
//  Solar_neutrino
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#ifndef visiible_anti_hpp
#define visiible_anti_hpp

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


class visible_anti_prob
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
    int prob_inside_sun();
    int prod_data_n_row;
    double Propagation(double E);

public:
    double Tan_Th12,Th13,Dm21,Delta,Tau;
    vec E_vec;
    int n_of_params;
    double *osc_params;
    std::string file_path;
    std::string outfile;
    double Energy;
    double P_ij;
    double pday[8],p_after_decay[8],pnight[8][9],pmed[8];
    vec Prob;
    double L;
    int Init_prob(std::string type_name,double,int);
    ~visible_anti_prob();
    int interpolate_data();
    int Calculate_decayed_flux(vec);
    double Calculate_decayed_flux(double);
    int Decayed_Flux_curve(double,double,double,bool);
    int free_data();
    int Print_probability();
    int Wrap_oscparams();
    
    
    
};

#endif /* visiible_anti_hpp */
