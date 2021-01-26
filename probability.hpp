//
//  probability.hpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 21/11/20.
//

#ifndef probability_hpp
#define probability_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include "read_files.hpp"

typedef std::vector<double> vec;

class Probability
{
private:
    regeneration *exp;
    int read_regen();
    int prod_data_n_row;
    double *pmed;
    
    
public:
    vec E_vec,E_p_day[8],E_p_night[8][9];
    vec osc_params;
    int fin_flav;
    int print_prob;
    std::string file_path;
    std::string outfile;
    double Energy;
    double *rhosol;
    vec *ff;
    vec Prob;
    double *pday,**pnight;
    vec P_day;
    ~Probability();
    int Init_probability_engine();
    int interpolate_data();
    int Calculate_probability(double,int,int);
    int Calculate_probability(vec,int);
    int Probability_curve(double,double,double,int,bool);
    int prob_inside_sun();
    int regeneration_earth();
    int free_data();
    int Print_probability();

};



#endif /* probability_hpp */
