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
    double *rhosol
    
public:
    vec E_vec;
    vec osc_params;
    std::string file_path;
    std::string outfile;
    double Energy;
    double P_ij;
    double pday[8],p_after_decay[8],pnight[8][9];
    double integrand(double);
    vec Prob;
    double Prob;
    double L;
    visible_anti_prob(std::string type_name);
    ~visible_anti_prob();
    int interpolate_data();
    int Calculate_probability(vec);
    double Calculate_probability(double);
    double Propagation(double)
    int Probability_curve(double,double,double,bool);
    int prob_inside_sun();
    int free_data();
    int Print_probability();
    
    
    
};

#endif /* visiible_anti_hpp */
