//
//  chisqmin.hpp
//  
//
//  Created by Dipyaman Pramanik on 24/01/21.
//

#ifndef chisqmin_hpp
#define chisqmin_hpp


#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include "read_files.hpp"
#include "probability.hpp"
#include "numerical.hpp"
#include "visible_anti.hpp"

class chisq_real_one
{
private:
    visible_anti_prob Proba_engine;
    file_reader exp_data,bkg_data;
    Event_generator _event;
public:
    int Init_chisq_calculator(Event_generator,std::string,std::string);
    int Read_data(std::string);
    int Read_background(strd::string);
  //  int Set_Systematic_func(sys_defn);
    double Calc_chi2_nosys(vec params);
    double Calc_chi2_sys(vec params);
}


double poiss_likelihood(double,double);

/*
class sys_defn
{
private:
    
public:
    double sys;
    
}
 */

#endif /* chisqmin_hpp */
