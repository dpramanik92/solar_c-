//
//  interactions.hpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#ifndef interactions_hpp
#define interactions_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

typedef std::vector<double> vec;

double sqr(double );

class weighted_differential
{
private:
    std::string which_type;
    
public:
    weighted_differential(std::string);
    double weighted_rate(double,double,double);
    
};

class Branching_ratio
{
private:
    std::string which_type;
    
public:
    Branching_ratio(std::string);
    double calc_branching(double);
    
};

#endif /* interactions_hpp */
