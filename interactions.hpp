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

typedef vector<double> vec;

double sqr(double );

class weighted_differential
{
private:
    std::string which_type;
    
public:
    weighted_differential(std::string);
    double weighted_rate(double,double,double);
    
}

#endif /* interactions_hpp */
