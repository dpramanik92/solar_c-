//
//  main.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 20/11/20.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>

#include "read_files.hpp"
#include "probability.hpp"
#include "numerical.hpp"
#include "KamLAND_anti.hpp"
#include "visible_anti.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

    
    double E_max = 16.0;
    
    visible_anti_prob scalar;
    scalar.Init_prob("Scalar",E_max,7);
    scalar.Tan_Th12 = 0.5;
    scalar.Th13 = 0.0;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = 0.05;
    scalar.Tau = 1e-5;
    scalar.L = 1.0;  /* Baseline in A.U. */
    
    Event_generator _event;
    
    _event.efficiency = 0.9;
    _event.resolution[0] = 0.0;
    _event.resolution[1] = 0.35;
    _event.resolution[2] = 0.0;
    
    _event.e_min = 1.8;
    _event.e_max = E_max;
    _event.n_bins = 50;
    
    
    _event.Set_probability_engine(scalar);
    _event.Init_evgen();
    _event.generate_events();
    
    
    

    return 0;
}
