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

#include "KamLAND_anti.hpp"
#include "visible_anti.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

    
    double E_max = 10.0;
    
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
    _event.resolution[1] = 0.25;
    _event.resolution[2] = 0.0;
    
    _event.e_min = 1.8;
    _event.e_max = E_max;
    _event.n_bins = 100;
    
    
    _event.Set_probability_engine(scalar);
    _event.Init_evgen();
    _event.Set_fast_event_generator(SOL_YES,1.8,10,30);
    _event.Init_fast_generator();
    _event.generate_events();
    

    
    ofstream ofl;
    ofl.open("event_test1.dat");
    
    for(int i=0;i<_event.n_bins;i++)
    {
        ofl<<_event.bin_center[i]<<"\t"<<_event.Events[i]<<endl;
        cout<<_event.bin_center[i]<<"\t"<<_event.Events[i]<<endl;

        
    }
    
    ofl.close();
    
    
    

    return 0;
}
