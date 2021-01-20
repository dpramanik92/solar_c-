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

using namespace std;


int main(int argc, const char * argv[]) {

    
  /*  vec true_params;

    double th12,th13,th23,del13,dm21,dm31;
    
    
    
    th12 = atan(sqrt(0.44));//33.46*(M_PI/180.0);
    th13 = asin(sqrt(0.023));//8.5*(M_PI/180.0);
    th23 = 45.0*(M_PI/180.0);
    del13 = -M_PI/2.0;
    dm21 = 8.0e-5;
    dm31 = 2.52e-3;
    
    
    true_params.push_back(th12);
    true_params.push_back(th13);
    true_params.push_back(th23);
    true_params.push_back(del13);
    true_params.push_back(dm21);
    true_params.push_back(dm31);
    
    double E_l = -2;
    double E_h = log10(200);
    double E_step = 0.01;
    
  
    Probability Prob;
    Prob.osc_params = true_params;
    Prob.file_path = "Probability_data/";
    Prob.outfile = "test_prob2";
    Prob.print_prob = 1;
    
    Prob.Probability_curve(E_l, E_h, E_step,0,true);

*/
    
    Event_generator Kamland;
    Kamland.efficiency = 0.9;
    Kamland.resolution[0] - 0.0;
    Kamland.resolution[1] - 0.35;
    Kamland.resolution[0] - 0.0;
    
    Kamland.e_min = 1.8;
    Kamland.e_max = 10.0;
    Kamland.n_bins = 50;
    
    Kamland.Init_evgen();
    Kamland.generate_events();
    
    ofstream ofl;
    ofl.open("test.dat");
    
    for(int i=0;i<Kamland.n_bins;i++)
    {
        ofl<<Kamland.bin_center[i]<<"\t"<<Kamland.Events[i]<<endl;
        
    }
    
    ofl.close();
    

    return 0;
}
