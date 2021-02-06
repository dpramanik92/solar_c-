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
#include <string>

#include "visible_anti.hpp"
//~ #include "invisible_e.hpp"
//~ #include "convers_prob.hpp"
#include "numerical.hpp"
#include "read_files.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

	string which_type = string(argv[1]);
	double tau = atof(argv[2]);
	string tau_s = string(argv[2]);
	double delta = atof(argv[3]);
	string delta_s = string(argv[3]);
	
    
    visible_anti_prob scalar;
    scalar.Init_prob(which_type,7);
    scalar.Tan_Th12 = 0.5;
    scalar.Th13 = 0.0;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = delta;
    scalar.Tau = tau;
    scalar.L = 1.0;  /* Baseline in A.U. */
    //~ scalar.Init_prob("Scalar",E_max,7);

	string out_file = "prob_"+which_type+"_"+tau_s+"_"+delta_s+".dat";

	cout<<"Writing output to "<<out_file<<endl;

    ofstream ofl;
    ofl.open(out_file);
    
    
    for(double E=0.1;E<15.5;E=E+0.01)
    {
        double P = scalar.Calculate_decayed_flux(E);
        ofl<<E<<"\t"<<P<<endl;
      //  cout<<E<<"\t"<<P<<endl;

        
        
    }
    
    ofl.close();
    
    

    return 0;
}
