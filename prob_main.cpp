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

#include "dec_probability.hpp"
// #include "visible_anti.hpp"
// #include "invisible_e.hpp"
// #include "convers_prob.hpp"
#include "numerical.hpp"
#include "read_files.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

	string which_type = string(argv[1]);
	double tau = atof(argv[2]);
	string tau_s = string(argv[2]);
	double delta = atof(argv[3]);
	string delta_s = string(argv[3]);
	
    
    
    dec_prob scalar;
    scalar.Init_prob(which_type,8);
    scalar.Tan_Th12 = square(tan(asin(sqrt(0.31))));
    scalar.Th13 = 0.0223;
    scalar.Dm21 = 7.39e-5;
    scalar.Delta = delta;
    scalar.Tau1 = tau;
    scalar.Tau2 = tau;
    scalar.L = 1.0;  /* Baseline in A.U. */
    //~ scalar.Init_prob("Scalar",E_max,7);

	string out_file1 = "flux_test_elec_"+which_type+"_"+tau_s+"_"+delta_s+".dat";

	string out_file2 = "flux_test_anti_"+which_type+"_"+tau_s+"_"+delta_s+".dat";
//	cout<<"Writing output to "<<out_file<<endl;

    ofstream ofl;
    ofl.open(out_file1);
    
    ofstream ofl1;
    ofl1.open(out_file2);

    
    for(double E = 0.0;E<15.0;E=E+0.01)
    {
        double P1 = scalar.Calculate_decayed_flux(E,0,1,2);
        double P2 = scalar.Calculate_decayed_flux(E,0,-1,2);
        ofl<<E<<"\t"<<P1<<endl;          
        ofl1<<E<<"\t"<<P2<<endl;
        
        
        
    }
    

    /*

    for(double E = -2;E<2;E=E+0.01)
    {
        double P = scalar.Calculate_decayed_flux(pow(10,E),0,1,1);

        ofl<<pow(10,E)<<"\t"<<P<<"\n";
    }

    ofl.close();
    
    */
    

    return 0;
}
