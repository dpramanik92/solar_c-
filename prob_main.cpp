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

/*
    string file_prefix = string(argv[1]);
	string which_type = string(argv[2]);
	double tau = atof(argv[3]);
	string tau_s = string(argv[3]);
	double delta = atof(argv[4]);
	string delta_s = string(argv[4]);
*/	
    
    
    dec_prob scalar;
    scalar.Init_prob("Scalar",8);
    scalar.Tan_Th12 = square(tan(asin(sqrt(0.51))));
    scalar.Th13 = 0.0223;
    scalar.Dm21 = 7.39e-5;
    scalar.Delta = 0.9;
    scalar.Tau1 = 1e-4;
    scalar.Tau2 = 1e-4;
    scalar.L = 1.0;  /* Baseline in A.U. */
    //~ scalar.Init_prob("Scalar",E_max,7);
    
    double Integral = scalar.integrate_flux();

    //cout<<Integral<<endl;

/*
    

	string out_file1 =  file_prefix +"_elec_"+which_type+"_"+tau_s+"_"+delta_s+".dat";

	string out_file2 = file_prefix + "_anti_"+which_type+"_"+tau_s+"_"+delta_s+".dat";
//	cout<<"Writing output to "<<out_file<<endl;
*/
/*    ofstream ofl;
    ofl.open("Prob_nu3nu1_1e-4.dat");
    
  //  ofstream ofl1;
  //  ofl1.open(out_file2);

    
    for(double E = 0.0;E<15.0;E=E+0.1)
    {
        double P1 = scalar.Calculate_decayed_flux(E,0,1,3);
        double P2 = scalar.Calculate_decayed_flux(E,0,-1,3);
        cout<<E<<"\t"<<P1<<"\t"<<P2<<endl;          
        ofl<<E<<"\t"<<P1<<"\t"<<P2<<endl;
        
        
        
    }
    
*/

    for(double tau=-6;tau<1;tau=tau+0.1)
    {
        scalar.Tau1 = pow(10,tau);
        scalar.Tau2 = pow(10,tau);

        double P1 = scalar.Calculate_decayed_flux(7.6,0,-1,4);
        double P2 = scalar.Calculate_decayed_flux(7.6,0,-1,5);

        cout<<pow(10,tau)<<"\t"<<P1<<"\t"<<P2<<"\n";
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
