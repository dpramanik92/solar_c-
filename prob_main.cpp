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

    
    double E_max = 16.0;
    
    visible_anti_prob scalar;
    scalar.Init_prob("Scalar",E_max,7);
    scalar.Tan_Th12 = 0.5;
    scalar.Th13 = 0.0;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = 0.05;
    scalar.Tau = 1e-5;
    scalar.L = 1.0;  /* Baseline in A.U. */
    
    ofstream ofl;
    ofl.open("prob_test.dat");
    
    
    for(double E= 1.8;E<16.0;E=E+0.01)
    {
        double P = scalar.Calculate_decayed_flux(E);
        ofl<<E<<"\t"<<P<<endl;
        
        
    }
    
    ofl.close();
    
    

    return 0;
}
