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

#include "visible_anti.hpp"
#include "invisible_e.hpp"
#include "convers_prob.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

 
    double E_max = 16.0;
    
    converse scalar;
    scalar.Tan_th12 = 0.5;
    scalar.Th13 = 0.0;
    scalar.Dm21 = 7.5e-5;
   // scalar.Delta = 0.05;
    scalar.conv = 1.0;
    scalar.L = 1.0;  /* Baseline in A.U. */
    scalar.Init_prob("Scalar",7);

    ofstream ofl;
    ofl.open("prob_test.dat");
    
    
    for(double E= -2;E<1.0;E=E+0.01)
    {
        double P = scalar.Calculate_probability(pow(10,E));
        ofl<<pow(10,E)<<"\t"<<P<<endl;
        
        
    }
    
    ofl.close();
    
    

    return 0;
}
