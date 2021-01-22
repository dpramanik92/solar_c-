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

    
    visible_anti_prob scalar("Scalar",10.0,7);
    scalar.Tan_Th12 = 0.5;
    scalar.Th13 = 0.0;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = 0.05;
    scalar.Tau = 1e-5;
    scalar.L = 1e6;
    
    for(double E = -2;E<log10(2);E=E+0.01)
    {
        double P = scalar.Calculate_decayed_flux(pow(10,E));
        
        cout<<pow(10,E)<<"\t"<<P<<endl;

    }
    
    
    

    return 0;
}
