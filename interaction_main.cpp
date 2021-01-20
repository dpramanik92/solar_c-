//
//  test1.cpp
//  
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>

#include "interactions.hpp"

using namespace std;

int calculate_weighted(string,string,double,double);
int calculate_branching(string,string);


int main(int argc, char* argv[])
{
    string Type[3];
    Type[0] = "Scalar";
    Type[1] = "Pseudo";
    Type[2] = "Mixed";
    
    string Delta[2]={"05","95"};
    double delta_value[2] = {0.05,0.95};
    
    for(int i=0;i<3;i++)
    {
        string file_name2 = Type[i]+"_Br.dat";
        calculate_branching(file_name2,Type[i]);
        
        for(int j=0;j<2;j++)
        {
            string file_name1 = Type[i]+"_"+Delta[j]+"_w.dat";

            calculate_weighted(file_name1,Type[i],delta_value[j],10.0);
        }
        
    }
    
    
    return 0;
}

int calculate_weighted(string file_name,string which_type,double delta,double E_alpha)
{
    ofstream ofl;
    ofl.open(file_name);
    
    weighted_differential test(which_type);
    
    for(double E_be = -1.8;E_be<=1.0;E_be=E_be+0.01)
    {
        double E_beta = pow(10,E_be);
        
        double w = test.weighted_rate(delta,E_alpha,E_beta);
        
        
        ofl<<E_beta/E_alpha<<"\t"<<w<<endl;
    }
    
    ofl.close();
    
    return 0;
    
}

int calculate_branching(string file_name,string which_type)
{
    ofstream ofl;
    ofl.open(file_name);
    
    Branching_ratio test(which_type);
    
    for(double delta=-2.0;delta<=0.0;delta=delta+0.01)
    {
        double Delta = pow(10,delta);
        
        double Br = test.calc_branching(Delta);
        
        ofl<<Delta<<"\t"<<Br<<endl;
    }
    
    ofl.close();
    
    return 0;
}
