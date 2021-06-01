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
#include "event.hpp"
#include "dec_probability.hpp"
#include "numerical.hpp"



using namespace std;

typedef vector<double> vec;


int main(int argc, const char * argv[]) {

    
    
    double E_max = 10.0;

    string Experiment = "SK-IV";
    
    dec_prob scalar;
    scalar.Tan_Th12 = tan(33.56*(M_PI/180.0));
    scalar.Th13 = 0.02;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = 0.0;
    scalar.Tau1 = 1.2e-3;
    scalar.Tau2 = 1.2e-3;
    scalar.L = 1.0;  // Baseline in A.U. 
    scalar.Init_prob("Scalar",8);   
    

     Event_generator _event;
    
     _event.efficiency = 0.05;
     _event.resolution[0] = 0.0;
     _event.resolution[1] = 0.06;
     _event.resolution[2] = 0.0;
    
     _event.e_min = 1.8;
     _event.e_max = E_max;
     _event.n_bins = 100;
    

    
     _event.Set_probability_engine(scalar);
    
     file_reader exp_data,bkg_data;


      

     if(Experiment=="SK-IV")
     {
        string data_file = "exp_data/SuperK-IV_data.dat";
        string bkg_file = "exp_data/SuperK-IV_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2970.1; /* days */
        Time = Time*86400;
        double fiducial = 22.5; /* ktons */

        _event.efficiency = 0.05;

        _event.Normalization = 20.0*1e-41*1e31;

        _event.Exposure = Time*fiducial;

     }
     else if(Experiment=="Kamland")
     {

        string data_file = "exp_data/KamLAND_data.dat";
        string bkg_file = "exp_data/KamLAND_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2343; /* days */
        Time = Time*86400;
        double fiducial = 0.705; /* ktons */

        _event.efficiency = 0.85;

        _event.Normalization = 26.5*1e-41*1e31;

        _event.Exposure = Time*fiducial;
     }
     else if(Experiment=="Borexino")
     {
         

        string data_file = "exp_data/Borexino_data.dat";
        string bkg_file = "exp_data/Borexino_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2485; /* days */
        Time = Time*86400;
        double Np = 1.32e30; /* ktons */

        _event.efficiency = 0.92;

        _event.Normalization = 12.5*1e-41;

        _event.Exposure = Time*Np;
     }
     else
     {
         cerr<<"ERROR!!Invalid Experiment name!!"<<endl;
         exit(-1);
     }



    
    for(int i=0;i<int(exp_data.data[0].size());i++)
    {
        _event.manual_bins.push_back(exp_data.data[0][i]);
		_event.manual_f.push_back(exp_data.data[1][i]);
    }
    
    
    
    _event.Man_bins = SOL_YES;
    
    _event.Init_evgen(0,-1,1);
    _event.res_stat = SOL_NO;
    _event.generate_events();  
    
    
    ofstream ofl;
    ofl.open("test_event.dat");
    
    double count = 0;
     
    for(int i=0;i<_event.n_bins-2;i++)
    {
        ofl<<_event.bin_i[i]<<"\t"<<_event.Events[i]+bkg_data.data[2][i]<<"\t"<<bkg_data.data[2][i]<<endl;
        

        count = count + _event.Events[i]+bkg_data.data[2][i];
    }
    
    ofl.close();
    
    cout<<count<<endl;
    

    return 0;
}
