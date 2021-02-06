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
#include "visible_anti.hpp"
#include "numerical.hpp"

using namespace std;


int main(int argc, const char * argv[]) {

    
    //double val = atof(argv[1]);
    //string file_name = argv[2];
    
    double E_max = 10.0;
    
    visible_anti_prob scalar;
    scalar.Tan_Th12 = tan(33.56*(M_PI/180.0));
    scalar.Th13 = 0.02;
    scalar.Dm21 = 7.5e-5;
    scalar.Delta = 0.96;
    scalar.Tau = 1e-6;
    scalar.L = 1.0;  /* Baseline in A.U. */
    scalar.Init_prob("Scalar",7);
    
    //~ conv_event _event;
    
    //~ _event.efficiency = 0.9;
    //~ _event.resolution[0] = 0.0;
    //~ _event.resolution[1] = 0.06;
    //~ _event.resolution[2] = 0.0;
    
    //~ _event.e_min = 1.8;
    //~ _event.e_max = E_max;
    //~ _event.n_bins = 100;
    
    //~ _event.Normalization = 1.35e-3;
    //~ _event.Exposure = 4.53;
    
    //~ _event.Set_probability_engine(scalar);
    
    //~ file_reader exp_data,bkg_data;


	//~ string data_file = "exp_data/KamLAND_data.dat";
	//~ string bkg_file = "exp_data/KamLAND_bkg.dat";
	
    
    //~ exp_data.read_file(data_file);
    //~ bkg_data.read_file(bkg_file);
    
    
    double Time = 2485;  /* days */
    Time = Time*86400;    /* sec */
    double Np = 1.32e31;  /* number of protons */
    
    
    
	Event_generator _event;
    
    _event.efficiency = 0.85;
    _event.resolution[0] = 0.0;
    _event.resolution[1] = 0.06;
    _event.resolution[2] = 0.0;
    
    _event.e_min = 1.8;
    _event.e_max = E_max;
    _event.n_bins = 100;
    
    _event.Normalization = 5.46e4*1e-47*1.53;
    _event.Exposure = Time*Np;
    
    _event.Set_probability_engine(scalar);
    
    file_reader exp_data,bkg_data;


	string data_file = "exp_data/Borexino_data.dat";
	string bkg_file = "exp_data/Borexino_bkg.dat";
	
    
    exp_data.read_file(data_file);
    bkg_data.read_file(bkg_file);
    
    
    
    for(int i=0;i<int(exp_data.data[0].size());i++)
    {
        _event.manual_bins.push_back(exp_data.data[0][i]);
		_event.manual_f.push_back(exp_data.data[1][i]);
    }
    
    
    
    _event.Man_bins = SOL_YES;
    
    _event.Init_evgen();
    _event.res_stat = SOL_YES;
    _event.generate_events();
    
//    cout<<"Writing output to "<<file_name<<endl;
    
    ofstream ofl;
    ofl.open("borexino_anti_e_scalar.dat");
    
    double count = 0;
     
    for(int i=0;i<_event.n_bins-2;i++)
    {
        ofl<<_event.bin_center[i]<<"\t"<<_event.Events[i]+bkg_data.data[2][i]<<"\t"<<bkg_data.data[2][i]<<endl;

        count = count + _event.Events[i];
    }
    
    ofl.close();
    
    cout<<count<<endl;
    

    return 0;
}
