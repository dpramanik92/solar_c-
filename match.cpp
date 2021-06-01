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
#include "convers_prob.hpp"
#include "event.hpp"
//#include "dec_probability.hpp"
#include "numerical.hpp"



using namespace std;

typedef vector<double> vec;


int main(int argc, const char * argv[]) {

    
    //double val = atof(argv[1]);
    //string file_name = argv[2];
    
    double E_max = 10.0;
   
    string Experiment = "SK-IV";

    converse scalar;
    scalar.Tan_th12 = tan(33.56*(M_PI/180.0));
    scalar.Th13 = 0.02;
    scalar.Dm21 = 7.5e-5;
    scalar.conv = 3.6e-4;
    scalar.L = 1.0;  // Baseline in A.U. 
    scalar.Init_prob("Scalar",8);   
    scalar.oscillation = SOL_YES;    


    ofstream proba;

    proba.open("test_flux.dat");

    for(double E = 1.0;E<=16.0;E=E+0.01)
    {
        proba<<E<<"\t"<<scalar.Calculate_flux(E)<<endl;
    }

    proba.close();




     conv_event _event;
    
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

        _event.Normalization = 15.0*1e-41*1e31;

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
        double fiducial = 1.0; /* ktons */

        _event.efficiency = 0.85;

        _event.Normalization = 8.5*1e-41*1e32;

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
        double Np = 1.32e31; /* number of protons */

        _event.efficiency = 0.85;

        _event.Normalization = 30.0*1e-41;

        _event.Exposure = Time*Np;
     }
     else
     {
         cerr<<"ERROR!!Invalid Experiment name!!"<<endl;
         exit(-1);
     }


     for( int i=0; i<int(exp_data.data[0].size()); i++)
     {
         _event.manual_bins.push_back(exp_data.data[0][i]);
         _event.manual_f.push_back(exp_data.data[1][i]);
     }
    
    _event.Man_bins = SOL_YES;
    
    _event.Init_evgen();
    _event.res_stat = SOL_NO;
    _event.generate_events();  
    
//    cout<<"Writing output to "<<file_name<<endl;
    
    ofstream ofl;
    ofl.open("test_event_SK.dat");
    
    double count = 0;
    double bkg = 0;
     
    for(int i=0;i<_event.n_bins-2;i++)
    {
        ofl<<_event.bin_i[i]<<"\t"<<_event.Events[i]+bkg_data.data[2][i]<<"\t"<<bkg_data.data[2][i]<<endl;
        cout<<_event.bin_i[i]<<"\t"<<_event.Events[i]+bkg_data.data[2][i]<<"\t"<<bkg_data.data[2][i]<<endl;
        

        count = count + _event.Events[i];
        bkg = bkg + bkg_data.data[2][i];
    }
    
    ofl.close();
    
    cout<<count<<"\t"<<bkg<<endl;
    

    return 0;
}
