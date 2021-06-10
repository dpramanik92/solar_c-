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
#include <stdlib.h>
#include <string.h>
#include <string>
#include "event.hpp"
#include "dec_probability.hpp"
#include "chisqmin.hpp"
#include "numerical.hpp"
#include "minimizer.hpp"


using namespace std;

typedef vector<double> vec;

class chi_sq
{
    private:
        dec_prob Prob;
        Event_generator _event;
        string Experiment;
        int Channel;
        string which_type;
        chisq::SK_IV chi2;

        vec start_values;
        chisq::sys_minimizer<chisq::SK_IV> _minimizer;


    public:
        void solGetExpParams(string,string ,int);
        void InitChi2Engine();
        double solCalculateChi2(double*);


    private:
};

void chi_sq::solGetExpParams(string _expname,string _which_type,int _chan)
{

    Experiment = _expname;
    which_type = _which_type;
    Channel = _chan;

}

void chi_sq::InitChi2Engine()
{
    double E_max = 10.0;

    Prob.Init_prob(which_type,8);


     _event.resolution[0] = 0.0;
     _event.resolution[1] = 0.06;
     _event.resolution[2] = 0.0;
    
     _event.e_min = 1.8;
     _event.e_max = E_max;
     _event.n_bins = 100;
    

    
     _event.Set_probability_engine(Prob);
    
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
        double Np = 1.32e31; /* ktons */

        _event.efficiency = 0.85;

        _event.Normalization = 30.0*1e-41;

        _event.Exposure = Time*Np;
     }
     else
     {
         cerr<<"ERROR!!Invalid Experiment name!!"<<endl;
         exit(-1);
     }



     _event.res_stat = SOL_NO;


     for(int i=0;i<16;i++)
     {
         start_values.push_back(0.0);
     }

     chi2.fin_flav = 0;
     chi2.channel = Channel;
     chi2.particle = -1;

     chi2.Set_sys(0.2,0.2,0.05,0.05);
     chi2.Init(_event,exp_data,bkg_data,SOL_NO);
     chi2.Set_binned_systematics(SOL_NO);
     chi2.statistics(start_values);
    


    



}

double chi_sq::solCalculateChi2(double* params)
{  
     
    vec fit_params;

    double tan_th12 = pow(tan(params[0]),2);

    fit_params.push_back(tan_th12);
    fit_params.push_back(params[1]);
    fit_params.push_back(params[2]);
    fit_params.push_back(params[3]);
    fit_params.push_back(params[4]);
    fit_params.push_back(params[5]);
    fit_params.push_back(params[6]);
    fit_params.push_back(params[7]);

    double res = _minimizer.Minimize(chi2,start_values,fit_params);

    return res;

}

extern "C"
{
    chi_sq chi2;
    void solInitExp(char* exp,char* which,int chan)
    {
        
        string _exp = "Kamland"; 
        string _which = "Scalar";
        chi2.solGetExpParams(exp,which,chan);
        chi2.InitChi2Engine();

    }
    double solChi2(double* fit_params)
    {
        double f = chi2.solCalculateChi2(fit_params);
        return f;
    }
    

    
}

int main()
{
  


}

