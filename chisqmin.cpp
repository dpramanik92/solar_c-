//
//  chisqmin.cpp
//  
//
//  Created by Dipyaman Pramanik on 24/01/21.
//

#include "chisqmin.hpp"

int chisq_real_one::Init_chisq_calculator(Event_generator event,std::string data_file,std::string,bkg_file,visible_anti_prob Proba)
{
    _event = event;
    
    exp_data.read_file(data_file);
    bkg_data.read_file(bkg_file);
    
    Proba_engine = Proba;
    _event.Set_probability_engine(Proba);
    _event.Init_evgen();
    
    for(int i=0;i<int(exp_data.data[0].size());i++)
    {
        _event.manual_bins.(exp_data.data[0][i]);
        
    }
    
    
    
    _event.Set_fast_event_generator(SOL_YES,exp_data.data[0][0],exp_data.data[0][int(exp_data.data[0].size())-1],30);
    
    _event.Man_bins = SOL_YES;
    _event.Init_fast_generator();
    
    
    std::cout<<"The chisq calculator is initialized."
    
    return 0;
}

double chisq_real_one::Calc_chi2_nosys(vec params);
{
    if(int(params.size())!=Proba_engine.n_of_params)
    {
        std::cerr<<"ERROR!Number of parameter mismatched\n"
        exit(-1);
    }
    
    Proba_engine.Tan_Th12 = params[0];
    Proba_engine.Th13 = params[1];
    Proba_engine.Dm21 = params[4];
    Proba_engine.Delta = params[3];
    Proba_engine.Tau = params[6];
    
    _event.generate_events();
    
    double res=0.0;
    
    for(int i = 0;i<_event.n_bins;i++)
    {
        double expected_rate = _event.Events[i]+bkg_data.data[1];
        
        res = res + poiss_likelihood(expected_rate,exp_data.data[1]);
    }
    
    return 0;
}

double poiss_likelihood(double n_test,double n_true)
{
    double res = 2.0*(n_true*log(n_true/n_test)+n_test-n_true);
    
    return res;
}

