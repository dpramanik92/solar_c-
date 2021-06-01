//
//  cchisqmin.cpp
//
//
//
//
//
//  Created by Dipyaman Pramanik on 28/01/21
//

#include "chisqmin2.hpp"


int chisq::chisq_real::Init_chisq_calculator(Event_generator event,file_reader exp_data,file_reader bkg_data)
{
    _event = event;
    
   _bkg_data = bkg_data;
   _exp_data = exp_data;
   
    sys_stat = _event.res_stat;

//    _event.Init_evgen(fin_flav,particle,channel);
    
    for(int i=0;i<int(exp_data.data[0].size());i++)
    {
        _event.manual_bins.push_back(exp_data.data[0][i]);
        _event.manual_f.push_back(exp_data.data[1][i]);  
    }
    
    
    
//    _event.Set_fast_event_generator(SOL_YES,exp_data.data[0][0],exp_data.data[0][int(exp_data.data[0].size())-1],30);
    
    _event.Man_bins = SOL_YES;
//    _event.Init_fast_generator();
//
    _event.Init_evgen(fin_flav,particle,channel);
//    
    
    
    std::cout<<"The chisq calculator is initialized.\n";
    
    return 0;
}


double chisq::No_sys::Calc_chi2_nosys(vec params)
{
    /*if(int(params.size())!=Proba_engine.n_of_params)
    {
        std::cerr<<"ERROR!Number of parameter mismatched\n";
        exit(-1);
    }
    */
    _event.Proba_engine.Tan_Th12 = params[0];
    _event.Proba_engine.Th13 = params[1];
    _event.Proba_engine.Dm21 = params[4];
    _event.Proba_engine.Delta = params[3];
    _event.Proba_engine.Tau1 = params[6];
    _event.Proba_engine.Tau2 = params[6];

    


    _event.Proba_engine.L = 1.0;


    
    _event.generate_events();
    
    double res=0.0;
    
        

    for(int i = 0;i<_event.n_bins-2;i++)
    {
        double expected_rate = _event.Events[i]+_bkg_data.data[2][i];

//        std::cout<<_event.Events[i]<<"\t"<<_exp_data.data[2][i]<<std::endl;
        
        res = res + gauss_likelihood(expected_rate,_exp_data.data[2][i]);
    }
    
    return res;
}



double chisq::SK_IV::Calc_chi2_nosys(vec params)
{
    /*if(int(params.size())!=Proba_engine.n_of_params)
    {
        std::cerr<<"ERROR!Number of parameter mismatched\n";
        exit(-1);
    }
    */
    _event.Proba_engine.Tan_Th12 = params[0];
    _event.Proba_engine.Th13 = params[1];
    _event.Proba_engine.Dm21 = params[4];
    _event.Proba_engine.Delta = params[3];
    _event.Proba_engine.Tau1 = params[6];
    _event.Proba_engine.Tau2 = params[6];

    


    _event.Proba_engine.L = 1.0;


    
    _event.generate_events();
    
    double res=0.0;
    
        

    for(int i = 0;i<_event.n_bins-2;i++)
    {
        double expected_rate = _event.Events[i]+_bkg_data.data[2][i];

//        std::cout<<_event.Events[i]<<"\t"<<_exp_data.data[2][i]<<std::endl;
        
        res = res + gauss_likelihood(expected_rate,_exp_data.data[2][i]);
    }
    
    return res;
}




double chisq::SK_IV::Function(vec sys_params)
{
    /*
    double Reactor[4] = {25.2,0.0,0.0,0.0};
    double Li[4] = {24.2,10.97,5.8,0.0};
    double NCQE[4] = {6.5,6.3,4.92,2.52};
    double nonNCQE[4] = {0.8,0.66,0.82,0.78};
    double Accidental[4] ={43.2,19.06,9.2,3.2};
*/

    double alpha = sys_params[0];
    double beta = sys_params[1];

    double sigma_alp = sigma[0];
    double sigma_bet = sigma[1];

    vec omega,omega_p;

    for(int i=0;i<4;i++)
    {
        omega.push_back(sys_params[i+2]);
        omega_p.push_back(sys_params[6+i]);
    }

    
    double res = 0.0;

    for(int i=0;i<4;i++)
    {
        double mu1 = 0.0;

       
        
            mu1 = mu1 + (1+alpha)*(1+omega[i])*_bkg_data.data[2][i];
       

        double mu2 = (1+beta)*(1+omega_p[i])*_event.Events[i];

        double mu = mu1 + mu2;

        res = res + poiss_likelihood(mu,_exp_data.data[2][i]) + pull_term(alpha,sigma_alp)+pull_term(beta,sigma_bet);
    }


    return res;
}


double chisq::SK_IV::pull_term(double alpha,double sigma)
{
    double res=0.0;

    res = (alpha/sigma)*(alpha/sigma);

    return res;
}

double chisq::SK_IV::pull_term(vec alpha,vec sigma)
{
    double res=0.0;

    for(int i=0;i<alpha.size();i++)
    {
        res = res + (alpha[i]/sigma[i])*(alpha[i]/sigma[i]);

    }


    return res;
}







double chisq::gauss_likelihood(double n_test,double n_true)
{
    double res=0;
//    std::cout<<n_test<<"\t"<<n_true<<"\n";
    if(n_test!=0)
    {
        res = (n_true-n_test)*(n_true-n_test)/n_test;
    }
    return res;
}


double chisq::poiss_likelihood(double n_test,double n_true)
{

    double res = 2*(n_test-n_true+n_true*log(n_true/n_test)); 

    return res;
}



template<class detector>
double chisq::sys_minimizer<detector>::Minimize(detector _detector,vec start_values,vec params)
{
    double res=0.0;

    if(_detector.sys_stat==0)
    {
        return _detector.Calc_chi2_nosys(params);

    }
    else
    {
    
        multi_minimizer<detector> _minimizer;

        _minimizer.init_minimizer(int(start_values.size()),start_values,_detector);

        _minimizer.powell();

        res = _minimizer.curr_val;
    }


    return res;
}


