//
//  cchisqmin.cpp
//
//  Created by Dipyaman Pramanik on 28/01/21
//

#include "chisqmin.hpp"


int chisq::chisq_real::Init_chisq_calculator(Event_generator event,file_reader exp_data,file_reader bkg_data,int _sys_stat,int _res_stat)
{
    _event = event;
    
   _bkg_data = bkg_data;
   _exp_data = exp_data;
   
    sys_stat = _sys_stat;
    res_stat = _res_stat;


    
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
    if(res_stat==SOL_YES)
    {
        _event.Set_fast_event_generator(SOL_YES,sampl_min,sampl_max,numsamp);
        _event.Init_fast_generator();
    }
//    
    
    
    std::cout<<"The chisq calculator is initialized.\n";
    
    return 0;
}

int chisq::chisq_real::Set_sampling_points(double _samp_min,double _samp_max,int num)
{
    sampl_min = _samp_min;
    sampl_max = _samp_max;
    numsamp = num;
    
    return 0;
}

int chisq::No_sys::Init(Event_generator event,file_reader exp_data,file_reader bkg_data,int _sys_stat,int _res_stat )
{
    sys_stat = _sys_stat;
    res_stat = _res_stat;
    
    

    Init_chisq_calculator(event,exp_data,bkg_data,sys_stat,res_stat);


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
    
        

    for(int i = 0;i<_event.n_bins;i++)
    {
        double expected_rate = _event.Events[i]+_bkg_data.data[2][i];

//        std::cout<<_event.Events[i]<<"\t"<<_exp_data.data[2][i]<<std::endl;
        
        res = res + poiss_likelihood(expected_rate,_exp_data.data[2][i]);
    }
    
    return res;
}


int chisq::SK_IV::expected_rates(vec params)
{
    _event.Proba_engine.Tan_Th12 = params[0];
    _event.Proba_engine.Th13 = params[1];
    _event.Proba_engine.Dm21 = params[4];
    _event.Proba_engine.Delta = params[3];
    _event.Proba_engine.Tau1 = params[6];
    _event.Proba_engine.Tau2 = params[6];

    _event.Proba_engine.L = 1.0;
    
    _event.generate_events();

    return 0;
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
    
        

    for(int i = 0;i<_event.n_bins;i++)
    {
        double expected_rate = _event.Events[i]+_bkg_data.data[2][i];
        
        
        res = res + poiss_likelihood(expected_rate,_exp_data.data[2][i]);
    }
    

    return res;
}



int chisq::SK_IV::statistics(vec sys_params)
{

     double res = 0.0;

    if(binned_sys!=0)
    {
    
        vec alpha,beta,sigma_alp,sigma_bet;
        vec omega,omega_p,sigma_om,sigma_bom;

        for(int i=0;i<4;i++)
        {
            alpha.push_back(sys_params[i]);
            beta.push_back(sys_params[i+4]);
            omega.push_back(sys_params[i+8]);
            omega_p.push_back(sys_params[i+12]);

            sigma_alp.push_back(norm_sig);
            sigma_bet.push_back(norm_bkg);
            sigma_om.push_back(calib_sig);
            sigma_bom.push_back(calib_bkg);

        

        }

    

        for(int i=0;i<4;i++)
        {
            double mu1 = 0.0;

       
        
            mu1 = mu1 + (1+alpha[i])*(1+omega[i])*_bkg_data.data[2][i];
       

            double mu2 = 0.0;

            mu2 = mu2 +  (1+beta[i])*(1+omega_p[i])*_event.Events[i];

            double mu = mu1 + mu2;


            res = res + poiss_likelihood(mu,_exp_data.data[2][i]);
        }

        double temp = res;

        double pulls = pull_term(beta,sigma_bet)+pull_term(omega,sigma_om)+pull_term(omega_p,sigma_bom);


        res = res+pulls;
    //    std::cout<<temp<<"\t"<<pulls<<"\t"<<res<<"\n";
        
    }
    else
    {
        double alpha,beta,omega,omega_p;
        double sigma_alp,sigma_bet,sigma_om,sigma_bom;

        alpha = sys_params[0];
        beta = sys_params[1];
        omega = sys_params[2];
        omega_p = sys_params[3];

        sigma_alp = norm_sig;
        sigma_bet = norm_bkg;
        sigma_om = calib_sig;
        sigma_bom = calib_bkg;


        for(int i=0;i<4;i++)
        {
            double mu1 = 0.0;

            mu1 = mu1 + (1+alpha)*(1+omega)*_bkg_data.data[2][i];

            double mu2 = 0.0;

            mu2 = mu2 + (1+beta)*(1+omega_p)*_event.Events[i];

            double mu = mu1 + mu2;

      //      std::cout<<mu1<<"\t"<<mu2<<"\t"<<mu<<"\t"<<poiss_likelihood(mu,_exp_data.data[2][i])<<"\n";
            
            res = res + poiss_likelihood(mu,_exp_data.data[2][i]);
        }

        double temp = res;

        double pulls = pull_term(alpha,sigma_alp) + pull_term(beta,sigma_bet) + pull_term(omega,sigma_om) + pull_term(omega_p,sigma_bom);




        res = temp + pulls;

  //      std::cout<<temp<<"\t"<<pulls<<"\t"<<res<<"\n";



    }

    return 0;
    

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


    double res = 0.0;

    if(binned_sys!=0)
    {
    
        vec alpha,beta,sigma_alp,sigma_bet;
        vec omega,omega_p,sigma_om,sigma_bom;

        for(int i=0;i<4;i++)
        {
            alpha.push_back(sys_params[i]);
            beta.push_back(sys_params[i+4]);
            omega.push_back(sys_params[i+8]);
            omega_p.push_back(sys_params[i+12]);

            sigma_alp.push_back(norm_sig);
            sigma_bet.push_back(norm_bkg);
            sigma_om.push_back(calib_sig);
            sigma_bom.push_back(calib_bkg);

        

        }

    

        for(int i=0;i<4;i++)
        {
            double mu1 = 0.0;

       
        
            mu1 = mu1 + (1+alpha[i])*(1+omega[i])*_bkg_data.data[2][i];
       

            double mu2 = 0.0;

            mu2 = mu2 +  (1+beta[i])*(1+omega_p[i])*_event.Events[i];

            double mu = mu1 + mu2;


            res = res + poiss_likelihood(mu,_exp_data.data[2][i]);
        }

        res = res+pull_term(alpha,sigma_alp)+pull_term(beta,sigma_bet)+pull_term(omega,sigma_om)+pull_term(omega_p,sigma_bom);
        
    }
    else
    {
        double alpha,beta,omega,omega_p;
        double sigma_alp,sigma_bet,sigma_om,sigma_bom;

        alpha = sys_params[0];
        beta = sys_params[1];
        omega = sys_params[2];
        omega_p = sys_params[3];

        sigma_alp = norm_sig;
        sigma_bet = norm_bkg;
        sigma_om = calib_sig;
        sigma_bom = calib_bkg;


        for(int i=0;i<4;i++)
        {
            double mu1 = 0.0;

            mu1 = mu1 + (1+alpha)*(1+omega)*_bkg_data.data[2][i];

            double mu2 = 0.0;

            mu2 = mu2 + (1+beta)*(1+omega_p)*_event.Events[i];

            double mu = mu1 + mu2;


            res = res + poiss_likelihood(mu,_exp_data.data[2][i]);
        }



        res = res + pull_term(alpha,sigma_alp) + pull_term(beta,sigma_bet) + pull_term(omega,sigma_om) + pull_term(omega_p,sigma_bom);



    }
   
   
    
//    res = res + pull_term(alpha,sigma_a;

    return res;
}

int chisq::SK_IV::Init(Event_generator event,file_reader exp_data,file_reader bkg_data,int _sys_stat,int _res_stat )
{
    sys_stat = _sys_stat;
    res_stat = _res_stat;
    binned_sys = 0;

    Init_chisq_calculator(event,exp_data,bkg_data,sys_stat,res_stat);
/*
    _event.Proba_engine.Tan_Th12 = params[0];
    _event.Proba_engine.Th13 = params[1];
    _event.Proba_engine.Dm21 = params[4];
    _event.Proba_engine.Delta = params[3];
    _event.Proba_engine.Tau1 = params[6];
    _event.Proba_engine.Tau2 = params[6];

    


    _event.Proba_engine.L = 1.0;

*/


    _event.generate_events();


    return 0;
}


int chisq::SK_IV::Set_sys(double _norm_sig,double _norm_bkg,double _calib_sig,double _calib_bkg)
{
    norm_sig = _norm_sig;
    norm_bkg = _norm_bkg;

    calib_sig = _calib_sig;
    calib_bkg = _calib_bkg;





    return 0;
}


int chisq::SK_IV::Set_binned_systematics(int _binned_sys)
{
    binned_sys = _binned_sys;
    return 0;
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
    
   // std::cout<<res<<"\n";

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
    double res;
    if(n_true>0)
    {
        res = 2*(n_test-n_true+n_true*log(n_true/n_test)); 
    }   
    else
    {
        res = 2*n_test;
    }
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

        int dimension = int(start_values.size());

        _detector.expected_rates(params);
    
        multi_minimizer<detector> _minimizer;

        

        _minimizer.init_minimizer(dimension,start_values,_detector);

        _minimizer.verbosity_m = 0;

        _minimizer.set_tolerance(1e-4);

        _minimizer.powell();

        for(int i=0;i<dimension;i++)
        {
            position.push_back(_minimizer.pos[i]);
        }

        res = _minimizer.curr_val;

        _minimizer.free_minimizer();
    }


    return res;
}



template class chisq::sys_minimizer<chisq::SK_IV>;
