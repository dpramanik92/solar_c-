//
//
//  chisqmin.hpp
//
//  Created by Dipyaman Pramanik on 28/03/21
//
//

#ifndef chisqmin_hpp
#define chisqmin_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include "read_files.hpp"
#include "probability.hpp"
#include "numerical.hpp"
#include "event.hpp"
#include "minimizer.hpp"

namespace chisq
{
    class chisq_real
    {
        private:
          
        public:
            int sys_stat,res_stat;
            double sampl_min,sampl_max,numsamp;

            dec_prob Proba_engine;
            file_reader _exp_data,_bkg_data;
            Event_generator _event;
            int fin_flav,channel,particle;
            int Init_chisq_calculator(Event_generator,file_reader,file_reader,int,
                                      int);
            int Set_sampling_points(double,double,int);
    };

    class No_sys: public chisq_real
    {
        public:
        int Init(Event_generator,file_reader,file_reader,int,int);
        double Calc_chi2_nosys(vec);

    };

    class SK_IV: public chisq_real
    {
        private:
            double norm_sig,norm_bkg,calib_sig,calib_bkg;
            int binned_sys;

        public:
            double Calc_chi2_nosys(vec);
            

            vec sigma;
            int Set_binned_systematics(int);
            int Set_sys(double,double,double,double);
            int Init(Event_generator,file_reader,file_reader,int,int);
            int expected_rates(vec);
            double Function(vec);
            int statistics(vec);
            double pull_term(vec,vec);
            double pull_term(double,double);

    };



    template<class detector>
    class sys_minimizer
    {
        public:
            double Minimize(detector,vec,vec);

            vec position;

    };


    
    double gauss_likelihood(double,double);

    double poiss_likelihood(double,double);


}

#endif /* chisqmin_hpp */

