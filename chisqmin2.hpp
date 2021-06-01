//
//
//  chisqmin.hpp
//
//  Created by Dipyaman Pramanik on 28/03/21
//
//

#ifndef chisqmin2_hpp
#define chisqmin2_hpp

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
            int sys_stat;
            dec_prob Proba_engine;
            file_reader _exp_data,_bkg_data;
            Event_generator _event;
            int fin_flav,channel,particle;
            int Init_chisq_calculator(Event_generator,file_reader,file_reader);
    };

    class No_sys: public chisq_real
    {
        public:
        double Calc_chi2_nosys(vec);

    };

    class SK_IV: public chisq_real
    {

        public:
            double Calc_chi2_nosys(vec);
            vec sigma;
            double Function(vec);
            double pull_term(vec,vec);
            double pull_term(double,double);

    };



    template<class detector>
    class sys_minimizer
    {
        public:

            double Minimize(detector,vec,vec);

    };


    
    double gauss_likelihood(double,double);

    double poiss_likelihood(double,double);


}

#endif /* chisqmin_hpp */

