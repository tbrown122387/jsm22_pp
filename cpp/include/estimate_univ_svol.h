#ifndef ESTIMATE_UNIV_SVOL_H
#define ESTIMATE_UNIV_SVOL_H

#include <string>

#include <ssme/ada_pmmh_mvn.h>
#include <ssme/parameters.h>
#include <pf/rv_eval.h>
#include <pf/resamplers.h>

#include "svol_leverage_mod.h"


#define DIMSTATE 1
#define DIMOBS 1
#define NUMPARAMS 4


using namespace pf::resamplers;



template<size_t numparts, typename float_t>
class svol_leverage_estimator : public ada_pmmh_mvn<NUMPARAMS,DIMOBS,numparts,float_t>
{
public:

    using psv = Eigen::Matrix<float_t,NUMPARAMS,1>;
    using psm = Eigen::Matrix<float_t,NUMPARAMS,NUMPARAMS>;
    using osv = Eigen::Matrix<float_t,DIMOBS,1>;
    
    svol_leverage_estimator(
            const psv &start_trans_theta,
            std::vector<std::string> tts,
            const unsigned int &num_mcmc_iters,
            const unsigned int &num_pfilters, 
            const std::string &data_file, 
            const std::string &samples_base_name, 
            const std::string &messages_base_name,
            const bool &mc,
            const unsigned int &t0, 
            const unsigned int &t1,
            const psm &C0,
            bool print_to_console,
            unsigned int print_every_k);
        
    float_t log_prior_eval(const param::pack<float_t,NUMPARAMS>& theta);

    float_t log_like_eval(const param::pack<float_t,NUMPARAMS>& theta, const std::vector<osv> &data);

};


template<size_t numparts, typename float_t>
svol_leverage_estimator<numparts,float_t>::svol_leverage_estimator(
                                                    const psv &start_trans_theta,
                                                    std::vector<std::string> tts,
                                                    const unsigned int &num_mcmc_iters,
                                                    const unsigned int &num_pfilters, 
                                                    const std::string &data_file, 
                                                    const std::string &samples_base_name, 
                                                    const std::string &messages_base_name,
                                                    const bool &mc,
                                                    const unsigned int &t0,
                                                    const unsigned int &t1,
                                                    const psm &C0,
                                                    bool print_to_console,
                                                    unsigned int print_every_k) 
    : ada_pmmh_mvn<NUMPARAMS,DIMOBS,numparts,float_t>(start_trans_theta, 
                                                      tts, 
                                                      num_mcmc_iters, 
                                                      num_pfilters,
                                                      data_file, 
                                                      samples_base_name, 
                                                      messages_base_name, 
                                                      mc, 
                                                      t0, 
                                                      t1, 
                                                      C0, 
                                                      print_to_console,
                                                      print_every_k)
{
}


template<size_t numparts, typename float_t>
float_t svol_leverage_estimator<numparts,float_t>::log_prior_eval(const param::pack<float_t,NUMPARAMS>& theta)
{
    // value to be returned
    float_t returnThis(0.0);
   
    // phi, mu, sigmaSq, rho
    // unpack parameters    
    float_t phi  = theta.get_untrans_params(0,0)(0);
    float_t mu   = theta.get_untrans_params(1,1)(0);
    float_t sigmaSq = theta.get_untrans_params(2,2)(0);
    float_t rho   = theta.get_untrans_params(3,3)(0);

    // phi ~ uniform(0, .99)
    returnThis += rveval::evalUniform<float_t>(phi, 0, .99, true);

    // mu ~ Normal(0,1)
    returnThis += rveval::evalUnivNorm<float_t>(mu, 0.0, 10.0, true);

    // ss ~ InverseGamma(.001, .001)
    returnThis += rveval::evalUnivInvGamma<float_t>(sigmaSq, .001, .001, true);

    // rho ~ Uniform(-1, 0)
    returnThis += rveval::evalUniform<float_t>(rho, -1.0, 0, true);


    return returnThis;    
}


template<size_t numparts, typename float_t>
float_t svol_leverage_estimator<numparts,float_t>::log_like_eval(const param::pack<float_t,NUMPARAMS>& theta, const std::vector<osv> &data)
{

    // jump out if there's a problem with the data
    if(data.empty())
        throw std::length_error("can't read in data\n");

    // the value to be returned
    float_t logLike(0.0);
    
    // instantiate model (needs phi, mu, sigma, rho)
    svol_leverage<numparts, mn_resampler<numparts,DIMSTATE,float_t>, float_t> 
        mod(theta.get_untrans_params(0,0)(0), 
            theta.get_untrans_params(1,1)(0),
            std::sqrt(theta.get_untrans_params(2,2)(0)),
            theta.get_untrans_params(3,3)(0),
            1); // set DTE to 1, but this doesn't matter at all because it isn't used
    
    // iterate through data and calculate sum of log p(y_t | y_{1:t-1}) s
    unsigned int row (0);
    while(row < data.size()){
    
    	if(row > 0){
	        mod.filter(data[row], data[row-1]);      	
    	}else{ 
            // row == 0
    		// covariates are just ignored, so we can plug in anything that's the right dimension
            mod.filter(data[row], data[row]);      	
    	}
        logLike += mod.getLogCondLike();
        row++;
    }

    // return the estimate of the log likelihood
    return logLike;    
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// do_ada_pmmh_svol_leverage ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


template<size_t numparts, typename float_t>
void do_ada_pmmh_svol_leverage(const std::string &datafile, 
                           const std::string &samples_base_name, 
                           const std::string &messages_base_name, 
                           unsigned int num_mcmc_iters, 
                           unsigned int num_pfilters,
                           bool multicore)
{

    // the chain's starting parameters
    using psv = Eigen::Matrix<float_t,NUMPARAMS,1>;
    using psm = Eigen::Matrix<float_t,NUMPARAMS,NUMPARAMS>;

    // phi, mu, sigmaSq, rho
    std::vector<std::string> tts {"logit", "null", "log", "twice_fisher"}; 
    psv start_trans_theta;
    start_trans_theta << rveval::logit<float_t>(.5), 0.0, std::log(2.0e-4), rveval::twiceFisher<float_t>(-.5);

    
    // the chain's initial covariance matrix 
    psm C0 = psm::Identity()*.15;
    unsigned int t0 = 150;  // start adapting the covariance at this iteration
    unsigned int t1 = 1000; // end adapting the covariance at this iteration
    svol_leverage_estimator<numparts,float_t> mcmcobj(
                                                    	start_trans_theta,
                                                    	tts,
                                                    	num_mcmc_iters, 
                                                        num_pfilters,
                                                    	datafile, 
                                                    	samples_base_name, 
                                                    	messages_base_name, 
                                                    	multicore, 
                                                    	t0,
                                                    	t1,
                                                    	C0,
                                                    	false, // print console
                                                        1); // print every
    mcmcobj.commence_sampling();

}



#endif //ESTIMATE_MSL1_H
