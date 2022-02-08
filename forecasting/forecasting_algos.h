#ifndef SVOL_MODS_H
#define SVOL_MODS_H

#include <limits>
#include <Eigen/Dense>

#include <pf/rv_samp.h> // for sampling random numbers
#include <pf/rv_eval.h> // for evaluating densities and pmfs
#include <pf/resamplers.h>

#include <ssme/liu_west_filter.h>
#include <ssme/pswarm_filter.h>

#include "../shared_cpp/svol_leverage_mod.h"


#define dimx 1
#define dimy 1
#define dimparam 4
#define dimcov 1

using namespace pf::filters;
using namespace pf;
using namespace pf::bases;
using namespace pf::resamplers;



/**
 * @brief Liu-West filter (sampling parameters from a parameterized 
 * parameter distribution with a prng, as opposed to drawing
 * unformly from samples that have been saved to a file)
 * parameter order: phi, mu, sigma, rho
 */
template<size_t nparts, typename float_t>
class svol_lw_1 : public LWFilterWithCovs<nparts, dimx, dimy, dimcov, dimparam, float_t>
{
public:

    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using csv = Eigen::Matrix<float_t, dimcov,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // used for sampling states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

public:
    // ctor
    svol_lw_1(const float_t& delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u);

    // functions that we need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_p1);    
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_p1); 
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    psv paramPriorSamp();

};


template<size_t nparts, typename float_t>
svol_lw_1<nparts,float_t>::svol_lw_1(const float_t& delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u)
    : LWFilterWithCovs<nparts,dimx,dimy,dimcov,dimparam,float_t>(
    		std::vector<std::string> tts {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            	delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_phi_sampler(phi_l, phi_u) 
    , m_mu_sampler(mu_l, mu_u)
    , m_sigma_sampler(sig_l, sig_u)      
    , m_rho_sampler(rho_l, rho_u)
{
}


template<size_t nparts, typename float_t>
float_t svol_lw_1<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_1<nparts,float_t>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    return untrans_pt(1) + untrans(0)*(xtm1 - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1); 
}


template<size_t nparts, typename float_t>
auto svol_lw_1<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t nparts, typename float_t>
auto svol_lw_1<nparts,float_t>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t nparts, typename float_t>
float_t svol_lw_1<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_1<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
auto svol_lw_1<nparts,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}


/**
 * @brief Liu-West filter (sampling parameters from a parameterized 
 * parameter distribution with a prng, as opposed to drawing
 * unformly from samples that have been saved to a file)
 * parameter order: phi, mu, sigma, rho
 */
template<size_t nparts, typename float_t>
class svol_lw_1_from_csv : public LWFilterWithCovs<nparts, dimx, dimy, dimcov, dimparam, float_t>
{
public:

    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using csv = Eigen::Matrix<float_t, dimcov,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // used for sampling states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from mcmc samples
    csv_param_sampler<dimparam, float_t> m_param_sampler;


public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_1_from_csv(const float_t& delta, const std::string &csv_of_samples);

    // functions that we need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_p1);    
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_p1); 
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    psv paramPriorSamp();

};


template<size_t nparts, typename float_t>
svol_lw_1_from_csv<nparts,float_t>::svol_lw_1(const float_t& delta, const std::string &csv_of_samples)
    : LWFilterWithCovs<nparts,dimx,dimy,dimcov,dimparam,float_t>(
    		std::vector<std::string> tts {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            	delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_of_samples)  
{
}


template<size_t nparts, typename float_t>
float_t svol_lw_1_from_csv<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_1_from_csv<nparts,float_t>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    return untrans_pt(1) + untrans(0)*(xtm1 - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1); 
}


template<size_t nparts, typename float_t>
auto svol_lw_1_from_csv<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t nparts, typename float_t>
auto svol_lw_1_from_csv<nparts,float_t>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t nparts, typename float_t>
float_t svol_lw_1_from_csv<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_1_from_csv<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
auto svol_lw_1_from_csv<nparts,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    return m_param_sampler.sample();  
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t nparts, typename float_t>
class svol_lw_2 : public LWFilter2WithCovs<nparts, dimx, dimy, dimcov, dimparam, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using csv = Eigen::Matrix<float_t, dimcov,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

public:
    // ctor
    svol_lw_2(const float_t &delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    float_t logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt);
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt);
    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt);    
    psv paramPriorSamp();

};


template<size_t nparts, typename float_t>
svol_lw_2<nparts,float_t>::svol_lw2(const float_t& delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u)
    : LWFilter2WithCovs<nparts,dimx,dimy,dimparam,float_t>(
    				std::vector<std::string> tts {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            			delta)           // PRIORS    // REMINDER: does output appear to be extremely sensitive to these?
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_phi_sampler(phi_l, phi_u) 
    , m_mu_sampler(mu_l, mu_u)
    , m_sigma_sampler(sig_l, sig_u)      
    , m_rho_sampler(rho_l, rho_u)
{
}


template<size_t nparts, typename float_t>
float_t svol_lw_2<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t nparts, typename float_t>
float_t svol_lw_2<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_2<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_2<nparts,float_t>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2<nparts,float_t>::qSamp(const ssv &xtm1, const osv &yt, const psv& untrans_pt) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t nparts, typename float_t>
float_t svol_lw_2<nparts,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2<nparts,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t nparts, typename float_t>
class svol_lw_2_from_csv : public LWFilter2WithCovs<nparts, dimx, dimy, dimcov, dimparam, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using csv = Eigen::Matrix<float_t, dimcov,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from mcmc samples
    csv_param_sampler<dimparam, float_t> m_param_sampler;

public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_2_from_csv(const float_t &delta, const std::string &csv_param_samples);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    float_t logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt);
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt);
    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt);    
    psv paramPriorSamp();

};

template<size_t nparts, typename float_t>
svol_lw_2_from_csv<nparts,float_t>::svol_lw2(const float_t& delta, const std::string &csv_param_samples)
    : LWFilter2WithCovs<nparts,dimx,dimy,dimparam,float_t>(
    				std::vector<std::string> tts {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            			delta)           // PRIORS    // REMINDER: does output appear to be extremely sensitive to these?
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_param_samples) 
{
}


template<size_t nparts, typename float_t>
float_t svol_lw_2_from_csv<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2_from_csv<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t nparts, typename float_t>
float_t svol_lw_2_from_csv<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_2_from_csv<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
float_t svol_lw_2_from_csv<nparts,float_t>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2_from_csv<nparts,float_t>::qSamp(const ssv &xtm1, const osv &yt, const psv& untrans_pt) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t nparts, typename float_t>
float_t svol_lw_2_from_csv<nparts,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans(0)*(xtm1(0) - untrans_pt(1)) + untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t nparts, typename float_t>
auto svol_lw_2_from_csv<nparts,float_t>::paramPriorSamp() -> psv
{
    return m_param_sampler.sample();
}



/**
 * @brief particle swarm filter (many bootstrap filters)
 * this samples parameters from parameterized distribution
 */

template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm_1 : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,dimx,float_t>, float_t>, 
                                1, 
                                n_state_parts, 
                                n_param_parts, 
                                dimy, dimx, dimcov, dimparam>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,dimx,float_t>, float_t>;
    using SwarmBase = Swarm<ModType, 1, n_state_parts, n_param_parts, dimy, dimx, dimparam>;
    using ssv = Eigen::Matrix<float_t,dimx,1>;
    using csv = Eigen::Matrix<float_t,dimcov,1>;
    using osv = Eigen::Matrix<float_t,dimy,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;
    using state_parm_func = typename SwarmBase::state_parm_func;

private:

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

public:

    svol_swarm_1() = delete;
    
    // default ctor
    svol_swarm_1(float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u) 
        : SwarmBase() 
    , m_phi_sampler(phi_l, phi_u) 
    , m_mu_sampler(mu_l, mu_u)
    , m_sigma_sampler(sig_l, sig_u)      
    , m_rho_sampler(rho_l, rho_u)
    {}

    // ctor
    svol_swarm_1(const std::vector<state_parm_func>& fs, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u)
        : SwarmBase(fs)
    , m_phi_sampler(phi_l, phi_u) 
    , m_mu_sampler(mu_l, mu_u)
    , m_sigma_sampler(sig_l, sig_u)      
    , m_rho_sampler(rho_l, rho_u)
    {
    }


    // functions tha twe need to define
    psv samp_untrans_params() override {
        psv untrans_param; //order: phi,mu,sigma,rho
        untrans_param(0) = m_phi_sampler.sample();
        untrans_param(1) = m_mu_sampler.sample();
        untrans_param(2) = m_sigma_sampler.sample();
        untrans_param(3) = m_rho_sampler.sample();
        return untrans_param; 
    }

    ModType instantiate_mod(const psv& untrans_params) {
        // order: phi, beta, sigma 
        return ModType(untrans_params(0), 
                       untrans_params(1),
                       untrans_params(2),
                       untrans_params(3));
    }
};


/**
 * @brief particle swarm filter (many bootstrap filters)
 * this samples parameters from csv file
 */

template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm_2 : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,dimx,float_t>, float_t>, 
                                1, 
                                n_state_parts, 
                                n_param_parts, 
                                dimy, dimx, dimcov, dimparam>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,dimx,float_t>, float_t>;
    using SwarmBase = Swarm<ModType, 1, n_state_parts, n_param_parts, dimy, dimx, dimparam>;
    using ssv = Eigen::Matrix<float_t,dimx,1>;
    using csv = Eigen::Matrix<float_t,dimcov,1>;
    using osv = Eigen::Matrix<float_t,dimy,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;
    using state_parm_func = typename SwarmBase::state_parm_func;

private:

    // for sampling from mcmc samples
    csv_param_sampler<dimparam, float_t> m_param_sampler;

public:

    // default ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_swarm_2(const std::string &param_csv_filename) 
    : SwarmBase() 
    , m_param_sampler(param_csv_filename)
    {    
    }

    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_swarm_2(const std::vector<state_parm_func>& fs)
        : SwarmBase(fs)
    	, m_param_sampler(param_csv_filename)
    {
    }


    // functions tha twe need to define
    psv samp_untrans_params() override {
    	return m_param_sampler.sample();
    }

    ModType instantiate_mod(const psv& untrans_params) {
        // order: phi, mu, sigma, rho
        auto param = m_param_sampler.sample();
        return ModType(param(0), 
                       param(1),
                       param(2),
                       param(3));
    }
};


/**
 * @brief particle swarm filter (many auxiliary particle filters)
 * TODO
 * parameter order: phi, beta, sigma
 */

#endif // SVOL_MODS_H
