#ifndef SVOL_MODS_H
#define SVOL_MODS_H

#include <limits>
#include <Eigen/Dense>
#include <pf/rv_samp.h>
#include <pf/rv_eval.h>
#include <pf/resamplers.h>
#include <ssme/liu_west_filter.h>
#include <ssme/pswarm_filter.h>
#include "svol_leverage_mod.h"

#define DIMX 1
#define DIMY 1
#define DIMPARAM 4
#define DIMCOV 1

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
template<size_t NUMPARTS, typename float_t>
class svol_lw_1_par
        : public LWFilterWithCovsFutureSimulator<NUMPARTS, DIMX, DIMY, DIMCOV, DIMPARAM, float_t>
{
public:

    using ssv = Eigen::Matrix<float_t, DIMX, 1>;
    using osv = Eigen::Matrix<float_t, DIMY, 1>;
    using csv = Eigen::Matrix<float_t, DIMCOV,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;

private:

    // used for sampling states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

    // how many days to expiration (aka how many days into the future you are simulating)
    unsigned int m_dte;

public:
    // ctor
    svol_lw_1_par(float_t delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l,
                  float_t sig_u, float_t rho_l, float_t rho_u, unsigned dte);

    // pure virtual functions that we need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_new_param) override;
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

    // these probably need to be refactored so that the class just asks for the number of time steps as a cl arg
    std::vector<std::array<osv,nparts>> sim_future(const osv &yt){
        return this->sim_future_obs(m_dte, yt);
    }
    sim_future_obs
};


template<size_t NUMPARTS, typename float_t>
svol_lw_1_par<NUMPARTS,float_t>::svol_lw_1_par(float_t delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u, unsigned dte)
        : LWFilterWithCovsFutureSimulator<NUMPARTS,DIMX,DIMY,DIMCOV,DIMPARAM,float_t>(
        std::vector<std::string>{"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
        delta)           // PRIORS
        // REMINDER: does output appear to be extremely sensitive to these?
        , m_phi_sampler(phi_l, phi_u)
        , m_mu_sampler(mu_l, mu_u)
        , m_sigma_sampler(sig_l, sig_u)
        , m_rho_sampler(rho_l, rho_u)
        , m_dte(dte)
{
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_par<NUMPARTS,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_par<NUMPARTS,float_t>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    xt(0) = untrans_old_param(1) + untrans_old_param(0)*(xtm1(0) - untrans_old_param(1));
    xt(0) += cov_data(0)*untrans_old_param(3)*untrans_old_param(2)*std::exp(-.5*xtm1(0));
    return xt;
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_par<NUMPARTS,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_par<NUMPARTS,float_t>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_new_param) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    xt(0) = untrans_new_param(1) + untrans_new_param(0)*(xtm1(0) - untrans_new_param(1)) + zt(0)*untrans_new_param(3)*untrans_new_param(2)*std::exp(-.5*xtm1(0));
    xt(0) += m_stdNormSampler.sample() * untrans_new_param(2) * std::sqrt( 1.0 - untrans_new_param(3) * untrans_new_param(3));
    return xt;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_par<NUMPARTS,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_par<NUMPARTS,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_par<NUMPARTS,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}

template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_par<NUMPARTS,float_t>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    osv ytsamp;
    ytsamp(0) = m_stdNormSampler.sample() * std::exp(.5*xt(0));
    return ytsamp;
}


/**
 * @brief Liu-West filter (sampling parameters from a parameterized 
 * parameter distribution with a prng, as opposed to drawing
 * unformly from samples that have been saved to a file)
 * parameter order: phi, mu, sigma, rho
 */
template<size_t NUMPARTS, typename float_t>
class svol_lw_1_from_csv
        : public LWFilterWithCovsFutureSimulator<NUMPARTS, DIMX, DIMY, DIMCOV, DIMPARAM, float_t>
{
public:

    using ssv = Eigen::Matrix<float_t, DIMX, 1>;
    using osv = Eigen::Matrix<float_t, DIMY, 1>;
    using csv = Eigen::Matrix<float_t, DIMCOV,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;

private:

    // used for sampling states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, float_t> m_param_sampler;

    // days to expiration (aka how many days into the future you want to simulate)
    unsigned int m_dte;

public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_1_from_csv(float_t delta, const std::string &csv_of_samples, unsigned dte);

    // functions that we need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_new_param) override;
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

};


template<size_t NUMPARTS, typename float_t>
svol_lw_1_from_csv<NUMPARTS,float_t>::svol_lw_1_from_csv(float_t delta, const std::string &csv_of_samples, unsigned dte)
    : LWFilterWithCovsFutureSimulator<NUMPARTS,DIMX,DIMY,DIMCOV,DIMPARAM,float_t>(
    		std::vector<std::string> {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            	delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_of_samples)
    , m_dte(dte)
{
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_from_csv<NUMPARTS,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_from_csv<NUMPARTS,float_t>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) -> ssv
{
    // phi, mu, sigma, rho
    ssv ans;
    ans(0) = untrans_old_param(1) + untrans_old_param(0)*(xtm1(0) - untrans_old_param(1));
    ans(0) += cov_data(0)*untrans_old_param(3)*untrans_old_param(2)*std::exp(-.5*xtm1(0));
    return ans;
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_from_csv<NUMPARTS,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_from_csv<NUMPARTS,float_t>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_new_param) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_new_param(1) + untrans_new_param(0)*(xtm1(0) - untrans_new_param(1)) + zt(0)*untrans_new_param(3)*untrans_new_param(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_new_param(2) * std::sqrt( 1.0 - untrans_new_param(3) * untrans_new_param(3));
    return xt;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_from_csv<NUMPARTS,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_1_from_csv<NUMPARTS,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_from_csv<NUMPARTS,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    return m_param_sampler.samp();
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_1_from_csv<NUMPARTS,float_t>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    osv ytsamp;
    ytsamp(0) = m_stdNormSampler.sample() * std::exp(.5*xt(0));
    return ytsamp;
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, mu, sigma, rho
 */
template<size_t NUMPARTS, typename float_t>
class svol_lw_2_par : public LWFilter2WithCovsFutureSimulator<NUMPARTS, DIMX, DIMY, DIMCOV, DIMPARAM, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, DIMX, 1>;
    using osv = Eigen::Matrix<float_t, DIMY, 1>;
    using csv = Eigen::Matrix<float_t, DIMCOV,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

    // days to expiration (aka how many days into the future you want to simulate)
    unsigned int m_dte;

public:
    // ctor
    svol_lw_2_par() = delete;
    svol_lw_2_par(float_t delta, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u, unsigned m_dte);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    float_t logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt) override;
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

};


template<size_t NUMPARTS, typename float_t>
svol_lw_2_par<NUMPARTS,float_t>::svol_lw_2_par(float_t delta, float_t phi_l, float_t phi_u, float_t mu_l,
                                             float_t mu_u, float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u,
                                             unsigned dte)
        : LWFilter2WithCovsFutureSimulator<NUMPARTS,DIMX,DIMY,DIMCOV,DIMPARAM,float_t>(
        std::vector<std::string> {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
        delta)           // PRIORS    // REMINDER: does output appear to be extremely sensitive to these?
        // REMINDER: does output appear to be extremely sensitive to these?
        , m_phi_sampler(phi_l, phi_u)
        , m_mu_sampler(mu_l, mu_u)
        , m_sigma_sampler(sig_l, sig_u)
        , m_rho_sampler(rho_l, rho_u)
        , m_dte(dte)
{
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_par<NUMPARTS,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_par<NUMPARTS,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_par<NUMPARTS,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_par<NUMPARTS,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_par<NUMPARTS,float_t>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_par<NUMPARTS,float_t>::qSamp(const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_par<NUMPARTS,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_par<NUMPARTS,float_t>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}

template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_par<NUMPARTS,float_t>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    osv ytsamp;
    ytsamp(0) = m_stdNormSampler.sample() * std::exp(.5*xt(0));
    return ytsamp;
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t NUMPARTS, typename float_t>
class svol_lw_2_from_csv
        : public LWFilter2WithCovsFutureSimulator<NUMPARTS, DIMX, DIMY, DIMCOV, DIMPARAM, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, DIMX, 1>;
    using osv = Eigen::Matrix<float_t, DIMY, 1>;
    using csv = Eigen::Matrix<float_t, DIMCOV,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, float_t> m_param_sampler;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_2_from_csv(float_t delta, const std::string &csv_param_samples, unsigned dte);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    float_t logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt) override;
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;
};

template<size_t NUMPARTS, typename float_t>
svol_lw_2_from_csv<NUMPARTS,float_t>::svol_lw_2_from_csv(float_t delta, const std::string &csv_param_samples,
                                                       unsigned dte)
    : LWFilter2WithCovsFutureSimulator<NUMPARTS,DIMX,DIMY,DIMCOV,DIMPARAM,float_t>(
    				std::vector<std::string>{"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            			delta)           // PRIORS    // REMINDER: does output appear to be extremely sensitive to these?
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_param_samples)
    , m_dte(dte)
{
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_from_csv<NUMPARTS,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_from_csv<NUMPARTS,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    // phi, mu, sigma, rho
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return x1samp;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_from_csv<NUMPARTS,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, sd, true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_from_csv<NUMPARTS,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, std::exp(.5*xt(0)), true);
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_from_csv<NUMPARTS,float_t>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_from_csv<NUMPARTS,float_t>::qSamp(const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) -> ssv
{
    // phi, mu, sigma, rho
    ssv xt;
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    xt(0) = mean + m_stdNormSampler.sample() * untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return xt;
}


template<size_t NUMPARTS, typename float_t>
float_t svol_lw_2_from_csv<NUMPARTS,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt)
{
    // phi, mu, sigma, rho
    float_t mean = untrans_pt(1) + untrans_pt(0)*(xtm1(0) - untrans_pt(1)) + cov_data(0)*untrans_pt(3)*untrans_pt(2)*std::exp(-.5*xtm1(0));
    float_t sd = untrans_pt(2) * std::sqrt( 1.0 - untrans_pt(3) * untrans_pt(3));
    return rveval::evalUnivNorm<float_t>(xt(0), mean, sd, true);
}


template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_from_csv<NUMPARTS,float_t>::paramPriorSamp() -> psv
{
    return m_param_sampler.samp();
}

template<size_t NUMPARTS, typename float_t>
auto svol_lw_2_from_csv<NUMPARTS,float_t>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    osv ytsamp;
    ytsamp(0) = m_stdNormSampler.sample() * std::exp(.5*xt(0));
    return ytsamp;
}


/**
 * @brief particle swarm filter (many bootstrap filters)
 * this samples parameters from parameterized distribution
 */
template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm_1 : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMX,float_t>, float_t>,
        1,
        n_state_parts,
        n_param_parts,
        DIMY, DIMX, DIMCOV, DIMPARAM>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMX,float_t>, float_t>;
    using SwarmBase = SwarmWithCovs<ModType, 1, n_state_parts, n_param_parts, DIMY, DIMX, DIMCOV, DIMPARAM>;
    using ssv = Eigen::Matrix<float_t,DIMX,1>;
    using csv = Eigen::Matrix<float_t,DIMCOV,1>;
    using osv = Eigen::Matrix<float_t,DIMY,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;
    using state_cov_parm_func = typename SwarmBase::state_cov_parm_func;

private:

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_mu_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_rho_sampler;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

public:

    svol_swarm_1() = delete;

    // ctor
    svol_swarm_1(const std::vector<state_cov_parm_func>& fs, float_t phi_l, float_t phi_u, float_t mu_l, float_t mu_u,
                 float_t sig_l, float_t sig_u, float_t rho_l, float_t rho_u, unsigned int dte)
            : SwarmBase(fs)
            , m_phi_sampler(phi_l, phi_u)
            , m_mu_sampler(mu_l, mu_u)
            , m_sigma_sampler(sig_l, sig_u)
            , m_rho_sampler(rho_l, rho_u)
            , m_dte(dte)
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
                       untrans_params(3),
                       m_dte);
    }
};



/**
 * @brief particle swarm filter (many bootstrap filters)
 * this samples parameters from csv file
 */

template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm_2
        : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMX,float_t>,float_t>,
                1,
                n_state_parts,
                n_param_parts,
                DIMY, DIMX, DIMCOV, DIMPARAM>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMX,float_t>, float_t>;
    using SwarmBase = SwarmWithCovs<ModType, 1, n_state_parts, n_param_parts, DIMY, DIMX, DIMCOV, DIMPARAM>;
    using ssv = Eigen::Matrix<float_t,DIMX,1>;
    using csv = Eigen::Matrix<float_t,DIMCOV,1>;
    using osv = Eigen::Matrix<float_t,DIMY,1>;
    using psv = Eigen::Matrix<float_t, DIMPARAM,1>;
    using state_cov_parm_func = typename SwarmBase::state_cov_parm_func;

private:

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, float_t> m_param_sampler;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

public:

    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_swarm_2(const std::string &param_csv_filename, const std::vector<state_cov_parm_func>& fs, unsigned int dte)
            : SwarmBase(fs)
            , m_param_sampler(param_csv_filename)
            , m_dte(dte)
    {
    }


    // functions tha twe need to define
    psv samp_untrans_params() override {
        return m_param_sampler.samp();
    }

    ModType instantiate_mod(const psv& untrans_params) {
        // order: phi, mu, sigma, rho
        auto param = m_param_sampler.samp();
        return ModType(param(0),
                       param(1),
                       param(2),
                       param(3),
                       m_dte);
    }
};

#endif // SVOL_MODS_H
