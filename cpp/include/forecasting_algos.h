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

// define directives are included from svol_leverage_mod.h!

using namespace pf::filters;
using namespace pf;
using namespace pf::bases;
using namespace pf::resamplers;



/**
* @brief Liu-West filter (sa
 * mpling parameters from a parameterized
* parameter distribution with a prng, as opposed to drawing
* unformly from samples that have been saved to a file)
* parameter order: phi, mu, sigma, rho
*/
template<size_t NUMPARTS>
class svol_lw_1_par
        : public LWFilterWithCovsFutureSimulator<NUMPARTS, DIMSTATE, DIMOBS, DIMCOV, DIMPARAM, FLOATTYPE>
{
private:

    // used for sampling states
    rvsamp::UnivNormSampler<FLOATTYPE> m_stdNormSampler;

    // for sampling from the parameter prior
    rvsamp::UniformSampler<FLOATTYPE> m_phi_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_mu_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_sigma_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_rho_sampler;

    // how many days to expiration (aka how many days into the future you are simulating)
    unsigned int m_dte;

public:
    // ctor
    svol_lw_1_par(FLOATTYPE delta, FLOATTYPE phi_l, FLOATTYPE phi_u, FLOATTYPE mu_l, FLOATTYPE mu_u, FLOATTYPE sig_l,
                  FLOATTYPE sig_u, FLOATTYPE rho_l, FLOATTYPE rho_u, unsigned dte);

    // pure virtual functions that we need to define
    FLOATTYPE logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_new_param) override;
    FLOATTYPE logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

    // these probably need to be refactored so that the class just asks for the number of time steps as a cl arg
    std::vector<std::array<osv,NUMPARTS>> sim_future(const osv &yt){
        return this->sim_future_obs(m_dte, yt);
    }

};


template<size_t NUMPARTS>
svol_lw_1_par<NUMPARTS>::svol_lw_1_par(FLOATTYPE delta, FLOATTYPE phi_l, FLOATTYPE phi_u, FLOATTYPE mu_l, FLOATTYPE mu_u, FLOATTYPE sig_l, FLOATTYPE sig_u, FLOATTYPE rho_l, FLOATTYPE rho_u, unsigned dte)
        : LWFilterWithCovsFutureSimulator<NUMPARTS,DIMSTATE,DIMOBS,DIMCOV,DIMPARAM,FLOATTYPE>(
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


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_par<NUMPARTS>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
auto svol_lw_1_par<NUMPARTS>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) -> ssv
{
    return sl::propMu(xtm1, cov_data, untrans_old_param);
}


template<size_t NUMPARTS>
auto svol_lw_1_par<NUMPARTS>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    return sl::muSamp(untrans_p1, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
auto svol_lw_1_par<NUMPARTS>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_new_param) -> ssv
{
    return sl::fSamp(xtm1, zt, untrans_new_param, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_par<NUMPARTS>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    // phi, mu, sigma, rho
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_par<NUMPARTS>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return sl::logGEv(yt, xt);
}


template<size_t NUMPARTS>
auto svol_lw_1_par<NUMPARTS>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}


template<size_t NUMPARTS>
auto svol_lw_1_par<NUMPARTS>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    return sl::gSamp(xt, m_stdNormSampler.sample());
}


/**
 * @brief Liu-West filter (sampling parameters from a parameterized 
 * parameter distribution with a prng, as opposed to drawing
 * unformly from samples that have been saved to a file)
 * parameter order: phi, mu, sigma, rho
 */
template<size_t NUMPARTS>
class svol_lw_1_from_csv
        : public LWFilterWithCovsFutureSimulator<NUMPARTS, DIMSTATE, DIMOBS, DIMCOV, DIMPARAM, FLOATTYPE>
{
private:

    // used for sampling states
    rvsamp::UnivNormSampler<FLOATTYPE> m_stdNormSampler;  

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, FLOATTYPE> m_param_sampler;

    // days to expiration (aka how many days into the future you want to simulate)
    unsigned int m_dte;

public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_1_from_csv(FLOATTYPE delta, const std::string &csv_of_samples, unsigned dte);

    // functions that we need to define
    FLOATTYPE logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv propMu  (const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    ssv fSamp (const ssv &xtm1, const csv &zt, const psv& untrans_new_param) override;
    FLOATTYPE logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

    // these probably need to be refactored so that the class just asks for the number of time steps as a cl arg
    std::vector<std::array<osv,NUMPARTS>> sim_future(const osv &yt){
        return this->sim_future_obs(m_dte, yt);
    }

};


template<size_t NUMPARTS>
svol_lw_1_from_csv<NUMPARTS>::svol_lw_1_from_csv(FLOATTYPE delta, const std::string &csv_of_samples, unsigned dte)
    : LWFilterWithCovsFutureSimulator<NUMPARTS,DIMSTATE,DIMOBS,DIMCOV,DIMPARAM,FLOATTYPE>(
    		std::vector<std::string> {"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            	delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_of_samples)
    , m_dte(dte)
{
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_from_csv<NUMPARTS>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
auto svol_lw_1_from_csv<NUMPARTS>::propMu(const ssv &xtm1, const csv &cov_data, const psv& untrans_old_param) -> ssv
{
    return sl::propMu(xtm1, cov_data, untrans_old_param);
}


template<size_t NUMPARTS>
auto svol_lw_1_from_csv<NUMPARTS>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    return sl::muSamp(untrans_p1, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
auto svol_lw_1_from_csv<NUMPARTS>::fSamp(const ssv &xtm1, const csv &zt, const psv& untrans_new_param) -> ssv
{
    return sl::fSamp(xtm1, zt, untrans_new_param, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_from_csv<NUMPARTS>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_1_from_csv<NUMPARTS>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return sl::logGEv(yt, xt);
}


template<size_t NUMPARTS>
auto svol_lw_1_from_csv<NUMPARTS>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    return m_param_sampler.samp();
}


template<size_t NUMPARTS>
auto svol_lw_1_from_csv<NUMPARTS>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    return sl::gSamp(xt, m_stdNormSampler.sample());
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, mu, sigma, rho
 */
template<size_t NUMPARTS>
class svol_lw_2_par : public LWFilter2WithCovsFutureSimulator<NUMPARTS, DIMSTATE, DIMOBS, DIMCOV, DIMPARAM, FLOATTYPE>
{
private:

    // use this for samplign states
    rvsamp::UnivNormSampler<FLOATTYPE> m_stdNormSampler;

    // for sampling from the parameter prior
    rvsamp::UniformSampler<FLOATTYPE> m_phi_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_mu_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_sigma_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_rho_sampler;

    // days to expiration (aka how many days into the future you want to simulate)
    unsigned int m_dte;

public:
    // ctor
    svol_lw_2_par() = delete;
    svol_lw_2_par(FLOATTYPE delta, FLOATTYPE phi_l, FLOATTYPE phi_u, FLOATTYPE mu_l, FLOATTYPE mu_u, FLOATTYPE sig_l, FLOATTYPE sig_u, FLOATTYPE rho_l, FLOATTYPE rho_u, unsigned m_dte);

    // functions tha twe need to define
    FLOATTYPE logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    FLOATTYPE logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt) override;
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    FLOATTYPE logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

    // these probably need to be refactored so that the class just asks for the number of time steps as a cl arg
    std::vector<std::array<osv,NUMPARTS>> sim_future(const osv &yt){
        return this->sim_future_obs(m_dte, yt);
    }


};


template<size_t NUMPARTS>
svol_lw_2_par<NUMPARTS>::svol_lw_2_par(FLOATTYPE delta, FLOATTYPE phi_l, FLOATTYPE phi_u, FLOATTYPE mu_l,
                                             FLOATTYPE mu_u, FLOATTYPE sig_l, FLOATTYPE sig_u, FLOATTYPE rho_l, FLOATTYPE rho_u,
                                             unsigned dte)
        : LWFilter2WithCovsFutureSimulator<NUMPARTS,DIMSTATE,DIMOBS,DIMCOV,DIMPARAM,FLOATTYPE>(
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


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_par<NUMPARTS>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
auto svol_lw_2_par<NUMPARTS>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    return sl::muSamp(untrans_p1, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_par<NUMPARTS>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_par<NUMPARTS>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return sl::logGEv(yt, xt);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_par<NUMPARTS>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    return sl::logFEv(xt, xtm1, cov_data, untrans_pt);
}


template<size_t NUMPARTS>
auto svol_lw_2_par<NUMPARTS>::qSamp(const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) -> ssv
{
    return sl::fSamp(xtm1, cov_data, untrans_pt, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_par<NUMPARTS>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt)
{
    return sl::logFEv(xt, xtm1, cov_data, untrans_pt);
}


template<size_t NUMPARTS>
auto svol_lw_2_par<NUMPARTS>::paramPriorSamp() -> psv
{
    // phi, mu, sigma, rho
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_mu_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    untrans_samp(3) = m_rho_sampler.sample();
    return untrans_samp;
}

template<size_t NUMPARTS>
auto svol_lw_2_par<NUMPARTS>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    return sl::gSamp(xt, m_stdNormSampler.sample());
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t NUMPARTS>
class svol_lw_2_from_csv
        : public LWFilter2WithCovsFutureSimulator<NUMPARTS, DIMSTATE, DIMOBS, DIMCOV, DIMPARAM, FLOATTYPE>
{
//public:
//    using ssv = Eigen::Matrix<FLOATTYPE, DIMSTATE, 1>;
//    using osv = Eigen::Matrix<FLOATTYPE, DIMOBS, 1>;
//    using csv = Eigen::Matrix<FLOATTYPE, DIMCOV,1>;
//    using psv = Eigen::Matrix<FLOATTYPE, DIMPARAM,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<FLOATTYPE> m_stdNormSampler;  

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, FLOATTYPE> m_param_sampler;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

public:
    // ctor
    // make sure parameter samples in csv are for phi,mu,sigma,rho
    svol_lw_2_from_csv(FLOATTYPE delta, const std::string &csv_param_samples, unsigned dte);

    // functions tha twe need to define
    FLOATTYPE logMuEv (const ssv &x1, const psv& untrans_p1) override;
    ssv q1Samp (const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1) override;
    FLOATTYPE logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt) override;
    FLOATTYPE logFEv (const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt) override;
    ssv qSamp (const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    FLOATTYPE logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) override;
    psv paramPriorSamp() override;
    osv gSamp(const ssv &xt, const psv &untrans_pt) override;

    // these probably need to be refactored so that the class just asks for the number of time steps as a cl arg
    std::vector<std::array<osv,NUMPARTS>> sim_future(const osv &yt){
        return this->sim_future_obs(m_dte, yt);
    }

};

template<size_t NUMPARTS>
svol_lw_2_from_csv<NUMPARTS>::svol_lw_2_from_csv(FLOATTYPE delta, const std::string &csv_param_samples,
                                                       unsigned dte)
    : LWFilter2WithCovsFutureSimulator<NUMPARTS,DIMSTATE,DIMOBS,DIMCOV,DIMPARAM,FLOATTYPE>(
    				std::vector<std::string>{"logit", "null", "log", "twice_fisher"}, // phi, mu, sigma, rho
            			delta)           // PRIORS    // REMINDER: does output appear to be extremely sensitive to these?
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_param_sampler(csv_param_samples)
    , m_dte(dte)
{
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_from_csv<NUMPARTS>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
auto svol_lw_2_from_csv<NUMPARTS>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    return sl::muSamp(untrans_p1, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_from_csv<NUMPARTS>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    return sl::logMuEv(x1, untrans_p1);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_from_csv<NUMPARTS>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    return sl::logGEv(yt, xt);
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_from_csv<NUMPARTS>::logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv& untrans_pt)
{
    return sl::logFEv(xt, xtm1, cov_data, untrans_pt);
}


template<size_t NUMPARTS>
auto svol_lw_2_from_csv<NUMPARTS>::qSamp(const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt) -> ssv
{
    return sl::fSamp(xtm1, cov_data, untrans_pt, m_stdNormSampler.sample());
}


template<size_t NUMPARTS>
FLOATTYPE svol_lw_2_from_csv<NUMPARTS>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt, const csv &cov_data, const psv& untrans_pt)
{
    return sl::logFEv(xt, xtm1, cov_data, untrans_pt);
}


template<size_t NUMPARTS>
auto svol_lw_2_from_csv<NUMPARTS>::paramPriorSamp() -> psv
{
    return m_param_sampler.samp();
}

template<size_t NUMPARTS>
auto svol_lw_2_from_csv<NUMPARTS>::gSamp(const ssv &xt, const psv &untrans_pt) -> osv
{
    return sl::gSamp(xt, m_stdNormSampler.sample());
}


/**
 * @brief particle swarm filter (many bootstrap filters)
 * this samples parameters from parameterized distribution
 */
template<size_t n_state_parts, size_t n_param_parts>
class svol_swarm_1 : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMSTATE,FLOATTYPE>>,
        1,
        n_state_parts,
        n_param_parts,
        DIMOBS, DIMSTATE, DIMCOV, DIMPARAM>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMSTATE,FLOATTYPE>>;
    using SwarmBase = SwarmWithCovs<ModType, 1, n_state_parts, n_param_parts, DIMOBS, DIMSTATE, DIMCOV, DIMPARAM>;
    using state_cov_parm_func = typename SwarmBase::state_cov_parm_func;

private:

    // for sampling from the parameter prior
    rvsamp::UniformSampler<FLOATTYPE> m_phi_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_mu_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_sigma_sampler;
    rvsamp::UniformSampler<FLOATTYPE> m_rho_sampler;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

public:

    svol_swarm_1() = delete;

    // ctor
    svol_swarm_1(const std::vector<state_cov_parm_func>& fs, FLOATTYPE phi_l, FLOATTYPE phi_u, FLOATTYPE mu_l, FLOATTYPE mu_u,
                 FLOATTYPE sig_l, FLOATTYPE sig_u, FLOATTYPE rho_l, FLOATTYPE rho_u, unsigned int dte)
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

template<size_t n_state_parts, size_t n_param_parts>
class svol_swarm_2
        : public SwarmWithCovs<svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMSTATE,FLOATTYPE>>,
                1,
                n_state_parts,
                n_param_parts,
                DIMOBS, DIMSTATE, DIMCOV, DIMPARAM>
{
public:

    using ModType = svol_leverage<n_state_parts,mn_resamp_fast1<n_state_parts,DIMSTATE,FLOATTYPE>>;
    using SwarmBase = SwarmWithCovs<ModType, 1, n_state_parts, n_param_parts, DIMOBS, DIMSTATE, DIMCOV, DIMPARAM>;
    using state_cov_parm_func = typename SwarmBase::state_cov_parm_func;

private:

    // for sampling from mcmc samples
    utils::csv_param_sampler<DIMPARAM, FLOATTYPE> m_param_sampler;

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
