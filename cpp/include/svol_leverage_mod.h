#ifndef SVOL_LEVERAGE_MOD_H
#define SVOL_LEVERAGE_MOD_H


#include <Eigen/Dense>

#include <pf/bootstrap_filter_with_covariates.h> // the boostrap particle filter
#include <pf/rv_samp.h> 			 // for sampling random numbers
#include <pf/rv_eval.h> 			 // for evaluating densities and pmfs
#include <pf/pf_base.h>                          // for GenFutureSimulator


using namespace pf;
using namespace pf::filters;
using pf::bases::GenFutureSimulator;
using FLOATTYPE = double;

#define DIMSTATE 1
#define DIMOBS 1
#define DIMCOV 1
#define DIMPARAM 4

/**
* inline functions defined for this class and all the swarm and liu/west filter stuff
*/
using ssv  = Eigen::Matrix<FLOATTYPE,1,1>;
using osv  = Eigen::Matrix<FLOATTYPE,1,1>;
using csv = Eigen::Matrix<FLOATTYPE,1,1>;
using psv  = Eigen::Matrix<FLOATTYPE,DIMPARAM,1>;
using psm  = Eigen::Matrix<FLOATTYPE,DIMPARAM,DIMPARAM>;

namespace sl {

    inline FLOATTYPE logMuEv(const ssv &x1, const psv &untrans_p1) {
        // parameter order phi, mu, sigma, rho
        return rveval::evalUnivNorm<FLOATTYPE>(
                x1(0),
                0.0,
                untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0) * untrans_p1(0)),
                true);
    }


    inline FLOATTYPE logGEv(const osv &yt, const ssv &xt) {
        return rveval::evalUnivNorm<FLOATTYPE>(
                yt(0),
                0.0,
                std::exp(.5 * xt(0)),
                true);
    }


    inline ssv propMu(const ssv &xtm1, const csv &cov_data, const psv &untrans_param) {
        // phi, mu, sigma, then rho
        ssv ans;
        ans(0) = untrans_param(1) + untrans_param(0) * (xtm1(0) - untrans_param(1));
        ans(0) += cov_data(0) * untrans_param(3) * untrans_param(2) * std::exp(-.5 * xtm1(0));
        return ans;
    }


    inline FLOATTYPE logFEv(const ssv &xt, const ssv &xtm1, const csv &cov_data, const psv &untrans_param) {
        // phi, mu, sigma, rho
        FLOATTYPE mean = untrans_param(1) + untrans_param(0) * (xtm1(0) - untrans_param(1)) +
                         cov_data(0) * untrans_param(3) * untrans_param(2) * std::exp(-.5 * xtm1(0));
        FLOATTYPE sd = untrans_param(2) * std::sqrt(1.0 - untrans_param(3) * untrans_param(3));
        return rveval::evalUnivNorm<FLOATTYPE>(xt(0), mean, sd, true);
    }


    inline ssv fSamp(const ssv &xtm1, const csv &ytm1, const psv &untrans_params, FLOATTYPE z_sample) {
        // phi, mu, sigma, rho
        ssv xt;
        FLOATTYPE mean = untrans_params(1) + untrans_params(0) * (xtm1(0) - untrans_params(1)) +
                         ytm1(0) * untrans_params(3) * untrans_params(2) * std::exp(-.5 * xtm1(0));
        xt(0) = mean + z_sample * untrans_params(2) * std::sqrt(1.0 - untrans_params(3) * untrans_params(3));
        return xt;
    }


    inline ssv muSamp(const psv &untrans_param, FLOATTYPE z_sample) {
        // phi, mu, sigma, rho
        ssv x1samp;
        x1samp(0) = z_sample * untrans_param(2) / std::sqrt(1.0 - untrans_param(0) * untrans_param(0));
        return x1samp;
    }


    inline osv gSamp(const ssv &xt, FLOATTYPE z_sample) {
        osv yt;
        yt(0) = z_sample * std::exp(.5 * xt(0));
        return yt;
    }

} // namespace sl

/**
 * @brief a particle filter class template for a Hull-White stochastic volatility model
 *
 */
template<size_t nparts, typename resampT>
class svol_leverage : public BSFilterWC<nparts, 1, 1, 1, resampT, FLOATTYPE>
        , public GenFutureSimulator<1,1,FLOATTYPE,nparts>
{
public:
    // parameters
    FLOATTYPE m_phi;
    FLOATTYPE m_mu;
    FLOATTYPE m_sigma;
    FLOATTYPE m_rho;
    psv m_untrans_params;

    // days to expiration (aka how many days into future you're simulating)
    unsigned int m_dte;

    // use this for sampling
    rvsamp::UnivNormSampler<FLOATTYPE> m_stdNormSampler; // for sampling

    // ctor
    svol_leverage() = default;
    svol_leverage(const FLOATTYPE &phi, const FLOATTYPE &mu, const FLOATTYPE &sigma, const FLOATTYPE& rho, unsigned int dte);

    // required by bootstrap filter base class
    FLOATTYPE logQ1Ev(const ssv &x1, const osv &y1, const csv &z1);
    FLOATTYPE logMuEv(const ssv &x1, const csv &z1);
    FLOATTYPE logGEv(const osv &yt, const ssv &xt, const csv& zt);
    auto stateTransSamp(const ssv &xtm1, const csv& zt) -> ssv;
    auto q1Samp(const osv &y1, const csv& z1) -> ssv;

    // required by FutureSimulator base class
    std::array<ssv,nparts> get_uwtd_samps() const;
    auto gSamp(const ssv &xt) -> osv;
    auto fSamp(const ssv &xtm1, const osv &ytm1) -> ssv;

};


template<size_t nparts, typename resampT>
svol_leverage<nparts,resampT>::svol_leverage(const FLOATTYPE &phi, const FLOATTYPE &mu, const FLOATTYPE &sigma,
                                                      const FLOATTYPE &rho, unsigned int dte)
        : m_phi(phi), m_mu(mu), m_sigma(sigma), m_rho(rho), m_dte(dte)
{
    // parameter order phi, mu, sigma, rho
    m_untrans_params << m_phi, m_mu, m_sigma, m_rho;
}


template<size_t nparts, typename resampT>
auto svol_leverage<nparts, resampT>::q1Samp(const osv &y1, const csv& z1) -> ssv
{
    return sl::muSamp(m_untrans_params, m_stdNormSampler.sample());
}


template<size_t nparts, typename resampT>
auto svol_leverage<nparts, resampT>::fSamp(const ssv &xtm1, const csv& zt) -> ssv
{
    return sl::fSamp(xtm1, zt, m_untrans_params, m_stdNormSampler.sample());
}


template<size_t nparts, typename resampT>
FLOATTYPE svol_leverage<nparts, resampT>::logGEv(const osv &yt, const ssv &xt, const csv& zt)
{
    return sl::logGEv(yt,xt);
}


template<size_t nparts, typename resampT>
auto svol_leverage<nparts, resampT>::gSamp(const ssv &xt) -> osv {
    return sl::gSamp(xt, m_stdNormSampler.sample());
}


template<size_t nparts, typename resampT>
FLOATTYPE svol_leverage<nparts, resampT>::logMuEv(const ssv &x1, const csv& z1)
{
    return sl::logMuEv(x1, m_untrans_params);
}


template<size_t nparts, typename resampT>
FLOATTYPE svol_leverage<nparts, resampT>::logQ1Ev(const ssv &x1samp, const osv &y1, const csv& z1)
{
    return sl::logMuEv(x1samp, m_untrans_params);
}


template<size_t nparts, typename resampT>
auto svol_leverage<nparts, resampT>::get_uwtd_samps() const -> std::array<ssv,nparts>
{
    return this->m_particles;
}



#endif //SVOL_LEVERAGE_MOD_H
