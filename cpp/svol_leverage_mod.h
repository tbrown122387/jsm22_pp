#ifndef SVOL_LEVERAGE_MOD_H
#define SVOL_LEVERAGE_MOD_H


#include <Eigen/Dense>

#include <pf/bootstrap_filter_with_covariates.h> // the boostrap particle filter
#include <pf/rv_samp.h> 			 // for sampling random numbers
#include <pf/rv_eval.h> 			 // for evaluating densities and pmfs


using namespace pf;
using namespace pf::filters;
using namespace pf::bases;


// questions:
// 1. Does FutureSimulator need to be refactored for pfs with covariates?

/**
 * NB: fSamp is reserved name for ForwardMod, so we call non-Markovian transition density stateTransSamp()
 *
 */
template<size_t nparts, typename resampT, typename float_t>
class svol_leverage : public BSFilterWC<nparts, 1, 1, 1, resampT, float_t>
                    , public GenFutureSimulator<1,1,float_t,nparts>
{
public:
    using ssv = Eigen::Matrix<float_t, 1, 1>;
    using osv = Eigen::Matrix<float_t, 1, 1>;
    using cvsv= Eigen::Matrix<float_t,1,1>; 
 
    // parameters
    float_t m_phi;
    float_t m_mu;
    float_t m_sigma;
    float_t m_rho;

    // use this for sampling
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler; // for sampling 

    // ctor
    svol_leverage(const float_t &phi, const float_t &mu, const float_t &sigma, const float_t& rho);
    
    // required by bootstrap filter base class
    float_t logQ1Ev(const ssv &x1, const osv &y1, const cvsv &z1);
    float_t logMuEv(const ssv &x1, const cvsv &z1);
    float_t logGEv(const osv &yt, const ssv &xt, const cvsv& zt);
    auto stateTransSamp(const ssv &xtm1, const cvsv& zt) -> ssv;
    auto q1Samp(const osv &y1, const cvsv& z1) -> ssv;

    // required by FutureSimulator base class
    std::array<ssv,nparts> get_uwtd_samps() const;
    auto gSamp(const ssv &xt) -> osv;
    auto fSamp(const ssv &xtm1, const osv &ytm1) -> ssv;

};


template<size_t nparts, typename resampT, typename float_t>
svol_leverage<nparts, 1, 1, resampT, float_t>::svol_leverage(const float_t &phi, const float_t &mu, const float_t &sigma, const float_t &rho) 
    : m_phi(phi), m_mu(beta), m_sigma(sigma), m_rho(rho)
{
}


template<size_t nparts, typename resampT, typename float_t>
auto svol_leverage<nparts, 1, 1, resampT, float_t>::q1Samp(const osv &y1, const cvsv& z1) -> ssv
{
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * m_sigma / std::sqrt(1.-m_phi*m_phi);
    return x1samp;
}


template<size_t nparts, typename resampT, typename float_t>
auto svol_leverage<nparts, 1, 1, resampT, float_t>::fSamp(const ssv &xtm1, const cvsv& zt) -> ssv
{
    // assuming zt is ytm1
    ssv xtsamp;
    float_t mean =  m_mu + m_phi * (xtm1(0) - m_mu) + m_rho*m_sigma*zt(0)*std::exp(-.5*xtm1(0)); 
    xtsamp(0) = mean + m_stdNormSampler.sample() * m_sigma * std::sqrt( 1.0 - m_phi*m_phi );
    return xtsamp;
}


// TODO this signature needs to work with the above one?
template<size_t nparts, typename resampT, typename float_t>
auto svol_leverage<nparts, 1, 1, resampT, float_t>::fSamp(const ssv &xtm1, const osv &ytm1) -> ssv
{
    ssv xtsamp;
    xtsamp(0) = m_phi * xtm1(0) + m_stdNormSampler.sample() * m_sigma;
    return xtsamp;
}


template<size_t nparts, typename resampT, typename float_t>
float_t svol_leverage<nparts, 1, 1, resampT, float_t>::logGEv(const osv &yt, const ssv &xt, const cvsv& zt)
{
    return rveval::evalUnivNorm<float_t>(
				   yt(0),
                                   0.0,
                                   std::exp(.5*xt(0)),
                                   true);
}


template<size_t nparts, typename resampT, typename float_t>
auto svol_leverage<nparts, 1, 1, resampT, float_t>::gSamp(const ssv &xt) -> osv {
    osv yt;
    yt(0) = m_stdNormSampler.sample() * std::exp(.5*xt(0));
    return yt;
}


template<size_t nparts, typename resampT, typename float_t>
float_t svol_leverage<nparts, 1, 1, resampT, float_t>::logMuEv(const ssv &x1, const cvsv& z1)
{
    return rveval::evalUnivNorm<float_t>(x1(0),
                                   0.0,
                                   m_sigma/std::sqrt(1.0 - m_phi*m_phi),
                                   true);
}


template<size_t nparts, typename resampT, typename float_t>
auto svol_leverage<nparts, 1, 1, resampT, float_t>::muSamp() -> ssv {

    ssv x1; 
    x1(0) = m_stdNormSampler.sample() * m_sigma/std::sqrt(1.0 - m_phi*m_phi);
    return x1;
}


template<size_t nparts, typename resampT, typename float_t>
float_t svol_leverage<nparts, 1, 1, resampT, float_t>::logQ1Ev(const ssv &x1samp, const osv &y1, const cvsv& z1)
{
    return rveval::evalUnivNorm<float_t>(x1samp(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
}


template<size_t nparts, typename resampT, typename float_t>    
auto svol_leverage<nparts, 1, 1, resampT, float_t>::get_uwtd_samps() const -> std::array<ssv,nparts> 
{
    return this->m_particles;
}

#endif //SVOL_LEVERAGE_MOD_H
