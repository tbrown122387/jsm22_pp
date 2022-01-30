#ifndef SVOL_MODS_H
#define SVOL_MODS_H

#include <limits>
#include <Eigen/Dense>

#include <pf/bootstrap_filter.h> 
#include <pf/auxiliary_pf.h> 
#include <pf/sisr_filter.h>
#include <pf/rv_samp.h> // for sampling random numbers
#include <pf/rv_eval.h> // for evaluating densities and pmfs
#include <pf/resamplers.h>

#include <ssme/liu_west_filter.h>
#include <ssme/pswarm_filter.h>

#define dimx 1
#define dimy 1
#define dimparam 3

#define phi_uniform_lb .94 
#define phi_uniform_ub .98
#define beta_uniform_lb .8
#define beta_uniform_ub 1.1
#define sigma_uniform_lb .07
#define sigma_uniform_ub .17

using namespace pf::filters;
using namespace pf;
using namespace pf::bases;
using namespace pf::resamplers;


/**
 * @brief single auxiliary particle filter for stochastic vol mod
 */
template<size_t nparts, typename float_t>
class svol_apf : public APF<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>,float_t>
               , public FutureSimulator<dimx,dimy,float_t,nparts>
{
public:
    using ModType = APF<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>,float_t>;
    using ssv = typename ModType::state_sized_vec;//Eigen::Matrix<float_t, dimx, 1>;
    using osv = typename ModType::obs_sized_vec;//Eigen::Matrix<float_t, dimy, 1>;

private:
    // parameters
    float_t m_phi;
    float_t m_beta;
    float_t m_sigma;

    // use this for samplign
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler; // for sampling 

public:
    // ctors
    svol_apf();
    svol_apf(const float_t &phi, const float_t &beta, const float_t &sigma);

    // functions tha twe need to define
    float_t logMuEv (const ssv& x1                );
    ssv     propMu  (const ssv& xtm1              );
    ssv     q1Samp  (const osv& y1                );  
    ssv     fSamp   (const ssv& xtm1              );  
    float_t logQ1Ev (const ssv& x1, const osv &y1 );
    float_t logGEv  (const osv& yt, const ssv &xt );

    // extra functions needed for FutureSimulator
    std::array<ssv,nparts> get_uwtd_samps() const;
    osv gSamp(const ssv& xt);
};

template<size_t nparts, typename float_t>
svol_apf<nparts, float_t>::svol_apf()  
    : m_phi(-2.0), m_beta(1.0), m_sigma(1.0)
{
}

template<size_t nparts, typename float_t>
svol_apf<nparts,float_t>::svol_apf(const float_t &phi, const float_t &beta, const float_t &sigma)
    : APF<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>,float_t>()
    , m_phi(phi), m_beta(beta), m_sigma(sigma) 
{
}


template<size_t nparts, typename float_t>
float_t svol_apf<nparts,float_t>::logMuEv(const ssv &x1 )
{
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
}


template<size_t nparts, typename float_t>
auto svol_apf<nparts,float_t>::propMu(const ssv &xtm1) -> ssv
{
    return m_phi * xtm1;
}

template<size_t nparts, typename float_t>
auto svol_apf<nparts,float_t>::q1Samp(const osv &y1) -> ssv
{
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * m_sigma / std::sqrt(1.-m_phi*m_phi);
    return x1samp;
}

template<size_t nparts, typename float_t>
auto svol_apf<nparts,float_t>::fSamp(const ssv &xtm1) -> ssv
{
    ssv xtsamp;
    xtsamp(0) = m_phi * xtm1(0) + m_stdNormSampler.sample()*m_sigma;
    return xtsamp;
}

template<size_t nparts, typename float_t>
float_t svol_apf<nparts,float_t>::logQ1Ev (const ssv &x1, const osv &y1)
{
    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
}


template<size_t nparts, typename float_t>
float_t svol_apf<nparts,float_t>::logGEv (const osv &yt, const ssv &xt)
{
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, m_beta * std::exp(.5*xt(0)), true);
}

template<size_t nparts, typename float_t>
auto svol_apf<nparts,float_t>::get_uwtd_samps() const -> std::array<ssv,nparts>
{
    return this->m_particles;
}


template<size_t nparts, typename float_t>
auto svol_apf<nparts,float_t>::gSamp(const ssv& xt) -> osv
{
    osv ytsamp;
    ytsamp(0) = m_stdNormSampler.sample() * m_beta * std::exp(.5*xt(0));
    return ytsamp;
}


//template<size_t nparts, typename float_t>
//class svol_sisr : public SISRFilter<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>,float_t>
//                , public FutureSimulator<dimx,dimy,float_t,nparts>
//{
//public:
//    using ModType = SISRFilter<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>,float_t>;
//    using ssv = typename ModType::state_sized_vec;//Eigen::Matrix<float_t, dimx, 1>;
//    using osv = typename ModType::obs_sized_vec;//Eigen::Matrix<float_t, dimy, 1>;
//
//private:
//    // parameters
//    float_t m_phi;
//    float_t m_beta;
//    float_t m_sigma;
//
//    // use this for samplign
//    rvsamp::UnivNormSampler<float_t> m_stdNormSampler; // for sampling 
//
//public:
//    // ctor
//    svol_sisr(const float_t &phi, const float_t &beta, const float_t &sigma);
//
//    // functions tha twe need to define
//    float_t logMuEv (const ssv &x1);
//    ssv q1Samp (const osv &y1);    
//    float_t logQ1Ev (const ssv &x1, const osv &y1 );
//    float_t logGEv (const osv &yt, const ssv &xt );
//    float_t logFEv (const ssv &xt, const ssv &xtm1 );
//    ssv qSamp (const ssv &xtm1, const osv &yt );
//    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt );    
//
//    // functions needed for FutureSimulator
//    std::array<ssv,nparts> get_uwtd_samps() const;
//    ssv fSamp(const ssv& xtm1);
//    osv gSamp(const ssv& xt);
//
//    float_t secondPredMoment(const ssv& xt) const { 
//        return std::exp(2.0 * std::log(m_beta) + m_phi*xt(0) + .5*m_sigma*m_sigma); 
//    }
//
//};
//
//
//template<size_t nparts, typename float_t>
//svol_sisr<nparts,float_t>::svol_sisr(const float_t &phi, const float_t &beta, const float_t &sigma)
//    : SISRFilter<nparts, dimx, dimy,mn_resamp_fast1<nparts,dimx,float_t>,float_t>()
//    , m_phi(phi), m_beta(beta), m_sigma(sigma) 
//{
//}
//
//
//template<size_t nparts, typename float_t>
//float_t svol_sisr<nparts,float_t>::logMuEv (const ssv &x1 )
//{
//    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
//}
//
//
//template<size_t nparts, typename float_t>
//auto svol_sisr<nparts,float_t>::q1Samp(const osv &y1) -> ssv
//{
//    ssv x1samp;
//    x1samp(0) = m_stdNormSampler.sample() * m_sigma / std::sqrt(1.-m_phi*m_phi);
//    return x1samp;
//}
//
//template<size_t nparts, typename float_t>
//float_t svol_sisr<nparts,float_t>::logQ1Ev (const ssv &x1, const osv &y1)
//{
//    return rveval::evalUnivNorm<float_t>(x1(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
//}
//
//
//template<size_t nparts, typename float_t>
//float_t svol_sisr<nparts,float_t>::logGEv (const osv &yt, const ssv &xt)
//{
//    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, m_beta * std::exp(.5*xt(0)), true);
//}
//
//
//template<size_t nparts, typename float_t>
//float_t svol_sisr<nparts,float_t>::logFEv(const ssv &xt, const ssv &xtm1)
//{
//    return rveval::evalUnivNorm<float_t>(xt(0), m_phi * xtm1(0), m_sigma, true);
//}
//
//
//template<size_t nparts, typename float_t>
//auto svol_sisr<nparts,float_t>::qSamp(const ssv &xtm1, const osv &yt) -> ssv
//{
//    ssv xtsamp;
//    xtsamp(0) = m_phi * xtm1(0) + m_stdNormSampler.sample()*m_sigma;
//    return xtsamp;
//}
//
//
//template<size_t nparts, typename float_t>
//float_t svol_sisr<nparts,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt)
//{
//    return rveval::evalUnivNorm<float_t>(xt(0), m_phi * xtm1(0), m_sigma, true);
//}
//
//
//template<size_t nparts, typename float_t>
//auto svol_sisr<nparts,float_t>::get_uwtd_samps() const -> std::array<ssv,nparts>
//{
//    return this->m_particles;
//}
//
//
//template<size_t nparts, typename float_t>
//auto svol_sisr<nparts,float_t>::fSamp(const ssv& xtm1) -> ssv
//{
//    ssv xtsamp;
//    xtsamp(0) = m_phi * xtm1(0) + m_stdNormSampler.sample()*m_sigma;
//    return xtsamp;
//}
//
//
//template<size_t nparts, typename float_t>
//auto svol_sisr<nparts,float_t>::gSamp(const ssv& xt) -> osv
//{
//    osv ytsamp;
//    ytsamp(0) = m_stdNormSampler.sample() * m_beta * std::exp(.5*xt(0));
//    return ytsamp;
//}




/**
 * @brief single bootstrap filter for svol mod
 */
template<size_t nparts, typename float_t>
class svol_bs : public BSFilter<nparts, dimx, dimy, mn_resamp_fast1<nparts,dimx,float_t>, float_t>
              , public ForwardMod<dimx,dimy,float_t>
              , public FutureSimulator<dimx,dimy,float_t,nparts>
{
public:
    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
  
    // parameters
    float_t m_phi;
    float_t m_beta;
    float_t m_sigma;

    // use this for samplign
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler; // for sampling 

    // ctors
    svol_bs();
    svol_bs(const float_t &phi, const float_t &beta, const float_t &sigma);
    
    // required by bootstrap filter base class
    float_t logQ1Ev(const ssv &x1, const osv &y1);
    float_t logMuEv(const ssv &x1);
    float_t logGEv(const osv &yt, const ssv &xt);
    auto fSamp(const ssv &xtm1) -> ssv;
    auto q1Samp(const osv &y1) -> ssv;

    // required by ForwardMod<> base class and/or FutureSimulator
    auto muSamp() -> ssv;
    auto gSamp(const ssv &xt) -> osv;
    std::array<ssv,nparts> get_uwtd_samps() const;

    // for expectation approximations
//    float_t secondPredMoment(const ssv& xt) const { 
//        return std::exp(2.0 * std::log(m_beta) + m_phi*xt(0) + .5*m_sigma*m_sigma); 
//    }
   
};

template<size_t nparts, typename float_t>
svol_bs<nparts, float_t>::svol_bs()  
    : m_phi(-2.0), m_beta(1.0), m_sigma(1.0)
{
}

template<size_t nparts, typename float_t>
svol_bs<nparts, float_t>::svol_bs(const float_t &phi, const float_t &beta, const float_t &sigma)  
    : m_phi(phi), m_beta(beta), m_sigma(sigma)
{
}


template<size_t nparts, typename float_t>
auto svol_bs<nparts, float_t>::q1Samp(const osv &y1) -> ssv
{
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * m_sigma / std::sqrt(1.-m_phi*m_phi);
    return x1samp;
}


template<size_t nparts, typename float_t>
auto svol_bs<nparts, float_t>::fSamp(const ssv &xtm1) -> ssv
{
    ssv xtsamp;
    xtsamp(0) = m_phi * xtm1(0) + m_stdNormSampler.sample() * m_sigma;
    return xtsamp;
}


template<size_t nparts, typename float_t>
float_t svol_bs<nparts, float_t>::logGEv(const osv &yt, const ssv &xt)
{
    return rveval::evalUnivNorm<float_t>(
                                   yt(0),
                                   0.0,
                                   m_beta * std::exp(.5*xt(0)),
                                   true);
}


template<size_t nparts, typename float_t>
auto svol_bs<nparts, float_t>::gSamp(const ssv &xt) -> osv {
    osv yt;
    yt(0) = m_stdNormSampler.sample() * m_beta * std::exp(.5*xt(0));
    return yt;
}


template<size_t nparts, typename float_t>
float_t svol_bs<nparts, float_t>::logMuEv(const ssv &x1)
{
    return rveval::evalUnivNorm<float_t>(x1(0),
                                   0.0,
                                   m_sigma/std::sqrt(1.0 - m_phi*m_phi),
                                   true);
}


template<size_t nparts, typename float_t>
auto svol_bs<nparts,float_t>::muSamp() -> ssv {

    ssv x1; 
    x1(0) = m_stdNormSampler.sample() * m_sigma/std::sqrt(1.0 - m_phi*m_phi);
    return x1;
}


template<size_t nparts, typename float_t>
float_t svol_bs<nparts,float_t>::logQ1Ev(const ssv &x1samp, const osv &y1)
{
    return rveval::evalUnivNorm<float_t>(x1samp(0), 0.0, m_sigma/std::sqrt(1.0 - m_phi*m_phi), true);
}


template<size_t nparts, typename float_t>
auto svol_bs<nparts,float_t>::get_uwtd_samps() const -> std::array<ssv,nparts>
{
    return this->m_particles;
}


/**
 * @brief Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t nparts, typename float_t>
class svol_lw : public LWFilter<nparts, dimx, dimy, dimparam, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_beta_sampler;
    //rvsamp::TruncUnivNormSampler<float_t> m_beta_sampler;
    //rvsamp::UnivInvGammaSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
public:
    // ctor
    svol_lw(const float_t &delta);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv propMu  (const ssv &xtm1, const psv& untrans_p1);   //** 
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    ssv fSamp (const ssv &xtm1, const psv& untrans_pt); //**
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    psv paramPriorSamp();

};


template<size_t nparts, typename float_t>
svol_lw<nparts,float_t>::svol_lw(const float_t& delta)
    : LWFilter<nparts,dimx,dimy,dimparam,float_t>(
            std::vector<std::string>{"twice_fisher","log","log"}, // phi, beta, sigma
            delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_phi_sampler(phi_uniform_lb, phi_uniform_ub) 
    , m_beta_sampler(beta_uniform_lb, beta_uniform_ub)  
    , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub)      //  PRIORS
{
}


template<size_t nparts, typename float_t>
float_t svol_lw<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 
                                         0.0,
                                         sd, 
                                         true);
}

template<size_t nparts, typename float_t>
auto svol_lw<nparts,float_t>::propMu(const ssv &xtm1, const psv& untrans_pt) -> ssv
{
    //parameter order: phi, beta, sigma
    return xtm1 * untrans_pt(0);
}

template<size_t nparts, typename float_t>
auto svol_lw<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * sd;
    return x1samp;
}


template<size_t nparts, typename float_t>
auto svol_lw<nparts,float_t>::fSamp(const ssv &xtm1, const psv& untrans_pt) -> ssv
{
    //parameter order: phi, beta, sigma
    ssv xt;
    xt(0) = untrans_pt(0)* xtm1(0) + m_stdNormSampler.sample() * untrans_pt(2);
    return xt;
}

template<size_t nparts, typename float_t>
float_t svol_lw<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 
                                         0.0,
                                         sd, 
                                         true);
}


template<size_t nparts, typename float_t>
float_t svol_lw<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    //parameter order: phi, beta, sigma
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, untrans_pt(1) * std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
auto svol_lw<nparts,float_t>::paramPriorSamp() -> psv
{
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_beta_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    return untrans_samp;
}


/**
 * @brief "alternative" Liu-West filter
 * parameter order: phi, beta, sigma
 */
template<size_t nparts, typename float_t>
class svol_lw2 : public LWFilter2<nparts, dimx, dimy, dimparam, float_t>
{
public:
    using ssv = Eigen::Matrix<float_t, dimx, 1>;
    using osv = Eigen::Matrix<float_t, dimy, 1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;

private:

    // use this for samplign states
    rvsamp::UnivNormSampler<float_t> m_stdNormSampler;  

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_beta_sampler;
    //rvsamp::TruncUnivNormSampler<float_t> m_beta_sampler;
    //rvsamp::UnivInvGammaSampler<float_t> m_sigma_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;
public:
    // ctor
    svol_lw2(const float_t &delta);

    // functions tha twe need to define
    float_t logMuEv (const ssv &x1, const psv& untrans_p1);
    ssv q1Samp (const osv &y1, const psv& untrans_p1);    
    float_t logQ1Ev (const ssv &x1, const osv &y1, const psv& untrans_p1);
    float_t logGEv (const osv &yt, const ssv &xt, const psv& untrans_pt);
    float_t logFEv (const ssv &xt, const ssv &xtm1, const psv& untrans_pt);
    ssv qSamp (const ssv &xtm1, const osv &yt, const psv& untrans_pt);
    float_t logQEv (const ssv &xt, const ssv &xtm1, const osv &yt, const psv& untrans_pt);    
    psv paramPriorSamp();

};


template<size_t nparts, typename float_t>
svol_lw2<nparts,float_t>::svol_lw2(const float_t& delta)
    : LWFilter2<nparts,dimx,dimy,dimparam,float_t>(
            std::vector<std::string>{"twice_fisher","log","log"}, // phi, beta, sigma
            delta)           // PRIORS
    // REMINDER: does output appear to be extremely sensitive to these?
    , m_phi_sampler(phi_uniform_lb, phi_uniform_ub) 
    , m_beta_sampler(beta_uniform_lb, beta_uniform_ub) 
    , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub) 
{
}


template<size_t nparts, typename float_t>
float_t svol_lw2<nparts,float_t>::logMuEv(const ssv &x1, const psv& untrans_p1)
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 
                                         0.0,
                                         sd, 
                                         true);
}


template<size_t nparts, typename float_t>
auto svol_lw2<nparts,float_t>::q1Samp(const osv &y1, const psv& untrans_p1) -> ssv
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    ssv x1samp;
    x1samp(0) = m_stdNormSampler.sample() * sd;
    return x1samp;
}


template<size_t nparts, typename float_t>
float_t svol_lw2<nparts,float_t>::logQ1Ev(const ssv &x1, const osv &y1, const psv& untrans_p1)
{
    //parameter order: phi, beta, sigma
    float_t sd = untrans_p1(2) / std::sqrt(1.0 - untrans_p1(0)*untrans_p1(0));
    return rveval::evalUnivNorm<float_t>(x1(0), 
                                         0.0,
                                         sd, 
                                         true);
}


template<size_t nparts, typename float_t>
float_t svol_lw2<nparts,float_t>::logGEv(const osv &yt, const ssv &xt, const psv& untrans_pt)
{
    //parameter order: phi, beta, sigma
    return rveval::evalUnivNorm<float_t>(yt(0), 0.0, untrans_pt(1) * std::exp(.5*xt(0)), true);
}


template<size_t nparts, typename float_t>
float_t svol_lw2<nparts,float_t>::logFEv(const ssv &xt, const ssv &xtm1, const psv& untrans_pt)
{
    //parameter order: phi, beta, sigma
    return rveval::evalUnivNorm<float_t>(xt(0), untrans_pt(0) * xtm1(0), untrans_pt(2), true);
}


template<size_t nparts, typename float_t>
auto svol_lw2<nparts,float_t>::qSamp(const ssv &xtm1, const osv &yt, const psv& untrans_pt) -> ssv
{
    //parameter order: phi, beta, sigma
    ssv xtsamp;
    xtsamp(0) = untrans_pt(0) * xtm1(0) + m_stdNormSampler.sample()*untrans_pt(2);
    return xtsamp;
}


template<size_t nparts, typename float_t>
float_t svol_lw2<nparts,float_t>::logQEv(const ssv &xt, const ssv &xtm1, const osv &yt, const psv& untrans_pt)
{
    //parameter order: phi, beta, sigma
    return rveval::evalUnivNorm<float_t>(xt(0), untrans_pt(0) * xtm1(0), untrans_pt(2), true);
}

template<size_t nparts, typename float_t>
auto svol_lw2<nparts,float_t>::paramPriorSamp() -> psv
{
    psv untrans_samp;
    untrans_samp(0) = m_phi_sampler.sample();
    untrans_samp(1) = m_beta_sampler.sample();
    untrans_samp(2) = m_sigma_sampler.sample();
    return untrans_samp;
}



/**
 * @brief particle swarm filter (many bootstrap filters)
 * parameter order: phi, beta, sigma
 */
template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm : public Swarm<svol_bs<n_state_parts,float_t>, 
                                1, 
                                n_state_parts, 
                                n_param_parts, 
                                1, 1, 3>
{
public:

    using ModType = svol_bs<n_state_parts,float_t>;
    using SwarmBase = Swarm<ModType, 1, n_state_parts, n_param_parts, 1, 1, 3>;
    using ssv = Eigen::Matrix<float_t,1,1>;
    using osv = Eigen::Matrix<float_t,1,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;
    using state_parm_func = typename SwarmBase::state_parm_func;

private:

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_beta_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;

public:

    // default ctor
    svol_swarm() 
        : SwarmBase() 
        , m_phi_sampler(phi_uniform_lb, phi_uniform_ub)
        , m_beta_sampler(beta_uniform_lb, beta_uniform_ub)
        , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub)
    {}

    // ctor
    svol_swarm(const std::vector<state_parm_func>& fs)
        : SwarmBase(fs)
        , m_phi_sampler(phi_uniform_lb, phi_uniform_ub)
        , m_beta_sampler(beta_uniform_lb, beta_uniform_ub)
        , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub)
    {
    }


    // functions tha twe need to define
    psv samp_untrans_params() override {
        psv untrans_param; //order: phi,beta,sigma
        untrans_param(0) = m_phi_sampler.sample();
        untrans_param(1) = m_beta_sampler.sample();
        untrans_param(2) = m_sigma_sampler.sample();
        return untrans_param; 
    }

    ModType instantiate_mod(const psv& untrans_params) {
        // order: phi, beta, sigma 
        return ModType(untrans_params(0), 
                       untrans_params(1),
                       untrans_params(2));
    }
};




/**
 * @brief particle swarm filter (many auxiliary particle filters)
 * parameter order: phi, beta, sigma
 */
template<size_t n_state_parts, size_t n_param_parts, typename float_t>
class svol_swarm2 : public Swarm<svol_apf<n_state_parts,float_t>, 
                                1, 
                                n_state_parts, 
                                n_param_parts, 
                                1, 1, 3>
{
public:

    using ModType = svol_apf<n_state_parts,float_t>;
    using SwarmBase = Swarm<ModType, 1, n_state_parts, n_param_parts, 1, 1, 3>;
    using ssv = Eigen::Matrix<float_t,1,1>;
    using osv = Eigen::Matrix<float_t,1,1>;
    using psv = Eigen::Matrix<float_t, dimparam,1>;
    using state_parm_func = typename SwarmBase::state_parm_func;

private:

    // for sampling from the parameter prior
    rvsamp::UniformSampler<float_t> m_phi_sampler;
    rvsamp::UniformSampler<float_t> m_beta_sampler;
    rvsamp::UniformSampler<float_t> m_sigma_sampler;

public:

    // default ctor
    svol_swarm2() 
        : SwarmBase() 
        , m_phi_sampler(phi_uniform_lb, phi_uniform_ub)
        , m_beta_sampler(beta_uniform_lb, beta_uniform_ub)
        , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub)
    {}

    // ctor
    svol_swarm2(const std::vector<state_parm_func>& fs)
        : SwarmBase(fs)
        , m_phi_sampler(phi_uniform_lb, phi_uniform_ub)
        , m_beta_sampler(beta_uniform_lb, beta_uniform_ub)
        , m_sigma_sampler(sigma_uniform_lb, sigma_uniform_ub)
    {
    }


    // functions tha twe need to define
    psv samp_untrans_params() override {
        psv untrans_param; //order: phi,beta,sigma
        untrans_param(0) = m_phi_sampler.sample();
        untrans_param(1) = m_beta_sampler.sample();
        untrans_param(2) = m_sigma_sampler.sample();
        return untrans_param; 
    }

    ModType instantiate_mod(const psv& untrans_params) {
        // order: phi, beta, sigma 
        return ModType(untrans_params(0), 
                       untrans_params(1),
                       untrans_params(2));
    }
};

#endif // SVOL_MODS_H
