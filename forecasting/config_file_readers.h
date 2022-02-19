#ifndef CONFIG_FILE_READERS_H
#define CONFIG_FILE_READERS_H

#include <string>
#include <map>



    // all config files are comma separated
    

/**
 * @brief option 1: (for run modes 1-3,7-9)
 * FILE FORMAT:
 * delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte
 *
 */
template<typename float_t>
class ConfigType1 {
public:
    ConfigType1() = delete;
    ConfigType1(const std::string& file);
    float_t m_delta, m_phi_l, m_phi_u, m_mu_l, m_mu_u, m_sig_l, m_sig_u, m_rho_l, m_rho_u;
    unsigned m_dte;
    void set_config_params(float_t& delta, float_t& phi_l, float_t& phi_u, float_t& mu_l, float_t& mu_u, float_t& sig_l, float_t& sig_u, float_t& rho_l, float_t& rho_u, unsigned& dte) const;
};


/**
 * @brief option 2: (for run modes 4-6, 10-12)
 * FILE FORMAT:
 * delta, param_samples_filename, dte
 */
template<typename float_t>
class ConfigType2 {
public:
    ConfigType2() = delete;
    ConfigType2(const std::string& file);
    float_t m_delta;
    unsigned m_dte;
    std::string m_param_samples_filename;
    void set_config_params(float_t& delta, std::string& param_samples_filename, unsigned& dte) const;
};


/**
 * @brief option 3: (for run modes 13-15)
 * FILE FORMAT:
 * phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte;
 */
template<typename float_t>
class ConfigType3 {
public:
    ConfigType3() = delete;
    ConfigType3(const std::string& file);
    float_t m_delta, m_phi_l, m_phi_u, m_mu_l, m_mu_u, m_sig_l, m_sig_u, m_rho_l, m_rho_u;
    unsigned m_dte;
    void set_config_params(float_t& phi_l, float_t& phi_u, float_t& mu_l, float_t& mu_u, float_t& sig_l, float_t& sig_u, float_t& rho_l, float_t& rho_u, unsigned& dte) const;
};


/**
 * @brief option 4: (for run modes 16-18)
 * FILE FORMAT:
 * param_samples_filename, dte;
 */
template<typename float_t>
class ConfigType4 {
public:
    ConfigType4() = delete;
    ConfigType4(const std::string& file);
    std::string m_param_samples_filename;
    unsigned m_dte;
    void set_config_params(std::string& param_samples_filename, unsigned& dte) const;
};



/**
 * @brief option 5: (for run modes 19-21)
 * FILE FORMAT:
 * phi, mu, sigma, rho, dte
 */
template<typename float_t>
class ConfigType5 {
public:
    ConfigType5() = delete;
    ConfigType5(const std::string& file);
    float_t m_phi, m_mu, m_sig, m_rho;
    unsigned m_dte;
    void set_config_params(float_t& phi, float_t& mu, float_t& sig, float_t& rho, unsigned& dte) const;
};




#endif // CONFIG_FILE_READERS_H


