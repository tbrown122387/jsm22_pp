#ifndef CONFIG_FILE_READERS_H
#define CONFIG_FILE_READERS_H

#include <string>
#include <fstream>
#include <sstream>


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

template<typename float_t>
ConfigType1<float_t>::ConfigType1(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _delta, _phi_l, _phi_u, _mu_l, _mu_u, _sig_l, _sig_u, _rho_l, _rho_u, _dte, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _delta, ',');
            std::getline(stream, _phi_l, ',');
            std::getline(stream, _phi_u, ',');
            std::getline(stream, _mu_l, ',');
            std::getline(stream, _mu_u, ',');
            std::getline(stream, _sig_l, ',');
            std::getline(stream, _sig_u, ',');
            std::getline(stream, _rho_l, ',');
            std::getline(stream, _rho_u, ',');
            std::getline(stream, _dte, ',');

            // store the info
            m_delta = std::stod(_delta);
            m_phi_l = std::stod(_phi_l);
            m_phi_u = std::stod(_phi_u);
            m_mu_l = std::stod(_mu_l);
            m_mu_u = std::stod(_mu_u);
            m_sig_l = std::stod(_sig_l);
            m_sig_u = std::stod(_sig_u);
            m_rho_l = std::stod(_rho_l);
            m_rho_u = std::stod(_rho_u);
            m_dte = std::stoi(_dte);
        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType1<float_t>::set_config_params(float_t& delta, float_t& phi_l, float_t& phi_u, float_t& mu_l, float_t& mu_u, float_t& sig_l, float_t& sig_u, float_t& rho_l, float_t& rho_u, unsigned& dte) const
{
    delta = m_delta;
    phi_l = m_phi_l;
    phi_u = m_phi_u;
    mu_l = m_mu_l;
    mu_u = m_mu_u;
    sig_l = m_sig_l;
    sig_u = m_sig_u;
    rho_l = m_rho_l;
    rho_u = m_rho_u;
    dte = m_dte;
}

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

template<typename float_t>
ConfigType2<float_t>::ConfigType2(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _delta, _param_samples_filename, _dte, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _delta, ',');
            std::getline(stream, _param_samples_filename, ',');
            std::getline(stream, _dte, ',');

            // store the info
            m_delta = std::stod(_delta);
            m_param_samples_filename = _param_samples_filename;
            m_dte = std::stoi(_dte);

        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType2<float_t>::set_config_params(float_t& delta, std::string& param_samples_filename, unsigned& dte) const
{
    delta = m_delta;
    param_samples_filename = m_param_samples_filename;
    dte = m_dte;
}


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

template<typename float_t>
ConfigType3<float_t>::ConfigType3(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _phi_l, _phi_u, _mu_l, _mu_u, _sig_l, _sig_u, _rho_l, _rho_u, _dte, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _phi_l, ',');
            std::getline(stream, _phi_u, ',');
            std::getline(stream, _mu_l, ',');
            std::getline(stream, _mu_u, ',');
            std::getline(stream, _sig_l, ',');
            std::getline(stream, _sig_u, ',');
            std::getline(stream, _rho_l, ',');
            std::getline(stream, _rho_u, ',');
            std::getline(stream, _dte, ',');

            // store the info
            m_phi_l = std::stod(_phi_l);
            m_phi_u = std::stod(_phi_u);
            m_mu_l = std::stod(_mu_l);
            m_mu_u = std::stod(_mu_u);
            m_sig_l = std::stod(_sig_l);
            m_sig_u = std::stod(_sig_u);
            m_rho_l = std::stod(_rho_l);
            m_rho_u = std::stod(_rho_u);
            m_dte = std::stoi(_dte);
        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType3<float_t>::set_config_params(float_t& phi_l, float_t& phi_u, float_t& mu_l, float_t& mu_u, float_t& sig_l, float_t& sig_u, float_t& rho_l, float_t& rho_u, unsigned& dte) const
{
    phi_l = m_phi_l;
    phi_u = m_phi_u;
    mu_l = m_mu_l;
    mu_u = m_mu_u;
    sig_l = m_sig_l;
    sig_u = m_sig_u;
    rho_l = m_rho_l;
    rho_u = m_rho_u;
    dte = m_dte;
}

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


template<typename float_t>
ConfigType4<float_t>::ConfigType4(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _param_samples_filename, _dte, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _param_samples_filename, ',');
            std::getline(stream, _dte, ',');

            // store the info
            m_param_samples_filename = _param_samples_filename;
            m_dte = std::stoi(_dte);

        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType4<float_t>::set_config_params(std::string& param_samples_filename, unsigned& dte) const
{
    param_samples_filename = m_param_samples_filename;
    dte = m_dte;
}

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

template<typename float_t>
ConfigType5<float_t>::ConfigType5(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _phi, _mu, _sig, _rho, _dte, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _phi, ',');
            std::getline(stream, _mu, ',');
            std::getline(stream, _sig, ',');
            std::getline(stream, _rho, ',');
            std::getline(stream, _dte, ',');

            // store the info
            m_phi = std::stod(_phi);
            m_mu = std::stod(_mu);
            m_sig = std::stod(_sig);
            m_rho = std::stod(_rho);
            m_dte = std::stoi(_dte);
        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType5<float_t>::set_config_params(float_t& phi, float_t& mu, float_t& sig, float_t& rho, unsigned& dte) const
{
    phi = m_phi;
    mu = m_mu;
    sig = m_sig;
    rho = m_rho;
    dte = m_dte;
}

/**
 * @brief option 6: (for run mode 22 which is estimation)
 * FILE FORMAT:
 * num_mcmc_iters, num_particle_filters, estimation_data_filename
 */
template<typename float_t>
class ConfigType6 {
public:
    ConfigType6() = delete;
    ConfigType6(const std::string& file);
    unsigned m_num_mcmc_iters, m_num_particle_filters;
    std::string m_est_data_filename;
    void set_config_params(unsigned int &num_mcmc_iters, unsigned int &num_particle_filters, std::string& est_data_filename) const;
};

template<typename float_t>
ConfigType6<float_t>::ConfigType6(const std::string& file)
{

    std::ifstream fs(file);
    if( fs.good() ){

        std::string _num_mcmc_iters, _num_particle_filters, _est_data_filename, line;
        while( std::getline(fs, line) ){

            // split up data
            std::istringstream stream(line);
            std::getline(stream, _num_mcmc_iters, ',');
            std::getline(stream, _num_particle_filters, ',');
            std::getline(stream, _est_data_filename, ',');

            // store the info
            m_num_mcmc_iters = std::stoi(_num_mcmc_iters);
            m_num_particle_filters = std::stoi(_num_particle_filters);
            m_est_data_filename = _est_data_filename;
        }

    }else{
        throw std::runtime_error("could not open symbol information file\n");
    }
}

template<typename float_t>
void ConfigType6<float_t>::set_config_params(unsigned int &num_mcmc_iters, unsigned int &num_particle_filters,
                                             std::string& est_data_filename) const
{
    num_mcmc_iters = m_num_mcmc_iters;
    num_particle_filters = m_num_particle_filters;
    est_data_filename = m_est_data_filename;
}
#endif // CONFIG_FILE_READERS_H


