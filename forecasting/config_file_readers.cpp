#include "config_file_readers.h"

#include <fstream>
#include <sstream>
#include <cmath> // round



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

    if( std::abs(mid_spread/min_tick - std::round(mid_spread/ min_tick)) > 0.0)
        throw std::runtime_error("spread must be a multiple of theminimum tick\n");
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

    if( std::abs(mid_spread/min_tick - std::round(mid_spread/ min_tick)) > 0.0)
        throw std::runtime_error("spread must be a multiple of theminimum tick\n");
}

template<typename float_t>
void ConfigType2<float_t>::set_config_params(float_t& delta, std::string& param_samples_filename, unsigned& dte) const
{
	delta = m_delta;
	param_samples_filename = m_param_samples_filename;
	dte = m_dte;
}


















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

    if( std::abs(mid_spread/min_tick - std::round(mid_spread/ min_tick)) > 0.0)
        throw std::runtime_error("spread must be a multiple of theminimum tick\n");
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

    if( std::abs(mid_spread/min_tick - std::round(mid_spread/ min_tick)) > 0.0)
        throw std::runtime_error("spread must be a multiple of theminimum tick\n");
}

template<typename float_t>
void ConfigType4<float_t>::set_config_params(std::string& param_samples_filename, unsigned& dte) const
{
	param_samples_filename = m_param_samples_filename;
	dte = m_dte;
}














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

    if( std::abs(mid_spread/min_tick - std::round(mid_spread/ min_tick)) > 0.0)
        throw std::runtime_error("spread must be a multiple of theminimum tick\n");
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


