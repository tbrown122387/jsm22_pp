#include <string>   // std::string
#include "estimate_univ_svol.h"
#include "data_reader.h"
#include "config_file_readers.h"
#include "forecasting_algos.h"


// template parameters
#define NUMPARTS 500
#define FLOATTYPE double
#define DIMSTATE 1
#define DIMOBS 1
#define DIMCOV 1
#define DIMPARAM 4
#define NUMXPARTS 100
#define NUMTHETAPARTS 100


/**
 * A few options to run with a Hull-White stochastic volatility model:
 *
 * 1. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from prior)
 * 2. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from prior)
 * 3. Run the Liu-West1 (original auxiliary style) filter for simulating future observations (sampling parameters from prior)
 * 4. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from csv)
 * 5. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from csv)
 * 6. Run the Liu-West1 (original auxiliary style) filter for simulating future observations (sampling parameters from csv)
 *
 * 7. Run the Liu-West2 filter for state output (sampling parameters from prior)
 * 8. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from prior)
 * 9. Run the Liu-West2 filter for simulating future observations (sampling parameters from prior)
 * 10. Run the Liu-West2 filter for state output (sampling parameters from csv)
 * 11. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from csv)
 * 12. Run the Liu-West2 filter for simulating future observations (sampling parameters from csv)
 *
 * 13. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from prior)
 * 14. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from prior)
 * 15. Run the Particle Swarm (bootstrap filters) algorithm for simulating future observations (sampling parameters from prior)
 * 16. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from csv)
 * 17. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from csv)
 * 18. Run the Particle Swarm (bootstrap filters) algorithm for simulating future observations (sampling parameters from csv)
 *
 * 19. Run standard auxiliary particle filter for state output (with given parameter estimates).
 * 20. Run standard auxiliary particle filter for conditional likelihoods (with given parameter estimates).
 * 21. Run standard auxiliary particle filter for simulating future observations (with given parameter estimates).
 *
 * 22. Run Pseudo-Marginal Metropolis-Hastings to estimate the parameters of the model on a historical stretch of data.

 * Note: for run modes 1-21, all models and variables are instantiated regardless of run_mode. This is done to ensure
 * that the difference in program run-time is only attributable to the run-time of algorithm. However, if estimation is
 * being performed, then it won't instantiate anything extra.
*/
int main(int argc, char* argv[]){

    ///////////////////////
    // some type aliases //
    ///////////////////////
    using floatT  = double;
    using DynMat  = Eigen::Matrix<floatT,Eigen::Dynamic,Eigen::Dynamic>;
    using osv     = Eigen::Matrix<floatT,DIMOBS,1>;
    using csv     = Eigen::Matrix<floatT,DIMCOV,1>;
    using ssv     = Eigen::Matrix<floatT,DIMSTATE,1>;
    using psv     = Eigen::Matrix<floatT,DIMPARAM,1>;
    using func   = std::function<const DynMat(const ssv&, const csv&, const psv&)>;
    using resampT = mn_resamp_fast1<NUMXPARTS,DIMSTATE,floatT>;

    ///////////////////////////////////////
    // command line arguments            //
    // 1. run_mode                       //
    // 2. name of file with returns data //
    ///////////////////////////////////////
    unsigned int run_mode;
    std::string data_filename;
    if( argc == 3 ){
        run_mode = std::stoi(argv[1]);
        data_filename = argv[2];
    }else{
        std::cout << "invalid command line arguments...need to pass in 1.) run mode, and 2.) name of data file \n";
        return 1;
    }



    ///////////////////////////////////
    // perform real-time forecasting //
    // parse config files            //
    // instantiate variables         //
    // get returns data              //
    ///////////////////////////////////

    // get data
    auto data = readInData<floatT>(data_filename, ',');

    // define filtering function that will be used for all algorithms
    std::vector<func> fs;
    fs.push_back(
            [](const ssv& xt, const csv& zt, const psv& pt)-> DynMat{
                return xt;
            }
    );

    // temporary variables that are rewritten several times
    floatT delta;
    floatT phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u;
    floatT phi, mu, sig, rho;
    std::string param_samples_filename;
    unsigned dte;

    // option 1: (for run modes 1-3,7-9)
    // delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u
    ConfigType1<floatT> cfg1("configs/config1.csv");
    cfg1.set_config_params(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte);
    svol_lw_1_par<NUMXPARTS,floatT> lw1p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); //run mode 1-3
    svol_lw_2_par<NUMXPARTS,floatT> lw2p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); //run mode 7-9

    // option 2: (for run modes 4-6, 10-12)
    // delta, param_samples_filename
    ConfigType2<floatT> cfg2("configs/config2.csv");
    cfg2.set_config_params(delta, param_samples_filename, dte);
    svol_lw_1_from_csv<NUMXPARTS,floatT> lw1c(delta, param_samples_filename, dte); //run mode 4-6
    svol_lw_2_from_csv<NUMXPARTS,floatT> lw2c(delta, param_samples_filename, dte); //run mode 10-12

    // option 3: (for run modes 13-15)
    // phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte
    ConfigType3<floatT> cfg3("configs/config3.csv");
    cfg3.set_config_params(phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte);
    svol_swarm_1<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_par(fs, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); // run mode 13-15

    // option 4: (for run modes 16-18)
    // param_samples_filename, dte
    ConfigType4<floatT> cfg4("configs/config4.csv");
    cfg4.set_config_params(param_samples_filename, dte);
    svol_swarm_2<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_csv(param_samples_filename, fs, dte); // run mode 16-18

    // option 5: (for run modes 19-21)
    // phi, mu, sigma, rho, dte
    ConfigType5<floatT> cfg5("configs/config5.csv");
    cfg5.set_config_params(phi, mu, sig, rho, dte);
    svol_leverage<NUMXPARTS,resampT,floatT> single_pf(phi, mu, sig, rho, dte); // run mode 19-21

    std::vector<int> liu_west_filter_run_modes  = {1,2,3,4,5,6,7,8,9,10,11,12};
    std::vector<int> swarm_run_modes  = {13,14,15,16,17,18};
    std::vector<int> single_pf_run_modes = {19,20,21};
    std::vector<int> estimation_run_mode = {22};
    bool is_liu_west_filter = std::find(liu_west_filter_run_modes.begin(), liu_west_filter_run_modes.end(), run_mode) != liu_west_filter_run_modes.end();
    bool is_swarm_filter = std::find(swarm_run_modes.begin(), swarm_run_modes.end(), run_mode) != swarm_run_modes.end();
    bool is_single_pfilter = std::find(single_pf_run_modes.begin(), single_pf_run_modes.end(), run_mode) != single_pf_run_modes.end();
    bool is_estimation = std::find(estimation_run_mode.begin(), estimation_run_mode.end(), run_mode) != estimation_run_mode.end();
    if( is_liu_west_filter ){
	std::cout << "filtering through data with liu west filter\n";

	for( const auto &yt : data ){
		std::cout << yt << "\n";
	}


    }else if( is_swarm_filter ){
	std::cout << "swarm\n";
    }else if( is_single_pfilter ){
	std::cout << "single\n";
    }else if( is_estimation ){

	std::cout << "estimating with particle Metropolis-Hastings\n";

	// get arguments and call function
        ConfigType6<floatT> cfg("configs/config6.csv");
        unsigned int num_mcmc_iters, num_pfilters;
        cfg.set_config_params(num_mcmc_iters, num_pfilters);
        do_ada_pmmh_svol_leverage<NUMPARTS,FLOATTYPE>(
                data_filename,
                "samples",
                "messages",
                num_mcmc_iters,
                num_pfilters,
                false); // use multicore?
    }


    return 0;
}
