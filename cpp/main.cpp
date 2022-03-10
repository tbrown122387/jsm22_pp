#include <string>   // std::string

#include "estimate_univ_svol.h"
#include "data_reader.h"
#include "config_file_readers.h"
#include "forecasting_algos.h"


// template parameters
#define NUMBOTHPARTS 500
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
    using DynMat  = Eigen::Matrix<FLOATTYPE,Eigen::Dynamic,Eigen::Dynamic>;
    using func   = std::function<const DynMat(const ssv&, const csv&, const psv&)>;
    using resampT = mn_resamp_fast1<NUMXPARTS,DIMSTATE,FLOATTYPE>;

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

    ////////////////////////////
    // instantiate everything //
    // parse config files     //
    // get returns data       //
    ////////////////////////////

    // get data
    auto data = readInData<FLOATTYPE>(data_filename, ',');

    // define filtering function that will be used for all algorithms
    std::vector<func> fs;
    fs.emplace_back([](const ssv& xt, const csv& zt, const psv& pt)-> DynMat{
                return xt;
            }
    );

    // temporary variables that are rewritten several times
    FLOATTYPE delta;
    FLOATTYPE phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u;
    FLOATTYPE phi, mu, sig, rho;
    std::string param_samples_filename;
    unsigned dte;

    // option 1: (for run modes 1-3,7-9)
    // delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u
    ConfigType1<FLOATTYPE> cfg1("configs/config1.csv");
    cfg1.set_config_params(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte);
    svol_lw_1_par<NUMBOTHPARTS> lw1p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); //run mode 1-3
    svol_lw_2_par<NUMBOTHPARTS> lw2p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); //run mode 7-9

    // option 2: (for run modes 4-6, 10-12)
    // delta, param_samples_filename
    ConfigType2<FLOATTYPE> cfg2("configs/config2.csv");
    cfg2.set_config_params(delta, param_samples_filename, dte);
    svol_lw_1_from_csv<NUMBOTHPARTS> lw1c(delta, param_samples_filename, dte); //run mode 4-6
    svol_lw_2_from_csv<NUMBOTHPARTS> lw2c(delta, param_samples_filename, dte); //run mode 10-12

    // option 3: (for run modes 13-15)
    // phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte
    ConfigType3<FLOATTYPE> cfg3("configs/config3.csv");
    cfg3.set_config_params(phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte);
    svol_swarm_1<NUMXPARTS,NUMTHETAPARTS> swarm_par(fs, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte); // run mode 13-15

    // option 4: (for run modes 16-18)
    // param_samples_filename, dte
    ConfigType4<FLOATTYPE> cfg4("configs/config4.csv");
    cfg4.set_config_params(param_samples_filename, dte);
    svol_swarm_2<NUMXPARTS,NUMTHETAPARTS> swarm_csv(param_samples_filename, fs, dte); // run mode 16-18

    // option 5: (for run modes 19-21)
    // phi, mu, sigma, rho, dte
    ConfigType5<FLOATTYPE> cfg5("configs/config5.csv");
    cfg5.set_config_params(phi, mu, sig, rho, dte);
    svol_leverage<NUMXPARTS,resampT> single_pf(phi, mu, sig, rho, dte); // run mode 19-21

    /////////////////////////////////
    // filter through returns data //
    // use appropriate algorithm   //
    // which depends on run_mode   //
    /////////////////////////////////
    if( run_mode == 1 ){
        // * 1. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1p.filter(data[i], data[i-1], fs);
            std::cout << lw1p.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 2){
        // * 2. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1p.filter(data[i], data[i-1], fs);
            std::cout << lw1p.getLogCondLike() << "\n";
        }
    }else if( run_mode == 3){
        //  * 3. Run the Liu-West1 (original auxiliary style) filter for simulating future observations (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1p.filter(data[i], data[i-1], fs);
            // TODO
        }
    }else if( run_mode == 4){
        //  * 4. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1c.filter(data[i], data[i-1], fs);
            std::cout << lw1c.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 5){
        //  * 5. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1c.filter(data[i], data[i-1], fs);
            std::cout << lw1c.getLogCondLike() << "\n";
        }
    }else if(run_mode == 6){
        //  * 6. Run the Liu-West1 (original auxiliary style) filter for simulating future observations (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw1c.filter(data[i], data[i-1], fs);
            // TODO std::cout << lw1c.sim_fu
        }
    }else if( run_mode == 7){
        //  * 7. Run the Liu-West2 filter for state output (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2p.filter(data[i], data[i-1], fs);
            std::cout << lw2p.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 8){
        //  * 8. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2p.filter(data[i], data[i-1], fs);
            std::cout << lw2p.getLogCondLike() << "\n";
        }
    }else if( run_mode == 9){
        //  * 9. Run the Liu-West2 filter for simulating future observations (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2p.filter(data[i], data[i-1], fs);
            // TODO
        }
    }else if( run_mode == 10){
        //  * 10. Run the Liu-West2 filter for state output (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2c.filter(data[i], data[i-1], fs);
            std::cout << lw2c.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 11){
        //  * 11. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2c.filter(data[i], data[i-1], fs);
            std::cout << lw2c.getLogCondLike() << "\n";
        }
    }else if( run_mode == 12){
        //  * 12. Run the Liu-West2 filter for simulating future observations (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            lw2c.filter(data[i], data[i-1], fs);
            //TODO
        }
    }else if( run_mode == 13){
        //  * 13. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_par.update(data[i], data[i-1]);
            std::cout << swarm_par.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 14){
        //  * 14. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_par.update(data[i], data[i-1]);
            std::cout << swarm_par.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 15){
        //  * 15. Run the Particle Swarm (bootstrap filters) algorithm for simulating future observations (sampling parameters from prior)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_par.update(data[i], data[i-1]);
            std::cout << swarm_par.getLogCondLike() << "\n";
        }
    }else if( run_mode == 16){
        //  * 16. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_csv.update(data[i], data[i-1]);
            std::cout << swarm_csv.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 17){
        //  * 17. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_csv.update(data[i], data[i-1]);
            std::cout << swarm_csv.getLogCondLike() << "\n";
        }
    }else if( run_mode == 18){
        //  * 18. Run the Particle Swarm (bootstrap filters) algorithm for simulating future observations (sampling parameters from csv)
        for( unsigned i = 1; i < data.size(); ++i){
            swarm_csv.update(data[i], data[i-1]);
            // TODO
        }
    }else if( run_mode == 19){
        //  * 19. Run standard auxiliary particle filter for state output (with given parameter estimates).
        std::vector<std::function<const DynMat(const ssv&, const csv&)>> pf_fs;
        pf_fs.emplace_back([](const ssv& xt, const csv& zt)-> DynMat{
                    return xt;
                }
        );
        for( unsigned i = 1; i < data.size(); ++i ) {
            single_pf.filter(data[i], data[i-1], pf_fs);
            std::cout << single_pf.getExpectations()[0](0,0) << "\n";
        }
    }else if( run_mode == 20){
        //  * 20. Run standard auxiliary particle filter for conditional likelihoods (with given parameter estimates).
        std::vector<std::function<const DynMat(const ssv&, const csv&)>> pf_fs;
        pf_fs.emplace_back([](const ssv& xt, const csv& zt)-> DynMat{
                    return xt;
                }
        );
        for( unsigned i = 1; i < data.size(); ++i ) {
            single_pf.filter(data[i], data[i-1], pf_fs);
            std::cout << single_pf.getLogCondLike() << "\n";
        }
    }else if( run_mode == 21){
        //  * 21. Run standard auxiliary particle filter for simulating future observations (with given parameter estimates).
        std::vector<std::function<const DynMat(const ssv&, const csv&)>> pf_fs;
        pf_fs.emplace_back([](const ssv& xt, const csv& zt)-> DynMat{
                    return xt;
                }
        );
        for( unsigned i = 1; i < data.size(); ++i ) {
            single_pf.filter(data[i], data[i-1], pf_fs);
        }
    }else if( run_mode == 22){
        //  * 22. Run Pseudo-Marginal Metropolis-Hastings to estimate the parameters of the model on a historical stretch of data.

        // get arguments and call function
        ConfigType6<FLOATTYPE> cfg("configs/config6.csv");
        unsigned int num_mcmc_iters, num_pfilters;
        cfg.set_config_params(num_mcmc_iters, num_pfilters);
        do_ada_pmmh_svol_leverage<NUMXPARTS>(
                data_filename,
                "samples",
                "messages",
                num_mcmc_iters,
                num_pfilters,
                false); // use multicore?
    }









    return 0;
}
