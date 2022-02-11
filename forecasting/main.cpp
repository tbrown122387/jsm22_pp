#include <fstream> // for ifstream

#include "data_reader.h"
#include "../shared_cpp/svol_leverage_mod.h"
#include "forecasting_algos.h"


// template parameters
#define DIMSTATE 1 
#define DIMOBS 1
#define DIMCOV 1
#define DIMPARAM 4
#define NUMXPARTS 100 
#define NUMTHETAPARTS 100



/**
 * A few options to run with a Hull-White svol model: 
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
*/
int main(int argc, char* argv[])
{

    ///////////////////////
    // some type aliases //
    ///////////////////////

    using floatT  = double;
    using DynMat  = Eigen::Matrix<floatT,Eigen::Dynamic,Eigen::Dynamic>;
    using osv     = Eigen::Matrix<floatT,DIMOBS,1>;
    using csv     = Eigen::Matrix<floatT,DIMCOV,1>;
    using ssv     = Eigen::Matrix<floatT,DIMSTATE,1>;
    using psv     = Eigen::Matrix<floatT,DIMPARAM,1>;
    using func1   = std::function<const DynMat(const ssv&)>;
    using func2   = std::function<const DynMat(const ssv&, const psv&)>;
    using resampT = mn_resamp_fast1<NUMXPARTS,DIMSTATE,floatT>;

    ////////////////////////////
    // command line arguments //
    ////////////////////////////
    unsigned int run_mode;
    std::string config_file, data_filename;
    if( argc == 4 ){
	run_mode = std::stoi(argv[1]);
	config_file = argv[2];
	data_filanem = argv[3];
    }else{
        std::cout << "invalid command line arguments...need to pass in TODO\n";
        return 1;
    }

    ///////////////////////
    // parse config file //
    /////////////////////// 

    // all config files are comma separated
    // option 1: (run modes 1-3,7-9) 
    // delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u); 
    // option 2: (run modes )
    // delta, param_filename

        date_filename = argv[1];
        daysToExpiration = std::stoi(argv[2]);
        run_mode = std::stoi(argv[3]);
        delta = std::stod(argv[4]);
	param_filename = argv[5];

    std::string data_filename, param_filename;
    unsigned int daysToExpiration;
    unsigned int run_mode;
    floatT delta;

    ////////////////////
    // do actual work //
    ////////////////////

    // get data
    auto data = readInData<floatT>(filename, ',');

    // initialize all models
    svol_lw_1_par<NUMXPARTS,floatT> lw1p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u); //run mode 1-3
    svol_lw_1_csv<NUMXPARTS,floatT> lw1c(delta, param_filename); //run mode 4-6
    svol_lw_2_par<NUMXPARTS,floatT> lw2p(delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u); //run mode 7-9
    svol_lw_2_csv<NUMXPARTS,floatT> lw2c(delta, param_filename); //run mode 10-12
    svol_swarm_1<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_par(phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u); // run mode 13-15
    svol_swarm_2<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_csv(param_filename); // run mode 16-18
    svol_leverage<NUMXPARTS,resampT,floatT> single_pf(phi, mu, sigma, rho); // run mode 19-21

//    // make a function for both the swarm and the liu-west filter
//    auto my_lambda3 = [](const ssv& xt, const psv& pt) -> const DynMat {
//        return xt;
//    };
//    std::vector<func2> swarm_filt_funcs;
//    swarm_filt_funcs.push_back( my_lambda3 );

    
//    // 1.
//    // create bootstrap swarm 
//    svol_swarm<NUMXPARTS,NUMTHETAPARTS,floatT> swarm(swarm_filt_funcs);
//
//    // 2. Create lw filter
//    svol_lw2<NUMXPARTS,floatT> myLWSVOL(delta);
//
//    // 3. create auxiliary swarm
//    svol_swarm2<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_apf(swarm_filt_funcs);
//
//    // 4. create regular/auxiliary lw filter
//    svol_lw<NUMXPARTS,floatT> myLWSVOLapf(delta);
//
//    // 55555. iterate through data with swarm, updating all models each time point
    for(size_t t = 0; t < data.size(); ++t) {

	    std::cout << data[t] << "\n";
    }
//
//        // first column is x, second column is y
//        if(run_mode == 3 || run_mode == 4 || run_mode == 5){
//            swarm.update(data[t].block(0,0,1,1));
//        }else if(run_mode == 1 || run_mode == 2){
//            myLWSVOL.filter(data[t].block(0,0,1,1), swarm_filt_funcs);
//        }else if( run_mode == 6 || run_mode == 7 || run_mode == 8){
//            swarm_apf.update(data[t].block(0,0,1,1));
//        }else if(run_mode == 9 || run_mode == 10){
//            myLWSVOLapf.filter(data[t].block(0,0,1,1), swarm_filt_funcs);
//        }
//
//        if(run_mode == 1){
//            std::cout << myLWSVOL.getExpectations()[0](0,0) << "\n"; 
//        }else if(run_mode == 2){
//            std::cout << myLWSVOL.getLogCondLike() << "\n";
//        }else if( run_mode == 3){
//            std::cout << swarm.getExpectations()[0](0,0) << "\n";
//        }else if( run_mode == 4){
//            std::cout << swarm.getLogCondLike() << "\n";
//        }else if( run_mode == 5){
//            // pass...does stuff after the filtering is done
//        }else if( run_mode == 6){
//            std::cout << swarm_apf.getExpectations()[0](0,0) << "\n";
//        }else if( run_mode == 7){
//            std::cout << swarm_apf.getLogCondLike() << "\n";
//        }else if( run_mode == 8){
//            // pass...does stuff after the filtering is done
//        }else if(run_mode == 9){
//            std::cout << myLWSVOLapf.getExpectations()[0](0,0) << "\n"; 
//        }else if(run_mode == 10){
//            std::cout << myLWSVOLapf.getLogCondLike() << "\n";
//
//
//        }else{
//            throw std::invalid_argument("invalid run_mode value");
//        }
//
//    }
//   
//
//    // simulate future stuff in the case of run_mode == 5
//    if( run_mode == 5){
//        // TODO: fix slow indexing
//        auto future_obs = swarm.simFutureObs(daysToExpiration);
//        for(unsigned int time = 0; time < daysToExpiration; ++time){
//            for(unsigned int prm = 0; prm < NUMTHETAPARTS; ++prm){
//                for(unsigned int idx = 0; idx < NUMXPARTS; ++idx){
//                    std::cout << future_obs[prm][time][idx] << ", ";
//                }
//            }
//            std::cout << "\n";
//        } 
//    }
//
//    // simulate future stuff in the case of run_mode == 8
//    if( run_mode == 8){
//        auto future_obs = swarm_apf.simFutureObs(daysToExpiration);
//        for(unsigned int time = 0; time < daysToExpiration; ++time){
//            for(unsigned int prm = 0; prm < NUMTHETAPARTS; ++prm){
//                for(unsigned int idx = 0; idx < NUMXPARTS; ++idx){
//                    std::cout << future_obs[prm][time][idx] << ", ";
//                }
//            }
//            std::cout << "\n";
//        } 
//    }
//

    return 0;

}
