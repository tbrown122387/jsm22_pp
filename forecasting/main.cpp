#include <fstream> // for ifstream

#include "data_reader.h"
#include "../shared_cpp/svol_leverage_mod.h"
#include "forecasting_algos.h"


// template parameters
#define DIMSTATE 1 
#define DIMOBS 1
#define DIMPARAM 4
#define NUMXPARTS 100 
#define NUMTHETAPARTS 100



/**
 * Implementation of the particle swarm algorithm on 
 * a simple stochastic volatility model.
 *
 * Different Run configurations:
 *
 * 1. Run the Liu-West2 filter for state output
 * 2. Run the Liu-West2 filter for conditional likelihoods
 * 3. Run the Particle Swarm (bootstrap filters) algorithm for state output
 * 4. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods
 * 5. Run the Particle Swarm (bootstrap filters) algorithm for simulating distant future  
 * 6. Run the Particle Swarm (auxiliary pfilters) algorithm for state output
 * 7. Run the Particle Swarm (auxiliary pfilters) algorithm for conditional likelihoods
 * 8. Run the Particle Swarm (auxiliary pfilters) algorithm for simulating distant future  
 * 9. Run the Liu-West (auxiliary style) filter for state output
 * 10.Run the Liu-West (auxiliary style) filter for conditional likelihoods.

*/
int main(int argc, char* argv[])
{

    ///////////////////////
    // some type aliases //
    ///////////////////////

    using floatT  = double;
//    using resampT = mn_resampler<NUMXPARTS,DIMSTATE,floatT>;
    //using SVOL1   = svol_sisr<NUMXPARTS, DIMSTATE, DIMOBS, resampT, floatT>;
    using SVOL2   = svol_bs<NUMXPARTS, floatT>;
    using obs_sv  = Eigen::Matrix<floatT,DIMOBS,1>;
    using DynMat  = Eigen::Matrix<floatT,Eigen::Dynamic,Eigen::Dynamic>;
    using ssv     = Eigen::Matrix<floatT,DIMSTATE,1>;
    using psv     = Eigen::Matrix<floatT,DIMPARAM,1>;
    using func1   = std::function<const DynMat(const ssv&)>;
    using func2   = std::function<const DynMat(const ssv&, const psv&)>;

    //////////////////////////////////////////
    // get file name from command line args //
    //////////////////////////////////////////
    std::string filename;
    unsigned int daysToExpiration;
    unsigned int run_mode;
    floatT delta;
    if( argc == 5){
        filename = argv[1];
        daysToExpiration = std::stoi(argv[2]);
        run_mode = std::stoi(argv[3]);
        delta = std::stod(argv[4]);
    }else{
        std::cout << "invalid command line arguments...need to pass in teh filename string, and the days to expiration, the run mode and delta (see dox for details)\n";
        return 1;
    }

    ////////////////////
    // do actual work //
    ////////////////////

    // get data
    auto data = readInData<floatT>(filename, ',');

    // make a function for both the swarm and the liu-west filter
    auto my_lambda3 = [](const ssv& xt, const psv& pt) -> const DynMat {
        return xt;
    };
    std::vector<func2> swarm_filt_funcs;
    swarm_filt_funcs.push_back( my_lambda3 );

    // 1.
    // create bootstrap swarm 
    svol_swarm<NUMXPARTS,NUMTHETAPARTS,floatT> swarm(swarm_filt_funcs);

    // 2. Create lw filter
    svol_lw2<NUMXPARTS,floatT> myLWSVOL(delta);

    // 3. create auxiliary swarm
    svol_swarm2<NUMXPARTS,NUMTHETAPARTS,floatT> swarm_apf(swarm_filt_funcs);

    // 4. create regular/auxiliary lw filter
    svol_lw<NUMXPARTS,floatT> myLWSVOLapf(delta);

    // 55555. iterate through data with swarm, updating all models each time point
    for(size_t t = 0; t < data.size(); ++t) {

        // first column is x, second column is y
        if(run_mode == 3 || run_mode == 4 || run_mode == 5){
            swarm.update(data[t].block(0,0,1,1));
        }else if(run_mode == 1 || run_mode == 2){
            myLWSVOL.filter(data[t].block(0,0,1,1), swarm_filt_funcs);
        }else if( run_mode == 6 || run_mode == 7 || run_mode == 8){
            swarm_apf.update(data[t].block(0,0,1,1));
        }else if(run_mode == 9 || run_mode == 10){
            myLWSVOLapf.filter(data[t].block(0,0,1,1), swarm_filt_funcs);
        }

        if(run_mode == 1){
            std::cout << myLWSVOL.getExpectations()[0](0,0) << "\n"; 
        }else if(run_mode == 2){
            std::cout << myLWSVOL.getLogCondLike() << "\n";
        }else if( run_mode == 3){
            std::cout << swarm.getExpectations()[0](0,0) << "\n";
        }else if( run_mode == 4){
            std::cout << swarm.getLogCondLike() << "\n";
        }else if( run_mode == 5){
            // pass...does stuff after the filtering is done
        }else if( run_mode == 6){
            std::cout << swarm_apf.getExpectations()[0](0,0) << "\n";
        }else if( run_mode == 7){
            std::cout << swarm_apf.getLogCondLike() << "\n";
        }else if( run_mode == 8){
            // pass...does stuff after the filtering is done
        }else if(run_mode == 9){
            std::cout << myLWSVOLapf.getExpectations()[0](0,0) << "\n"; 
        }else if(run_mode == 10){
            std::cout << myLWSVOLapf.getLogCondLike() << "\n";


        }else{
            throw std::invalid_argument("invalid run_mode value");
        }

    }
   

    // simulate future stuff in the case of run_mode == 5
    if( run_mode == 5){
        // TODO: fix slow indexing
        auto future_obs = swarm.simFutureObs(daysToExpiration);
        for(unsigned int time = 0; time < daysToExpiration; ++time){
            for(unsigned int prm = 0; prm < NUMTHETAPARTS; ++prm){
                for(unsigned int idx = 0; idx < NUMXPARTS; ++idx){
                    std::cout << future_obs[prm][time][idx] << ", ";
                }
            }
            std::cout << "\n";
        } 
    }

    // simulate future stuff in the case of run_mode == 8
    if( run_mode == 8){
        auto future_obs = swarm_apf.simFutureObs(daysToExpiration);
        for(unsigned int time = 0; time < daysToExpiration; ++time){
            for(unsigned int prm = 0; prm < NUMTHETAPARTS; ++prm){
                for(unsigned int idx = 0; idx < NUMXPARTS; ++idx){
                    std::cout << future_obs[prm][time][idx] << ", ";
                }
            }
            std::cout << "\n";
        } 
    }


    return 0;

}
