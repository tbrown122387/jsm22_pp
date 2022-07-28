##############################
## HIGH LEVEL PHONY TARGETS ##
##############################

.PHONY: all
all: refresh_data run_filters run_vis

## Refresh and prepare all data
.PHONY: refresh_data
refresh_data: data/SPY_returns* 


## Run all of the filters on all returns data
# TODO: add forecasting bit when it becomes available
STATE_FILTER_TARGETS = data/state_estimates/lw_aux_prior.txt \
					   data/state_estimates/lw_aux_csv.txt \
					   data/state_estimates/lw2_prior.txt\
					   data/state_estimates/lw2_csv.txt \
					   data/state_estimates/swarm_prior.txt \
					   data/state_estimates/swarm_csv.txt \
					   data/state_estimates/pf_est.txt 
	       
CLIKE_FILTER_TARGETS = data/cond_likes/lw_aux_prior.txt \
					   data/cond_likes/lw_aux_csv.txt \
					   data/cond_likes/lw2_prior.txt\
					   data/cond_likes/lw2_csv.txt \
					   data/cond_likes/swarm_prior.txt \
					   data/cond_likes/swarm_csv.txt \
					   data/cond_likes/pf_est.txt 

LW_POSTERIOR_TARGETS = data/posterior_samps/lw_aux_posterior.txt \
					   data/posterior_samps/lw2_prior_posterior.txt

.PHONY: run_filters
run_filters: $(STATE_FILTER_TARGETS) $(CLIKE_FILTER_TARGETS) $(LW_POSTERIOR_TARGETS)


## Run all of the visualization scripts, checking for up-to-date dep.s first
## TODO write out all MCMC targets
.PHONY: run_vis
run_vis: plots/mcmc_vis/* plots/state_vis/filter_vis.pdf plots/cond_likes_vis/clike_vis.pdf

## remove all target, output and extraneous files
.PHONY: clean
clean: 	
	rm -f *~ *.Rout *.RData *.docx *.pdf *.html *.RData ./messages_* ./samples_* $(STATE_FILTER_TARGETS) $(CLIKE_FILTER_TARGETS)

## run the shiny applet 
.PHONY: run_shiny
run_shiny: 	
	@echo "\n####### running shiny app ###########"
	R --quiet -e "shiny::runApp(launch.browser=TRUE)"

##############################
## DATA VISUALIZATION STEPS ##
##############################

# instructions in R/vis_mcmc_samps.R
plots/mcmc_vis/* data/mcmc_numerical_diagnostics.txt : data/posterior_samps/param_samples.csv R/vis_mcmc_samps.R $(LW_POSTERIOR_TARGETS)
	echo "\n## running MCMC visualization code"
	Rscript R/vis_mcmc_samps.R

# instructions in R/vis_state_ests.R
plots/state_vis/filter_vis.pdf: $(STATE_FILTER_TARGETS) R/vis_state_ests.R
	echo "\n## running state vis code"
	Rscript R/vis_state_ests.R

# R/vis_conditional_likes.R
plots/cond_likes_vis/clike_vis.pdf: $(CLIKE_FILTER_TARGETS) R/vis_conditional_likes.R
	echo "\n## running conditional likes vis code"
	Rscript R/vis_conditional_likes.R


#########################
## DATA ANALYSIS STEPS ##
#########################


# instructions in R/run_mcmc.R
data/posterior_samps/param_samples.csv data/posterior_samps/param_samples.csv.bck data/posterior_samps/mcmc_time_taken.txt: R/run_mcmc.R $(cppprog) data/SPY_returns_estimation.csv
	@echo "\n## running MCMC"
	Rscript R/run_mcmc.R

# instructions in R/run_conditional_likes.R
$(CLIKE_FILTER_TARGETS): R/run_conditional_likes.R data/SPY_returns.csv configs/* $(cppprog)
	@echo "\n## running conditional likelihoods"
	Rscript R/run_conditional_likes.R

# instructions in R/run_state_ests.R
$(STATE_FILTER_TARGETS): R/run_state_ests.R data/SPY_returns.csv configs/* $(cppprog)
	@echo "\n### running state estimates"
	Rscript R/run_state_ests.R

# instructions in R/run_lw_posterior_samples.R
$(LW_POSTERIOR_TARGETS): R/run_lw_posterior_samples.R data/SPY_returns_estimation.csv $(cppprog)
	@echo "\n### running lw filter to get posterior samples"
	Rscript R/run_lw_posterior_samples.R

#data/forecast_samps.csv: SPY.csv SPY_returns.csv param_samples.csv 
#	./forecasting/build_all_cpp.sh

#########################
## DATA CLEANING STEPS ##
#########################

# R/prices_to_returns.R
data/SPY_returns_estimation.csv data/SPY_returns.csv: R/prices_to_returns.R data/SPY.csv
	@echo "\n## updating SPY returns data for estimation and filtering"
	Rscript R/prices_to_returns.R

# R/create_configs.R
configs/*.csv: R/create_configs.R data/posterior_samps/param_samples.csv
	@echo "\n## updating config files"
	Rscript R/create_configs.R

# update data if it has been 24 hours
data/SPY.csv: last_updated.txt
	@echo "\n######## downloading fresh data and updating SPY.csv #########"
	Rscript R/update_csv.R

$(cppprog): cpp/main.cpp
	$(cmake) --build cpp/cmake-build-release --target jsmpp_v2 -- -j 6

# signal that it's been >24 hours since data download
last_updated.txt: $(TS24)
	@echo "\n##### updating last_updated.txt#####"
	touch "$@"


##################
## Extra Things ##
##################

TS24 := .timestamp24
DUMMY := $(shell touch -d 'yesterday' "$(TS24)")
cppprog := cpp/cmake-build-release/jsmpp_v2
cmake := /home/trb5me_admin/Downloads/clion-2021.3.3/bin/cmake/linux/bin/cmake 
