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
.PHONY: run_filters
run_filters: data/cond_likes/* data/state_estimates/* 

## Run all of the visualization scripts, checking for up-to-date dep.s first
.PHONY: run_vis
run_vis: plots/*

## remove all target, output and extraneous files
.PHONY: clean
clean: 	rm -f *~ *.Rout *.RData *.docx *.pdf *.html *.RData estimation/samps_* estimation/messages_* data/messages_* data/samps_*

## run the shiny applet 
.PHONY: run_shiny
run_shiny: 	
	@echo "\n\n####### running shiny app ###########\n\n"
	R --quiet -e "shiny::runApp(launch.browser=TRUE)"

##############################
## DATA VISUALIZATION STEPS ##
##############################

#R/vis_mcmc_samps.R
plots/mcmc_vis/* data/mcmc_numerical_diagnostics.txt : data/param_samples.csv R/vis_mcmc_samps.R
	echo "## running MCMC visualization code\n"
	Rscript R/vis_mcmc_samps.R

# R/vis_state_ests.R
plots/state_vis/*: data/state_estimates/* R/vis_state_ests.R
	echo "## running state vis code \n"
	Rscript R/vis_state_ests.R

# R/vis_conditional_likes.R
plots/cond_likes_vis/*: data/cond_likes/* R/vis_conditional_likes.R
	echo "## running conditional likes vis code\n"
	Rscript R/vis_conditional_likes.R


#########################
## DATA ANALYSIS STEPS ##
#########################

# R/run_mcmc.R
data/param_samples.csv data/param_samples.csv.bck data/mcmc_time_taken.txt: R/run_mcmc.R $(cppprog) data/SPY_returns_estimation.csv
	@echo "## running MCMC\n"
	Rscript R/run_mcmc.R

# R/run_conditional_likes.R
data/cond_likes/*: R/run_conditional_likes.R data/SPY_returns.csv configs/* $(cppprog)
	@echo "## running conditional likelihoods\n"
	Rscript R/run_conditional_likes.R

# R/run_state_ests.R
data/state_estimates/*: R/run_state_ests.R data/SPY_returns.csv configs/* $(cppprog)
	@echo "### running state estimates  \n"
	Rscript R/run_state_ests.R


#data/forecast_samps.csv: SPY.csv SPY_returns.csv param_samples.csv 
#	./forecasting/build_all_cpp.sh

#########################
## DATA CLEANING STEPS ##
#########################

# R/prices_to_returns.R
data/SPY_returns_estimation.csv data/SPY_returns.csv: R/prices_to_returns.R data/SPY.csv
	@echo "## updating SPY returns data for estimation and filtering \n"
	Rscript R/prices_to_returns.R

# R/create_configs.R
configs/*.csv: R/create_configs.R data/param_samples.csv
	@echo "## updating config files \n"
	Rscript R/create_configs.R

# update data if it has been 24 hours
data/SPY.csv: last_updated.txt
	@echo "\n######## downloading fresh data and updating SPY.csv #########\n"
	Rscript R/update_csv.R

$(cppprog): cpp/*
	$(cmake) --build cpp/cmake-build-release --target jsmpp_v2 -- -j 6

# signal that it's been >24 hours since data download
last_updated.txt: $(TS24)
	@echo "\n##### updating last_updated.txt#####\n"
	touch "$@"


##################
## Extra Things ##
##################

TS24 := .timestamp24
DUMMY := $(shell touch -d 'yesterday' "$(TS24)")
cppprog := cpp/cmake-build-release/jsmpp_v2
cmake := /home/trb5me_admin/Downloads/clion-2021.3.3/bin/cmake/linux/bin/cmake 
