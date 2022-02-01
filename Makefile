TS24 := .timestamp24
DUMMY := $(shell touch -d 'yesterday' "$(TS24)")

# update data if it has been 24 hours
data/SPY.csv: last_updated.txt
	@echo "\n######## downloading fresh data and updating SPY.csv #########\n"
	Rscript R/update_csv.R

# signal that it's been >24 hours since data download
last_updated.txt: $(TS24)
	@echo "\n##### updating last_updated.txt#####\n"
	touch "$@"

# TODO
estimation/estimation_up_to_date.txt: SPY.csv SPY_returns.csv
	./estimation/build_all_cpp.sh

forecasting/forecasting_up_to_date.txt: SPY.csv SPY_returns.csv param_samples.csv 
	./forecasting/build_all_cpp.sh

.PHONY: run_shiny
run_shiny: 
	@echo "\n\n####### running shiny app ###########\n\n"
	R --quiet -e "shiny::runApp(launch.browser=TRUE)"

## remove all target, output and extraneous files
.PHONY: clean
clean:
	rm -f *~ *.Rout *.RData *.docx *.pdf *.html *-syntax.R *.RData