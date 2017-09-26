OBJS = main.o calc_fru_avg_local.o calc_fru_flux_local.o Delta_V.o fcn.o fcn_fru.o fru_fire_rxn.o fru_rates_local.o initialize_ran.o initialize_state.o mpi_master.o mpi_slave.o mt19937_local.o read_parameters.o ReleaseUnit.o rk_pd54m.o simfru_local.o spirit_wrapper.o StochModel.o StochModel_IO.o 

###################################################
#   Compilers    #
###################################################
COMP=mpic++
CFLAGS=-O3 
LDFLAGS=-fopenmp
#Set to -openmp if Intel or empty if GNU
OMPFLAGS=
NAME=stoch3d
BOOST=/home/mwalker/boost_1_48_0

##########################################
#               main program             #
##########################################
all:	$(OBJS)
		$(COMP) $(OBJS) $(CFLAGS) $(OMPFLAGS) $(LDFLAGS) -o $(NAME)  

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cc  parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

calc_fru_avg_local.o: calc_fru_avg_local.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

calc_fru_flux_local.o: calc_fru_flux_local.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

Delta_V.o: Delta_V.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

fcn.o: fcn.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

fcn_fru.o: fcn_fru.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

fru_fire_rxn.o: fru_fire_rxn.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

fru_rates_local.o: fru_rates_local.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

initialize_ran.o: initialize_ran.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

initialize_state.o: initialize_state.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

mpi_master.o: mpi_master.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

mpi_slave.o: mpi_slave.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

mt19937_local.o: mt19937_local.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

read_parameters.o: read_parameters.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

ReleaseUnit.o: ReleaseUnit.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

rk_pd54m.o: rk_pd54m.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

simfru_local.o: simfru_local.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

spirit_wrapper.o: spirit_wrapper.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  -I$(BOOST) $<

StochModel.o: StochModel.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<

StochModel_IO.o: StochModel_IO.cc parameters.h indices.h optionparser.h ReleaseUnit.h spirit_wrapper.h StochModel.h ExceptionBase.h
	$(COMP) -c $(CFLAGS) $(OMPFLAGS)  $<
