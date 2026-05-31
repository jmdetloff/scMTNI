LFLAG = -lgsl -lgslcblas 
SRC = Framework.C MetaLearner.C MetaMove.C SpeciesDataManager.C	SpeciesDistance.C common/Error.C common/Potential.C	common/SlimFactor.C	common/VariableManager.C common/Matrix.C common/PotentialManager.C common/Utils.C common/Variable.C

GSL_PREFIX = $(shell brew --prefix gsl)
LIBPATH = $(GSL_PREFIX)/lib
INCLPATH1 = $(GSL_PREFIX)/include

INCLPATH2 = common/
INCLPATH3 = orthomapsparser/

CC=g++
CFLAGS = -g -std=gnu++14 -gdwarf-2

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2) -I $(INCLPATH3) -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o scMTNI



