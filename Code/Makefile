LFLAG = -lgsl -lgslcblas 
SRC = Framework.C MetaLearner.C MetaMove.C SpeciesDataManager.C	SpeciesDistance.C common/Error.C common/Potential.C	common/SlimFactor.C	common/VariableManager.C common/FactorGraph.C common/Matrix.C common/PotentialManager.C common/Utils.C common/Variable.C

INCLPATH2 = common/
INCLPATH3 = orthomapsparser/

CC=g++
CFLAGS = -g -std=gnu++14 -gdwarf-2

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH2) -I $(INCLPATH3)  $(LFLAG) $(CFLAGS) -o scMTNI



