LFLAG = -lgsl -lgslcblas 
CC=g++
CFLAGS = -g

BIN=NCALearner
SRC=EdgeList.C Error.C EvidenceManager.C Framework.C Graph.C main.C Matrix.C MemoryCheck.C NCA.C NCALearner.C Variable.C VariableManager.C VariableSelection.C

$(BIN): $(SRC)
	$(CC) $(SRC) $(LFLAG) $(CFLAGS) -o $(BIN)

clean:
	rm -f $(BIN) *~
