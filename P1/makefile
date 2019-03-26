# Definición de macros para definir las carpetas de trabajo
BIN=./bin
OBJ=./obj
SRC=./src
INC=./include
LIB=./lib
ZIP=./zip
IN=./input
OUT=./output

# Opciones de compilación
# -Wall muestra todas las advertencias
# -O2 para optimizacion
OPT=-Wall -O2

# Nombre de la práctica
PRJ=practica1

# Macros
all: execute

execute: $(BIN)/greedy $(BIN)/localSearch
	./execute.sh

measureLS: $(BIN)/localSearch $(BIN)/stats
	./measure.sh $(BIN)/localSearch $(OUT)/measureLS.dat 50

measureLSD: $(BIN)/localSearchDet $(BIN)/stats
	./measure.sh $(BIN)/localSearchDet $(OUT)/measureLSD.dat 10

stats: $(BIN)/stats
	${BIN}/stats

# Compilamos
$(BIN)/greedy :
	g++ $(OPT) $(SRC)/greedy.cpp -o $(BIN)/greedy

$(BIN)/localSearch :
	g++ $(OPT) $(SRC)/localSearch.cpp -o $(BIN)/localSearch

$(BIN)/localSearchDet :
	g++ $(OPT) $(SRC)/localSearchDet.cpp -o $(BIN)/localSearchDet

$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats


examplesGreedy: $(BIN)/greedy
	$(BIN)/greedy < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/greedy < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/greedy < $(IN)/MDG-a_21_n2000_m200.txt

examplesLS: $(BIN)/localSearch
	$(BIN)/localSearch < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/localSearch < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/localSearch < $(IN)/MDG-a_21_n2000_m200.txt

examplesLSD: $(BIN)/localSearchDet
	$(BIN)/localSearchDet < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/localSearchDet < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/localSearchDet < $(IN)/MDG-a_21_n2000_m200.txt

clean:
	rm -f $(BIN)/*

zip: clean
	mkdir $(ZIP)
	zip -r $(ZIP)/$(PRJ).zip *
