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
all: examplesGreedy examplesLS examplesLSD examplesGS

execute: $(BIN)/greedy $(BIN)/localSearch
	./execute.sh

# Compilamos
$(BIN)/greedy :
	g++ $(OPT) $(SRC)/greedy.cpp -o $(BIN)/greedy

$(BIN)/localSearch :
	g++ $(OPT) $(SRC)/localSearch.cpp -o $(BIN)/localSearch

$(BIN)/localSearchDet :
	g++ $(OPT) $(SRC)/localSearchDet.cpp -o $(BIN)/localSearchDet

$(BIN)/greedySearch :
	g++ $(OPT) $(SRC)/greedySearch.cpp -o $(BIN)/greedySearch

$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats

# Ejemplos

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

examplesGS: $(BIN)/greedySearch
		$(BIN)/greedySearch < $(IN)/GKD-c_1_n500_m50.txt
		$(BIN)/greedySearch < $(IN)/SOM-b_16_n400_m160.txt
		$(BIN)/greedySearch < $(IN)/MDG-a_21_n2000_m200.txt

evolution: $(BIN)/localSearch $(BIN)/localSearchDet $(BIN)/greedySearch
	./evolution.sh

# Estadísticas

measureAll: measureLS measureLSD measureGS

measureLS: $(BIN)/localSearch $(BIN)/stats
	./measure.sh $(BIN)/localSearch $(OUT)/measureLS.dat 500

measureLSD: $(BIN)/localSearchDet $(BIN)/stats
	./measure.sh $(BIN)/localSearchDet $(OUT)/measureLSD.dat 500

measureGS: $(BIN)/greedySearch $(BIN)/stats
	./measure.sh $(BIN)/greedySearch $(OUT)/measureGS.dat 2

stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*

zip: clean
	mkdir $(ZIP)
	zip -r $(ZIP)/$(PRJ).zip *
