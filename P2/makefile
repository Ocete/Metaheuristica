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
$(BIN)/AGGp :
	g++ $(OPT) $(SRC)/AGGp.cpp -o $(BIN)/AGGp

$(BIN)/AGGu :
	g++ $(OPT) $(SRC)/AGGu.cpp -o $(BIN)/AGGu

$(BIN)/AGEu :
	g++ $(OPT) $(SRC)/AGEu.cpp -o $(BIN)/AGEu

$(BIN)/AGEp :
	g++ $(OPT) $(SRC)/AGEp.cpp -o $(BIN)/AGEp

$(BIN)/AM :
	g++ $(OPT) $(SRC)/AM.cpp -o $(BIN)/AM



$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats

# Ejemplos

examplesAGGu: $(BIN)/AGGu
	$(BIN)/AGGu < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AGGu < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AGGu < $(IN)/MDG-a_21_n2000_m200.txt

examplesAGGp: $(BIN)/AGGp
	$(BIN)/AGGp < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AGGp < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AGGp < $(IN)/MDG-a_21_n2000_m200.txt

examplesAGEu: $(BIN)/AGEu
	$(BIN)/AGEu < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AGEu < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AGEu < $(IN)/MDG-a_21_n2000_m200.txt

examplesAGEp: $(BIN)/AGEp
	$(BIN)/AGEp < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AGEp < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AGEp < $(IN)/MDG-a_21_n2000_m200.txt

examplesAM: $(BIN)/AM
	$(BIN)/AM < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AM < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AM < $(IN)/MDG-a_21_n2000_m200.txt



evolution: $(BIN)/localSearch $(BIN)/localSearchDet $(BIN)/greedySearch
	./evolution.sh

# Estadísticas

measureAll: measureLS measureLSD measureGS

measureLS: $(BIN)/localSearch $(BIN)/stats
	./measure.sh $(BIN)/localSearch $(OUT)/measureLS.dat 200



stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*
