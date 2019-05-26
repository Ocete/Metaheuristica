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
PRJ=practica3

# Macros
all:

# Compilamos
$(BIN)/BMB :
	g++ $(OPT) $(SRC)/BMB.cpp -o $(BIN)/BMB

$(BIN)/GRASP :
	g++ $(OPT) $(SRC)/GRASP.cpp -o $(BIN)/GRASP

$(BIN)/ILS :
	g++ $(OPT) $(SRC)/ILS.cpp -o $(BIN)/ILS

$(BIN)/ES :
	g++ $(OPT) $(SRC)/ES.cpp -o $(BIN)/ES

$(BIN)/ILS-ES :
	g++ $(OPT) $(SRC)/ILS-ES.cpp -o $(BIN)/ILS-ES


$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats

# Ejemplos

examplesBMB: $(BIN)/BMB
	$(BIN)/BMB < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/BMB < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/BMB < $(IN)/MDG-a_21_n2000_m200.txt

examplesGRASP: $(BIN)/GRASP
	$(BIN)/GRASP < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/GRASP < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/GRASP < $(IN)/MDG-a_21_n2000_m200.txt

examplesILS: $(BIN)/ILS
	$(BIN)/ILS < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/ILS < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/ILS < $(IN)/MDG-a_21_n2000_m200.txt

examplesES: $(BIN)/ES
	$(BIN)/ES < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/ES < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/ES < $(IN)/MDG-a_21_n2000_m200.txt

examplesILS-ES: $(BIN)/ILS-ES
	$(BIN)/ILS-ES < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/ILS-ES < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/ILS-ES < $(IN)/MDG-a_21_n2000_m200.txt


# Experimento 3
# evolution: clean $(BIN)/AM $(BIN)/AMM $(BIN)/AGEu
# 	./evolution.sh

# Estadísticas

measureAll: clean measureAM measureAMM measureLSD measureAGGu measureAGEu measureAGGp measureAGEp

measureBMB: $(BIN)/BMB $(BIN)/stats
	./measure.sh $(BIN)/BMB $(OUT)/BMB.dat 30

measureGRASP: $(BIN)/GRASP $(BIN)/stats
	./measure.sh $(BIN)/GRASP $(OUT)/GRASP.dat 30

measureILS: $(BIN)/ILS $(BIN)/stats
	./measure.sh $(BIN)/ILS $(OUT)/ILS.dat 30

measureES: $(BIN)/ES $(BIN)/stats
	./measure.sh $(BIN)/ES $(OUT)/ES.dat 30

measureILS-ES: $(BIN)/ILS-ES $(BIN)/stats
	./measure.sh $(BIN)/ILS-ES $(OUT)/ILS-ES.dat 30


stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*
