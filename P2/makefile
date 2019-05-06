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
PRJ=practica2

# Macros
all:

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

$(BIN)/AMM :
	g++ $(OPT) $(SRC)/AMM.cpp -o $(BIN)/AMM

$(BIN)/LSD :
	g++ $(OPT) $(SRC)/LSD.cpp -o $(BIN)/LSD

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

examplesAMM: $(BIN)/AMM
	$(BIN)/AMM < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/AMM < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/AMM < $(IN)/MDG-a_21_n2000_m200.txt

examplesLSD: $(BIN)/LSD
	$(BIN)/LSD < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/LSD < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/LSD < $(IN)/MDG-a_21_n2000_m200.txt

# Practica 1
# evolution: $(BIN)/localSearch $(BIN)/localSearchDet $(BIN)/greedySearch
# 	./evolution.sh

# Estadísticas

measureAll: clean measureAM measureAMM measureLSD measureAGGu measureAGEu measureAGGp measureAGEp

measureAGGu: $(BIN)/AGGu $(BIN)/stats
	./measure.sh $(BIN)/AGGu $(OUT)/AGGu.dat 30

measureAGGp: $(BIN)/AGGp $(BIN)/stats
	./measure.sh $(BIN)/AGGp $(OUT)/AGGp.dat 30

measureAGEp: $(BIN)/AGEp $(BIN)/stats
	./measure.sh $(BIN)/AGEp $(OUT)/AGEp.dat 30

measureAGEu: $(BIN)/AGEu $(BIN)/stats
	./measure.sh $(BIN)/AGEu $(OUT)/AGEu.dat 30

measureAM: $(BIN)/AM $(BIN)/stats
	./measure.sh $(BIN)/AM $(OUT)/AM2.dat 30

measureAMM: $(BIN)/AMM $(BIN)/stats
	./measure.sh $(BIN)/AMM $(OUT)/AMM.dat 30

measureLSD: $(BIN)/LSD $(BIN)/stats
	./measure.sh $(BIN)/AMM $(OUT)/AMM.dat 30

stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*
