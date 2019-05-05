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


# Practica 1
# evolution: $(BIN)/localSearch $(BIN)/localSearchDet $(BIN)/greedySearch
# 	./evolution.sh

# Estadísticas

$repts="2"

measureAll: clean measureAGGu measureAGEu  # measureAGGp measureAGEp measureAM

measureAGGu: $(BIN)/AGGu $(BIN)/stats
	./measure.sh $(BIN)/AGGu $(OUT)/AGGu.dat $repts

measureAGGp: $(BIN)/AGGp $(BIN)/stats
	./measure.sh $(BIN)/AGGp $(OUT)/AGGp.dat $repts

measureAGEp: $(BIN)/AGEp $(BIN)/stats
	./measure.sh $(BIN)/AGEp $(OUT)/AGEp.dat $repts

measureAGEu: $(BIN)/AGEu $(BIN)/stats
	./measure.sh $(BIN)/AGEu $(OUT)/AGEu.dat $repts

measureAM: $(BIN)/AM $(BIN)/stats
	./measure.sh $(BIN)/AM $(OUT)/AM2.dat $repts



stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*
