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

$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats

# Ejemplos

examplesBMB: $(BIN)/BMB
	$(BIN)/BMB < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/BMB < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/BMB < $(IN)/MDG-a_21_n2000_m200.txt

# Experimento 3
# evolution: clean $(BIN)/AM $(BIN)/AMM $(BIN)/AGEu
# 	./evolution.sh

# Estadísticas

measureAll: clean measureAM measureAMM measureLSD measureAGGu measureAGEu measureAGGp measureAGEp

measureBMB: $(BIN)/BMB $(BIN)/stats
	./measure.sh $(BIN)/BMB $(OUT)/BMB.dat 30

stats: $(BIN)/stats
	${BIN}/stats

clean:
	rm -f $(BIN)/*
