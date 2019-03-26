# Definición de macros para definir las carpetas de trabajo
BIN=./bin
OBJ=./obj
SRC=./src
INC=./include
LIB=./lib
ZIP=./zip
IN=./input

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
	./measureLS.sh

stats: $(BIN)/stats
	${BIN}/stats

# Compilamos
$(BIN)/greedy :
	g++ $(OPT) $(SRC)/greedy.cpp -o $(BIN)/greedy

$(BIN)/localSearch :
	g++ $(OPT) $(SRC)/localSearch.cpp -o $(BIN)/localSearch

$(BIN)/stats :
	g++ $(OPT) $(SRC)/stats.cpp -o $(BIN)/stats


examples: $(BIN)/p1
	$(BIN)/greedy < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/localSearch < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/greedy < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/localSearch < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/greedy < $(IN)/MDG-a_21_n2000_m200.txt
	$(BIN)/localSearch < $(IN)/MDG-a_21_n2000_m200.txt

clean:
	rm -f $(BIN)/*

zip: clean
	mkdir $(ZIP)
	zip -r $(ZIP)/$(PRJ).zip *
