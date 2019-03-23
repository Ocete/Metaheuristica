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

execute: $(BIN)/p1
	./execute.sh

# Compilamos p1.cpp
$(BIN)/p1 :
	g++ $(OPT) $(SRC)/p1.cpp -o $(BIN)/p1

examples: $(BIN)/p1
	$(BIN)/p1 < $(IN)/GKD-c_1_n500_m50.txt
	$(BIN)/p1 < $(IN)/SOM-b_16_n400_m160.txt
	$(BIN)/p1 < $(IN)/MDG-a_21_n2000_m200.txt

clean:
	rm -f $(BIN)/*

zip: clean
	mkdir $(ZIP)
	zip -r $(ZIP)/$(PRJ).zip *
