# Makefile - Partie II, version sequentielle
#
# Compilation : make
# Execution   : ./sequentiel parametres.txt
# Nettoyage   : make clean

CC     = gcc
CFLAGS = -O2 -Wall
LDLIBS = -lm

EXE = sequentiel
SRC = src/main.c src/fonctions.c src/prod_matvec.c src/solveur.c
HDR = src/fonctions.h src/prod_matvec.h src/solveur.h

# On recompile tout des qu un .c ou un .h change (pas de fichiers .o intermediaires)
$(EXE): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $(EXE) $(SRC) $(LDLIBS)

clean:
	rm -f $(EXE) solution.dat

.PHONY: clean
