# Compilateur utilisé
CC=g++

# Options en mode optimisé - La variable NDEBUG est définie
OPTIM_FLAG = -O3 -DNDEBUG -I ~/Info/TPC++/EigenLibrary/Eigen -std=c++11 -Wall
# Options en mode debug
DEBUG_FLAG = -g -I ~/Info/TPC++/EigenLibrary/Eigen -std=c++11 -Wall

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.cc Laplacian2D.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)
	mv $(PROG) .

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
