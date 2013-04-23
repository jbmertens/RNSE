# Generic Makefile

CC = gcc

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c, %.o, $(SRC))

EXEC = rnse

CC_OPTS = -fopenmp
CC_LINKS = -lm -lfftw3 -fopenmp -lhdf5

# enable debug mode
ifeq ($(debug), 1)
	CC_OPTIMIZE = -O0 -g -Wall -DDEBUG -std=c99
else
	CC_OPTIMIZE = -O3 -std=c99
endif

ifeq ($(fast), 1)
	CC_OPTIMIZE += -ffast-math
endif

# Require all object files and then link
all: $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(CC_OPTS) $(CC_OPTIMIZE) $(CC_LINKS)

# Just compile every .c file
%.o: %.c
	$(CC) -c $< $(CC_OPTIMIZE) $(CC_OPTS)

clean:
	rm $(EXEC) $(OBJ)
