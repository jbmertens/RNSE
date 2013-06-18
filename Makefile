# Generic Makefile

CC = gcc

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c, %.o, $(SRC))

EXEC = rnse

CC_OPTS = -fopenmp
CC_LINKS = -lm -lfftw3 -fopenmp -lhdf5


# Special settings for running on cwru's cluster
ifeq ($(cluster), 1)
	CC = h5cc
	CC_OPTS = -fopenmp $(FFTW_CFLAGS) $(FFTW_LIBS)
	CC_LINKS = -lm -fopenmp $(FFTW_CFLAGS) $(FFTW_LIBS)
endif

# enable debug mode
ifeq ($(debug), 1)
	CC_OPTIMIZE = -O0 -g -Wall -DDEBUG -std=c99
else
	CC_OPTIMIZE = -O3 -std=c99
endif

# fast math
ifeq ($(fast), 1)
	CC_OPTIMIZE += -ffast-math
endif

# profiling
ifeq ($(prof), 1)
	CC_OPTS += -pg
endif

# Require all object files and then link
all: $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(CC_OPTS) $(CC_OPTIMIZE) $(CC_LINKS)

# Just compile every .c file
%.o: %.c
	$(CC) -c $< $(CC_OPTIMIZE) $(CC_OPTS)

clean:
	rm $(EXEC) $(OBJ)
