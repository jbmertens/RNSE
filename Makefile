# Generic Makefile

CC= gcc

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c, %.o, $(SRC))

EXEC = rnse

CC_OPTS = -lm -fopenmp -lhdf5

# enable debug mode
ifeq ($(debug), 1)
	CC_OPTIMIZE = -O0 -g -Wall -DDEBUG
else
	CC_OPTIMIZE = -O2
endif

# Require all object files and then link
all: $(OBJ)
	$(CC) -o $(EXEC) $(CC_OPTS) $(CC_OPTIMIZE) $(OBJ)

# Just compile every .c file
%.o: %.c
	$(CC) -c $(CC_OPTS) $(CC_OPTIMIZE) -c $<

clean:
	rm $(EXEC) *.o
