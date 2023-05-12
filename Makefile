# Makefile for C libraries
# by Mark A. Post, November 2010

# Sources, Binaries, Libraries
SRCS	= fixmath.c gentables.c
OBJS	= $(SRCS:.c=.o)
LIBS	= libumath.a
BIN	= gentables testmath
LDLIBS	= -lumath -lm

# Compiler, Linker Defines
CC	= gcc
AR	= ar
LDFLAGS	= -L.
ifdef DEBUG
	CFLAGS	= -Wall -g -DDEBUG
else
	CFLAGS	= -Wall -O2
endif

# Declarations
.PHONY: all clean

# Compile and Assemble C Source Files into Object Files
all: $(LIBS) $(BIN)

# Clean Up Objects, Binaries out of source directory
clean:
	rm -f $(OBJS) $(LIBS) $(BIN)

libumath.a: $(OBJS)
	$(AR) rcs $@ $?

gentables: gentables.o
	$(CC) -Wall $? -o $@ -lm

