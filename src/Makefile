TARGET = eos
LIBS = -lm -lgsl -lgslcblas
CC = gcc
CFLAGS = -g -Wall

# Add path for gsl from Homebrew in OSX
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	INCLUDE = /usr/local/include
	LIB_DIR = /usr/local/lib
endif

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@ -I$(INCLUDE)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall -L$(LIB_DIR) $(LIBS) -o $@
	mv $(TARGET) ..

clean:
	-rm -f *.o
