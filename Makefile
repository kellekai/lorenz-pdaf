FC := mpifort
FCFLAGS := -g -c -fdefault-real-8

OBJECTS = $(patsubst %.F90, %.o, $(wildcard *.F90))
$(info $$OBJECTS is [${OBJECTS}])

PROGRAM := model

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FC) $(FLFLAGS) -o $@ $^

%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod

cleanall:
	rm -f *.o *.mod $(PROGRAM)
