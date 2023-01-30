include $(ALTIM)/include/config.mk

EXEC	= $(BIN)/boxcar $(BIN)/qfit

all: $(EXEC)

clean:
	$(RM) -f $(EXEC) *.o

$(BIN)/boxcar:	boxcar.o $(RSSUBS)
	$(FC) -o $@ boxcar.o $(RSSUBS)

$(BIN)/qfit:	qfit.o $(RSSUBS)
	$(FC) -o $@ qfit.o $(RSSUBS) $(MATHLIB)
