include $(ALTIM)/include/config.mk

EXEC	= $(BIN)/flux2asc $(BIN)/flux2bin \
	$(BIN)/ephem2asc $(BIN)/ephem2bin \
	$(BIN)/slr2asc $(BIN)/slr2bin \
	$(BIN)/res2asc $(BIN)/res2bin \
	$(BIN)/orb2asc $(BIN)/orb2bin

all:	$(EXEC)

clean:
	$(RM) -f *.o $(EXEC)

$(BIN)/flux2asc:	flux2asc.o
	$(FC) -o $@ flux2asc.o

$(BIN)/flux2bin:	flux2bin.o
	$(FC) -o $@ flux2bin.o

$(BIN)/ephem2asc:	ephem2asc.o
	$(FC) -o $@ ephem2asc.o

$(BIN)/ephem2bin:	ephem2bin.o
	$(FC) -o $@ ephem2bin.o

$(BIN)/slr2asc:		slr2asc.o
	$(FC) -o $@ slr2asc.o

$(BIN)/slr2bin:		slr2bin.o
	$(FC) -o $@ slr2bin.o

$(BIN)/res2asc:		res2asc.o
	$(FC) -o $@ res2asc.o

$(BIN)/res2bin:		res2bin.o
	$(FC) -o $@ res2bin.o

$(BIN)/orb2asc:		orb2asc.o
	$(FC) -o $@ orb2asc.o

$(BIN)/orb2bin:		orb2bin.o
	$(FC) -o $@ orb2bin.o
