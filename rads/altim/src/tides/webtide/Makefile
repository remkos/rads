include $(ALTIM)/include/config.mk

# Makefile for WebTide

all:	webtide_c webtide_f

test_%:	webtide_%
	webtide_$* $(ALTIM)/data/WebTide/data/sshelf < input_test.txt | diff output_test.txt -

webtide.o webtide_c.o:	webtide.h

webtide_c : webtide.o webtide_c.o
	$(CC) webtide_c.o webtide.o -o $@

webtide_ele : webtide.o webtide_ele.o
	$(CC) webtide_ele.o webtide.o -o $@

webtide_f : webtide.o webtide_f.o
	$(FC) webtide_f.o webtide.o -o $@ $(RSSUBS)
