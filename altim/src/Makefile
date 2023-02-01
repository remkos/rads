all:
	$(MAKE) -C tools
	$(MAKE) -C rssubs
#	$(MAKE) -C pmplot
	$(MAKE) -C splib
	$(MAKE) -C numlib
	$(MAKE) -C gridlib
	$(MAKE) -C tides
	$(MAKE) -C iono
	$(MAKE) -C iono/iri95
	$(MAKE) -C iono/iri2007
	$(MAKE) -C xover
	$(MAKE) -C pim
	$(MAKE) -C grid
	$(MAKE) -C ingrid
	$(MAKE) -C inimini
	$(MAKE) -C max
	$(MAKE) -C orbits
	$(MAKE) -C geograph
	$(MAKE) -C xrgb
	$(MAKE) -C gdyn
	$(MAKE) -C gdyn/tle2elem
	$(MAKE) -C gdyn/convert
	$(MAKE) -C filter

rads:
	$(MAKE) -C $(RADSROOT)/src install

man:
	$(MAKE) -C rssubs man
	$(MAKE) -C pmplot man
	$(MAKE) -C splib man
	$(MAKE) -C gridlib man
	$(MAKE) -C tides man
	$(MAKE) -C xover man
	$(MAKE) -C pim man

expired:
	$(MAKE) -C gdrselect

update:
	cvs update

clean:
	find . -name "*.o" -o -name "*.a" -o -name "*.mod" -o -name "*genmod*" | xargs rm -f
	$(MAKE) -C tides clean
	$(MAKE) -C pmplot spotless
