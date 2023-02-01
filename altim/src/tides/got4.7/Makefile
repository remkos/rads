include $(RADSROOT)/src/config.mk

test_gottide.out:	test_gottide
	test_gottide > $@

test_perth.out:	test_perth fort.30
	test_perth > $@

test_perth:	perth3.f test_perth.f
	$(FC) -O -o $@ test_perth.f perth3.f

test_gottide:	test_gottide.f90 lpeqmt.f $(TIDES)
	$(FC) -O -o $@ test_gottide.f90 lpeqmt.f $(TIDES) $(NETCDF)

fort.30:
	zcat $(RADSROOT)/altim/data/GOT4.7/grids_oceantide/{q1,o1,p1,s1,k1,n2,m2,s2,k2,m4}.d.gz > fort.30
