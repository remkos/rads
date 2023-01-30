      integer*4 maxpass
      parameter (maxpass=200000)
     
      integer*4	p_sat(maxpass),p_sta(maxpass),p_nr(maxpass),npass
      logical   p_ql(maxpass)
      real*8	p_t0(maxpass),p_t1(maxpass)
      
      common /ccatalog/ p_t0,p_t1,p_ql,p_sat,p_sta,p_nr,npass
     
