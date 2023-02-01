* gbfsort.inc

      integer*4 maxpass,maxmeas,list,mode,dt,uniq
      parameter (maxpass=200000,maxmeas=2000000,dt=1800)
      integer*4 p_idsat(maxpass),p_mtype(maxpass),p_istat(maxpass),
     |		p_nr(maxpass)
      real*8	p_t0(maxpass),p_t1(maxpass)
      integer*4 p_dble(maxpass)
      integer*4	s_dble(9999),s_dbnr(9999),s_used(9999),s_usnr(9999)
      real*8	s_t0(9999),s_t1(9999)
      character*68 m_gbf(maxmeas)
      integer*4	m_pass(maxmeas),pnt(maxmeas)
      real*8	m_time(maxmeas),t0,t1,t_t0,t_t1
      integer*4 nmeas,npass,t_dble,t_dbnr,t_used,t_usnr
      common /cpass/ p_idsat,p_mtype,p_istat,p_t0,p_t1,p_dble,p_nr
      common /cstat/ s_t0,s_t1,s_dble,s_dbnr,s_used,s_usnr
      common /cmeas/ m_pass,m_time,m_gbf
      common /cglob/ t0,t1,nmeas,npass,pnt,
     | t_t0,t_t1,t_dble,t_dbnr,t_used,t_usnr,uniq,list,mode
