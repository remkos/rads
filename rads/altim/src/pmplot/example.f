* file: example.f
*
      program example
      real x(2),y(2),a(100),b(100),scale
      character*80 text

      call getarg(1,text)
      if (text.eq.' ') text='?'
      call pgopen(text)
      call pgvstd
*
* Select projection and window (conic, Europe).
*
      call PMDEF(22,0.,100.,100.)
      call PMSWIN(-10.,40.,30.,75.)
*
* Fill land.
* Draw coastlines (thick) and boundaries (thin, dashed) (WDB 1; all ranks)
*
      call pgscr(2,0.75,0.75,0.75)
      call PMWDB('1.lnd',2,0)
      call pgslw(2)
      call PMWDB('1.cil',1,0)
      call pgslw(1)
      call pgsls(2)
      call PMWDB('1.bdy',1,0)
      call pgsls(1)
*
* Draw two markers (Madrid and Moscow) and connect them by a great circle
*
      x(1)=-3.74
      y(1)=40.48
      x(2)=37.4
      y(2)=56.0
      call grtcir(x(1),y(1),x(2),y(2),100,a,b)
      call PMCONV(100,a,b)
      call pgline(100,a,b)
      a(1)=90.0
      a(2)=90.0
      call PMCVEC(2,x,y,a)
      call pgpt(2,x,y,850)
      call pgptxt(x(1),y(1),a(1),1.2,'Madrid')
      call pgptxt(x(2),y(2),a(2),-.2,'Moscow')
*
* Draw meridians and parallels.
* Draw scale bar, etc.
*
      call PMBOX('BCNGST',10.,5,'BCNGST',10.,5)
      call pgsch(.6)
      call PMBAR('CT',1.,0.,5,4)
      call PMQDEF(text,scale)
      call pgmtxt('R',1.,0.,0.,text)
      call pgend
      end
