      subroutine repack
      include "ingrid.inc"
      integer number,kbi(MBUF),id1,id2,i,nspace,k,kb
      number=0
      idnext=1
      do kb=1,MBUF
         if (used(kb)) then
            number=number+1
            kbi(number)=kb
         endif
      enddo
      call bubble(kbi,id,number)
      do k=1,number
         kb=kbi(k)
         id1=id(kb)-1
         id2=idnext-1
         nspace=nx(kb)*ny(kb)
         if (id2.lt.id1) then
            do i=1,nspace
               a(id2+i)=a(id1+i)
	    enddo
            id(kb)=idnext
         endif
         idnext=idnext+nspace
      enddo
      end
