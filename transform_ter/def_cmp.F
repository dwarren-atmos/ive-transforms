
#include "pointer.icl"

      ive_ptr function def_cmp(in1,in2,ndims,dims,stag,min,max,miss,
     >                               dunits,disunits,dim_names)
*
*  def_cmp takes two indices (typed as characters) as input and
*    returns a ive_ptr to the appropriate component of the 
*    deformation tensor
*
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) in1,in2,dunits,disunits,dim_names(4)

      ive_ptr duptr,dvptr,dwptr,deriv,getmem
      character*5 sub1,sub2,tmp
      real del
      integer strbeg,strend

      def_cmp = 0

      call setstr(sub1,in1)
      call setstr(sub2,in2)

      if ((sub2(1:1).eq.'1').or.
     >          ((sub2(1:1).eq.'2').and.(sub1(1:1).eq.'3'))) then
        tmp = sub1
        sub1 = sub2
        sub2 = tmp
      end if

      if (sub1(1:1).eq.sub2(1:1)) then

        dwptr = deriv('w','z',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dvptr = deriv('v','y',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        duptr = deriv('u','x',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)

        if (sub1(1:1).eq.'1') then
          call def_aa(%val(duptr),%val(dvptr),%val(dwptr),dims(1),
     >                          dims(2),dims(3),dims(4))
          def_cmp = duptr
          call freemem(dvptr)
          call freemem(dwptr)
        else if (sub1(1:1).eq.'2') then
          call def_aa(%val(dvptr),%val(duptr),%val(dwptr),dims(1),
     >                           dims(2),dims(3),dims(4))
          def_cmp = dvptr
          call freemem(duptr)
          call freemem(dwptr)
        else if (sub1(1:1).eq.'3') then
          call def_aa(%val(dwptr),%val(duptr),%val(dwptr),dims(1),
     >                           dims(2),dims(3),dims(4))
          def_cmp = dwptr
          call freemem(duptr)
          call freemem(dvptr)
        end if

      else if ((sub1(1:1).eq.'1').and.(sub2(1:1).eq.'2')) then

        duptr = deriv('u','y',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dvptr = deriv('v','x',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dims(2) = dims(2) - 2
        def_cmp = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (def_cmp.eq.0) then
          write(*,*) 'Error:  could not allocate space for def_cmp.'
          return
        end if
        
        call def_12(%val(def_cmp),%val(duptr),%val(dvptr),dims(1),
     >                                dims(2),dims(3),dims(4))
        call freemem(duptr)
        call freemem(dvptr)

        del = (max(2)-min(2))/(dims(2)+1)
        min(2) = min(2) + del
        max(2) = max(2) - del

      else if ((sub1(1:1).eq.'1').and.(sub2(1:1).eq.'3')) then
   
        dwptr = deriv('w','x',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        duptr = deriv('u','z',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dims(1) = dims(1) - 2
        def_cmp =  getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (def_cmp.eq.0) then
          write(*,*) 'Error:  could not allocate space for def_cmp.'
          return
        end if
 
        call def_13(%val(def_cmp),%val(duptr),%val(dwptr),dims(1),
     >                                 dims(2),dims(3),dims(4))
        call freemem(duptr)
        call freemem(dwptr)

        del = (max(1)-min(1))/(dims(1)+1)
        min(1) = min(1) + del
        max(1) = max(1) - del

      else if ((sub1(1:1).eq.'2').and.(sub2(1:1).eq.'3')) then

        dwptr = deriv('w','y',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dvptr = deriv('v','z',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        dims(2) = dims(2) - 2
        def_cmp = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (def_cmp.eq.0) then
          write(*,*) 'Error:  could not allocate space for def_cmp.'
          return
        end if 

        call def_23(%val(def_cmp),%val(dvptr),%val(dwptr),dims(1),
     >                                    dims(2),dims(3),dims(4))
        call freemem(dwptr)
        call freemem(dvptr)

        del = (max(2) - min(2))/(dims(2)+1)
        min(2) = min(2) + del
        max(2) = max(2) - del

      else 

        write(*,*) 'Error:  cannot compute deformation component ',
     >            in1(strbeg(in1):strend(in1)),
     >            in2(strbeg(in2):strend(in2))
        def_cmp = 0

       end if

       return
       end                

*************************************************************************

      subroutine def_aa(fld1,fld2,fld3,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real fld1(nx,ny,nz,nt),fld2(nx,ny,nz,nt),fld3(nx,ny,nz,nt)

      integer i,j,k,n

      do i=1,nx
      do j=1,ny
      do k=1,nz
      do n=1,nt

        fld1(i,j,k,n) = 4.*fld1(i,j,k,n)/3. 
     >                    - 2.*(fld2(i,j,k,n)+fld3(i,j,k,n))/3.

      enddo
      enddo
      enddo
      enddo

      return
      end

*********************************************************************
      
      subroutine def_12(def,dudy,dvdx,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real def(nx,ny,nz,nt),dudy(nx+2,ny,nz,nt),dvdx(nx,ny+2,nz,nt)   
 
      integer i,j,k,n

      do i=1,nx
      do j=1,ny
      do k=1,nz
      do n=1,nt

        def(i,j,k,n) = dudy(i+1,j,k,n) + dvdx(i,j+1,k,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

*************************************************************************

      subroutine def_13(def,dudz,dwdx,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real def(nx,ny,nz,nt),dudz(nx+2,ny,nz,nt),dwdx(nx,ny,nz+2,nt)  

      integer i,j,k,n

      do i=1,nx
      do j=1,ny
      do k=1,nz
      do n=1,nt

        def(i,j,k,n) = dudz(i+1,j,k,n) + dwdx(i,j,k+1,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

***************************************************************************

      subroutine def_23(def,dvdz,dwdy,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real def(nx,ny,nz,nt),dvdz(nx,ny+2,nz,nt),dwdy(nx,ny,nz+2,nt)

      integer i,j,k,n

      do i=1,nx
      do j=1,ny
      do k=1,nz
      do n=1,nt

        def(i,j,k,n) = dvdz(i,j+1,k,n) + dwdy(i,j,k+1,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

