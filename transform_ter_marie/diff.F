

      subroutine diff(der,var,ter,gz,ztop,dx,dy,dt,nx,ny,nz,nt,idir)

      integer nx,ny,nz,nt,idir
      real der(nx,ny,nz,nt),var(nx,ny,nz,nt),ter(nx,ny),gz(nz)
      real ztop,dx,dy,dt

      integer i,j,k,n,MAXNZ,kp1,km1
      parameter (MAXNZ=201)

      real zlbnd(MAXNZ),zrbnd(MAXNZ),dzet(MAXNZ),dvdzet,dzetdx,
     >                   dzetdy,dzetdz

      if (nz.gt.201) then
        write(*,*) 'Error:  nz too big for derivative calculation.'
        write(*,*) '        Adjust MAXNZ in diff.F'
      end if

      do k=1,nz-1
        dzet(k) = gz(k+1)-gz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(nz) = 2. 
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      if (idir.eq.1) then

        do n=1,nt
        do k=1,nz
          kp1 = min(k+1,nz)
          km1 = max(k-1,1)
        do j = 1,ny

          do i = 2,nx-1

            dvdzet = zrbnd(k)*(var(i,j,kp1,n)-var(i,j,k,n))
     >                               /(2.*dzet(k))
     >        + zlbnd(k)*(var(i,j,k,n)-var(i,j,km1,n))/(2.*dzet(km1)) 
            dzetdx = (gz(k)-ztop)/(ztop-ter(i,j))
     >                        *0.5*(ter(i+1,j)-ter(i-1,j))/dx
            der(i,j,k,n) = 0.5*(var(i+1,j,k,n)-var(i-1,j,k,n))/dx
     >                      + dvdzet*dzetdx

          enddo

          dvdzet = zrbnd(k)*(var(1,j,kp1,n)-var(1,j,k,n))
     >                                /(2.*dzet(k))
     >      + zlbnd(k)*(var(1,j,k,n)-var(1,j,km1,n))/(2.*dzet(km1))
          dzetdx = (gz(k)-ztop)/(ztop-ter(1,j))
     >                            *(ter(2,j)-ter(1,j))/dx
          der(1,j,k,n) = (var(2,j,k,n)-var(1,j,k,n))/dx
     >                      + dvdzet*dzetdx 

          dvdzet = zrbnd(k)*(var(nx,j,kp1,n)-var(nx,j,k,n))
     >                                /(2.*dzet(k))
     >      + zlbnd(k)*(var(nx,j,k,n)-var(nx,j,km1,n))/(2.*dzet(km1))
          dzetdx = (gz(k)-ztop)/(ztop-ter(nx,j))
     >                            *(ter(nx,j)-ter(nx-1,j))/dx
          der(nx,j,k,n) = (var(nx,j,k,n)-var(nx-1,j,k,n))/dx
     >                      + dvdzet*dzetdx 

        enddo
        enddo
        enddo

      else if (idir.eq.2) then

        do n=1,nt
        do k=1,nz
          kp1 = min(k+1,nz)
          km1 = max(k-1,1)
        do i=1,nx

          do j=2,ny-1

            dvdzet = zrbnd(k)*(var(i,j,kp1,n)-var(i,j,k,n))
     >                             /(2.*dzet(k))
     >        + zlbnd(k)*(var(i,j,k,n)-var(i,j,km1,n))/(2.*dzet(km1))
            dzetdy = (gz(k)-ztop)/(ztop-ter(i,j))
     >                       *0.5*(ter(i,j+1)-ter(i,j-1))/dy
            der(i,j,k,n) = 0.5*(var(i,j+1,k,n)-var(i,j-1,k,n))/dy
     >                             + dvdzet*dzetdy

          enddo

          dvdzet = zrbnd(k)*(var(i,1,kp1,n)-var(i,1,k,n))
     >                                /(2.*dzet(k))
     >       + zlbnd(k)*(var(i,1,k,n)-var(i,1,km1,n))/(2.*dzet(km1))
          dzetdy = (gz(k)-ztop)/(ztop-ter(i,1))
     >                          *(ter(i,2)-ter(i,1))/dy
          der(i,1,k,n) = (var(i,2,k,n)-var(i,1,k,n))/dy
     >                        + dvdzet*dzetdy
          dvdzet = zrbnd(k)*(var(i,ny,kp1,n)-var(i,ny,k,n))
     >                                /(2.*dzet(k))
     >       + zlbnd(k)*(var(i,ny,k,n)-var(i,ny,km1,n))/(2.*dzet(km1))
          dzetdy = (gz(k)-ztop)/(ztop-ter(i,ny))
     >                          *(ter(i,ny)-ter(i,ny-1))/dy
          der(i,ny,k,n) = (var(i,ny,k,n)-var(i,ny-1,k,n))/dy
     >                        + dvdzet*dzetdy

        enddo
        enddo
        enddo

      else if (idir.eq.3) then

        do n=1,nt
        do i=1,nx
        do j=1,ny
 
          dzetdz = ztop/(ztop-ter(i,j))

          do k=2,nz-1

            der(i,j,k,n) = (0.5*(var(i,j,k+1,n)-var(i,j,k,n))/dzet(k)
     >                + 0.5*(var(i,j,k,n)-var(i,j,k-1,n))/dzet(k-1))
     >                             *dzetdz

          enddo

          der(i,j,1,n) = (var(i,j,2,n)-var(i,j,1,n))*dzetdz/dzet(1)
          der(i,j,nz,n) = (var(i,j,nz,n)-var(i,j,nz-1,n))*dzetdz
     >                           /dzet(nz-1)

        enddo
        enddo
        enddo

      else if (idir.eq.4) then

        do i=1,nx
        do j=1,ny
        do k=1,nz

          do n=2,nt-1

            der(i,j,k,n) = 0.5*(var(i,j,k,n+1)-var(i,j,k,n-1))/dt

          enddo

          der(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))/dt
          der(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))/dt

        enddo
        enddo
        enddo

      end if

      return
      end
