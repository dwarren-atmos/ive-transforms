

      subroutine buoy_freq(dtdz,thet,nx,ny,nz,nt,dunits,disunits,
     >                             missing)
*
*  buoy_freq calculates the Brunt-Vaisala frequency.  The call
*    parameters and code should be self-evident.
*
      integer nx,ny,nz,nt
      real dtdz(nx,ny,nz,nt),thet(nx,ny,nz+1,nt),missing
      character*(*) dunits,disunits

      real g
      PARAMETER (g = 9.806)
      integer strbeg,strend,i,j,k,n

      do i=1,nx
      do j=1,ny
      do k=1,nz
      do n=1,nt


        if (dtdz(i,j,k,n).gt.0.) then  
          dtdz(i,j,k,n) = sqrt( 2.*g*dtdz(i,j,k,n) / 
     >                   (thet(i,j,k+1,n)+thet(i,j,k,n)))
        else
          dtdz(i,j,k,n) = missing 
        end if

      enddo
      enddo
      enddo
      enddo

      dunits = '1/s'
      disunits = '1/s'

      return
      end

