      subroutine traj_index
!
!   traj_index reads in traj_dump.dat and finds the array index of each point
!   to output to traj_index.dat
!
      integer inmax,i,npts,flag(4)
      parameter (inmax = 500)
      real phys(4,inmax),index(4,inmax)
      character*13 infile
      character*14 outfile

      data flag/1,1,1,1/

      infile  = 'traj_dump.dat'
      outfile = 'traj_index.dat'

      npts = 0
      open(UNIT=5,FILE=infile)
      do i=1,inmax
        read(UNIT=5,FMT=*,END=10) 
     &  phys(1,i),phys(2,i),phys(3,i),phys(4,i)
        npts = npts + 1
      enddo
10    close(5)

      call phys_2_index(phys,index,flag,4,npts)

      open(UNIT=6,FILE=outfile)
      do i = 1,npts
        write(UNIT=6,FMT=*) 
     &  index(1,i),index(2,i),index(3,i),index(4,i)
      enddo
      close(6)

      return
      end
