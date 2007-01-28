
#include "pointer.icl"

      subroutine make_out(outname,idcdf)

#include "default.icl"

      include '/usr/local/include/netcdf.inc'

      character*(*) outname
      integer idcdf

      integer strbeg,strend,ierr,id_dum
      integer iunit,istrt(2),ilen(2),istrt1,ilen1,fptr,k
      real r_val
      logical new,err
      character*(80) cdl_name,cdf_name,tmp,runname,command,command2
      character*(80) grid

      common/flcount/numfls
      integer numfls

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)
*
* check to see if we have space for another file
*
      fptr = numfls+1

      if (fptr.gt.10) then
        write(*,*) 'Error:  too many files open!'
        return
      end if
*
* check to see if file already open
*
      do k=1,numfls
        if (outname(strbeg(outname):strend(outname)).eq.
     >       flnam(k)(strbeg(flnam(k)):strend(flnam(k)))) then
          write(*,*) 'Error:  file ',flnam(k)(strbeg(flnam(k)):
     >            strend(flnam(k))), ' already open.'
          return
        end if
      enddo
*
* open cdl file and set up dimensions
*
      cdl_name(1:len(cdl_name)) = ' '
      cdl_name = outname(strbeg(outname):strend(outname))//'.cdl'

      write(6,*) 
      write(6,*) 'TRANS_OUT:  creating output file'
      write(6,*) 'TRANS_OUT:  cdl file is ',
     >                cdl_name(strbeg(cdl_name):strend(cdl_name))

      iunit = 454
      open(unit=iunit,file=cdl_name(strbeg(cdl_name):
     >                      strend(cdl_name)),status='new')

      write(iunit,*) 'netcdf TransOut{'
      write(iunit,*)

      write(iunit,*) 'dimensions:'
      write(iunit,*) 'nx=',nx,';'
      write(iunit,*) 'ny=',ny,';'
      write(iunit,*) 'nz=',nz,';'
      write(iunit,*) 'nxp1=',nx+1,';'
      if (twod.eq.1) then
        write(iunit,*) 'nyp1=',1,';'
      else
        write(iunit,*) 'nyp1=',ny+1,';'
      end if
      write(iunit,*) 'nzp1=',nz+1,';'
      write(iunit,*) 'one=',1,';'
      write(iunit,*) 'time=UNLIMITED;'
      write(iunit,*)
*
* set up standard variables
*
      write(iunit,*) 'variables:'
      write(iunit,*) 
      write(iunit,*) 'float f_cor(one);'
      write(iunit,*) 'f_cor:units="1/s";'
      write(iunit,*) 'f_cor:def="coriolis parameter";'
      write(iunit,*) 'f_cor:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float ztop(one);'
      write(iunit,*) 'ztop:units="m";'
      write(iunit,*) 'ztop:def="domain depth";'
      write(iunit,*) 'ztop:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float time(time);'
      write(iunit,*) 'time:units="s";'
      write(iunit,*)
      write(iunit,*) 'float sgz(nz);'
      write(iunit,*) 'sgz:units="m";'
      write(iunit,*) 'sgz:def="thermodynamic grid heights";'
      write(iunit,*) 'sgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float wgz(nzp1);'
      write(iunit,*) 'wgz:units="m";'
      write(iunit,*) 'wgz:def="vertical velocity grid heights";'
      write(iunit,*) 'wgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_p(ny,nx);'
      write(iunit,*) 'zbot_p:units="m";'
      write(iunit,*) 'zbot_p:def="terrain at thermo pts";'
      write(iunit,*) 'zbot_p:x_min=',zbot_min(1),';'
      write(iunit,*) 'zbot_p:y_min=',zbot_min(2),';'
      write(iunit,*) 'zbot_p:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_u(ny,nxp1);'
      write(iunit,*) 'zbot_u:units="m";'
      write(iunit,*) 'zbot_u:def="terrain at u pts";'
      write(iunit,*) 'zbot_u:x_min=',0.,';'
      write(iunit,*) 'zbot_u:y_min=',zbot_min(2),';'
      write(iunit,*) 'zbot_u:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_v(nyp1,nx);'
      write(iunit,*) 'zbot_v:units="m";'
      write(iunit,*) 'zbot_v:def="terrain at v pts";'
      write(iunit,*) 'zbot_v:x_min=',zbot_min(1),';'
      write(iunit,*) 'zbot_v:y_min=',0.,';'
      write(iunit,*) 'zbot_v:no_button=1;'
      write(iunit,*)
*
* write out domain info
*
      write(iunit,*) '//global attributes:'
      write(iunit,*)
      write(iunit,*) ':x_min=',plmin(1),';'
      write(iunit,*) ':x_max=',plmax(1),';'
      write(iunit,*) ':x_delta=',dx,';'
      write(iunit,*) ':x_units="m";'
      write(iunit,*) ':x_label="x";'
      write(iunit,*) ':x_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':y_min=',plmin(2),';'
      write(iunit,*) ':y_max=',plmax(2),';'
      write(iunit,*) ':y_delta=',dy,';'
      write(iunit,*) ':y_units="m";'
      write(iunit,*) ':y_label="y";'
      write(iunit,*) ':y_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':z_min=',plmin(3),';'
      write(iunit,*) ':z_max=',plmax(3),';'
      write(iunit,*) ':z_delta=',dz,';'
      write(iunit,*) ':z_units="m";'
      write(iunit,*) ':z_label="z";'
      write(iunit,*) ':z_display_units="km";'
      write(iunit,*)
*
* write runname and grid (analysis or coarse) as comments
*
      runname(1:len(runname)) = ' '
      call ncagtc(idcdf,NCGLOBAL,'runname',runname,len(runname),ierr)
      if (ierr) then
        write(*,*) 'Warning:  unable to find attribute runname.'
      else
        write(iunit,*) ':runname="',runname(strbeg(runname):
     >                  strend(runname)),'";'
        write(iunit,*)
      end if

      write(iunit,*) '}'
*
* close cdl file
*
      close(iunit)
*
* generate netcdf file from cdl file and remove cdl file
*
      cdf_name = outname(strbeg(outname):strend(outname))//'.cdf'
      command(1:len(command)) = ' '
      write(*,*) 'TRANS_OUT:  netcdf file is ',
     >         cdf_name(strbeg(cdf_name):strend(cdf_name))
      write(*,*)
      write(command,*) 'ncgen -o',cdf_name(strbeg(cdf_name):
     >                 strend(cdf_name)),' ',
     >            cdl_name(strbeg(cdl_name):strend(cdl_name))
      call system(command)
      command2(1:len(command2)) = ' '
      write(command2,*) 'rm -f ',cdl_name(strbeg(cdl_name):
     >                      strend(cdl_name))
      call system(command2)
*
* open new netcdf file
*
      id_out(fptr) = ncopn(cdf_name(strbeg(cdf_name):strend(cdf_name)),
     >                       NCWRITE,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  cannot open netcdf file ',
     >        cdf_name(strbeg(cdf_name):strend(cdf_name))
        return
      end if
*
* get IDs for dimensions (used later when defining variables)
*
      x_dims(1,fptr,1) = nx
      x_dims(1,fptr,2) = ncdid(id_out(fptr),'nx',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nx.'
        return
      end if
      x_dims(2,fptr,1) = nx+1
      x_dims(2,fptr,2) = ncdid(id_out(fptr),'nxp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nxp1.'
        return
      end if
      x_dims(3,fptr,1) = 1
      x_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      y_dims(1,fptr,1) = ny
      y_dims(1,fptr,2) = ncdid(id_out(fptr),'ny',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension ny.'
        return
      end if
      if (twod.eq.1) then
        y_dims(2,fptr,1) = 1
      else
        y_dims(2,fptr,1) = ny+1
      end if
      y_dims(2,fptr,2) = ncdid(id_out(fptr),'nyp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nyp1.'
        return
      end if
      y_dims(3,fptr,1) = 1
      y_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      z_dims(1,fptr,1) = nz
      z_dims(1,fptr,2) = ncdid(id_out(fptr),'nz',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nz.'
        return
      end if
      z_dims(2,fptr,1) = nz+1
      z_dims(2,fptr,2) = ncdid(id_out(fptr),'nzp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nzp1.'
        return
      end if
      z_dims(3,fptr,1) = 1
      z_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      t_dims(1,fptr,1) = ntime
      t_dims(1,fptr,2) = ncdid(id_out(fptr),'time',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension time.'
        return
      end if
      t_dims(2,fptr,1) = 1
      t_dims(2,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if
*
* write standard variables to new file
*
      id_dum = ncvid(id_out(fptr),'f_cor',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable f_cor.'
        return
      end if
      call ncvpt1(id_out(fptr),id_dum,1,fcor,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write coriolis frequency.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'ztop',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable ztop.'
        return
      end if
      call ncvpt1(id_out(fptr),id_dum,1,ztop,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable ztop.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'time',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable time.'
        return
      end if
      istrt1 = 1
      ilen1 = ntime 
      call ncvpt(id_out(fptr),id_dum,istrt1,ilen1,%val(time),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable time.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'sgz',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable sgz.'
        return
      end if
      istrt1 = 1
      ilen1 = nz
      call ncvpt(id_out(fptr),id_dum,istrt1,ilen1,%val(sgz),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable sgz.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'wgz',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable wgz.'
        return
      end if
      istrt1 = 1
      ilen1 = nz+1
      call ncvpt(id_out(fptr),id_dum,istrt1,ilen1,%val(wgz),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable wgz.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'zbot_p',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable zbot_p.'
        return
      end if
      istrt(1) = 1
      istrt(2) = 1
      ilen(1) = nx
      ilen(2) = ny
      call ncvpt(id_out(fptr),id_dum,istrt,ilen,%val(zbot_p),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable zbot_p.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'zbot_u',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable zbot_u.'
        return
      end if
      istrt(1) = 1
      istrt(2) = 1
      ilen(1) = nx+1
      ilen(2) = ny
      call ncvpt(id_out(fptr),id_dum,istrt,ilen,%val(zbot_u),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable zbot_u.'
        return
      end if

      id_dum = ncvid(id_out(fptr),'zbot_v',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for variable zbot_v.'
        return
      end if
      istrt(1) = 1
      istrt(2) = 1
      ilen(1) = nx
      if (twod.eq.1) then
        ilen(2) = 1
      else
        ilen(2) = ny+1
      end if
      call ncvpt(id_out(fptr),id_dum,istrt,ilen,%val(zbot_v),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to write variable zbot_v.'
        return
      end if

      numout(fptr) = 8
*
* add new file to list of viable output files
*
      numfls = numfls+1
      flnam(numfls) = outname(strbeg(outname):strend(outname))

      return
      end

*---------------------------------------------------------------------------

      subroutine write_var(outfile,varnm,numvars)

#include "default.icl"

      include '/usr/local/include/netcdf.inc'
      

      integer numvars
      character*(*) outfile,varnm(numvars)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(80) dunits,disunits,dim_names(4)

      integer strbeg,strend,id_dims(4),istrt(4),ilen(4),varid
      integer ierr,ierr1,ierr2,ierr3,ierr4,ierr5,k,n,p,fptr,lgth
      logical new
      ive_ptr var

      common/flcount/numfls
      integer numfls

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2)
      integer t_dims(1,10,2),numout(10)
      real*8 xmin,ymin
      real delx,dely
      character*(80) flnam(10)

      k = 1
      fptr = 0
 345  if ((fptr.eq.0).and.(k.le.numfls)) then
        if (flnam(k)(strbeg(flnam(k)):strend(flnam(k))).eq.
     >         outfile(strbeg(outfile):strend(outfile))) then
          fptr = k
        else
          k = k+1
        end if
        goto 345
      end if

      if (fptr.eq.0) then
        write(*,*) 'Error:  file ',outfile(strbeg(outfile):
     >                 strend(outfile)),' is not a currently open'
        write(*,*) 'output file.'
        return
      end if

      call ncagt(id_out(fptr),NCGLOBAL,'x_min',xmin,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem getting x_min.'
        return
      end if
      call ncagt(id_out(fptr),NCGLOBAL,'y_min',ymin,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem getting y_min.'
        return
      end if
      delx = plmin(1)-xmin
      dely = plmin(2)-ymin

      WRITE(6,*) 'xmin = ',xmin
      WRITE(6,*) 'ymin = ',ymin
      WRITE(6,*) 'delx = ',delx
      WRITE(6,*) 'dely = ',dely

      write(*,*)

      do k=1,numvars

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
        min(1) = min(1)-delx
        min(2) = min(2)-dely
        
        if (var.ne.0) then

          do p=1,4
            id_dims(p) = -1
          enddo
          
          n = 1
 675      if ((id_dims(1).eq.-1).and.(n.le.3)) then
            if (dims(1).eq.x_dims(n,fptr,1)) then
              id_dims(1) = x_dims(n,fptr,2)
            else
              n = n+1
            end if
            goto 675
          end if

          n = 1
 821      if ((id_dims(2).eq.-1).and.(n.le.3)) then
            if (dims(2).eq.y_dims(n,fptr,1)) then
              id_dims(2) = y_dims(n,fptr,2)
            else
              n = n+1
            end if
            goto 821
          end if

          n = 1
 822      if ((id_dims(3).eq.-1).and.(n.le.3)) then
            if (dims(3).eq.z_dims(n,fptr,1)) then
              id_dims(3) = z_dims(n,fptr,2)
            else
              n = n+1
            end if
            goto 822
          end if

          n = 1
 823      if ((id_dims(4).eq.-1).and.(n.le.2)) then
            if (dims(4).eq.t_dims(n,fptr,1)) then
              id_dims(4) = t_dims(n,fptr,2)
            else
              n = n+1
            end if
            goto 823
          end if

          if ((id_dims(1).ge.0).and.(id_dims(2).ge.0).and.
     >          (id_dims(3).ge.0).and.(id_dims(4).ge.0)) then

            call ncredf(id_out(fptr),ierr)
            if (ierr.eq.0) then
              varid = ncvdef(id_out(fptr),
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >              NCFLOAT,4,id_dims,ierr)
              call ncapt(id_out(fptr),varid,'x_min',NCFLOAT,1,
     >                   min(1),ierr1) 
              call ncapt(id_out(fptr),varid,'y_min',NCFLOAT,1,
     >                   min(2),ierr2)
              call ncapt(id_out(fptr),varid,'z_min',NCFLOAT,1,
     >                   min(3),ierr3)
              lgth = strend(dunits) - strbeg(dunits) + 1
              call ncaptc(id_out(fptr),varid,'units',NCCHAR,
     >              lgth,dunits(strbeg(dunits):strend(dunits)),ierr4)
            else
              write(*,*) 'Error:  failed to put file ',
     >                outfile(strbeg(outfile):strend(outfile)),
     >                ' into definition mode.'
            end if
            call ncendf(id_out(fptr),ierr5)

            if ((ierr.eq.0).and.(ierr1.eq.0).and.(ierr2.eq.0)
     >                .and.(ierr3.eq.0).and.(ierr4.eq.0).and.
     >                    (ierr5.eq.0)) then
              do p=1,4
                istrt(p) = 1
                ilen(p) = dims(p)
              enddo
              write(*,*) 'TRANS_OUT:  writing variable ',
     >               varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
              write(*,*) '            to file ',
     >                outfile(strbeg(outfile):strend(outfile))
              call ncvpt(id_out(fptr),varid,istrt,ilen,%val(var),ierr)
              if (ierr.ne.0) then
                write(*,*) 'Error:  failed to write variable ',
     >               varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >               ' to output.'
              end if

            else

              write(*,*) 'Error:  problem either with definition ',
     >             'mode or variable ID for: ' 
              write(*,*) '  file:  ',
     >            outfile(strbeg(outfile):strend(outfile))
              write(*,*) '  variable:  ',
     >            varnm(k)(strbeg(varnm(k)):strend(varnm(k)))

            end if

          else

            write(*,*) 'Error:  failed to find appropriate ',
     >                   'dimensions for variable ',
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k)))

          end if

        else

          write(*,*) 'Error:  unable to access variable ',
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k)))

        end if

      enddo

      write(*,*)

      call ncsnc(id_out(fptr),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem synchronizing variables to disk.'
      end if

      numout(fptr) = numout(fptr)+1

      return
      end

*---------------------------------------------------------------------------

      subroutine close_out(outname)

      character*(*) outname

      integer ierr,i,j,k,n,fptr,strbeg,strend

      common/flcount/numfls
      integer numfls

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2)
      integer t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      k = 1
      fptr = 0
 432  if ((fptr.eq.0).and.(k.le.numfls)) then
        if (flnam(k)(strbeg(flnam(k)):strend(flnam(k))).eq.
     >            outname(strbeg(outname):strend(outname))) then
          fptr = k
        else
          k = k+1
        end if
        goto 432
      end if

      if (fptr.eq.0) then
        write(*,*)  'Error:  file ',
     >         outname(strbeg(outname):strend(outname)),
     >       ' is not a currently open output file.'
        return
      end if

      write(*,*)
      write(*,*) 'TRANS_OUT:  closing file ',
     >         outname(strbeg(outname):strend(outname))
      write(*,*)

      call ncclos(id_out(fptr),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to close file ',
     >         outname(strbeg(outname):strend(outname))
        return
      end if

      do k=fptr,numfls-1

        id_out(k) = id_out(k+1)
        do i=1,3
        do j=1,2
          x_dims(i,k,j) = x_dims(i,k+1,j)
          y_dims(i,k,j) = y_dims(i,k+1,j)
        enddo
        enddo
        do i=1,5
        do j=1,2
          z_dims(i,k,j) = z_dims(i,k+1,j)
        enddo
        enddo
        t_dims(1,k,1) = t_dims(1,k+1,1)
        t_dims(1,k,2) = t_dims(1,k+1,2)
        numout(k) = numout(k+1)
        flnam(k) = flnam(k+1)

      enddo

      write(*,*)

      numfls = numfls-1

      return
      end

*------------------------------------------------------------------------

#include "pointer.icl"

      subroutine open_out(outname,idcdf)

#include "default.icl"

      include '/usr/local/include/netcdf.inc'

      character*(*) outname
      integer idcdf

      integer strbeg,strend,ierr,id_dum,fptr,k
      integer nxf,nyf,nzf,ntf
      logical new,err
      character*(80) cdf_name,tmp

      common/flcount/numfls
      integer numfls

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      write(*,*)
      write(*,*) 'TRANS_OUT:  opening file ',
     >               outname(strbeg(outname):strend(outname))
      write(*,*)
*
* check to see if we have space for another file
*
      fptr = numfls+1

      if (fptr.gt.10) then
        write(*,*) 'Error:  too many files open!'
        return
      end if
*
* check to see if file already open
*
      do k=1,numfls
        if (outname(strbeg(outname):strend(outname)).eq.
     >       flnam(k)(strbeg(flnam(k)):strend(flnam(k)))) then
          write(*,*) 'Error:  file ',flnam(k)(strbeg(flnam(k)):
     >            strend(flnam(k))), ' already open.'
          return
        end if
      enddo
*
* open new netcdf file
*
      cdf_name = outname(strbeg(outname):strend(outname))//'.cdf'
      id_out(fptr) = ncopn(
     >      cdf_name(strbeg(cdf_name):strend(cdf_name)),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  cannot open netcdf file ',
     >        cdf_name(strbeg(cdf_name):strend(cdf_name))
        return
      end if
*
* get IDs for dimensions (used later when defining variables)
*
      x_dims(1,fptr,2) = ncdid(id_out(fptr),'nx',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nx.'
        return
      end if
      x_dims(2,fptr,2) = ncdid(id_out(fptr),'nxp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nxp1.'
        return
      end if
      x_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      y_dims(1,fptr,2) = ncdid(id_out(fptr),'ny',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension ny.'
        return
      end if
      y_dims(2,fptr,2) = ncdid(id_out(fptr),'nyp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nyp1.'
        return
      end if
      y_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      z_dims(1,fptr,2) = ncdid(id_out(fptr),'nz',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nz.'
        return
      end if
      z_dims(2,fptr,2) = ncdid(id_out(fptr),'nzp1',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension nzp1.'
        return
      end if
      z_dims(3,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if

      t_dims(1,fptr,2) = ncdid(id_out(fptr),'time',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension time.'
        return
      end if
      t_dims(2,fptr,2) = ncdid(id_out(fptr),'one',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find id for dimension one.'
        return
      end if
*
* check to see if dimensions match those of current data set
*
      call ncdinq(id_out(fptr),x_dims(1,fptr,2),tmp,nxf,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find dimension nx.'
        return
      end if
      if (nxf.eq.nx) then
        x_dims(1,fptr,1) = nx
        x_dims(2,fptr,1) = nx+1
        x_dims(3,fptr,1) = 1
      else
        write(*,*) 'Error:  dimensions of named file do not match ',
     >                       'those of data set.'
        return
      end if 

      call ncdinq(id_out(fptr),y_dims(1,fptr,2),tmp,nyf,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find dimension ny.'
        return
      end if
      if (nyf.eq.ny) then
        y_dims(1,fptr,1) = ny
        y_dims(2,fptr,1) = ny+1
        y_dims(3,fptr,1) = 1
      else
        write(*,*) 'Error:  dimensions of named file do not match ',
     >                        'those of data set.'
        return
      end if

      call ncdinq(id_out(fptr),z_dims(1,fptr,2),tmp,nzf,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find dimension nz.'
        return
      end if
      if (nzf.eq.nz) then
        z_dims(1,fptr,1) = nz
        z_dims(2,fptr,1) = nz+1
        z_dims(3,fptr,1) = 1
      else
        write(*,*) 'Error:  dimensions of named file do not match ',
     >                        'those of data set.'
        return
      end if

      call ncdinq(id_out(fptr),t_dims(1,fptr,2),tmp,ntf,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  failed to find dimension time.'
        return
      end if
      if ((ntime.eq.ntf).or.(ntime.eq.1)) then
        t_dims(1,fptr,1) = ntf
        t_dims(2,fptr,1) = 1
      else
        write(*,*) 'Error:  dimensions of named file do not match ',
     >                        'those of data set.'
        return
      end if
*
* add new file to list of viable output files
*
      numfls = numfls+1
      flnam(numfls) = outname(strbeg(outname):strend(outname))

      return
      end

*------------------------------------------------------------------------

      subroutine join(outfl,newfl,nfl)      

#include "default.icl"

      include '/usr/local/include/netcdf.inc'
 
      integer nfl 
      character*(*) outfl(nfl),newfl

      integer maxfl
      parameter (maxfl=10)
  
      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=50,maxdim=15,maxgatt=25,maxvatt=7)

      common/newfl_info/ndims,dsz,dnm,nvars,vnm,ngatts,gattnm,igatt,
     >          rgatt,cgatt,nvdims,vdims,nvatts,vattnm,ivatt,rvatt,
     >          cvatt
      integer ndims,dsz(maxdim),nvars,ngatts,igatt(maxgatt)
      integer nvdims(maxvar),vdims(4,maxvar),nvatts(maxvar)
      integer ivatt(maxvatt,maxvar)
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dnm(maxdim),vnm(maxvar),gattnm(maxgatt,2)
      character*(80) cgatt(maxgatt),vattnm(maxvatt,maxvar,2)
      character*(80) cvatt(maxvatt,maxvar)

      integer strbeg,strend,ierr,fptr(maxfl),k,n,q,p,ifnd
      integer idnew,itmp(5),numt,size,idtime
      real r_val
      character*(80) ctmp(2)
      logical err
      ive_ptr getmem,times,var
      
      common/flcount/numfls
      integer numfls
 
      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      if (nfl.gt.maxfl) then
        write(*,*) 'Error:  too many files.  Increase maxfl in ',
     >                         'routine join.'
        return
      end if

      numt = 0

      do k=1,nfl 
        fptr(k) = 0
        q = 1
 888    if ((fptr(k).eq.0).and.(q.le.numfls)) then
          if (outfl(k)(strbeg(outfl(k)):strend(outfl(k))).eq.
     >           flnam(q)(strbeg(flnam(q)):strend(flnam(q)))) then
            fptr(k) = q
          else
            q = q+1
          end if
          goto 888
        end if
        if (fptr(k).eq.0) then
          write(*,*) 'Error:  file ',
     >         outfl(k)(strbeg(outfl(k)):strend(outfl(k))),
     >             ' is not a currently open output file.'
          return
        end if
        numt = numt + t_dims(1,fptr(k),1)
      enddo

      times = getmem(numt)
      if (times.eq.0) then
        write(*,*) 'Error:  unable to allocate space for times.'
        return
      end if

      call set_times(%val(times),outfl,fptr,nfl,numt,err)
      if (err) then
        call freemem(times)
        return
      end if

      call get_stuff(fptr,nfl,err)
      if (err) then
        call freemem(times)
        return
      end if

      call set_newfl(newfl)
 
      newfl = newfl(strbeg(newfl):strend(newfl))//'.cdf'
      idnew = ncopn(newfl(strbeg(newfl):strend(newfl)),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem opening new data file.'
        call freemem(times)
        return
      end if

      size = x_dims(2,fptr(1),1)*y_dims(2,fptr(1),1)
     >                       *z_dims(2,fptr(1),1)
      var = getmem(size)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to allocate space for var.'
        call freemem(times)
        return
      end if 

      idtime = ncvid(id_out(fptr(1)),'time',ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem getting id for time variable.'
        call freemem(var)
        call freemem(times)
        return
      end if

      do k=1,nvars

        if (k.ne.idtime) then

          q = 1
          ifnd = 0
 87       if ((ifnd.eq.0).and.(q.le.nvdims(k))) then
            if (vdims(q,k).eq.t_dims(1,fptr(1),2)) then
              ifnd = 1
            else 
              q = q+1
            end if
            goto 87
          end if

          if (ifnd.eq.0) then
            call joinvar(%val(var),%val(times),size,numt,idnew,fptr,
     >                              nfl,k,0,err)
            if (err) then
              call freemem(var)
              call freemem(times)
              return
            end if
          end if

        end if

      enddo

      do n=1,numt
      
        call ncvpt1(idnew,idtime,n,r_val(%val(times),n),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem writing time.'
          call freemem(var)
          call freemem(times)
          return
        end if

        do k=1,nvars

          if (k.ne.idtime) then

            q = 1
            ifnd = 0
 897        if ((ifnd.eq.0).and.(q.le.nvdims(k))) then
              if (vdims(q,k).eq.t_dims(1,fptr(1),2)) then
                ifnd = 1
              else 
                q = q+1
              end if
              goto 897
            end if

            if (ifnd.eq.1) then
              call joinvar(%val(var),%val(times),size,numt,idnew,fptr,
     >                                     nfl,k,n,err)
              if (err) then
                call freemem(var)
                call freemem(times)
                return
              end if
            end if

          end if

        enddo

      enddo

      call ncclos(idnew,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem closing new file.'
        call freemem(var)
        call freemem(times)
        return
      end if 

      call freemem(var)
      call freemem(times)
        
      return
      end
      
*-----------------------------------------------------------------

      subroutine set_times(times,outfl,fptr,nfl,nt,err)

      include '/usr/local/include/netcdf.inc'

      integer nfl,nt,fptr(nfl) 
      real times(nt)
      character*(*) outfl(nfl)
      logical err

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      integer ierr,k,q,p,itmp,strbeg,strend,idtime
      real rtmp
      character*(80) ctmp

      err = .false.

      do k=1,nfl
        idtime = ncvid(id_out(fptr(k)),'time',ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem getting id for variable time.'
          err = .true.
          return
        end if
        call ncvgt1(id_out(fptr(k)),idtime,1,times(k),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem reading time.'
          err = .true.
          return
        end if
      enddo
         
      do k=1,nfl-1
      do q=nfl-1,k,-1
        if (times(q+1).lt.times(q)) then
          ctmp = outfl(q)
          outfl(q) = outfl(q+1)
          outfl(q+1) = ctmp
          rtmp = times(q)
          times(q) = times(q+1)
          times(q+1) = rtmp
          itmp = fptr(q)
          fptr(q) = fptr(q+1)
          fptr(q+1) = itmp
        end if
      enddo
      enddo

      itmp = 1
      do k=1,nfl
        idtime = ncvid(id_out(fptr(k)),'time',ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem getting id for variable time.'
          err = .true.
          return
        end if
        call ncvgt(id_out(fptr(k)),idtime,1,t_dims(1,fptr(k),1),
     >                            times(itmp),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem reading times.'
          err = .true.
          return
        end if
        if (k.gt.1) then
          if (times(itmp).le.times(itmp-1)) then
            write(*,*) 'Error:  times do not increase ',
     >                    'monotonically from ',
     >           outfl(k-1)(strbeg(outfl(k-1)):strend(outfl(k-1)))
            write(*,*) '        to ',
     >           outfl(k)(strbeg(outfl(k)):strend(outfl(k)))
            err = .true.
            return
          end if
        end if
        itmp = itmp + t_dims(1,fptr(k),1)
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine get_stuff(fptr,nfl,err)

      include '/usr/local/include/netcdf.inc'

      integer nfl,fptr(nfl)
      logical err

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=50,maxdim=15,maxgatt=25,maxvatt=7)

      common/newfl_info/ndims,dsz,dnm,nvars,vnm,ngatts,gattnm,igatt,
     >          rgatt,cgatt,nvdims,vdims,nvatts,vattnm,ivatt,rvatt,
     >          cvatt
      integer ndims,dsz(maxdim),nvars,ngatts,igatt(maxgatt)
      integer nvdims(maxvar),vdims(4,maxvar),nvatts(maxvar)
      integer ivatt(maxvatt,maxvar)
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dnm(maxdim),vnm(maxvar),gattnm(maxgatt,2)
      character*(80) cgatt(maxgatt),vattnm(maxvatt,maxvar,2)
      character*(80) cvatt(maxvatt,maxvar)
 
      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      integer strbeg,strend,k,p,q,idrec,vtyp,ierr,itmp(6),ifnd
      integer atyp,alen
      real*8 atmp
      character*(80) ctmp(2)

      err = .false.

      call ncinq(id_out(fptr(1)),ndims,nvars,ngatts,idrec,ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem with file inquiry.'
        err = .true.
        return
      end if

      do k=2,nfl
        call ncinq(id_out(fptr(k)),itmp(1),itmp(2),itmp(3),idrec,ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem with file inquiry.'
          err = .true.
          return
        else if (ndims.ne.itmp(1)) then
          write(*,*) 'Error:  files have different numbers of ',
     >                              'dimensions.'
          err = .true.
          return
        else if (nvars.ne.itmp(2)) then 
          write(*,*) 'Error:  files have different numbers of ',
     >                              'variables.'
          err = .true.
          return
        else if (ngatts.ne.itmp(3)) then
          write(*,*) 'Error:  files have different numbers of ',
     >                              'global attributes.'
          err = .true.
          return
        end if
      enddo

      if (nvars.gt.maxvar) then
        write(*,*) 'Error:  too many variables.  Increase ',
     >                       'maxvar in routine join.'
        err = .true.
        return
      else if (ndims.gt.maxdim) then
        write(*,*) 'Error:  too many dimensions.  Increase ',
     >                       'maxdim in routine join.'
        err = .true.
        return
      else if (ngatts.gt.maxgatt) then
        write(*,*) 'Error:  too many global attributes.  Increase ',
     >                       'maxgatt in routine join.'
        err = .true.
        return
      end if

      do p=1,ndims

        call ncdinq(id_out(fptr(1)),p,dnm(p),dsz(p),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem with dimension inquiry.'
          err = .true.
          return
        end if

        do k=2,nfl
          call ncdinq(id_out(fptr(k)),p,ctmp(1),itmp(1),ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem with dimension inquiry.'
            err = .true.
            return
          else if (ctmp(1)(strbeg(ctmp(1)):strend(ctmp(1))).ne.
     >          dnm(p)(strbeg(dnm(p)):strend(dnm(p)))) then
            write(*,*) 'Error:  files have different dimensions.'
            err = .true.
            return
          else if ((itmp(1).ne.dsz(p))
     >               .and.(p.ne.t_dims(1,fptr(1),2))) then
            write(*,*) 'Error:  files have different dimensions.'
            err = .true.
            return
          end if
       enddo

      enddo 

      do p=1,nvars

        call ncvinq(id_out(fptr(1)),p,vnm(p),vtyp,nvdims(p),
     >                       vdims(1,p),nvatts(p),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem with variable inquiry.'
          err = .true.
          return
        else if (nvatts(p).gt.maxvatt) then
          write(*,*) 'Error:  too many variable attributes. ',
     >                 'Increase maxvatts in routine join.'
          err = .true.
          return
        end if 

        do k=2,nfl

          ifnd = 0
          q = 1
 1234     if ((ifnd.eq.0).and.(q.le.nvars)) then
            call ncvinq(id_out(fptr(k)),q,ctmp(1),vtyp,itmp(5),
     >                              itmp(1),itmp(6),ierr)
            if (ierr.ne.0) then
              write(*,*) 'Error:  problem with variable inquiry.'
              err = .true.
              return
            else if (ctmp(1)(strbeg(ctmp(1)):strend(ctmp(1))).eq.
     >             vnm(p)(strbeg(vnm(p)):strend(vnm(p)))) then
              ifnd = 1
            else
              q = q+1
            end if
            goto 1234
          end if

          if ((ifnd.eq.0).or.(itmp(5).ne.nvdims(p))) then
            write(*,*) 'Error:  files have different variables.'
            err = .true.
            return
          else
            q = 1
 2233       if (q.le.nvdims(p)) then
              if (itmp(q).ne.vdims(q,p)) then
                write(*,*) 'Error:  files have different variables.'
                err = .true.
                return
              else
                q = q+1
              end if
              goto 2233
            end if
          end if

        enddo   

        do k=1,nvatts(p)

          call ncanam(id_out(fptr(1)),p,k,vattnm(k,p,1),ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem getting attribute name.'
            err = .true. 
            return
          end if

          call ncainq(id_out(fptr(1)),p,vattnm(k,p,1),atyp,alen,ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem with attribute inquiry.'
            err = .true.
            return
          end if

          if (atyp.eq.NCCHAR) then
            vattnm(k,p,2) = 'c'
            call ncagtc(id_out(fptr(1)),p,vattnm(k,p,1),cvatt(k,p),80,
     >                                      ierr)
            if (ierr.ne.0) then
              write(*,*) 'Error:  problem reading string attribute.'
              err = .true.
              return
            end if
          else if (atyp.eq.NCLONG) then
            vattnm(k,p,2) = 'i'
            call ncagt(id_out(fptr(1)),p,vattnm(k,p,1),ivatt(k,p),ierr)
            if (ierr.ne.0) then
              write(*,*) 'Error:  problem reading integer attribute.'
              err = .true.
              return
            end if
          else if (atyp.eq.NCFLOAT) then
            vattnm(k,p,2) = 'r'
            call ncagt(id_out(fptr(1)),p,vattnm(k,p,1),rvatt(k,p),ierr)
            if (ierr.ne.0) then
              write(*,*) 'Error:  problem reading real attribute.'
              err = .true.
              return
            end if
          else if (atyp.eq.NCDOUBLE) then
            vattnm(k,p,2) = 'r'
            call ncagt(id_out(fptr(1)),p,vattnm(k,p,1),atmp,ierr)
            if (ierr.ne.0) then
              write(*,*) 'Error:  problem reading real attribute.'
              err = .true.
              return
            end if
            rvatt(k,p) = real(atmp)
          else
            write(*,*) 'Error:  attribute type not recognized.' 
            err = .true.
            return
          end if

        enddo

      enddo

      do p=1,ngatts

        call ncanam(id_out(fptr(1)),NCGLOBAL,p,gattnm(p,1),ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem getting global attribute name.'
          err = .true.
          return
        end if 

        call ncainq(id_out(fptr(1)),NCGLOBAL,gattnm(p,1),atyp,
     >                                 alen,ierr)
        if (ierr.ne.0) then
          write(*,*) 'Error:  problem inquiring about global ',
     >                                'attribute.'
          err = .true.
          return
        end if

        if (atyp.eq.NCCHAR) then
          gattnm(p,2) = 'c'
          call ncagtc(id_out(fptr(1)),NCGLOBAL,gattnm(p,1),
     >                     cgatt(p),80,ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem reading string global ',
     >                              'attribute.'
            err = .true.
            return
          end if
        else if (atyp.eq.NCLONG) then
          gattnm(p,2) = 'i'
          call ncagt(id_out(fptr(1)),NCGLOBAL,gattnm(p,1),
     >                    igatt(p),ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem reading integer global ',
     >                              'attribute.'
            err = .true.
            return
          end if
        else if (atyp.eq.NCFLOAT) then
          gattnm(p,2) = 'r'
          call ncagt(id_out(fptr(1)),NCGLOBAL,gattnm(p,1),
     >                    rgatt(p),ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem reading real global ',
     >                              'attribute.'
            err = .true.
            return
          end if
        else if (atyp.eq.NCDOUBLE) then
          gattnm(p,2) = 'r'
          call ncagt(id_out(fptr(1)),NCGLOBAL,gattnm(p,1),
     >                    atmp,ierr)
          if (ierr.ne.0) then
            write(*,*) 'Error:  problem reading real global ',
     >                              'attribute.'
            err = .true.
            return
          end if
          rgatt(p) = real(atmp)
        else
          write(*,*) 'Error:  global attribute type not recognized.'
          err = .true.
          return
        end if

      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine set_newfl(newfl)

      include '/usr/local/include/netcdf.inc'

      character*(*) newfl

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=50,maxdim=15,maxgatt=25,maxvatt=7)

      common/newfl_info/ndims,dsz,dnm,nvars,vnm,ngatts,gattnm,igatt,
     >          rgatt,cgatt,nvdims,vdims,nvatts,vattnm,ivatt,rvatt,
     >          cvatt
      integer ndims,dsz(maxdim),nvars,ngatts,igatt(maxgatt)
      integer nvdims(maxvar),vdims(4,maxvar),nvatts(maxvar)
      integer ivatt(maxvatt,maxvar)
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dnm(maxdim),vnm(maxvar),gattnm(maxgatt,2)
      character*(80) cgatt(maxgatt),vattnm(maxvatt,maxvar,2)
      character*(80) cvatt(maxvatt,maxvar)

      integer strbeg,strend,k,q,p,slen,iunit
      character*(80) outstr,cdl_name,cdf_name,command,command2

      cdl_name(1:len(cdl_name)) = ' '
      cdl_name = newfl(strbeg(newfl):strend(newfl))//'.cdl'

      write(*,*)
      write(*,*) 'TRANS_OUT:  creating new file'
      write(*,*) 'TRANS_OUT:  cdl file is ',
     >                 cdl_name(strbeg(cdl_name):strend(cdl_name))

      iunit = 777
      open(unit=iunit,file=cdl_name(strbeg(cdl_name):
     >                      strend(cdl_name)),status='new')

      write(iunit,*) 'netcdf TransOut{'
      write(iunit,*) 

      write(iunit,*) 'dimensions:'
      write(iunit,*)
      do k=1,ndims
        if (dnm(k)(strbeg(dnm(k)):strend(dnm(k))).eq.'time') then
          write(iunit,*) 'time=UNLIMITED;'
        else
          write(iunit,*) dnm(k)(strbeg(dnm(k)):strend(dnm(k))),'=',
     >                         dsz(k),';'
        end if
      enddo
      write(iunit,*)

      write(iunit,*) 'variables:'
      write(iunit,*)
      do k=1,nvars

        outstr(1:len(outstr)) = ' '
        outstr = 'float '//vnm(k)(strbeg(vnm(k)):strend(vnm(k)))
        slen = 6 + strend(vnm(k)) - strbeg(vnm(k)) + 1
        outstr = outstr(1:slen)//'('
        slen = slen + 1
        do q=nvdims(k),1,-1
          outstr = outstr(1:slen)//
     >       dnm(vdims(q,k))(strbeg(dnm(vdims(q,k))):
     >                             strend(dnm(vdims(q,k))))
          slen = slen + strend(dnm(vdims(q,k)))
     >                        - strbeg(dnm(vdims(q,k))) + 1
          outstr = outstr(1:slen)//','
          slen = slen + 1
        enddo
        outstr = outstr(1:slen-1)//');'

        write(iunit,*) outstr(1:slen+1)

        do p=1,nvatts(k)

          outstr(1:len(outstr)) = ' '
          outstr = vnm(k)(strbeg(vnm(k)):strend(vnm(k)))//':'
          outstr = outstr(strbeg(outstr):strend(outstr))//
     >            vattnm(p,k,1)(strbeg(vattnm(p,k,1))
     >                               :strend(vattnm(p,k,1)))
          outstr = outstr(strbeg(outstr):strend(outstr))//'='
          if (vattnm(p,k,2)(1:1).eq.'c') then
            write(iunit,*) outstr(strbeg(outstr):strend(outstr)),'"',
     >           cvatt(p,k)(strbeg(cvatt(p,k)):strend(cvatt(p,k))),
     >                             '";'
          else if (vattnm(p,k,2)(1:1).eq.'i') then 
            write(iunit,*) outstr(strbeg(outstr):strend(outstr)),
     >                      ivatt(p,k),';'
          else if (vattnm(p,k,2)(1:1).eq.'r') then
            write(iunit,*) outstr(strbeg(outstr):strend(outstr)),
     >                      rvatt(p,k),';'
          end if

        enddo 

        write(iunit,*)

      enddo

      write(iunit,*) '//global attributes:'
      write(iunit,*)

      do p=1,ngatts

        outstr(1:len(outstr)) = ' '
        outstr = ':'//gattnm(p,1)(strbeg(gattnm(p,1)):
     >                                   strend(gattnm(p,1)))
        outstr = outstr(strbeg(outstr):strend(outstr))//'='
        if (gattnm(p,2)(1:1).eq.'c') then
          write(iunit,*) outstr(strbeg(outstr):strend(outstr)),'"',
     >         cgatt(p)(strbeg(cgatt(p)):strend(cgatt(p))),'";'
        else if (gattnm(p,2)(1:1).eq.'i') then
          write(iunit,*) outstr(strbeg(outstr):strend(outstr)),
     >                    igatt(p),';'
        else if (gattnm(p,2)(1:1).eq.'r') then
          write(iunit,*) outstr(strbeg(outstr):strend(outstr)),
     >                    rgatt(p),';'
        end if

      enddo

      write(iunit,*) 
      write(iunit,*) '}'

      close(iunit)

      cdf_name = newfl(strbeg(newfl):strend(newfl))//'.cdf'
      command(1:len(command)) = ' '
      write(*,*) 'TRANS_OUT:  netcdf file is ',
     >         cdf_name(strbeg(cdf_name):strend(cdf_name))
      write(*,*)
      write(command,*) 'ncgen -o',cdf_name(strbeg(cdf_name):
     >                 strend(cdf_name)),' ',
     >            cdl_name(strbeg(cdl_name):strend(cdl_name))
      call system(command)
      write(command2,*) 'rm -f ',cdl_name(strbeg(cdl_name):
     >                      strend(cdl_name))
      call system(command2)     

      return 
      end 

*----------------------------------------------------------------------

      subroutine joinvar(var,times,size,numt,idnew,fptr,nfl,k,n,err)

      include '/usr/local/include/netcdf.inc'

      integer size,numt,idnew,nfl,fptr(nfl),k,n
      real var(size),times(numt)
      logical err

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=50,maxdim=15,maxgatt=25,maxvatt=7)

      common/newfl_info/ndims,dsz,dnm,nvars,vnm,ngatts,gattnm,igatt,
     >          rgatt,cgatt,nvdims,vdims,nvatts,vattnm,ivatt,rvatt,
     >          cvatt
      integer ndims,dsz(maxdim),nvars,ngatts,igatt(maxgatt)
      integer nvdims(maxvar),vdims(4,maxvar),nvatts(maxvar)
      integer ivatt(maxvatt,maxvar)
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dnm(maxdim),vnm(maxvar),gattnm(maxgatt,2)
      character*(80) cgatt(maxgatt),vattnm(maxvatt,maxvar,2)
      character*(80) cvatt(maxvatt,maxvar)

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >              t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      integer idvar,ist(4),iln(4),ierr,q,p,ifnd,itmp

      err = .false.

      itmp = 1
      q = 1
      ifnd = 0
 26   if ((ifnd.eq.0).and.(q.le.nfl)) then
        if (n.le.itmp+t_dims(1,fptr(q),1)-1) then
          ifnd = 1
        else
          itmp = itmp + t_dims(1,fptr(q),1)
          q = q+1
        end if
        goto 26
      end if

      if (ifnd.eq.0) then
        write(*,*) 'Error:  time not found in joinvar.'
        err = .true.
        return
      end if

      idvar = ncvid(id_out(fptr(q)),vnm(k),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  variable not found in joinvar.'
        err = .true.
        return
      end if

      do p=1,nvdims(k)
        if (vdims(p,k).eq.t_dims(1,fptr(q),2)) then
          ist(p) = n-itmp+1
          iln(p) = 1
        else
          ist(p) = 1
          iln(p) = dsz(vdims(p,k))
        end if
      enddo

      call ncvgt(id_out(fptr(q)),idvar,ist,iln,var(1),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem reading variable.'
        err = .true.
        return
      end if

      do p=1,nvdims(k)
        if (vdims(p,k).eq.t_dims(1,fptr(q),2)) then
          ist(p) = n
          iln(p) = 1
        else
          ist(p) = 1
          iln(p) = dsz(vdims(p,k))
        end if
      enddo

      call ncvpt(idnew,k,ist,iln,var(1),ierr)
      if (ierr.ne.0) then
        write(*,*) 'Error:  problem writing variable.'
        err = .true.
        return
      end if

      return
      end

      

