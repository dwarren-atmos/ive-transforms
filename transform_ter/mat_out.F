
#include "pointer.icl"

      subroutine make_mat(filename,fld,nfld,z,n,ispc,jspc)

#include "default.icl"

      integer n,nfld,ispc,jspc
      real z
      character*(*) fld(nfld),filename

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(80) dunits,disunits,dim_names(4)
   
      integer*8 matOpen,matClose,fp,stat
      integer strbeg,strend,iu,iv,iw,k,nt,nxout,nyout,i0,j0
      character*(80) outfl
      logical requ,new
      ive_ptr var,var2d,getmem

      nt = ntime

      nxout = nx/ispc
      nyout = ny/jspc
      i0 = mod(nx,ispc)
      j0 = mod(ny,jspc)
      i0 = i0/2 + 1
      j0 = j0/2 + 1

      var2d = getmem(4*nxout*nyout)
      if (var2d.eq.0) then
        write(6,*) 'Error:  could not allocate space for var2d.'
        return
      end if

      outfl(1:len(outfl)) = ' '
      outfl = '/home/disk/misty/cepi/matlab/'
      outfl = outfl(strbeg(outfl):strend(outfl))//
     >              filename(strbeg(filename):strend(filename))
      outfl = outfl(strbeg(outfl):strend(outfl))//'.mat'
      fp = matOpen(outfl(strbeg(outfl):strend(outfl)),"w")

      do k=1,nfld

        var = getvar(fld(k),ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (var.eq.0) then
          write(6,*) 'Error:  could not access field ',
     >                fld(k)(strbeg(fld(k)):strend(fld(k)))
          return
        end if

        iu = 0
        iv = 0
        iw = 0
        if (requ(min(1),plmin(1))) iu = 1
        if (requ(min(2),plmin(2)).and.(twod.ne.1)) iv = 1
        if (requ(min(3),plmin(3))) iw = 1

        call v_z2d(%val(var2d),%val(var),%val(zbot_p),%val(wgz),
     >                  %val(sgz),ztop,z,n,trans_z,nx,ny,nz,nt,
     >                  nxout,nyout,i0,j0,ispc,jspc,iu,iv,iw)

        call write_mat(%val(var2d),fld(k),fp,nxout,nyout)

      enddo

      stat = matClose(fp)

      call freemem(var2d)

      return
      end

*------------------------------------------------------------------------

      subroutine v_z2d(var2d,var,zs,wgz,sgz,ztop,z,n,itrans,
     >                         nx,ny,nz,nt,nxout,nyout,i0,j0,
     >                         ispc,jspc,iu,iv,iw)

      integer n,itrans,nx,ny,nz,nt,nxout,nyout,i0,j0
      integer ispc,jspc,iu,iv,iw
      double precision var2d(nxout,nyout)
      real var(nx+iu,ny+iv,nz+iw,nt),zs(nx,ny)
      real wgz(nz+1),sgz(nz),ztop,z

      integer i,j,iout,jout,kplus,k
      real z1,z2,var1,var2,tfac,znode,z0

      do jout=1,nyout
      do iout=1,nxout

        j = j0 + (jout-1)*jspc
        i = i0 + (iout-1)*ispc

        if (itrans.ne.0) then 
          z0 = zs(i,j)
        else
          z0 = 0.
        end if
        tfac = (ztop-z0)/ztop

        kplus = 0
        k = 1

 565    if ((kplus.eq.0).and.(k.le.(nz+iw))) then
          if (iw.eq.1) then
            znode = z0 + tfac*wgz(k)
          else
            znode = z0 + tfac*sgz(k)
          end if
          if (znode.gt.z) then
            kplus = k 
          else
            k = k+1
          end if
          goto 565
        end if

        if (kplus.eq.0) kplus = nz+iw
        if (kplus.eq.1) kplus = 2

        if (iw.eq.1) then
          z1 = z0 + tfac*wgz(kplus-1)
          z2 = z0 + tfac*wgz(kplus)
        else
          z1 = z0 + tfac*sgz(kplus-1)
          z2 = z0 + tfac*sgz(kplus)
        end if
        var1 = 0.5*(var(i,j,kplus-1,n)+var(i+iu,j+iv,kplus-1,n))
        var2 = 0.5*(var(i,j,kplus,n)+var(i+iu,j+iv,kplus,n))
       
        var2d(iout,jout) = (z-z2)/(z1-z2)*var1 + (z-z1)/(z2-z1)*var2

      enddo
      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine write_mat(var2d,fld,fp,nxout,nyout)

      integer*8 fp
      integer nxout,nyout
      double precision var2d(nxout,nyout)
      character*(*) fld

      integer*8 mxCreateFull,matPutMatrix,mp,stat
      integer*8 mxGetPr
      integer strbeg,strend

      mp = mxCreateFull(nxout,nyout,0)
      call mxSetName(mp,fld(strbeg(fld):strend(fld)))
      call mxCopyReal8ToPtr(var2d,mxGetPr(mp),nxout*nyout)
      stat = matPutMatrix(fp,mp)

      return
      end 

