
# include "pointer.icl"

      ive_ptr function uflx(whch,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) whch,dunits,disunits,dim_names(4)

      integer ndims2,dims2(4)
      real stag2(4),min2(4),max2(4),miss2
      character*(80) dunits2,disunits2,dim_names2(4)

      integer twod,iu,iv,iw,k,f
      real dom_min(4),dom_max(4),r_val
      logical new,err
      ive_ptr getvar,getmem,var,fcor

      uflx = 0

      fcor = getvar('f_cor',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
      if (fcor.eq.0) then
        write(*,*) 'Warning:  could not find coriolis parameter.'
        f = 0
      else
        if (r_val(%val(fcor),1).eq.0.) then
          f = 0
        else
          f = 1
        end if
      end if

      call getrarr('plmin',dom_min(4),4,err) 
      if (err) then
        write(*,*) 'Error:  cannot access physical domain minimums.'
        return
      end if
      call getrarr('plmax',dom_max(4),4,err)
      if (err) then
        write(*,*) 'Error:  cannot access physical domain maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else
        twod = 0
      end if

      if (whch(1:1).eq.'X') then
        var = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                           dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access u.'
          return
        end if
      else if (whch(1:1).eq.'Y') then
        var = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access v.'
          return
        end if
      else if (whch(1:1).eq.'Z') then
        var = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access w.'
          return
        end if
      else
        write(*,*) 'Error:  do not know how to compute momentum flux ',
     >                        'for direction ',whch(1:1)
        return
      end if

      miss = 8675309

      uflx = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (uflx.eq.0) then
        write(*,*) 'Error:  unable to allocate memory for uflx.'
        return
      end if

      do k=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(uflx),k,r_val(%val(var),k))
      enddo

      if (whch(1:1).eq.'Y') then
        call driver('free = v',0)
      else if (whch(1:1).eq.'Z') then
        call driver('free = w',0)
      end if

      if (f.eq.0) then

        if (whch(1:1).eq.'X') then
          call var_2_upts(%val(uflx),%val(var),miss,dims(1)-1,dims(2),
     >                     dims(3),dims(4),1,0,0)
        else
          var = getvar('u',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >                   disunits2,dim_names2,new)
          if (var.eq.0) then
            write(*,*) 'Error:  unable to access u.'
            call freemem(uflx)
            return
          end if
          if (whch(1:1).eq.'Y') call var_2_vpts(%val(uflx),%val(var),
     >        miss,dims2(1)-1,dims2(2),dims2(3),dims2(4),1,0,0,twod)
          if (whch(1:1).eq.'Z') call var_2_wpts(%val(uflx),%val(var),
     >        miss,dims2(1)-1,dims2(2),dims2(3),dims2(4),1,0,0)
        end if

      else

        var = getvar('UABS',ndims2,dims2,stag2,min2,max2,miss2,
     >             dunits2,disunits2,dim_names2,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access uabs.'
          call freemem(uflx)
          return
        end if
        call driver('free = eta',0)
        if (whch(1:1).eq.'X') call var_2_upts(%val(uflx),%val(var),
     >        miss,dims2(1),dims2(2),dims2(3),dims2(4),0,0,0)   
        if (whch(1:1).eq.'Y') call var_2_vpts(%val(uflx),%val(var),
     >        miss,dims2(1),dims2(2),dims2(3),dims2(4),0,0,0,twod)   
        if (whch(1:1).eq.'Z') call var_2_wpts(%val(uflx),%val(var),
     >        miss,dims2(1),dims2(2),dims2(3),dims2(4),0,0,0)   
        call driver('free = uabs',0)

      end if

      var = getvar('RHO',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >                      disunits2,dim_names2,new)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access density.'
        call freemem(uflx)
        return
      end if
      if (whch(1:1).eq.'X') call var_2_upts(%val(uflx),%val(var),miss,
     >     dims2(1),dims2(2),dims2(3),dims2(4),0,0,0)   
      if (whch(1:1).eq.'Y') call var_2_vpts(%val(uflx),%val(var),miss,
     >     dims2(1),dims2(2),dims2(3),dims2(4),0,0,0,twod)   
      if (whch(1:1).eq.'Z') call var_2_wpts(%val(uflx),%val(var),miss,
     >     dims2(1),dims2(2),dims2(3),dims2(4),0,0,0)   

      dunits = 'kg/m/s/s'
      disunits = 'kg/m/s/s'

      return
      end 

*-----------------------------------------------------------------------

      subroutine var_2_wpts(var_w,var,miss,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_w(nx,ny,nz+1,nt),var(nx+iu,ny+iv,nz+iw,nt),miss

      integer i,j,k,n

      if (iw.eq.1) then
  
        do n=1,nt
        do k=1,nz+1
        do j=1,ny
        do i=1,nx
            var_w(i,j,k,n) = var_w(i,j,k,n)*var(i,j,k,n)
        enddo
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do j=1,ny
        do i=1,nx
          var_w(i,j,1,n) = var_w(i,j,1,n)
     >              *0.25*(3.*(var(i,j,1,n)+var(i+iu,j+iv,1,n))
     >                      -var(i,j,2,n)-var(i+iu,j+iv,2,n))
          var_w(i,j,nz+1,n) = var_w(i,j,nz+1,n)
     >              *0.25*(3.*(var(i,j,nz,n)+var(i+iu,j+iv,nz,n))
     >                 -var(i,j,nz-1,n)-var(i+iu,j+iv,nz-1,n))
        enddo
        enddo
        enddo

        do n=1,nt
        do k=2,nz
        do j=1,ny
        do i=1,nx
          var_w(i,j,k,n) = var_w(i,j,k,n)*0.25*(var(i,j,k,n)
     >           +var(i+iu,j+iv,k,n)+var(i,j,k-1,n)
     >                   +var(i+iu,j+iv,k-1,n))
        enddo
        enddo
        enddo
        enddo

      end if

      return
      end

*-------------------------------------------------------------------------

      subroutine var_2_upts(var_u,var,miss,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_u(nx+1,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt),miss

      integer i,j,k,n

      if (iu.eq.1) then

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx+1
          var_u(i,j,k,n) = var_u(i,j,k,n)*var(i,j,k,n)
        enddo
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do k=1,nz
        do j=1,ny
          var_u(1,j,k,n) = var_u(1,j,k,n)
     >              *0.25*(3.*(var(1,j,k,n)+var(1,j+iv,k+iw,n))
     >                      -var(2,j,k,n)-var(2,j+iv,k+iw,n))
          var_u(nx+1,j,k,n) = var_u(nx+1,j,k,n)
     >              *0.25*(3.*(var(nx,j,k,n)+var(nx,j+iv,k+iw,n))
     >                 -var(nx-1,j,k,n)-var(nx-1,j+iv,k+iw,n))
        enddo
        enddo
        enddo

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=2,nx
          var_u(i,j,k,n) = var_u(i,j,k,n)*0.25*(var(i,j,k,n)
     >           +var(i,j+iv,k+iw,n)+var(i-1,j,k,n)
     >                     +var(i-1,j+iv,k+iw,n))
        enddo
        enddo
        enddo
        enddo

      end if

      return
      end

*------------------------------------------------------------------------

      subroutine var_2_vpts(var_v,var,miss,nx,ny,nz,nt,iu,iv,iw,twod)

      integer nx,ny,nz,nt,iu,iv,iw,twod
      real var_v(nx,ny+1-twod,nz,nt),var(nx+iu,ny+iv,nz+iw,nt),miss

      integer i,j,k,n

      if (twod.eq.1) then

        do n=1,nt
        do k=1,nz
        do i=1,nx
          var_v(i,1,k,n) = var_v(i,1,k,n)
     >               *0.5*(var(i,1,k,n)+var(i+iu,1,k+iw,n))
        enddo
        enddo
        enddo
     
      else if (iv.eq.1) then

        do n=1,nt
        do k=1,nz
        do j=1,ny+1
        do i=1,nx
          var_v(i,j,k,n) = var_v(i,j,k,n)*var(i,j,k,n)
        enddo
        enddo
        enddo
        enddo

      else 

        do n=1,nt
        do k=1,nz
        do i=1,nx
          var_v(i,1,k,n) = var_v(i,1,k,n)
     >              *0.25*(3.*(var(i,1,k,n)+var(i+iu,1,k+iw,n))
     >                      -var(i,2,k,n)-var(i+iu,2,k+iw,n))
          var_v(i,ny+1,k,n) = var_v(i,ny+1,k,n)
     >              *0.25*(3.*(var(i,ny,k,n)+var(i+iu,ny,k+iw,n))
     >                 -var(i,ny-1,k,n)-var(i+iu,ny-1,k+iw,n))
        enddo
        enddo
        enddo

        do n=1,nt
        do k=1,nz
        do j=2,ny
        do i=1,nx
          var_v(i,j,k,n) = var_v(i,j,k,n)*0.25*(var(i,j,k,n)
     >           +var(i+iu,j,k+iw,n)+var(i,j-1,k,n)
     >                   +var(i+iu,j-1,k+iw,n))
        enddo
        enddo
        enddo
        enddo

      end if

      return
      end 
  
