c  ....  c_2d (i,j)  ... Concentration values in 2-D space
c  ....  u_2d (i,j)  ... X component of velocities in 2-D space
c  ....  v_2d (i,j)  ... Y component of velocities in 2-D space 
c
c  ....  cold (i)    ... Old Concentration values in 1-D space 
c  ....  cnew (i)    ... New Concentration values in 1-D space 
c  ....  c_mid (i)   ... Intermediate Concentration values in 1-D space 
c  ....  uu (i)      ... Velocity array in 1_D space  
c
c  ....  Nx     ........ Number of Global Columns 
c  ....  Ny     ........ Number of Global Rows 
c  ....  Nro    ........ Number of horizontal lines (Rows)
c  ....  Nco    ........ Number of Vertical lines (Columns)
c  ....  Nbg    ........ Number of Boundary Groups
c
c  ....  Ibr (i) ....... First column number of i_th row i=1,Nro           
c  ....  Ier (i) ....... Last column number of i_th row i=1,Nro
c  ....  Jr (i)  ....... Global Row number of i_th row i=1,Nro
c  ....  Kbr (i) ....... Upstream Boundary Group number of i_th row i=1,Nro
c  ....  Ker (i) ....... Downstream Boundary Group number of i_th row i=1,Nro
c  ....  Jbc (i) ....... First row number of i_th column i=1,Nco           
c  ....  Jec (i) ....... Last column number of i_th column i=1,Nco
c  ....  Ic (i)  ....... Global Row number of i_th column i=1,Nco
c  ....  Kbc (i) ....... Upstream Boundary Group number of i_th column i=1,Nco
c  ....  Kec (i) ....... Downstream Boundary Group number of i_th column i=1,Nco
c
c  ....  Rbg (i) ....... Boundary value of i_th Group i=1,Nbg
c  ....  Kbg (i) ....... Boundary code of i_th Group i=1,Nbg

      parameter ( mx = 500 )
      real, allocatable :: c_2d (:,:) , u_2d (:,:) , v_2D(:,:) , 
     $      Rbg(:) , cross (:,:), dif_2d (:,:) , zd (:,:), sos ( :,:) , 
     $      sos_c (:,:) , Pbg ( :) 
      real cold ( mx ) , cnew ( mx ) , c_mid ( mx ) , uu ( mx )
      integer, allocatable :: Ibr(:) , Ier (:) , Jr (:) , 
     $                        kbr (:) , Ker (:) 
      integer, allocatable :: Jbc(:) , Jec(:) , Ic(:) ,
     $    Kbc(:) , Kec (:) , Kbg(:) 
      real dif ( mx ) , ss ( mx ) , b0 ( mx ) , zd_1 ( mx )
      common / bl1 / nx , ny
      common / bl2 / x00 , y00 , sig , sig2
      common / xy/ Dx , Dy
      common /dd8/ k_init , v_init , k_bed , v_bed , 
     $     k_Tr_per , k_tr , v_tr , j_init , J_bed , j_tr      
      common / qq/ kq , nq , xq , yq , zq , zt 
      real zt ( 100 , 5 )
      real xq ( 100, 5 ) , yq ( 100,5 ) , zq ( 100,5 )
      integer nq ( 5 )
      common /ncr/  Nco , Nro
      character*40 fn



      ss = 0.
      open ( 4 , file = 'GRD_2d_project.txt' )
92    format ( a40 )      
      do k = 1 , 6 
          kk = k + 4
          read ( 4 , 92 ) fn

          open ( kk , file = fn )
      end do    
      close ( 4 )
      
      read ( 8 , * ) kq
      do k = 1 , kq
          read ( 8 , * ) nq ( k )
          do i = 1 , nq ( k )
              read ( 8 , * ) xq ( i,k ) , yq ( i,k ) , zq ( i,k )
          end do
      end do    
      close ( 8 )  
      call plane ( )  
      read ( 5 , * ) Dx , Dy
      read ( 5 , * ) nx , ny
      read ( 5 , * ) Nro
      allocate ( Ibr( Nro) , Ier (Nro) , Jr (Nro) , 
     $                        kbr (Nro) , Ker (Nro ) )
      allocate ( c_2d ( nx , ny ) , dif_2d(nx,ny) , cross ( nx,ny) ,
     $             U_2d ( nx,ny ) , V_2d ( nx , ny ) , zd (nx,ny) , 
     $             sos ( nx , ny ) , sos_c ( nx , ny ) )       
      do i = 1 , Nro
          read ( 5 , * ) Ibr ( i ) , Ier ( i ) , Jr ( i ) , 
     $                   Kbr ( i ) , Ker ( i )
      end do

      read ( 5 , * ) Nco
      
      allocate (  Jbc(Nco) , Jec(Nco) , Ic(Nco) ,
     $    Kbc(Nco) , Kec (Nco)  )

      do i = 1 , Nco
          read ( 5 , * ) Jbc(i) , Jec(i) , Ic(i) ,
     $    Kbc(i) , Kec (i)  

      end do 

      read ( 6 , * )
      read ( 6 , * ) Dt , Tend , npt
      read ( 6 , * )
      read ( 6 , * ) Dt_2 , Tend_2 , npt_2
      read ( 6 , * )
      
      read ( 6 , * ) Nbg
    
      allocate ( Kbg ( Nbg ) , Rbg ( Nbg) , Pbg ( Nbg )  )
      
      read ( 6 , * )
      do i = 1 , Nbg
          read ( 6 , * ) Kbg ( i ) , Rbg ( i ) , Pbg ( i ) 
         
      end do
      read ( 6 , * )
      read ( 6 , * ) k_init , v_init , j_init
      read ( 6 , * )
      read ( 6 , * ) k_bed , v_bed, j_bed
      read ( 6 , * )  
      read ( 6 , * ) k_Tr_per
      read ( 6 , * )
      read ( 6 , * ) k_tr , v_tr , j_tr      
c..............................................  Initial condition

      
C$$$$$$       do j = 1 , Nro
C$$$$$$         
C$$$$$$           i1 = Ibr ( j )
C$$$$$$           i2 = ier ( j )
C$$$$$$           j0 = Jr ( j )
C$$$$$$           y0 = ( j0 - .5 ) * Dy
C$$$$$$           do i = i1 , i2
C$$$$$$               x0 = ( i - .5 ) * Dx
C$$$$$$               rt = ( x0 - x00 ) ** 2 + ( y0 - y00 ) ** 2
C$$$$$$               c_2d ( i , j0 ) = exp ( - rt / sig2 )
C$$$$$$           end do
C$$$$$$       end do

      C_2D  = 2.
      U_2d = 0.
      V_2d = 0.       
      nt = int ( tend / Dt ) + 1
      cross = 0.
      sos = 0.
      write ( 7 , 19 )
19    format ( 'VARIABLES=X  Y  U  V  C D' )
      
      sc = Dx * Dy / Dt
      
c------------------------------------------------------------------------------------
c     Ground Water Solver  Start
c
17    format ( 'Time step no.',i5,'     Time =',f10.3 )
18    format ( 'step no.',i5,'          Time =',f10.3 )
      do k = 1 , nt
c        nt
        
          tim0 = ( k - 1 ) * Dt
          if ( mod ( k-1 , npt ) .eq. 0 ) then
              call plot ( tim0 , dif_2D ,  
     $             c_2d , Ibr , Ier , Jr , Jbc , Jec ,Ic , U_2d , 
     $             V_2d , zd , dif_2d )
              if ( k .eq. 1 ) call solve ( sos , sos_c )
          end if       
                    
          write(*,17) k , tim0
       
c          pause
          difmax = 0.
          do j = 1 , Nro
            
              i1 = Ibr ( j )
              i2 = Ier ( j )
              j0 = Jr ( j )
              
              kb = Kbr ( j )
              ke = Ker ( j )

              do i = i1 , i2
                  cold ( i ) = c_2d ( i , j0 )

                  dif ( i ) = dif_2d ( i , j0 )
                  b0 ( i ) = cross ( i , j0 ) + sos ( i , j0 ) / sc 
                 

              end do
              
              kbb = Kbg ( kb )
              kee = Kbg ( ke )
              rbb = Rbg ( kb )
              ree = Rbg ( ke )

c              call Fromm_o ( i1 , i2 , cold , cnew , uu , Dx , 
c     $          Dt , kbb , rbb , kee , ree )
              call diff_imp ( i1 , i2 , dt , dx , Dif , cold , cnew, ss, 
     $                        kbb , rbb , kee , ree , b0 ) 

              do i = i1 , i2
                  c_2D ( i ,j0 ) = cnew ( i )
                  
                  cross ( i , j0 ) = b0 ( i )
                  rt = abs ( cold ( i ) - cnew ( i ) )
                  if ( rt .gt. difmax ) difmax = rt                  
              end do
          end do

          
         
          do i = 1 , Nco
            
              j1 = Jbc ( i )
              j2 = Jec ( i )
              i0 = Ic ( i )
              
              kb = Kbc ( i )
              ke = Kec ( i )

              do j = j1 , j2
                  cold ( j ) = c_2d ( i0 , j )
                  uu ( j ) = v0
                  dif ( j ) = dif_2d ( i0 , j )
                  b0 ( j ) = cross ( i0 , j )
              end do
     
              kbb = Kbg ( kb )
              kee = Kbg ( ke )
              rbb =  Rbg ( kb )
              ree = Rbg ( ke )

c              call Fromm_o ( j1 , j2 , cold , cnew , uu , Dx , 
c     $          Dt , kbb , rbb , kee , ree )
              call diff_imp ( j1 , j2 , dt , dx , Dif , cold , cnew, ss, 
     $                        kbb , rbb , kee , ree , b0  )     
      
              do j = j1 , j2
                  c_2D ( i0 ,j ) = cnew ( j )
                  cross ( i0 , j ) = b0 ( j ) + sos ( i0 , j ) / sc
                  rt = abs ( cold ( j ) - cnew ( j ) )
                  if ( rt .gt. difmax ) difmax = rt
              end do
              
          end do
          
       end do
       umax = -10000.
       do i = 1 , nx
           do j = 1 , ny
               h = c_2d ( i , j ) - zd ( i , j )
               U_2d ( i , j ) = U_2d ( i , j ) / h
               V_2d ( i , j ) = V_2d ( i , j ) / h
               zd ( i , j ) = h
               u0 = abs ( u_2d ( i , j ) )
               v0 = abs ( v_2d ( i , j ) )
               if ( u0 .gt. umax ) then
                   umax = u0
                   hmax = h
                   imax = i
                   jmax = j
               end if    
                   
               if ( v0 .gt. umax ) then
                  
                   umax = v0
                   hmax = h
                   imax = i
                   jmax = j
               end if    
                   
           end do
       end do 
       close (7)
       
c------------------------------------------------------------------------------------
c     Ground Water Solver  End
c------------------------------------------------------------------------------------




      write ( 10 , 29 )
29    format ( 'VARIABLES=X  Y  U  V  C' )
      Dt = Dt_2
      Tend = Tend_2
      npt = npt_2
      
      Dif0 = 0     
      nt = int ( tend / Dt ) + 1
        
      sc = Dx * Dy / Dt         
c------------------------------------------------------------------------------------
c     Transpor Solver  Start
c

c      call plot_2 ( tim0 , dif_2D ,  
c     $             c_2d , Ibr , Ier , Jr , Jbc , Jec ,Ic , U_2d , 
c     $             V_2d  )      

      c_2d = 0.
      dd = min ( Dx , Dy )
      ddt = .8 * dd / umax
      write(*,*) '.....End of Hydrodynamic Part.........'
      write(*,*) '......................................'
      write(*,*)
      write(*,*) '..... Maximum velocity =',umax
      write(*,*) '..... Make sure that Dt has been set to less'
      write(*,8) ddt
8     format (   '.. ... than ',f4.1,' to fulfill CFL condition')      
      write(*,*)
      write(*,*) '..... Press RETURN to start Transport part'
c      write(*,*)  umax , hmax , imax , jmax , 
c     $           c_2d ( imax , jmax ) , zd ( imax , jmax )
      pause
         
      do k = 1 , nt
        
          tim0 = ( k - 1 ) * Dt
          write(*,18) k , tim0 
          if ( mod ( k-1 , npt ) .eq. 0 ) call plot_2 ( tim0 , dif_2D ,  
     $             c_2d , Ibr , Ier , Jr , Jbc , Jec ,Ic , U_2d , 
     $             V_2d  )
          su = 0.
          do i = 1 , nx
              do j = 1 , ny
                  su = su + c_2d ( i , j ) * zd ( i , j ) * Dx * Dy 
              end do
          end do        


          do j = 1 , Nro
            
              i1 = Ibr ( j )
              i2 = Ier ( j )
              j0 = Jr ( j )
              
              kb = Kbr ( j )
              ke = Ker ( j )

              do i = i1 , i2
                  cold ( i ) = c_2d ( i , j0 )
                  uu ( i ) = u_2d ( i , j0 )
c                  write(*,*) i , uu(i)
                  dif ( i ) = Dif0
                  q0 = sos ( i , j0 )
                  
                  if ( q0 .gt. 0. ) then
                      b0 ( i ) = .5 * sos_c ( i , j0 ) * q0 * Dt / Dy
    
                    else
                      b0 ( i ) = .5 * cold ( i ) * q0 * Dt / Dy  
                  end if
                      
                  zd_1 ( i ) = zd ( i , j0 )
                 
              end do
              
              kbb = Kbg ( kb )
              kee = Kbg ( ke )
              rbb =  Pbg ( kb )
              ree = Pbg ( ke )

              call Fromm_o ( i1 , i2 , cold , cnew , uu , Dx , 
     $          Dt , kbb , rbb , kee , ree , b0 , zd_1 )
c              call diff_imp ( i1 , i2 , dt , dx , Dif , cold , cnew, ss, 
c     $                        kbb , rbb , kee , ree , b0 ) 

              do i = i1 , i2
                  c_2D ( i ,j0 ) = cnew ( i )
                  
                 
              end do
          end do

          
         
          do i = 1 , Nco
            
              j1 = Jbc ( i )
              j2 = Jec ( i )
              i0 = Ic ( i )
              
              kb = Kbc ( i )
              ke = Kec ( i )

              do j = j1 , j2
                  cold ( j ) = c_2d ( i0 , j )
                  uu ( j ) = v_2d ( i0 , j )
c                  write(*,*) j , uu(j)
                  dif ( j ) = dif0
                  q0 = sos ( i0 , j )
                  if ( q0 .gt. 0. ) then
                    
                      b0 ( j ) = .5 * q0 * sos_c ( i0 , j ) * Dt / Dx
                    
                    else
                      b0 ( j ) = .5 * q0 * cold ( j ) * Dt / Dx
                  end if    
                  zd_1 ( j ) = zd ( i0 , j )
              end do
     
              kbb = Kbg ( kb )
              kee = Kbg ( ke )
              rbb =  Pbg ( kb )
              ree = pbg ( ke )

              call Fromm_o ( j1 , j2 , cold , cnew , uu , Dx , 
     $          Dt , kbb , rbb , kee , ree , b0 , zd_1 )
c              call diff_imp ( j1 , j2 , dt , dx , Dif , cold , cnew, ss, 
c     $                        kbb , rbb , kee , ree , b0  )     
      
              do j = j1 , j2
                  c_2D ( i0 ,j ) = cnew ( j )
                  
              end do
          end do
       end do

c------------------------------------------------------------------------------------
c     Transport Solver  End
c------------------------------------------------------------------------------------
       

2      do i = 5 , 10
           close ( i )
       end do
       end   


c-----------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine Fromm_o (  i1 , i2 , c_old , c_new , uu , Dx , 
     $          Dt , kbb , rbb , kee , ree , b0 , zd_1 )
     
      parameter ( mx = 500 )

      real c_old(mx) , c_new(mx) , mass ( mx )
      real slope ( mx ) , uu ( mx ) , b0 ( mx ) , zd_1 ( mx )
   
      do j = i1 , i2
          mass ( j ) = c_old ( j ) * Dx * zd_1 ( j ) + b0 ( j ) 
         
      end do
   
      do j = i1+1 , i2-1
          slope ( j ) = .5*( c_old ( j+1 ) - c_old ( j-1 ) ) / Dx
      end do
      slope ( i2 ) = 0.  
      slope ( i1 ) = 0.   
c      udt = u * Dt
     
      do j = i1+1 , i2
          udt = uu ( j ) * Dt 
          er = udt / Dx
          if ( abs ( er ) .gt. .9 ) then
              write (*,*) '  ....!!!!! WARNING ... CFD number=',er,'>1' 
              write(*,*) j , uu ( j ) , Dt
              pause
          end if    
          jm = j-1
          jp = j

          z0 = .5 * ( zd_1 ( jm ) + zd_1 ( jp ) )

c          rr = ( c_old ( jp ) - c_old ( jm ) ) / Dx
c          r2 = - Dt * Dif * rr
      
          if ( udt .gt. 0. ) then
              d0 = .5 * ( Dx - udt ) 
              c0 = c_old ( jm ) + slope ( jm )* d0
            else
              d0 = .5 * ( Dx + udt )
              c0 = c_old ( jp ) - slope ( jp ) * d0
              
              
          end if      
          flux = c0 * udt * z0
          mass ( jm ) = mass ( jm ) - flux
          mass ( jp ) = mass ( jp ) + flux
      end do

      if ( kbb .eq. 2 ) then
          
          mass ( i1 ) = mass ( i1 ) + rbb
      end if

      if ( kee .eq. 2 ) then
          mass ( i2 ) = mass ( i2 ) - ree
      end if      

      do j = i1 , i2
          c_new ( j ) = mass ( j ) / Dx / zd_1 ( j )
      end do

      if ( kbb .eq. 1 ) c_new ( i1 ) = rbb
      if ( kee .eq. 1 ) c_new ( i2 ) = ree  
      return
      end  
c-----------------------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine plot_2 ( tim0 , dif , cc , Ibr , Ier , Jr , Jbc , Jec ,
     $                 Ic , U_2d , V_2d  ) 
   
      integer  kele ( 30000 ,4)
      real xp ( 30000 ) , yp ( 30000 ) , vv ( 30000 )
      integer Ike ( 30000 ) , Jke ( 30000 ) 
      integer*1 nm ( 300000 )
      integer Ibr( Nro) , Ier (Nro) , Jr (Nro) 
      integer  Jbc(Nco) , Jec(Nco) , Ic(Nco) 
      real cc ( nx , ny ) , Dif ( nx , ny ) , U_u ( 30000 ) , 
     $         U_v ( 30000 ) ,zzd ( 30000 )
      common / bl2 / x00 , y00 , sig , sig2
      common / bl1 / nx , ny
      common / xy/ Dx , Dy

      common /ncr/  Nco , Nro

      real u1d ( 300 ) , U_2d ( nx , ny ) , V_2d ( nx , ny )
      common / a1 / kele , xp , yp , Ike , Jke , nm
      common / a2 /  nnode , nelem 
c            write(*,*) c_2d(nx-4 , ny-4 )



      
      do i = 1 , nnode
          vv ( i ) = 0
         
      end do    
      U_u = 0.
      U_v = 0.
      zzd = 0.
                    
      write ( 10 , 12 ) tim0 , nnode , nelem

      do i = 1 , nelem 
          i0 = Ike ( i )
          j0 = Jke ( i )
          do j = 1 , 4
              k0 = kele ( i , j )
              
              vv ( k0 ) = vv ( k0 ) + cc ( i0 , j0 ) 
c              if ( i0 .eq. 40 .and. j0 .eq. 80 ) then
cc                  write(*,*) cc(i0,j0) , k0 , vv(k0)
c                  pause
c              end if     
              U_u ( k0 ) = U_u ( k0 ) + U_2d ( i0 , j0 )
              U_v ( k0 ) = U_v ( k0 ) + V_2d ( i0 , j0 ) 
           
          end do    
      end do 

      do i = 1 , nnode
          rt = 1. / real ( nm ( i ) )
          c0 = vv ( i ) * rt
          u0 = U_u ( i ) * rt
          v0 = U_v ( i ) * rt 
   
       
          write ( 10 , 11 ) xp ( i ) , yp ( i ) , u0 , v0 , c0
      end do   

      do i = 1 , nelem
          write ( 10 , 13 ) (kele ( i , j ) , j = 1 , 4 )
      end do
11    format ( 8e15.7 )
12    format ('ZONE T="',f10.1,'" N=',i6,' E=',i6,
     $            ' F=FEPOINT ET=QUADRILATERAL')         
13    format ( 4i6 )
      return
      end
c--------------------------------------------------------------------------------------
      subroutine plot ( tim0 , dif , cc , Ibr , Ier , Jr , Jbc , Jec ,
     $                 Ic , U_2d , V_2d , zd , dif_2d ) 
c       call plot ( tim0 , dif_2D ,  
c     $             cc , Ibr , Ier , Jr , Jbc , Jec ,Ic , U_2d , 
c     $             V_2d  )     
      integer  kele ( 30000 ,4)
      real xp ( 30000 ) , yp ( 30000 ) , vv ( 30000 )
      integer Ike ( 30000 ) , Jke ( 30000 ) 
      integer*1 nm ( 300000 )
      integer Ibr( Nro) , Ier (Nro) , Jr (Nro) 
      integer  Jbc(Nco) , Jec(Nco) , Ic(Nco) 
      real cc ( nx , ny ) , Dif ( nx , ny ) , U_u ( 30000 ) , 
     $         U_v ( 30000 ) , zd ( nx , ny ) ,zzd ( 30000 )
      common / bl2 / x00 , y00 , sig , sig2
      common / bl1 / nx , ny
      common / xy/ Dx , Dy
      common / a1 / kele , xp , yp , Ike , Jke , nm
      common / a2 /  nnode , nelem      
      common /dd8/ k_init , v_init , k_bed , v_bed , 
     $     k_Tr_per , k_tr , v_tr , j_init , J_bed , j_tr
      common /ncr/  Nco , Nro
      common / min0 / x_min , y_min 

      real u1d ( 300 ) , U_2d ( nx , ny ) , V_2d ( nx , ny )
      real z_res ( 5 ) , dif_2d ( nx , ny ) 
      save kele , xp , yp , Ike , Jke , nm
      save nnode , nelem
c            write(*,*) c_2d(nx-4 , ny-4 )
      read ( 5 , * , end = 10 ) nnode , nelem
      nm = 0
      do i = 1 , nnode
          read ( 5 , * ) xp ( i ) , yp ( i )
      end do
      
      x_min = 1.E15
      y_min = x_min
      do i = 1 , nelem
          read ( 5 , * ) ( kele ( i , j ) , j = 1 , 4 )
     $                   , Ike ( i ) , Jke ( i )  

          xx = 0.
          yy = 0.

          ii = Ike ( i )
          jj = Jke ( i )
          
          do j = 1 , 4
              k0 = kele ( i , j )
              nm ( k0 ) = nm ( k0 ) + 1  
              xx = xx + xp ( k0 )
              yy = yy + yp ( k0 )
          end do    
          x0 = .25 * xx
          y0 = .25 * yy
          if ( x0 .lt. x_min ) x_min = x0
          if ( y0 .lt. y_min ) y_min = y0  
          call interpolate ( x0 , y0 , z_res )
          rt = ( x0 - x00 ) ** 2 + ( y0 - y00 ) ** 2
c          cc ( Ike ( i ) , Jke ( i )  ) = exp ( - rt / sig2 )
          if ( k_init .lt. 0 ) then
              zd (  Ike ( i ) , Jke ( i ) ) = v_init 
            else
              cc ( ii , jj ) = z_res ( j_init )
          end if  
                          
          if ( k_bed .lt. 0 ) then  
              zd (  ii , jj ) = v_bed
            else
              zd (  ii , jj ) = z_res ( j_bed )
          end if
          if ( k_tr .lt. 0. ) then
              q0 = v_tr
            else
              q0 = z_res ( j_tr )
          end if      
            
          if ( k_Tr_per .eq. 1 ) then
              dif_2d (  ii , jj ) = q0
            else
              h0 = cc ( ii  , jj ) - 
     $             zd (  ii , jj )
              dif_2d (  ii , jj ) = q0 * h0
          end if     
               

                
      end do

10    continue
      U_2d = 0.
      V_2d = 0. 
      do i = 1 , nnode
          vv ( i ) = 0
         
      end do    
      U_u = 0.
      U_v = 0.
      zzd = 0.
      do j = 1 , Nro
          
          i1 = Ibr ( j )
          i2 = Ier ( j )
          j0 = Jr ( j )
          do i = i1 , i2
              u1d ( i ) = 0.
          end do    

          do i = i1 + 1 , i2
              ip = i
              im = i-1 
              u0 = - dif ( i,j0 ) * ( cc ( ip , j0 ) - cc ( im , j0 ) ) 
     $         / dx / 2.
              u1d ( im ) = u1d ( im ) + u0
              u1d ( ip ) = u1d ( ip ) + u0
          end do     

          do i = i1 , i2
              U_2d ( i , j0 ) = u1d ( i )
          end do
      end do        
              
      do i = 1 , Nco
         
          j1 = Jbc ( i )
          j2 = Jec ( i )
          i0 = Ic ( i )          
          do j = j1+1 , j2
              u1d ( j ) = 0.
          end do    
      
          do j = j1+1 , j2
              jp = j
              jm = j-1 
              u0 = - dif ( i0,j ) * ( cc ( i0 , jp ) - cc ( i0 , jm ) ) 
     $         / dy / 2.
              u1d ( jm ) = u1d ( jm ) + u0
              u1d ( jp ) = u1d ( jp ) + u0
          end do    
          do j = j1 , j2
              V_2d ( i0 , j ) = u1d ( j )
          end do                       
      end do     
                    
      write ( 7 , 12 ) tim0 , nnode , nelem

      do i = 1 , nelem 
          i0 = Ike ( i )
          j0 = Jke ( i )
          do j = 1 , 4
              k0 = kele ( i , j )
              
              vv ( k0 ) = vv ( k0 ) + cc ( i0 , j0 )  
              U_u ( k0 ) = U_u ( k0 ) + U_2d ( i0 , j0 )
              U_v ( k0 ) = U_v ( k0 ) + V_2d ( i0 , j0 ) 
              zzd ( k0 ) = zzd ( k0 ) + zd ( i0 , j0 )
          end do    
      end do 

      do i = 1 , nnode
          rt = 1. / real ( nm ( i ) )
          c0 = vv ( i ) * rt
          u0 = U_u ( i ) * rt
          v0 = U_v ( i ) * rt 
          e0 = zzd ( i ) * rt
          write ( 7 , 11 ) xp ( i ) , yp ( i ) , u0 , v0 , c0, e0
      end do   

      do i = 1 , nelem
          write ( 7 , 13 ) (kele ( i , j ) , j = 1 , 4 )
      end do
11    format ( 8e15.7 )
12    format ('ZONE T="',f10.1,'" N=',i6,' E=',i6,
     $            ' F=FEPOINT ET=QUADRILATERAL')         
13    format ( 4i6 )
      return
      end




c--------------------------------------------------
      SUBROUTINE diff_imp ( i1 , i2 , dt , dx , Dif , co , cn, ss, 
     $   kb , rb , ke , re  , b0 )
              
      parameter ( mx = 500 )
      real co ( mx ) , cn ( mx ) , dif ( mx ) , ss ( mx ) , b0 ( mx )
      real a1 ( mx ) , a2 ( mx ) , a3 ( mx ) , a0 ( mx )
      real*8 x ( mx )
      
      q1 = .5 * Dt / Dx / dx
      
      do j = i1+1 , i2-1
        
          q0 = q1 * Dif ( j )
          a1 ( j) = -q0
          a2 ( j ) =1. + 2 * q0
          a3 ( j ) = -q0
          a0 ( j ) = q0 * ( co ( j-1 ) - 2. * co ( j ) + co ( j+1 ) ) +
     $            co ( j ) + ss ( j ) * Dt + b0 ( j )
  
          
      end do
      j = i1
      if ( kb .eq. 1 ) THEN
          a1 ( j ) = 0.
          a2 ( j ) = 1.
          a3 ( j ) = 0.
          a0 ( j ) = rb
        else
          q0 = q1 * Dif ( j )
          a1 ( j ) = 0.
          a3 ( j ) = - q0
          a2 ( j ) =1. + q0
          a0 ( j ) = q0 * ( - 1. * co ( j ) + co ( j+1 ) ) +
     $            co ( j ) + ss ( j ) * Dt + b0 ( j ) + rb * Dt / Dx
     
  
      end if
    
      j = i2   
      if ( ke .eq. 1 ) then
          
          a1 ( j )  = 0.
          a2 ( j ) = 1.
          a3 ( j ) = 0.
          a0 ( j ) = re
        else
          q0 = q1 * Dif ( j )
          a1 ( j ) = - q0
          a2 ( j ) = 1. + q0
          a3 ( j ) = 0.
          a0 ( j ) = q0 * ( co ( j-1 ) - 1. * co ( j )  ) +
     $            co ( j ) + ss ( j ) * Dt + b0 ( j ) + re* Dt / Dx
      end if

      call tomas (  i1 , i2 , a1 , a2 , a3 , a0 , x )

      do i = i1 , i2
          cn ( i ) = x ( i )
      end do

      do j = i1+1 , i2-1
          q0 = q1 * Dif ( j )

          b0(j) = 2 *q0 * ( cn ( j-1 ) - 2. * cn ( j ) + cn ( j+1 ) ) 
     $          
      end do    
      b0 ( i1 ) = - 2 * q0 * ( cn ( i1 ) - cn ( i1+1 ) )  
      b0 ( i2 ) = - 2.* q0 * ( cn ( i2 ) - cn ( i2-1 ) )
      return
      end    





c-----------------------------------------------
       subroutine tomas( j1 , j2 , a1 , a2 , a3 , a0 , x )
      parameter ( mx = 500 ) 
      real a1 ( mx ) , a2 ( mx ), a3 ( mx ) , a0 ( mx )
      Real*8 rt , E ( 0:mx ), F( 0:mx ), x ( mx )

      E(j1-1) = 0.
      F(j1-1) = 0.

      do j = j1 , j2

         rt = a1( j ) * E( j-1 ) + a2 ( j )
         E ( j )= -a3 ( j ) / rt
         F ( j )= ( a0( j )-a1 ( j ) * F( j-1 ) ) / rt
      end do

      x ( j2 ) = F ( j2 )
      do j = j2-1 , j1, -1
          x( j ) = E( j )* x ( j+1 ) + F( j )
      end do
      return
      end
c------------------------------------------------------------------
      subroutine interpolate ( x0 , y0 , z_res )
      common / qq/ kq , nq , xq , yq , zq , zt
      common / abc / aa , bb , cc
      real xq ( 100,5 ) , yq ( 100,5 ) , zq ( 100,5 )
      real zt ( 100,5 ) , z_res ( 5 )
      integer nq ( 5 ) 
      real aa ( 5 ) , bb ( 5 ) , cc ( 5 )
     
      do k = 1 , kq
          sum1 = 0.
          sum2 = 0.
          a = aa ( k )
          b = bb ( k )
          c = cc ( k )
          z2 = a * x0 + b * y0 + c
          do i = 1 , nq ( k )
              x1 = xq ( i , k )
              y1 = yq ( i , k )
              z1 = zt ( i , k )

              d0 = ( ( x0 - x1 ) ** 2 + ( y0 - y1 ) ** 2 ) + 1.E-6

              sum1 = sum1 + 1. / d0
              sum2 = sum2 + z1 / d0 

          end do
          z_res ( k ) = sum2 / sum1 + z2

      end do 
      return
      end    
c---------------------------------------------------------------------------------
      subroutine plane ( )        
      common / abc / aa , bb , cc
      common / qq/ kq , nq , xq , yq , zq , zt
      real xq ( 100,5 ) , yq ( 100,5 ) , zq ( 100,5 )
      real zt ( 100,5 )      
      integer nq ( 5 ) 
      real aa ( 5 ) , bb ( 5 ) , cc ( 5 )     
      do k = 1 , kq 
          xx = 0.
          xy = 0.
          xz = 0.
          yy = 0.
          yz = 0.
          zz = 0.
          nq0 = nq ( k ) 

          do i = 1 , nq0
        
              x = xq ( i,k )
              y = yq ( i,k )
              z = zq ( i,k )


              xx = xx + x * x
              xy = xy + x * y
              xz = xz + x * z
              yy = yy + y * y
              yz = yz + y * z
              zz = zz + z
        
          end do

          e1 = xy ** 2
          e2 = xx * yy
          e3 = xy * yz
          e4 = xz * yy
          e5 = xy * xz
          e6 = yz * xx
      
          rt = e1 - e2
   
          a = ( e3 - e4 ) / rt
          b = ( e5 - e6 ) / rt
          c = zz / real ( nq0 )

          do i = 1 , nq ( k )
        
              z = a * xq ( i ,k) + b * yq ( i,k ) + c
              zt ( i,k ) = zq ( i ,k) - z
              
          end do   
          aa ( k ) = a
          bb ( k ) = b
          cc ( k ) = c 
      end do    
      return
      end
c---------------------------------------------------------------
      subroutine solve ( sos , sos_c )      
      
      common / bl1 / nx , ny
      common / xy/ Dx , Dy  
      common / min0 / x_min , y_min        

      real sos ( nx , ny ) , sos_c ( nx , ny )
c      open ( 4 , file = 'source.dat' )

      read ( 9 , * ) n0
      do i = 1 , n0
          read ( 9 , * ) x1 , x2 , y1 , y2 , q0 , c0
          
          i1 = int ( (x1-x_min) / Dx ) + 1
          i2 = int( (x2-x_min) / Dx ) + 2
          j1 = int ( (y1-y_min) / Dy ) + 1
          j2 = int ( (y2-Y_min) / Dy ) + 2

          nn = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
          q00 = q0 / real ( nn )
c          c00 = c0 / real ( nn )
          do l1 = i1 , i2
              do l2 = j1 , j2
                  sos ( l1 , l2 ) = sos ( l1 , l2 ) + q00
                  sos_c ( l1 , l2 ) = sos_c ( l1 , l2 ) + c0
  
              end do
          end do  
      end do  
      close ( 9 )    
          return
          end    

      
      
      
      