      PROGRAM test
      implicit none 
      integer :: i,k
      integer, parameter :: N_sp=9, N_react=42
      real*8 :: xnH, T_K, y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, 
     &     y_He, y_Hep, y_Hepp, y(N_sp), xk(N_react), C(N_sp), D(N_sp),
     &     yHe, tiny, t_ff
      real*8 :: time, dt, t_chem, t_chem_1

      tiny = 1.0d-20
      yHe  = 8.33333d-2
      xnH = 1.0d8
      T_K = 1.0d3

!     initial abundance
      y_H2   = 1.0d-6
      y_Hm   = 1.0d-12
      y_H2p  = 1.0d-12
      y_Hp   = 1.0d-4
      y_H    = 1.0d0 - 2.0d0*y_H2 - 2.0d0*y_H2p - y_Hm - y_Hp
      y_He   = yHe-2.0d0*tiny
      y_Hep  = tiny
      y_Hepp = tiny             ! yHe-y_He(i,j)-y_Hep(i,j)
      y_e    = y_Hp + y_H2p - y_Hm + y_Hep + 2.0d0*y_Hepp

!     convert to y(N_sp)
      y(1) = y_H
      y(2) = y_H2
      y(3) = y_e
      y(4) = y_Hp
      y(5) = y_H2p
      y(6) = y_Hm
      y(7) = y_He
      y(8) = y_Hep
      y(9) = y_Hepp
  
!========== time evolution ==========!     
      time = 0.0d0
      dt = 1.0d-3
      do i = 1,1000000000
         time = time + dt
         t_chem = 1.0d100
         !if (time>1.0d19) exit
         if (time>8.d10) exit
         call react_coef(xnH, T_K, y(1), y(2), xk)
         !call solve_chem_simple(xk,xnH,T_K,y,dt,t_chem)
         !print*, time
         call solve_chem_implicit(xnH, T_K, y, dt, t_chem, i) 
!     set time step
c         dt = 0.01d0 * t_chem
         t_ff = 1.0d15/dsqrt(xnH) * 1.0d2
         if(dt<t_ff) then 
            dt = dt*1.01d0
         else
            dt = t_ff
         endif

         if(mod(i,1)==0) then
c            print*, i,time, dt, y(1), y(4), y(2),
c     &           y(1)+2.0d0*y(2)+y(4)+2.0d0*y(5)+y(6),
c     &           y(3)+y(6), y(4)+y(5)+y(8)+2.0d0*y(9)     
            write(10,*), time, y(1), y(2), y(3), y(4), y(5), y(6),
     &           y(7), y(8), y(9)
         endif

c         print*, y(1)+2.0d0*y(2)+y(4)+2.0d0*y(5)+y(6),
c     &        y(3)+y(6), y(4)+y(5)+y(8)+2.0d0*y(9)

c         write(10,*), xnH, T_K
c     &        C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), 
c     &        D(1), D(2), D(3), D(4), D(5), D(6), D(7), D(8), D(9)
c         write(10,*) T_K, 
c     &        xk(1), xk(2), xk(3), xk(4), xk(5), xk(6), 
c     &        xk(7), xk(8), xk(9), xk(10), xk(11), xk(12), xk(13),
c     &        xk(14), xk(15), xk(16), xk(17), xk(18), xk(19), xk(20)
c     &        xk(31), xk(32), xk(33), xk(34), xk(35), xk(36),
c     &        xk(37), xk(38), xk(39), xk(40), xk(41), xk(42)
      enddo


      END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE solve_chem_implicit(xnH, T_K, y, dt, t_chem, i)
!     Non Equilibrium (Low Density)
      implicit none
      integer,parameter :: N_sp=9, N_react=42
      real*8 :: eps=1.d-8, eps_y=1.d-10
      integer :: isp, jsp, itr, i
      real*8 :: xnH, T_K, y(N_sp), y_init(N_sp), dy(N_sp), xk(N_react), 
     &     dt, t_chem, C(N_sp), D(N_sp), r_f(N_sp), y_tmp(N_sp), 
     &     r_f_fw(N_sp), r_f_bw(N_sp), delta_y, ddy(N_sp), err_0,
     &     dr_fdy(N_sp,N_sp), A(N_sp,N_sp), r_f_big, dr_f, err_max, tch,
     &     yHe, y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp

      do isp=1,N_sp
         y_init(isp)=y(isp)
         dy(isp)=0.d0
      enddo

      call react_coef(xnH,T_K,y(1),y(2),xk)
      do 20 itr=1,20
         call react_rat(xk,xnH,T_K,y,C,D,r_f)
         
         do jsp=1,N_sp
            if(dabs(y(jsp)).le.1.d-100) then
               delta_y=eps_y
            else
               delta_y=eps*y(jsp)
            endif
            do isp=1,N_sp
               if(isp.eq.jsp) then
                  y_tmp(isp)=y(isp)+5.d-1*delta_y
               else
                  y_tmp(isp)=y(isp)
               endif
            enddo
            call react_rat(xk,xnH,T_K,y_tmp,C,D,r_f_fw)
            y_tmp(jsp)=y(jsp)-0.5d0*delta_y
            call react_rat(xk,xnH,T_K,y_tmp,C,D,r_f_bw)
            
            do isp=1,N_sp
               dr_fdy(isp,jsp)=(r_f_fw(isp)-r_f_bw(isp))
     &              /delta_y
               r_f_big=max(dabs(r_f_fw(isp)),dabs(r_f_bw(isp)))
               if(r_f_big.ne.0.d0) then
                  dr_f=dabs(r_f_fw(isp)-r_f_bw(isp))
                  if(dr_f/r_f_big.lt.1.d-15) then
                     dr_fdy(isp,jsp)=0.d0
                  endif
               endif
            enddo
         enddo
         
         ! print*, "**********************"
         ! do isp=1,5
         !    do jsp=1,5
         !       print*, dr_fdy(isp,jsp)
         !    enddo
         ! enddo
         ! print*, "**********************"

C ------ set the matrix A
         do isp=1,N_sp
            do jsp=1,N_sp
               if(isp.ne.jsp) then
                  A(isp,jsp)=-dt*dr_fdy(isp,jsp)
               else
                  A(isp,jsp)=1.d0-dt*dr_fdy(isp,jsp)
               endif
            enddo
         enddo
         
C ------- set the vector ddy; WLI:initially dy(isp)=0
         do isp=1,N_sp
            ddy(isp)=r_f(isp)*dt-dy(isp)
         enddo

C ------- solve linear equations 
         call gaussj(A,N_sp,N_sp,ddy,1,1,1)
         
         do isp=1,N_sp
            dy(isp)=dy(isp)+ddy(isp)
            y(isp)=y(isp)+ddy(isp)
            if(y(isp).lt.0.d0) y(isp)=y_init(isp)            
         enddo
         !WLI: err_max is largest relative change dy among all species 
         err_max=0.d0
         do isp=1,N_sp
            if(y(isp).ne.0.d0) then
               err_0=dabs(ddy(isp)/y(isp))
            else
               err_0=0.d0
            endif
            err_max=max(err_0,err_max)
         enddo
         !WLI: if err_max small, consider reset the reaction timescale: t_chem
         if(err_max.lt.1.d-8) then
            t_chem=1.d20
            do isp=1,N_sp
               if(y(isp)-y_init(isp).ne.0.d0) then
                  tch=dabs((y(isp)+y_init(isp))
     &                 /(2.d0*(y(isp)-y_init(isp))))*dt                  
                  t_chem=dmin1(tch,t_chem)
               endif               
            enddo
            go to 100
         endif
         
 20   enddo
 100  continue

      y_H    = y(1)
      y_H2   = y(2)
      y_e    = y(3)
      y_Hp   = y(4)
      y_H2p  = y(5)
      y_Hm   = y(6)
      y_He   = y(7)
      y_Hep  = y(8)
      y_Hepp = y(9)
      
      ! 为了保证质子数守恒
      if(y_H<y_Hp) y_Hp = 1.0d0-y_H -2.0d0*y_H2-2.0d0*y_H2p-y_Hm
      if(y_H>y_Hp) y_H  = 1.0d0-y_Hp-2.0d0*y_H2-2.0d0*y_H2p-y_Hm

      yHe = 8.3333333d-2
      if(y_He>y_Hep) then
         if(y_He>y_Hepp) then
            y_He = yHe - y_Hep - y_Hepp
         endif
      endif
      if(y_Hep>y_He) then
         if(y_Hep>y_Hepp) then
            y_Hep = yHe - y_He - y_Hepp
         endif
      endif
      if(y_Hepp>y_He) then
         if(y_Hepp>y_Hep) then
            y_Hepp = yHe - y_He - y_Hep
         endif
      endif

      y_e = y_Hp + y_Hep + 2.0d0*y_Hepp + y_H2p - y_Hm

      y(1) = y_H
      y(2) = y_H2
      y(3) = y_e
      y(4) = y_Hp
      y(5) = y_H2p
      y(6) = y_Hm
      y(7) = y_He
      y(8) = y_Hep
      y(9) = y_Hepp

      return
      END

      ! in solving reaction: call gaussj(a=A,n=N_sp,np=N_sp,b=ddy,m=1,mp=1,ind=1)
      SUBROUTINE gaussj(a,n,np,b,m,mp,ind)
      INTEGER m,mp,n,np,NMAX
      PARAMETER (NMAX=200)
      DOUBLE PRECISION a(np,np),b(np,mp),b_max(np,mp)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv
      if(n == 1) then
         a(1,1)=1.d0/a(1,1)
         b(1,1)=a(1,1)*b(1,1)
         return
      endif
      if(ind == 1) then
         eps=1.d-10
      else
         eps=0.d0
      endif

      do 11 j=1,n
         ipiv(j)=0
 11   continue
      
      do ll=1,np
         do l=1,mp
            b_max(ll,l)=0.d0
         enddo
      enddo
      
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
           if(ipiv(j).ne.1)then
              do 12 k=1,n
                 if (ipiv(k).eq.0) then
                    if (abs(a(j,k)).ge.big)then
                       big=abs(a(j,k))
                       irow=j
                       icol=k
                    endif
                 else if (ipiv(k).gt.1) then
                    print *,'singular matrix in gauss', ind
!     pause
                    go to 25 
                 endif
 12           continue
           endif
 13     continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
           do 14 l=1,n
              dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
 14        continue
           do 15 l=1,m
              dum=b(irow,l)
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
 15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
           print *, 'singular matrix in gaussj', ind
!     pause
           go to 25
        endif
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
           a(icol,l)=a(icol,l)*pivinv
 16     continue
        do 17 l=1,m
           b(icol,l)=b(icol,l)*pivinv
 17     continue 
        
        do 21 ll=1,n
           if(ll.ne.icol)then
              dum=a(ll,icol)
              a(ll,icol)=0.d0
              do 18 l=1,n
                 a(ll,l)=a(ll,l)-a(icol,l)*dum                 
                 if(a(icol,l)*dum+a(ll,l).ne.a(ll,l)) then
                    if(dabs(a(ll,l)/(a(icol,l)*dum)).lt.eps) then
                       a(ll,l)=0.d0
                    endif
                 endif
 18           continue
              do 19 l=1,m
                 b(ll,l)=b(ll,l)-b(icol,l)*dum
                 if(b(icol,l)*dum+b(ll,l).ne.b(ll,l)) then
                    if(dabs(b(ll,l)/(b(icol,l)*dum)).lt.eps) then
                       b(ll,l)=0.d0
                    endif
                 endif
                 b_max(ll,l)=max(b_max(ll,l),dabs(b(icol,l)*dum)) 
 19           continue
           endif
 21     continue
 22   continue
      
      do ll=1,np
         do l=1,mp
            if(b_max(ll,l).ne.0.d0) then
               if(dabs(b(ll,l)/b_max(ll,l)).le.eps) then
                  b(ll,l)=0.d0
               endif
            endif
         enddo
      enddo
      
      
      do 24 l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         continue
         endif
 24   continue
 25   continue
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE solve_chem_simple(xk,xnH,T_K,y,dt,t_chem)
      implicit none
      integer,parameter :: N_sp=9,N_react=42
      integer :: i
      real*8 :: xnH, T_K, y(N_sp), xk(N_react), dt, C(N_sp), D(N_sp), 
     &     yHe, y_H, y_H2, y_e, y_Hp, y_H2p, y_Hm, y_He, y_Hep, y_Hepp,
     &     t_chem, r_f(N_sp)

      yHe    = 8.33333d-2
      y_H    = y(1)
      y_H2   = y(2)
      y_e    = y(3)
      y_Hp   = y(4)
      y_H2p  = y(5)
      y_Hm   = y(6)
      y_He   = y(7)
      y_Hep  = y(8)
      y_Hepp = y(9)

!     Hm [6]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_Hm  = (C(6)*dt+y_Hm)/(1.0d0+D(6)*dt)
!     H2+ [5]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_H2p = (C(5)*dt+y_H2p)/(1.0d0+D(5)*dt)
!     H2 [2]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_H2  = (C(2)*dt+y_H2)/(1.0d0+D(2)*dt)
!     H [1]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_H   = (C(1)*dt+y_H)/(1.0d0+D(1)*dt)
!     H+ [4]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_Hp  = (C(4)*dt+y_Hp)/(1.0d0+D(4)*dt)
!     H nuclei conservation
      if(y_H<y_Hp) y_Hp = 1.0d0-y_H -2.0d0*y_H2-2.0d0*y_H2p-y_Hm
      if(y_H>y_Hp) y_H  = 1.0d0-y_Hp-2.0d0*y_H2-2.0d0*y_H2p-y_Hm

!     He [7]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_He  = (C(7)*dt+y_He)/(1.0d0+D(7)*dt)
!     He+ [8]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_Hep  = (C(8)*dt+y_Hep)/(1.0d0+D(8)*dt)
!     He++ [9]
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      y_Hepp  = (C(9)*dt+y_Hepp)/(1.0d0+D(9)*dt)
!     He nuclei conservation
      if(y_He>y_Hep) then
         if(y_He>y_Hepp) then
            y_He = yHe - y_Hep - y_Hepp
         endif
      endif
      if(y_Hep>y_He) then
         if(y_Hep>y_Hepp) then
            y_Hep = yHe - y_He - y_Hepp
         endif
      endif
      if(y_Hepp>y_He) then
         if(y_Hepp>y_Hep) then
            y_Hepp = yHe - y_He - y_Hep
         endif
      endif
!     charge neutrality
!     e [3]
      y_e = y_Hp + y_Hep + 2.0d0*y_Hepp + y_H2p - y_Hm

!     convert
      y(1) = y_H
      y(2) = y_H2
      y(3) = y_e
      y(4) = y_Hp
      y(5) = y_H2p
      y(6) = y_Hm
      y(7) = y_He
      y(8) = y_Hep
      y(9) = y_Hepp

!     chemical time
      call react_rat(xk,xnH,T_K,y,C,D,r_f)
      do i=1,3
         t_chem = dmin1(t_chem, y(i)/dabs(C(i)-D(i)*y(i)))
      enddo
      return
      END


      SUBROUTINE react_coef(xnH, T_K, y_H, y_H2, xk)
      implicit none

      integer, parameter :: N_react=42, n_max=5
      integer :: n
      real*8 :: xnH, T_K, T_eV, xlnT_eV, y_H, y_H2, xk(N_react)
      real*8 :: gamma_h1, gamma_l1, gamma_h2, gamma_l2, nc_1, nc_2,
     &     pp, gamma, xk_L, xk_H, lgT4, n_cr, a, n_cr_H, n_cr_H2, 
     &     zH, zH2, chi_H2, z_n, xk_rr_A, xk_rr_B, xk_di

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Chemical reaction rate coefficients    !
!     mostly from GLover & Abel (2008)       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      chi_H2  =  5.1965d4      
      T_eV = 8.61735d-5*T_K
      xlnT_eV = dlog(T_eV)
      lgT4 = dlog10(T_K*1.d-4)
      n_cr_H  = 10.0d0**(3.0d0-0.416d0*lgT4-0.327d0*lgT4**2)
      n_cr_H2 = 10.d0**(4.845d0-1.3d0*lgT4+1.62d0*lgT4**2)
      n_cr = 1.0d0/(y_H/n_cr_H + 2.0d0*y_H2/n_cr_H2)
      a=1.d0/(1.d0+xnH/n_cr)

      zH=0.0d0
      do n=1,n_max
         z_n=2.d0*dble(n**2)*dexp(-1.57798d5*(1.d0-1.d0/dble(n**2))/T_K)
         zH=zH+z_n
      enddo

      zH2=1.0d1**(2.20859d0  -1.8089d0*dlog10(T_K)
     &     + 0.451858d0*(dlog10(T_K)**2.0d0))
            
      ! print*, n_cr, zH, zH2

!================= Hydrogen (1-20) =================!
!
!  1)   H     +   e     ->   H+    + 2 e   
!     Abel et al. (1997): GA08
      xk(1)=dexp(-32.71396786d0+(13.536556d0
     &     +(-5.73932875d0+(1.56315498d0+(-0.2877056d0
     &     +(3.48255977d-2+(-2.63197617d-3
     &     +(1.11954395d-4-2.03914985d-6*xlnT_eV)*xlnT_eV)
     &     *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)

!  2)   H+    +   e     ->   H     +   ph.
!     Glover & Jappsen (2007): GA08
!     (Case A)
!       xk(2)=1.269d-13*((315614d0/T_K)**(1.503d0))
!     &      *((1.0d0+(604625d0/T_K)**(0.47d0))**(-1.923d0))
!     (Case B)
      xk(2)=2.753d-14*((315614d0/T_K)**(1.5d0))
     &     *((1.0d0+(115188d0/T_K)**(0.407d0))**(-2.242d0))

!  3)   H     +   e     ->   H-    +   ph.
!     Wishart (1979): GA08 (01)
      if(T_K<=6.0d3) then
         xk(3) = 10.0d0**(-17.845d0 + 0.762d0 * dlog10(T_K)
     &        + 0.1523d0*(dlog10(T_K)**2.0d0)
     &        - 0.03274d0*(dlog10(T_K)**3.0d0))
      else
         xk(3) = 10.0d0**(-16.4199d0
     &        + 0.1998d0*(dlog10(T_K)**2.0d0)
     &        - 5.447d-3*(dlog10(T_K)**4.0d0)
     &        + 4.0415d-5*(dlog10(T_K)**6.0d0))  
      endif

!  4)   H-    +   H     ->   H2    +   e
!     Kreckel et al. (2010): 
!     The most recent experiment
      xk(4) = 1.35d-9 * (T_K**9.8493d-2 
     &      + 0.32852d0*T_K**0.5561d0 
     &      + 2.771d-7*T_K**2.1826d0)
     &      /(1.0d0 + 6.1910d-3*T_K**1.0461d0 
     &      + 8.9712d-11*T_K**3.0424d0  
     &      + 3.2576d-14*T_K**3.7742d0)
      
!  5)   H     +   H+    ->   H2+   +   ph.
!     Coppola et al. (2011):
      if(T_K<=30.0d0) then
         xk(5) = 2.1d-20/(T_K/30.0d0)**0.15d0
      else
         xk(5) = 10.0d0**(-18.20d0
     &        - 3.194d0*dlog10(T_K) 
     &        + 1.786d0*dlog10(T_K)**2.0d0 
     &        - 0.2072d0*dlog10(T_K)**3.0d0)
      endif

!  6)   H2+   +   H     ->   H2    +   H+         
!     Galli & Palla (1998): 
      xk(6)=6.4d-10 
      
!  7)   H2    +   H     -> 3 H              
!     Martin et al. (1996)
   ! use Glover & Abel instead; wli0313
      gamma_h1=-1.784239d2-6.842243d1*dlog10(T_K)
     &     +4.320243d1*(dlog10(T_K)**2.0)-4.633167d0*(dlog10(T_K)**3.0)
     &     +6.970086d1*dlog10(1.0d0+4.087038d4/T_K)
      gamma_h2=-2.370570d4/T_K
      gamma_l1=1.288953d2-5.391334d1*dlog10(T_K)
     &     +5.315517d0*(dlog10(T_K)**2.0)
     &     -1.973427d1*dlog10(1.0d0+1.678095d4/T_K)
      gamma_l2=-2.578611d4/T_K
      nc_1=1.482123d1-4.890915d0*dlog10(T_K)
     &     +4.749030d-1*(dlog10(T_K)**2.0)-1.338283d2/T_K
      nc_2=-1.164408d0+nc_1
      pp=8.227443d-1+5.864073d-1*dexp(-T_K/1850d0)
     &     -2.056313d0*dexp(-T_K/440d0)
      nc_1=1.0d1**(nc_1)
      nc_2=1.0d1**(nc_2)
      gamma=gamma_h1
     &     -(gamma_h1-gamma_l1)/(1.0d0+(xnH/nc_1)**pp)
     &     +gamma_h2
     &     -(gamma_h2-gamma_l2)/(1.0d0+(xnH/nc_2)**pp)
      xk(7)=1.0d1**(gamma)

!  8)   H2    +   H+    ->   H2+   +   H        
!     Savin (2004) low density 
      xk_L=(-3.3232183d-7 + 3.3735382d-7*dlog(T_K)
     &     - 1.4491368d-7*  (dlog(T_K)**2.0d0)
     &     + 3.4172805d-8*  (dlog(T_K)**3.0d0)
     &     - 4.7813720d-9*  (dlog(T_K)**4.0d0)
     &     + 3.9731542d-10* (dlog(T_K)**5.0d0)
     &     - 1.8171411d-11* (dlog(T_K)**6.0d0)
     &     + 3.5311932d-13* (dlog(T_K)**7.0d0))
     &     *dexp(-21237.15d0/T_K)
!     Coppola et al.(2011) high density
      xk_H=dexp(-33.081d0+6.3173d-5*T_K-2.3478d4/T_K
     &     -1.8691d-9*(T_K**2.0d0))*1.0d6   
!     connect between v=0 and LTE
      if(a.eq.1.d0) then
         xk(8)=xk_L
      elseif(a.eq.0.d0) then
         xk(8)=xk_H
      else
         xk(8)=(xk_H**(1.d0-a))*(xk_L**a)
      endif    

!  9)   H2    +   e     -> 2 H     +   e
!     Trevisan & Tennyson (2002): GA08
      xk_L = 4.49d-9*(T_K**0.11d0) *dexp(-101858d0/T_K)
      xk_H = 1.91d-9*(T_K**0.136d0)*dexp(-53407.1d0/T_K)
!     connect between v=0 and LTE
      if(a.eq.1.d0) then
         xk(9)=xk_L
      elseif(a.eq.0.d0) then
         xk(9)=xk_H
      else
         xk(9)=(xk_H**(1.d0-a))*(xk_L**a)
      endif

!  10)   H-    +   e     ->   H     + 2 e              
!     Janev et al. (1987): GA08
      xk(10)=dexp(-18.01849334d0+(2.3608522d0
     &     +(-0.28274430d0+(1.62331664d-2+(-3.36501203d-2
     &     +(1.17832978d-2+(-1.65619470d-3+(1.06827520d-4
     &     -2.63128581d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
     &     *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)

!  11)   H-    +   H+    -> 2 H
!     Croft, Dickinson & Gadea (1999): GA08
      xk(11) = 2.4d-6/dsqrt(T_K)*(1.0d0+T_K/2.0d4)

!  12)   H-    +   H+    ->   H2+   +   e          
!     Poulaert et al. (1978): GA08
      if(T_K.le.8.d3) then
         xk(12)=6.9d-9/(T_K**3.5d-1)
      else
         xk(12)=9.6d-7/(T_K**9.d-1)
      endif

!  13)   H2+   +   e     -> 2 H                
!     Savin et al. (2004): GA08
      if(T_K<=617.0d0) then
         xk(13)=1.0d-8
      else
         xk(13)=1.32d-6*(T_K**(-0.76d0))
      endif

!  14)   H2+   +   H-    ->   H2    +   H         
!     Dalgarno & Lepp (1987): GA08
      xk(14)=1.40d-7/dsqrt(T_K/300.d0)

!  15) 3 H               ->   H2    +   H                
!     inverse rate of dissociation
      xk(15)=xk(7)*zH2/(zH**2)*1.493d-20
     &     /(T_K**1.5d0)*dexp(chi_H2/T_K)

!  16) 2 H2           -> 2 H    +   H2
!     Martin et al. (1998)
      xk_L=5.996d-30*(T_K**4.1881d0)
     &     /(1.0d0+6.761d-6*T_K)**5.6881d0
     &     *dexp(-54657.4d0/T_K)
      xk_H=1.3d-9*dexp(-53300.0d0/T_K)
!     connect between v=0 and LTE
      if(a.eq.1.d0) then
         xk(16)=xk_L
      elseif(a.eq.0.d0) then
         xk(16)=xk_H
      else
         xk(16)=xk_H**(1.d0-a)*xk_L**a
      endif

!  17) 2 H    +   H2  -> 2 H2
!     inverse rate of dissociation
      xk(17)=xk(16)*zH2/(zH**2)*1.493d-20
     &     /(T_K**1.5d0)*dexp(chi_H2/T_K)

!  18) H-   +   H  -> 2 H   +  e 
!     Janev et al. (1987): GA08
      if(T_eV<=0.1d0) then
         xk(18) = 2.5634d-9 * (T_eV**1.78186d0)
      else
         xk(18) = dexp(-2.0372609d1 
     &        +1.13944933d0*xlnT_eV
     &        -1.4210135d-1*(xlnT_eV**2.0d0)
     &        +8.4644554d-3*(xlnT_eV**3.0d0)
     &        -1.4327641d-3*(xlnT_eV**4.0d0)
     &        +2.0122503d-4*(xlnT_eV**5.0d0)
     &        +8.6639632d-5*(xlnT_eV**6.0d0)
     &        -2.5850097d-5*(xlnT_eV**7.0d0)
     &        +2.4555012d-6*(xlnT_eV**8.0d0)
     &        -8.0683825d-8*(xlnT_eV**9.0d0))
      endif

!  19) H-   +   H2+  -> 3 H
!     Dalgarno & Lepp (1987): GA08
      xk(19) = 1.4d-7/dsqrt(T_K/3.0d2)

!  20) H2   +   e   ->  H-  +  H
!     Schulz & Asundi (1967): GA08
      xk(20) = 2.7d-8/T_K**1.27d0*dexp(-43000.0d0/T_K)

!  21-30) blank:
      do n=21,30
         xk(n)=0.0d0
      enddo
!================= Helium (31-) =================! 

!  31) He   +   e   ->   He+   +   2 e
      xk(31) = dexp(-4.409864886d1 
     &     + 2.391596563d1*xlnT_eV
     &     - 1.07532302d1 *(xlnT_eV)**2.0d0
     &     + 3.05803875d0 *(xlnT_eV)**3.0d0
     &     - 5.68511890d-1*(xlnT_eV)**4.0d0
     &     + 6.79539123d-2*(xlnT_eV)**5.0d0
     &     - 5.00905610d-3*(xlnT_eV)**6.0d0
     &     + 2.06723616d-4*(xlnT_eV)**7.0d0
     &     - 3.64916141d-6*(xlnT_eV)**8.0d0)

!  32) He+  +   e   ->   He   +   ph.
      xk_rr_A = 1.0d-11/dsqrt(T_K)
     &     *( 12.72d0 - 1.615*dlog10(T_K)
     &     -  0.3162d0*(dlog10(T_K))**2.0d0
     &     +  0.0493d0*(dlog10(T_K))**3.0d0)
      xk_rr_B = 1.0d-11/dsqrt(T_K)
     &     *( 11.19d0 - 1.676d0*dlog10(T_K)
     &     -  0.2852d0*(dlog10(T_K))**2.0d0
     &     +  4.433d-2*(dlog10(T_K))**3.0d0)
      xk_di = 1.9d-3/(T_K**1.5d0)*dexp(-473421.d0/T_K)
     &     *(1.0d0 + 0.3d0*dexp(-94684d0/T_K))

      xk(32) = 0.68d0*xk_rr_A + 0.32d0*xk_rr_B + xk_di

!  33) He+  +   e   ->   He++  +   2 e
      xk(33) = dexp(-6.87104099d1 
     &     + 4.393347633d1 *xlnT_eV
     &     - 1.848066990d1 *(xlnT_eV)**2.0d0
     &     + 4.701626490d0 *(xlnT_eV)**3.0d0
     &     - 7.69246630d-1 *(xlnT_eV)**4.0d0
     &     + 8.11304200d-2 *(xlnT_eV)**5.0d0
     &     - 5.32402063d-3 *(xlnT_eV)**6.0d0
     &     + 1.97570531d-4 *(xlnT_eV)**7.0d0
     &     - 3.16558106d-6 *(xlnT_eV)**8.0d0)

!  34) He++  +   e   ->   He+  +   ph.    
      xk(34) = 5.506d-14*((1262456.0d0/T_K)**1.5d0)
     &     /((1.0d0+(460752.0d0/T_K)**0.407d0)**2.242d0)

!  35) H2   +   He   ->   He   +  H  +  H
      xk_L = 10.0d0**(-27.029d0+3.801d0*dlog10(T_K)-29487.0d0/T_K)
      xk_H = 10.0d0**(-2.729d0 -1.75d0 *dlog10(T_K)-23474.0d0/T_K)
      if(a.eq.1.d0) then
         xk(35)=xk_L
      elseif(a.eq.0.d0) then
         xk(35)=xk_H
      else
         xk(35)=xk_H**(1.d0-a)*xk_L**a
      endif

!  36) H2   +   He+   ->   He   +  H  +  H+
      xk(36) = 3.7d-14*dexp(35.0d0/T_K)

!  37) H2   +   He+   ->   H2+   +   He
      xk(37) = 7.2d-15

!  38) He+  +   H   ->   He   +   H+
      xk(38) = 1.2d-15*(T_K/3.0d2)**0.25d0

!  39) He   +   H+  ->   He+  +   H
      if(T_K<=1.0d4) then
         xk(39) = 1.26d-9/(T_K**0.75d0)*dexp(-127500.0d0/T_K)
      else
         xk(39) = 4.0d-37*(T_K**4.74d0)
      endif
      
!  40) He+  +   H-  ->   He   +   H
      xk(40) = 2.32d-7/(T_K/3.0d0)**0.52d0 * dexp(T_K/22400.0d0)  !wli //  should be T_K/3.0d2   see GA08 TableA1-28

!  41) He   +   H-  ->   He   +   H   +   e
      xk(41) = 4.1d-17*(T_K**2.0d0)*dexp(-19870.0d0/T_K)

!  42) H  +  H  +  He  ->  H2  +  He
      xk(42) = 6.9d-32/T_K**0.4d0
 
      return
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE react_rat(xk,xnH,T_K,y,C,D,r_f_tot)
      implicit none
****************************************************************
*     dy(i)/dt=r_f_tot(i)                                      * 
*     This subroutine returns reaction rate for each spieces,  *
*     r_f_tot(i).                                              *
****************************************************************
C     N_sp = number of spiecies
C     N_react = number of reactions
      integer,parameter :: N_sp=9,N_react=42
      integer :: ire, isp
      real*8 :: xnH, T_K, y(N_sp), r_f(N_react,N_sp), r_f_tot(N_sp),
     &     xk(N_react), r_react, C(N_sp), D(N_sp)
c
****************************************************
*     SPECIES                                      *
*     1 : H      2 : H2     3 : e      4 : H+      *
*     5 : H2+    6 : H-                            *
*     7 : He     8 : He+    9 : He++               *
****************************************************

      do isp=1,N_sp
         do ire=1,N_react
            r_f(ire,isp)=0.d0
         enddo
      enddo
!  1)   H     +   e     ->   H+    +  2 e
!       1         3          4        2 * 3
      r_react  = xk(1) * y(1) * y(3) * xnH
      r_f(1,1) =-r_react
      r_f(1,3) = r_react
      r_f(1,4) = r_react

!  2)   H+    +   e     ->   H     +   ph.
!       4         3          1
      r_react  = xk(2) * y(3) * y(4) * xnH
      r_f(2,1) = r_react
      r_f(2,3) =-r_react
      r_f(2,4) =-r_react

!  3)   H     +   e     ->   H-    +   ph.
!       1         3          6
      r_react  = xk(3) * y(1) * y(3) * xnH
      r_f(3,1) =-r_react
      r_f(3,3) =-r_react
      r_f(3,6) = r_react

!  4)   H-    +   H     ->   H2    +   e
!       6         1          2         3
      r_react  = xk(4) * y(1) * y(6) * xnH
      r_f(4,1) =-r_react
      r_f(4,2) = r_react
      r_f(4,3) = r_react
      r_f(4,6) =-r_react

!  5)   H     +   H+    ->   H2+   +   ph.
!       1         4          5
      r_react  = xk(5) * y(1) * y(4) * xnH
      r_f(5,1) =-r_react
      r_f(5,4) =-r_react
      r_f(5,5) = r_react
  
!  6)   H2+   +   H     ->   H2    +   H+
!       5         1          2         4
      r_react  = xk(6) * y(1) * y(5) * xnH
      r_f(6,1) =-r_react
      r_f(6,2) = r_react
      r_f(6,4) = r_react
      r_f(6,5) =-r_react

!  7)   H2    +   H     -> 3 H
!       2         1        3 * 1
      r_react  = xk(7) * y(1) * y(2) * xnH
      r_f(7,1) = r_react*2.0d0
      r_f(7,2) =-r_react
 
!  8)   H2    +   H+    ->   H2+   +   H
!       2         4          5         1
      r_react  = xk(8) * y(2) * y(4) * xnH
      r_f(8,1) = r_react
      r_f(8,2) =-r_react
      r_f(8,4) =-r_react
      r_f(8,5) = r_react

!  9)   H2    +   e     -> 2 H     +   e
!       2         3        2 * 1       3
      r_react  = xk(9) * y(2) * y(3) * xnH
      r_f(9,1) = r_react*2.0d0
      r_f(9,2) =-r_react

!  10)  H-    +   e     ->   H     +  2 e
!       6         3          1        2 * 3
      r_react  = xk(10) * y(3) * y(6) * xnH
      r_f(10,1) = r_react
      r_f(10,3) = r_react
      r_f(10,6) =-r_react
  
!  11)  H-    +   H+    -> 2 H
!       6         4        2 * 1
      r_react  = xk(11) * y(4) * y(6) * xnH
      r_f(11,1) = r_react*2.0d0
      r_f(11,4) =-r_react
      r_f(11,6) =-r_react
  
!  12)  H-    +   H+    ->   H2+   +   e
!       6         4          5         3
      r_react  = xk(12) * y(4) * y(6) * xnH
      r_f(12,3) = r_react
      r_f(12,4) =-r_react
      r_f(12,5) = r_react
      r_f(12,6) =-r_react

!  13)  H2+   +   e     -> 2 H
!       5         3        2 * 1
      r_react  = xk(13) * y(3) * y(5) * xnH
      r_f(13,1) = r_react*2.0d0
      r_f(13,3) =-r_react
      r_f(13,5) =-r_react
  
!  14)  H2+   +   H-    ->   H2    +   H
!       5         6          2         1
      r_react  = xk(14) * y(5) * y(6) * xnH
      r_f(14,1) = r_react
      r_f(14,2) = r_react
      r_f(14,5) =-r_react
      r_f(14,6) =-r_react

!  15) 3 H              ->   H2    +   H
!      3 * 1                 2         1
      r_react  = xk(15) * y(1) * y(1) * y(1) * xnH * xnH
      r_f(15,1) =-2.0d0*r_react
      r_f(15,2) = r_react

!  16) 2 H2             -> 2 H     +   H2
!      2 * 2               2 * 1       1
      r_react  = xk(16) * y(2) * y(2) * xnH
      r_f(16,1) = r_react*2.0d0
      r_f(16,2) =-r_react

!  17) 2 H    +   H2    -> 2 H2
!      2 * 1      2        2 * 2
      r_react  = xk(17) * y(1) * y(1) * y(2) * xnH * xnH
      r_f(17,1) =-2.0d0*r_react
      r_f(17,2) = r_react

!  18)  H-    +   H     -> 2 H     +   e
!       6         1        2 * 1       3
      r_react  = xk(18) * y(1) * y(6) * xnH
      r_f(18,1) = r_react
      r_f(18,3) = r_react
      r_f(18,6) =-r_react

!  19)  H-    +   H2+   -> 3 H
!       6         5        3 * 1
      r_react  = xk(19) * y(5) * y(6) * xnH
      r_f(19,1) = r_react*3.0d0
      r_f(19,5) =-r_react
      r_f(19,6) =-r_react

!  20)  H2    +   e     ->   H-    +   H
!       2         3          6         1
      r_react  = xk(20) * y(2) * y(3) * xnH
      r_f(20,1) = r_react
      r_f(20,2) =-r_react
      r_f(20,3) =-r_react
      r_f(20,6) = r_react

!  31)  He    +   e     ->   He+   +  2 e
!       7         3          8        2 * 3
      r_react  = xk(31) * y(3) * y(7) * xnH
      r_f(31,3) = r_react
      r_f(31,7) =-r_react
      r_f(31,8) = r_react

!  32)  He+   +   e     ->   He    +   ph.
!       8         3          7
      r_react  = xk(32) * y(3) * y(8) * xnH
      r_f(32,3) =-r_react
      r_f(32,7) = r_react
      r_f(32,8) =-r_react

!  33)  He+   +   e     ->   He++  +  2 e
!       8         3          9        2 * 3
      r_react  = xk(33) * y(3) * y(8) * xnH
      r_f(33,3) = r_react
      r_f(33,8) =-r_react
      r_f(33,9) = r_react

!  34)  He++  +   e     ->   He+   +   ph.
!       9         3          8
      r_react  = xk(34) * y(3) * y(9) * xnH
      r_f(34,3) =-r_react
      r_f(34,8) = r_react
      r_f(34,9) =-r_react

!  35)  H2    +   He    ->   He    +  2 H
!       2         7          7        2 * 1
      r_react  = xk(35) * y(2) * y(7) * xnH
      r_f(35,1) = r_react*2.0d0
      r_f(35,2) =-r_react
   
!  36)  H2    +   He+   ->   He    +   H   +   H+
!       2         8          7         1       4
      r_react  = xk(36) * y(2) * y(8) * xnH
      r_f(36,1) = r_react
      r_f(36,2) =-r_react
      r_f(36,4) = r_react
      r_f(36,7) = r_react
      r_f(36,8) =-r_react

!  37)  H2    +   He+   ->   H2+   +   He
!       2         8          5         7
      r_react  = xk(37) * y(2) * y(8) * xnH
      r_f(37,2) =-r_react
      r_f(37,5) = r_react
      r_f(37,7) = r_react
      r_f(37,8) =-r_react
 
!  38)  He+   +   H     ->   He    +   H+
!       8         1          7         4
      r_react  = xk(38) * y(1) * y(8) * xnH
      r_f(38,1) =-r_react
      r_f(38,4) = r_react
      r_f(38,7) = r_react
      r_f(38,8) =-r_react

!  39)  He    +   H+    ->   He+   +   H
!       7         4          8         1
      r_react  = xk(39) * y(4) * y(7) * xnH
      r_f(39,1) = r_react
      r_f(39,4) =-r_react
      r_f(39,7) =-r_react
      r_f(39,8) = r_react

!  40)  He+   +   H-    ->   He    +   H
!       8         6          7         1
      r_react  = xk(40) * y(6) * y(8) * xnH
      r_f(40,1) = r_react
      r_f(40,6) =-r_react
      r_f(40,7) = r_react
      r_f(40,8) =-r_react

!  41)  He    +   H-    ->   He    +   H   +   e
!       7         6          7         1       3
      r_react  = xk(41) * y(6) * y(7) * xnH
      r_f(41,1) = r_react
      r_f(41,3) = r_react
      r_f(41,6) =-r_react
 
!  42) 2 H    +   He    ->   H2    +   He
!      2 * 1      7          2         7
      r_react  = xk(42) * y(1) * y(1) * y(7) * xnH * xnH
      r_f(42,1) =-2.0d0*r_react
      r_f(42,2) = r_react


      do isp=1,N_sp      
         r_f_tot(isp)=0.d0
         C(isp)=0.0d0
         D(isp)=0.0d0
         do ire=1,N_react
            r_f_tot(isp)=r_f_tot(isp)+r_f(ire,isp)
            if(r_f(ire,isp)>0.0d0) then
               C(isp)=C(isp)+r_f(ire,isp)
            else
               D(isp)=D(isp)+r_f(ire,isp)
            endif
         enddo
         D(isp) = -D(isp)/y(isp)
      enddo

      return
      END

      
      
