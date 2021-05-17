   !    PROGRAM test
   !    implicit none 
   !    integer :: i,k
   !    real*8 :: xnH, T_K, xNc_CII, y_m, y_a, y_e,
   !   &     esc_10, tau_c, xLd_CII, a, beta_esc

   !    xnH = 1.0d3
   !    T_K = 1.0d3
   !    xNc_CII = 1.d10
   !    y_m = 1.d-3
   !    y_a = 1.0
   !    y_e = 1.d-3
   !    esc_10 = 0.5
   !    tau_c = 0.5
   !    xLd_CII = 0.0

   !    call CIIcool(xnH,T_K,xNc_CII,y_m,y_a,y_e,esc_10,tau_c,xLd_CII)

   !    END

c   (alias gg='gfortran cooling_rate.f CII.f GaussJordan.f -o cooling.o && ./cooling.o')
c   gfortran OI.f CII.f GaussJordan.f -o cooling.out && ./cooling.out

      subroutine CIIcool(xnH,T_K,xNc_CII,y_m,y_a,y_e,
     &     esc_10,tau_c,xLd_CII)
      implicit real*8(a-h,o-z)
      data xk_B/1.38d-16/,h_P/6.63d-27/,xm_p/1.67d-24/,
     &     pi/3.14159265358979d0/
      common /line2/Q_10,A_10,C_10,C_01,g_0,g_1,tau0,tau_cnt
      tau_cnt=tau_c
C     set coefficients
      g_0=2.d0
      g_1=4.d0
c
      DT_10=92.d0
      DE_10=DT_10*xk_B
      xnu_10=DE_10/h_P
c
      Q_10=Q_bg(DT_10) ! wli: = 0
c
      A_10=2.4d-6
      gamma_e=2.8d-7/dsqrt(T_K*1.d-2)
      gamma_H=8.0d-10*(T_K*1.d-2)**7.d-2
      gamma_H2=5.d-1*gamma_H
      C_10=xnH*(y_e*gamma_e+y_a*gamma_H+y_m*gamma_H2)
      C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)

      xm_C=12.d0*xm_p
      v_th=dsqrt(2.d0*xk_B*T_K/xm_C)
      tau0=(A_10/8.d0/pi)*(3.d10/xnu_10)**3*xNc_CII/v_th

      call CIIpop(f_0,f_1,esc_10)

      S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
      xLd_CII=DE_10*A_10*f_1*esc_10*(1.d0-Q_10/S_10)/xnH

      print*, xLd_CII
      return
      END

      subroutine CIIpop(f_0,f_1,esc_10)
      implicit real*8(a-h,o-z)
      parameter (err_eps=1.d-4)
      dimension esc(1),error(1),desc(1),     
     &     esc_f(1),error_f(1),A(1,1)

      esc(1)=esc_10

      do 10 itr=1,1000
         call twolevel(f_0,f_1,esc,error)

c     evaluate error
         err_max=0.d0
         do i=1,1
            err_max=max(dabs(error(i)),err_max)
         enddo
         if(err_max.lt.err_eps) go to 20

         do j=1,1
            if(esc(j).eq.0.d0) then
               delta_esc=1.d-10
            else
               delta_esc=1.d-5*esc(j)
            endif

            do jj=1,1
               if(jj.eq.j) then
                  esc_f(jj)=esc(jj)+delta_esc
               else
                  esc_f(jj)=esc(jj)
               endif
            enddo

            call twolevel(f_0,f_1,esc_f,error_f)
C ------ set the matrix A
            do i=1,1              
               A(i,j)=(error_f(i)-error(i))/(esc_f(j)-esc(j))
            enddo
         enddo

C ------- set the vector desc
         do i=1,1
            desc(i)=-error(i)
         enddo

C ------- solve linear equations
         call gaussj(A,1,1,desc,1,1,0)
c     !wli in solving reaction: call gaussj(a=A,n=N_sp,np=N_sp,b=ddy,m=1,mp=1,ind=1)

         fact=1.d0
         if(itr.gt.20) then
            if(err_max.gt.1.d0) fact=1.d-2
         endif
         do i=1,1
            if(esc(i)*desc(i).ne.0.d0) then
               fact=min(fact,4.d-1*dabs(esc(i)/desc(i)))
            endif 
         enddo

c     improve esc
         do i=1,1
            esc(i)=esc(i)+fact*desc(i)
         enddo 
 10   continue
      
 20   continue            
      esc_10=esc(1)
      return
      END

      subroutine twolevel(f_0,f_1,esc,error)
c     solves level population (f_0, f_1) of two-level system
c     for given escape probability (esc)
c     returns the difference of old and new esc as (error). 
      implicit real*8(a-h,o-z)
      dimension esc(1),error(1)
      common /line2/Q_10,A_10,C_10,C_01,g_0,g_1,tau0,tau_cnt
      esc_10=esc(1)
C     population of level 0 & 1
      R_01=(g_1/g_0)*A_10*esc_10*Q_10+C_01
      R_10=esc_10*A_10*(1.d0+Q_10)+C_10
      f_1=R_01/(R_10+R_01)
      f_0=R_10/(R_10+R_01)
      tau_10=tau0*(f_0*g_1/g_0-f_1)
      esc(1)=beta_esc(tau_10,tau_cnt)
      error(1)=esc_10-esc(1)
      return
      END


      FUNCTION Q_bg(T_nu)
      IMPLICIT REAL*8(a-h,o-z)
c      parameter(T_rad=27.3d0)
c      x=T_nu/T_rad
c      if(x.gt.1.d2) then
c         Q_bg_CMB=0.d0
c      else
c         Q_bg_CMB=1.d0/(dexp(x)-1.d0)
c      endif
c      Q_bg=Q_bg_CMB
      Q_bg=0.d0
      return
      END

      FUNCTION beta_esc(tau_L,tau_C)
      IMPLICIT REAL*8(a-h,o-z)      
      if(tau_L.lt.0.d0) then
         beta_esc=1.d0
      elseif(tau_L.lt.1.d-5) then 
         beta_esc=dexp(-tau_C)
      else
         beta_esc=dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
      endif
      return
      END

