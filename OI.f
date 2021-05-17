      PROGRAM test
      implicit none 
      integer :: i,k
      real*8 :: xnH, T_K, xNc_OI, y_m, y_a, y_e,
     &          tau_c, xLd_OI, a, beta_esc
      real*8 :: esc(1:3)

      xnH = 1.0d3
      T_K = 1.0d3
      xNc_OI = 1.d10
      y_m = 1.d-3
      y_a = 1.0
      y_e = 1.d-3
      esc(1) = 0.5
      esc(2) = 0.5
      esc(3) = 0.5
      tau_c = 0.
      xLd_OI = 0.0

      call OIcool(xnH,T_K,xNc_OI,y_m,y_a,y_e,esc,tau_c,xLd_OI)
      print*, xLd_OI
      END

c   gfortran OI.f CII.f GaussJordan.f -o cooling.out && ./cooling.out


      subroutine OIcool(xnH,T_K,xNc_OI,y_m,y_a,y_e,
     &     esc,tau_c,xLd_OI)
      implicit real*8(a-h,o-z)
      real*8 :: esc(1:3)
      data xk_B/1.38d-16/,h_P/6.63d-27/,xm_p/1.67d-24/
      common /lines3/xnu_10,xnu_20,xnu_21,Q_10,Q_20,Q_21,
     &     A_10,A_20,A_21,C_10,C_20,C_21,C_01,C_02,C_12,
     &     g_0,g_1,g_2,xNclmn,v_th,tau_cnt
      tau_cnt=tau_c
      if(T_K < 10.d0) then
         xLd_OI=0.d0
         return
      endif
C     set coefficients
      g_0=5.d0
      g_1=3.d0
      g_2=1.d0
c
      DT_10=2.3d2
      DT_20=3.28d2
      DT_21=9.8d1
c
      DE_10=DT_10*xk_B
      DE_20=DT_20*xk_B
      DE_21=DT_21*xk_B
c     
      xnu_10=DE_10/h_P
      xnu_20=DE_20/h_P
      xnu_21=DE_21/h_P
c
      A_10=9.0d-5
      A_20=1.0d-10
      A_21=1.7d-5
c
      Q_10=Q_bg(DT_10)
      Q_20=Q_bg(DT_20)
      Q_21=Q_bg(DT_21)
c
      gamma10_e=1.4d-8
      gamma20_e=1.4d-8
      gamma21_e=5.0d-9
c
      gamma10_H=9.2d-11*(T_K*1.d-2)**0.67d0
      gamma20_H=4.3d-11*(T_K*1.d-2)**0.80d0
      gamma21_H=1.1d-10*(T_K*1.d-2)**0.44d0
c
      gamma10_H2=5.d-1*gamma10_H
      gamma20_H2=5.d-1*gamma20_H
      gamma21_H2=5.d-2*gamma21_H
c
      C_10=xnH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2)
      C_20=xnH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2)
      C_21=xnH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2)
      C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)
      C_02=(g_2/g_0)*C_20*dexp(-DT_20/T_K)
      C_12=(g_2/g_1)*C_21*dexp(-DT_21/T_K)
c
      xNclmn=xNc_OI
      xm_O=16.d0*xm_p
      v_th=dsqrt(2.d0*xk_B*T_K/xm_O)           
      call pop3lev(f_0,f_1,f_2,esc)
c
      S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
      S_20=1.d0/(g_2*f_0/(g_0*f_2)-1.d0)
      S_21=1.d0/(g_2*f_1/(g_1*f_2)-1.d0)
c
      x_10=esc(1)*(1.d0-Q_10/S_10)
      x_20=esc(2)*(1.d0-Q_20/S_20)
      x_21=esc(3)*(1.d0-Q_21/S_21)
c
      xLd_OI=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20
     &     +DE_21*A_21*f_2*x_21)/xnH
      return
      END
      

      subroutine pop3lev(f_0,f_1,f_2,esc)
      implicit real*8(a-h,o-z)
      parameter (err_eps=1.d-5)
      dimension esc(3),error(3),desc(3),     
     &     esc_f(3),error_f(3),A(3,3)
c
      do 10 itr=1,100
c     determin population under given esc
         call thrlev(f_0,f_1,f_2,esc,error)
c     evaluate error
         err_max=0.d0
         do i=1,3
            err_max=max(dabs(error(i)),err_max)
         enddo
c     error small enough ?
         if(err_max.lt.err_eps) go to 20
c     if not, improve guess for esc
c     by Newton-Raphson scheme
         do j=1,3
            if(esc(j).eq.0.d0) then
               delta_esc=1.d-10
            else
               delta_esc=1.d-5*esc(j)
            endif
c 
            do jj=1,3
               if(jj.eq.j) then
                  esc_f(jj)=esc(jj)+delta_esc
               else
                  esc_f(jj)=esc(jj)
               endif
            enddo            
            call thrlev(f_0,f_1,f_2,esc_f,error_f)
C ------ set the matrix A
            do i=1,3               
               A(i,j)=(error_f(i)-error(i))/(esc_f(j)-esc(j))
            enddo
         enddo
C ------- set the vector desc
         do i=1,3
            desc(i)=-error(i)
         enddo
C ------- solve linear equations
         call gaussj(A,3,3,desc,1,1,0)
c     !wli in solving reaction: call gaussj(a=A,n=N_sp,np=N_sp,b=ddy,m=1,mp=1,ind=1)

         fact=1.d0
         if(itr.gt.20) then
            if(err_max.gt.1.d0) fact=1.d-2
         endif
         do i=1,3
            if(esc(i)*desc(i).ne.0.d0) then
               fact=min(fact,4.d-1*dabs(esc(i)/desc(i)))
            endif 
         enddo
c     next guess for esc
         do i=1,3
            esc(i)=esc(i)+fact*desc(i)
         enddo 
 10   continue
c      
 20   continue     
      return
      END


      subroutine thrlev(f_0,f_1,f_2,esc,error)
      implicit real*8(a-h,o-z)
      dimension esc(3),error(3)
      common /lines3/xnu_10,xnu_20,xnu_21,Q_10,Q_20,Q_21,
     &     A_10,A_20,A_21,C_10,C_20,C_21,C_01,C_02,C_12,
     &     g_0,g_1,g_2,xNclmn,v_th,tau_cnt
      data pi/3.14159265358979d0/

      esc_10=esc(1)
      esc_20=esc(2)
      esc_21=esc(3)

      R_10=esc_10*A_10*(1.d0+Q_10)+C_10
      R_20=esc_20*A_20*(1.d0+Q_20)+C_20
      R_21=esc_21*A_21*(1.d0+Q_21)+C_21
      R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
      R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
      R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12
      
      f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
     &     /( (R_01+R_02+R_20)*(R_10+R_12+R_21) 
     &     -(R_01-R_21)*(R_10-R_20) )
      f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
      f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

      xNc_0=xNclmn*f_0
      xNc_1=xNclmn*f_1
      xNc_2=xNclmn*f_2

      tau_10=(A_10/8.d0/pi)*(3.d10/xnu_10)**3
     &     *(xNc_0*g_1/g_0-xNc_1)/v_th 
      tau_20=(A_20/8.d0/pi)*(3.d10/xnu_20)**3
     &     *(xNc_0*g_2/g_0-xNc_2)/v_th 
      tau_21=(A_21/8.d0/pi)*(3.d10/xnu_21)**3
     &     *(xNc_1*g_2/g_1-xNc_2)/v_th 

      err_10=esc_10-beta_esc(tau_10,tau_cnt)
      err_20=esc_20-beta_esc(tau_20,tau_cnt)
      err_21=esc_21-beta_esc(tau_21,tau_cnt)
      
      error(1)=err_10
      error(2)=err_20
      error(3)=err_21
      
      return
      END



