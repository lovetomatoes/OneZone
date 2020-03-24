
cd      PROGRAM set_nu
!     gfortran set_freq.f -o set_freq.o && ./set_freq.o
      implicit none
      integer, parameter :: NF=60, N_nu_th=9
      integer :: i_sw_fr, k, k_fine, k_th(N_nu_th)
      real*8 :: nu_prev, nu(NF), h_p, eV, xnurat1, xnurat2, nu_min, 
     &     xnurat, eVhp09, xk_B
      real*8 :: dnu(NF), xk_pd_Hm, xk_pd_H2p, pi=4.d0*datan(1.d0),
     &     Cph=3.d10 !speed of light
      real*8,parameter :: nu_th(N_nu_th)=(/0.74d0, 2.65d0, 12.24d0,
     &     13.51d0, 13.6d0, 15.4d0, 24.586d0, 30.0d0, 54.42d0/)
      real*8 :: sigma_H_nu(NF), sigma_He_nu(NF), sigma_Hep_nu(NF), 
     &     sigma_Hm_nu(NF), sigma_H2p_nu(NF), Planck(NF), T_rad,
     &     Flux(NF)

      real*8 :: Ta(10),rnua(45),sigmaa(10,45),
     &     rlTaa(148),rlaa(148)
      common /op_H2p_bf/Ta,rnua,sigmaa,rlTaa,rlaa

      h_p = 6.63d-27
      eV = 1.60217657d-12
      xk_B = 1.38d-16
      T_rad = 2.d05
      print*, "T_rad=", T_rad

      xnurat1 = 1.1d0
      xnurat2 = 1.05d0
      nu_min  = 3.0d15
C wli changed: from hp nu = 0.5 eV
      nu_min = 0.45*eV/h_p

      k_fine=5

      i_sw_fr=5
      xnurat=xnurat1
      nu_prev=nu_min!/xnurat
      eVhp09=eV/h_p*0.9999

      call H2p_bf_data
      do k=1,NF
         if(i_sw_fr.gt.N_nu_th+1)then
            print *, "i_sw_fr > N_nu_th+1"
            stop
         endif
         nu(k)=nu_prev*xnurat
         if(nu(k).gt.nu_th(i_sw_fr)*eVhp09 .and.
     &        i_sw_fr .le. N_nu_th)then
            print *, nu_th(i_sw_fr)
            nu(k)=nu_th(i_sw_fr)*eV/h_p
            k_th(i_sw_fr)=k
            xnurat=xnurat2
            i_sw_fr=i_sw_fr+1
         endif
         if(i_sw_fr>1)then
            if(k>k_th(i_sw_fr-1)+k_fine) xnurat=xnurat1
            if(i_sw_fr.eq.4) xnurat=10.
         endif
         nu_prev=nu(k)
         Planck(k) = nu(k)**3/(dexp(h_p*nu(k)/xk_B/T_rad)-1.0d0)
         Flux(k) = (h_p*nu(k)/eV)**(-1.5d0)
         Flux(k) = Planck(k)
         
         call cross_section_H_He(h_p*nu(k)/eV,sigma_H_nu(k),
     &        sigma_He_nu(k),sigma_Hep_nu(k),sigma_Hm_nu(k),
     &        sigma_H2p_nu(k))
         print*, k, h_p*nu(k)/eV,xnurat ! wli
      enddo

! wli add, integration of pd rate.
      xk_pd_Hm = 0.d0
      xk_pd_H2p = 0.d0

      do k=1,NF-1
         dnu(k)=nu(k+1)-nu(k)
         if (nu(k).lt.13.6d0*eV/h_p) then
            xk_pd_Hm = xk_pd_Hm 
     &         + sigma_Hm_nu(k)*4*pi*Planck(k) 
     &         /h_p/nu(k)*dnu(k)
            xk_pd_H2p = xk_pd_H2p
     &         + sigma_H2p_nu(k)*4*pi*Planck(k) 
     &         /h_p/nu(k)*dnu(k)
         endif
         write(10,*) k, h_p*nu(k)/eV, Planck(k),
     &        sigma_H_nu(k), sigma_He_nu(k), sigma_Hep_nu(k),
     &        sigma_Hm_nu(k), sigma_H2p_nu(k),
     &        xk_pd_Hm, xk_pd_H2p

      enddo
      
      print*,''
      print*, NF,'            NF: frequency grid'
      print*, 'nu_min(eV, Hz)',h_p*nu(1) /eV,nu(1)
      print*, 'nu_max(eV, Hz)' , h_p*nu(NF)/eV,nu(NF)
      print*, ''
      print*, xk_pd_Hm/((12.646-12.4)*Planck(34)+
     &     (12.4-11.496)*Planck(35))*(12.646-11.496)*1.d-21,
     & xk_pd_H2p/((12.646-12.4)*Planck(34)+
     &     (12.4-11.496)*Planck(35))*(12.646-11.496)*1.d-21,
     & xk_pd_Hm/((12.646-12.4)*Planck(34)+
     &     (12.4-11.496)*Planck(35))*(12.646-11.496)*1.d-21
     &     /1.39d-12


      ! if(i_sw_frset.eq.0)then
      !   i_sw_frset=1
      !   i_sw_fr=0
      !   nu_prev=1e13/1.05
      !   do k=1,NF
      !      nu(k)=nu_prev*1.05
      !      if(i_sw_fr.eq.0 .and. nu(k).gt.3.2865e15)then
      !         nu(k)=13.6*eV/h_p
      !         k_Hion=k
      !         i_sw_fr=1
      !      else if(i_sw_fr.eq.1 .and. nu(k).gt.5.9413e15)then
      !         nu(k)=24.586*eV/h_p
      !         k_Heion=k
      !         i_sw_fr=2
      !      else if(i_sw_fr.eq.2 .and. nu(k).gt.1.315e16)then
      !         nu(k)=54.42*eV/h_p
      !         k_Hepion=k
      !         i_sw_fr=3
      !      endif
      !      nu_prev=nu(k)
      !      print*, k,nu(k)
      !   enddo
      ! print*,''
      ! print*, NF,'            NF: frequency grid'
      ! print*, h_p*nu(1) /eV,'nu_min(eV)'
      ! print*, h_p*nu(NF)/eV,'nu_max(eV)'
      ! print*, ''
      
      END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE cross_section_H_He(E,sigma_H,sigma_He,sigma_Hep,
     &     sigma_Hm,sigma_H2p)

      implicit none
      real*8 :: E, sigma_H, sigma_He, sigma_Hep, sigma_Hm, sigma_H2p
      real*8 :: E_H,x,y_w,y_0,y_1,y,P,y_a
      real*8 :: E_T,a1,s,a2,a3
      real*8 :: rlmbd_0,rlmbd, C1, C2, C3, C4, C5, C6, f
      real*8 :: T_H2, H2p_bf_CrossSec

!     H cross section
      E_H=13.6d0
      if(E<E_H) then
         sigma_H=0.0d0
      else
         y_0=0.0d0
         x=(E/4.298d-1)-y_0
         y_1=0.0d0
         y=dsqrt(x*x+y_1*y_1)
         y_w=0.0d0
         P=2.963d0
         y_a=3.288d1

         sigma_H = 5.475d4*1.0d-18
     &        *((x-1.0d0)**2.0d0 +y_w**2.0d0)
     &        *(y**(P/2.0d0-5.5d0))
     &        /((1.0d0+dsqrt(y/y_a))**P)
      endif

!     He cross section
!     Yan 1998
      E_T=24.586d0
      if(E<E_T) then
         sigma_He=0.0d0
      else
         a1=7.3861d0
         s=3.9119d0
         a2=-3.2491d0
         a3=1.1783d0
         x=E/24.58d0

         sigma_He=7.4d6*(a1/x**s + (1.0d0-a1)/x**(s+1.0d0))
     &        +733.0d0/((E/1.0d3)**3.5d0)
     &        *(1.0d0+a2/(x**0.5d0)*dexp(-a3/dsqrt(x)))
         sigma_He=sigma_He*1.0d-24
      endif

!     He+ cross section
!     VF 1996
      E_T = 54.42d0
      if(E<E_T) then
         sigma_Hep=0.0d0
      else
         y_0=0.0d0
         x=(E/1.72d0)-y_0
         y_1=0.0d0
         y=dsqrt(x*x+y_1*y_1)
         y_w=0.0d0
         P=2.963d0
         y_a=3.288d1

         sigma_Hep=1.369d4*1.0d-18
     &        *((x-1.0d0)**2.0d0 +y_w**2.0d0)
     &        *(y**(P/2.0d0-5.5d0))
     &        /((1.0d0+dsqrt(y/y_a))**P)
      endif

!     H- photodetachment
!     H-   +  ph.   ->   H    +    e
!     photon energy E in eV
!     ref. T.L.John (1988 A&A 193,189)
      rlmbd_0=1.6419d0
      C1=152.519d0
      C2=49.534d0
      C3=-118.858d0
      C4=92.536d0
      C5=-34.194d0
      C6=4.982d0
      rlmbd=1.23985d0/E
      if(rlmbd.le.rlmbd_0) then
         x=(1.d0/rlmbd-1.d0/rlmbd_0)**0.5d0
         f=C1+x*(C2+x*(C3+x*(C4+x*(C5+x*C6))))
         sigma_Hm=1.d-18*(rlmbd**3)*(x**3)*f
      else
         sigma_Hm=0.d0
      endif
      if(rlmbd.le.0.09d0) then
         sigma_Hm=0.0d0
      endif

!     H2+ photodetachment
      T_H2 = 8.0d3
      sigma_H2p = H2p_bf_CrossSec(E,T_H2)

      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION H2p_bf_CrossSec(rnu,T_K)
      implicit real*8(a-h,o-z)
C     cross section of H2+ b-f
C     H2+   +  ph.   ->   H   +   H+
C     frequency rnu in eV
C     ref. P.C.Stancil (1994 ApJ 430,360)
      common /op_H2p_bf/Ta(10),rnua(45),sigmaa(10,45),
     &     rlTaa(148),rlaa(148)
      dimension sa(10)

      if(T_K.le.Ta(2)) then
         rlT=dlog(T_K)
         if(rlT.gt.rlTaa(1)) then
            call linear(rlTaa,rlaa,148,rlT,rla)
            a=dexp(rla)
         else
            a=0.d0
         endif
         if(rnu.lt.rnua(1)) then
            s1=sigmaa(2,1)
            s2=sigmaa(2,2)
            rlnu=dlog(rnu)
            rlnu1=dlog(rnua(1))
            rlnu2=dlog(rnua(2))
            rls1=dlog(s1)
            rls2=dlog(s2)
            rls=((rlnu-rlnu1)/(rlnu2-rlnu1))*rls2
     &        +((rlnu2-rlnu)/(rlnu2-rlnu1))*rls1
            H2p_bf_CrossSec=a*dexp(rls)
         elseif(rnu.gt.rnua(45)) then
            H2p_bf_CrossSec=0.d0
         else
            call bilinear(Ta,rnua,sigmaa,10,45,Ta(1),rnu,s1)
            call bilinear(Ta,rnua,sigmaa,10,45,Ta(2),rnu,s2)
            H2p_bf_CrossSec=(1.d0-a)*s1+a*s2
         endif
      else
         T_KK=min(T_K,Ta(10)) !T_KK = 8000K 
         if(rnu.le.rnua(1)) then
            do i=1,10
               sa(i)=sigmaa(i,1)
            enddo
            call linear(Ta,sa,10,T_KK,s1)
            do i=1,10
               sa(i)=sigmaa(i,2)
            enddo
            call linear(Ta,sa,10,T_KK,s2)
            rlnu=dlog(rnu)
            rlnu1=dlog(rnua(1))
            rlnu2=dlog(rnua(2))
            rls1=dlog(s1)
            rls2=dlog(s2)
            rls=((rlnu-rlnu1)/(rlnu2-rlnu1))*rls2
     &           +((rlnu2-rlnu)/(rlnu2-rlnu1))*rls1
            H2p_bf_CrossSec=dexp(rls)
         elseif(rnu.gt.rnua(45)) then
            H2p_bf_CrossSec=0.d0
         else
            call bilinear(Ta,rnua,sigmaa,10,45,T_KK,rnu,sigma)
            H2p_bf_CrossSec=sigma
         endif
      endif
      
 11   continue
      return
      END


      SUBROUTINE linear(xa,ya,m,x,y)
      implicit real*8(a-h,o-z)
      dimension xa(m),ya(m)
      if((x.lt.xa(1)).or.(x.gt.xa(m))) then
         print *,xa(1),x,xa(m)
      endif
      do 11 i=1,m
         if(x-xa(i).le.0.d0) then
            ms=i
            go to 12
         endif
 11   continue
 12   continue
      if(ms.eq.1) ms=2
      y1=ya(ms-1)
      y2=ya(ms)
      t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
      y=(1.d0-t)*y1+t*y2
      return
      END


      SUBROUTINE bilinear(x1a,x2a,ya,m,n,x1,x2,y)
      implicit real*8(a-h,o-z)
      dimension x1a(m),x2a(n),ya(m,n)
      if((x1.lt.x1a(1)).or.(x1.gt.x1a(m))
     &     .or.(x2.lt.x2a(1)).or.(x2.gt.x2a(n))) then
         stop
      endif
      do 11 i=1,m
         if(x1-x1a(i).le.0.d0) then
            ms=i
            go to 12
         endif
 11   continue
 12   continue
      do 13 i=1,n
         if(x2-x2a(i).le.0.d0) then
            ns=i
            go to 14
         endif
 13   continue
 14   continue
      if(ms.eq.1) ms=2
      if(ns.eq.1) ns=2
      y1=ya(ms-1,ns-1)
      y2=ya(ms,ns-1)
      y3=ya(ms,ns)
      y4=ya(ms-1,ns)
      t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
      u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
      y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
      return
      END


      SUBROUTINE H2p_bf_data
      implicit real*8(a-h,o-z)
      common /op_H2p_bf/Ta(10),rnua(45),sigmaa(10,45),
     &     rlTaa(148),rlaa(148)
      Ta(1)=0.d0
      Ta(2)=2.52d3
      Ta(3)=3.15d3
      Ta(4)=4.2d3
      Ta(5)=5.04d3
      Ta(6)=6.3d3
      Ta(7)=8.4d3
      Ta(8)=1.26d4
      Ta(9)=1.68d4
      Ta(10)=2.52d4

      open(10,file='./H2p/sigma_0K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(1,i)
      !   print*, i, rnua(i), sigmaa(1,i) ! wli 
      enddo
      close(10)

      open(10,file='./H2p/sigma_2520K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(2,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_3150K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(3,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_4200K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(4,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_5040K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(5,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_6300K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(6,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_8400K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(7,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_12600K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(8,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_16800K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(9,i)
      enddo
      close(10)

      open(10,file='./H2p/sigma_25200K',status='old')
      do i=1,45
         read(10,*) rnua(i),sigmaa(10,i)
      enddo
      close(10)

      open(10,file='./H2p/a',status='old')
      do i=1,148
         read(10,*) Taa,aa
         rlTaa(i)=dlog(Taa)
         rlaa(i)=dlog(aa)
      enddo
      close(10)

      return
      END

