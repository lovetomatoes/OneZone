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

