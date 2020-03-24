PROGRAM analysis

implicit none

integer, parameter :: ntree = 10
integer, parameter :: N = 30000000, NN =30000, M = 100, P=1000
integer :: id_node_a(N), id_gen_a(N), Nchild_a(N),
&     id_child_a(N), id_sibling_a(N), id_parent_a(N)
integer :: id_node(NN), id_gen(NN), Nchild(NN), i, j, jmax, 
&     id_child(NN), id_sibling(NN), id_parent(NN), k, kmin, kmax,
&     N_tree_tot, N_tree(P)
real*8 :: mhalo_a(N), z_a(N), mhalo(NN), z(NN)

!     output (M should be the number of z bin)
real*8 :: Omega_m, Omega_L, H0, d, Omega_mz, z1, Delta_crit, pi
integer :: id_halo_m(M), a,b,c
real*8 :: mhalo_m(M), z_m(M), T_vir(M), mu, Mhalo_obs, V_circ_obs,
&     mhalo_ave(M), T_vir_ave(M), mhalo_var(M), T_vir_var(M)

pi = 3.141592d0
Omega_m = 0.3089d0
Omega_L = 1.0d0 - Omega_m
H0      = 0.6774d0
mu = 0.6d0

open(10, file='Ntrees.txt')
read (10,*) (N_tree(k), k=1,ntree )
close(10)
N_tree_tot = 0
do k = 1,ntree 
   N_tree_tot = N_tree_tot + N_tree(k)
enddo
print*,"N_tree_tot = " ,N_tree_tot

!     initialize
do j=1,N
   id_node_a(j)   = 0.0d0
   id_gen_a(j)    = 0.0d0
   mhalo_a(j)     = 0.0d0
   z_a(j)         = 0.0d0
   Nchild_a(j)    = 0.0d0
   id_child_a(j)  = 0.0d0
   id_sibling_a(j)= 0.0d0
   id_parent_a(j) = 0.0d0
enddo
do j=1,M
   mhalo_ave(j) = 0.0d0
   T_vir_ave(j) = 0.0d0
   mhalo_var(j) = 0.0d0
   T_vir_var(j) = 0.0d0
enddo
!      open(10, file='average.txt')
!      read (10,*) (z_m(k), mhalo_ave(k), T_vir_ave(k), k=1,99)
!      close(10)  

open(10, file='trees.txt')
read (10,*) (id_node_a(j), id_gen_a(j), mhalo_a(j), z_a(j), 
&     Nchild_a(j), id_child_a(j), id_sibling_a(j), id_parent_a(j),
&     j = 1, N_tree_tot)
close(10)

!     check read
c      do j = 1,N_tree_tot
c         print*, id_node_a(j), id_gen_a(j), mhalo_a(j), z_a(j), 
c     &        Nchild_a(j), id_child_a(j), id_sibling_a(j), 
c     &        id_parent_a(j)
c      enddo

N_tree_tot = 0.0d0
do i = 1,ntree 
c         print*, 'ntree = ',i
   do j=1,N_tree(i)
      id_node(j)    = id_node_a(j+N_tree_tot)
      id_gen(j)     = id_gen_a(j+N_tree_tot)
      mhalo(j)      = mhalo_a(j+N_tree_tot)
      z(j)          = z_a(j+N_tree_tot)
      Nchild(j)     = Nchild_a(j+N_tree_tot)
      id_child(j)   = id_child_a(j+N_tree_tot)
      id_sibling(j) = id_sibling_a(j+N_tree_tot)
      id_parent(j)  = id_parent_a(j+N_tree_tot)
c            print*, i, id_node(j), id_gen(j), mhalo(j), z(j),
c     &           Nchild(j), id_child(j), id_sibling(j),
c     &           id_parent(j)  
   enddo
   N_tree_tot = N_tree_tot + N_tree(i)
!--------------------------------------------------
!     a = ID of most massive ones in a generation
!     b = total number of progenitors 
   a = 1     
   id_halo_m(1) = id_node(1)
   mhalo_m(1) = mhalo(1)
   z_m(1) = z(1)
c         do j=1,1
c            print*, j,id_halo_m(j), mhalo_m(j)/1.0d11, z_m(j) 
c         enddo
   do j=2,M
      if(a==1) then
         b = Nchild(a)
         a = a + 1
         id_halo_m(j) = a
         mhalo_m(j) = mhalo(a)
         z_m(j) = z(a)
!     print*, a,b
      else
         kmin = a
         kmax = b + a - 1
         a = a + b
         b = 0
         do k = kmin, kmax
            b = b + Nchild(k)
         enddo
         id_halo_m(j) = a
         mhalo_m(j) = mhalo(a)
         z_m(j) = z(a)
c      print*, a,b
      endif
      if(b==0) exit
c            print*, j,id_halo_m(j), mhalo_m(j)/1.0d11, z_m(j)
      jmax = j
   enddo
   
   do j=1,jmax-1
      z1 = 1.0d0 + z_m(j)
      Omega_mz = Omega_m*z1**3.0d0/(Omega_m*z1**3.0d0 + Omega_L)
      d = Omega_mz-1.0d0
      Delta_crit = 18.0d0*pi*pi + 82.0d0*d - 39.0d0*d*d
      Delta_crit = Delta_crit/(18.0d0*pi*pi)

      T_vir(j) = 1.98d4 * (mu/0.6d0)
&           * (mhalo_m(j)/1.0d8*H0)**(2.0d0/3.0d0) * (z1/10.0d0)
&           * (Omega_m/Omega_mz*Delta_crit)**(1.0d0/3.0d0)

      mhalo_ave(j) = mhalo_ave(j) + mhalo_m(j)
      T_vir_ave(j) = T_vir_ave(j) + T_vir(j)
      mhalo_var(j) = mhalo_var(j) 
&           + (mhalo_m(j)-mhalo_ave(j))**2.0d0
      T_vir_var(j) = T_vir_var(j)
&           + (T_vir(j)-T_vir_ave(j))**2.0d0
      if (j==1) write(11,*) "j ", "id_m ", "mhalo_m/1.e11Ms "
&           ,"z_m ","T_{vir,m} "
      write(11,*) j, id_halo_m(j), mhalo_m(j)/1.0d11, z_m(j), 
&           T_vir(j)
   enddo
   write(11,*) ' '
enddo ! i
! write average.txt file acturally
do j=1,jmax-1
   if (j==1) write(12,*) "z_m ","mhalo_{ave} ", "T_{vir,ave} "
   write(12,*) z_m(j),mhalo_ave(j)/ntree, T_vir_ave(j)/ntree
enddo
c      do j=1,jmax-1
c         write(13,*) z_m(j), mhalo_ave(j), T_vir_ave(j),
c     &        dsqrt(mhalo_var(j)/ntree), dsqrt(T_vir_var(j)/ntree)
c      enddo

!     J1208−0200
!     V_circ =143 km/s
!     z = 6.1165

V_circ_obs = 143.0d0
z1 = 1.0d0 + 6.1165d0
Omega_mz = Omega_m*z1**3.0d0/(Omega_m*z1**3.0d0 + Omega_L)
d = Omega_mz-1.0d0
Delta_crit = 18.0d0*pi*pi + 82.0d0*d - 39.0d0*d*d
Delta_crit = Delta_crit/(18.0d0*pi*pi)

Mhalo_obs = 1.0d8/H0 *(V_circ_obs/23.4d0)**3.0d0
&     / (Omega_m/Omega_mz*Delta_crit)**0.5d0 
&     / (z1/10.0d0)**1.5d0
c      print*, Mhalo_obs/1.0d11, Omega_mz, Delta_crit

END
