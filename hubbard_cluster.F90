!--------------------------------------------------
!
! hubbard_cluster.F90
! module: hubbard_mod
! requirements: hop_mod, bas_mod, matrix_mod, underwood_mod
!
! created by Kun Fang
!
! This module contains definations and interface for Hubbard
! model creations and calculations. The Hubbard model is 
! represented by combination of basis, hopping matrix and 
! its Hamiltonian. The subroutines for solving eigenvalues
! and eigenvectors are also included.
!
! The module provides interfaces for extracting the whole 
! matrix of the Hamiltonian or only a part of the Hamiltonian
! for specific number of electrons and spins. The interfaces
! of solving the Hamiltonian using Lanczos method is also
! included. There is an option whether to calculate a given
! number of eigenstates or a to find eigenstates within a
! given energy range. It is default to choose the first choice.
!
! types:
! hubbard_cluster_type
!  |- hop
!  |- basis_type
!
!------------------------------------------------------

module hubbard_mod
  use hop_mod
  use bas_mod
  use matrix_mod
  use underwood_mod
  implicit none

  ! type for the Hubbard model
  ! type(hubbard_cluster_type)
  type hubbard_cluster_type
    integer::nsite=0,nbas=0
    type(hop),pointer::hop=>NULL()
    type(basis_type),pointer::bas=>NULL()
  end type hubbard_cluster_type

  type(complex_matrix),private,pointer::comm_matrix=>NULL() ! a temperary matrix for calculations
  logical,public::checkEnergy=.false. ! option whether to calcuate within an energy range
  real(8),private,parameter::ediff=3.5d0
  real(8),private,parameter::epsilon=1.d-5
  

contains

!-----------------hubbard cluster------------------------

  ! initialize the Hubbard model
  !
  ! input:
  ! cl - type(hop): contains cluster coordinations and 
  !                 the hopping matrix
  ! b - type(basis_type): contains all the basis of the
  !                       cluster
  !
  ! output:
  ! hc - type(hubbard_cluster_type): initialized Hubbard
  !                                  model 
  function hubbard_init(cl,b) result(hc)
    type(hubbard_cluster_type),pointer::hc
    type(hop),pointer::cl
    type(basis_type),pointer::b
    allocate(hc)
    hc%bas=>b
    hc%hop=>cl
    hc%nsite=hc%hop%site
    hc%nbas=basis_get_n_basis(hc%bas)
  end function
  
  subroutine hubbard_clean(hc)
    type(hubbard_cluster_type),pointer::hc
    call basis_clean(hc%bas)
    call cluster_delete(hc%hop)
    deallocate(hc);
  end subroutine
  
  ! extract a partition of the Hamiltonian matrix
  !
  ! input:
  ! cl - type(hop): cluster coordinations and the hopping 
  !                 matrix
  ! bas - type(basis_type): basis of the cluster
  ! si, sf - integer: start and end of the matrix partition
  !
  ! output:
  ! M - type(complex_matrix): complex matrix of the partition of
  !                           the Hamiltonian matrix
  function cluster_Hamiltonian_Partition(cl,bas,si,sf) result(M)
    type(hop),pointer::cl
    type(basis_type),pointer::bas
    integer::si,sf
    type(complex_matrix),pointer::M
    integer::i,j,r,s,n,spin,spin2,site,k1,k2,alpha,beta,mr,ms
    complex(8),pointer::t(:,:)
    complex(8)::x
    n=sf-si+1
    M=>complex_matrix_init(n)
    site=cl%site
    t=>cluster_update(cl)
    do r=1,n
      mr=r+si-1
      x=complex_matrix_show(M,r,r)+basis_double_occupy(mr,bas)*cl%U*(-1)**(cl%nambu)
      call complex_matrix_set(M,r,r,x)
      do spin=1,2
        do alpha=1,site
          i=basis_c(bas,alpha,spin,mr)
          if(i==0) cycle
          k1=i/abs(i)
          i=abs(i)
          do spin2=1,2
            do beta=1,site
              if(t(alpha+(spin-1)*site,beta+(spin2-1)*site)==0) cycle
              ms=basis_cplus(bas,beta,spin2,i)
              if(ms==0) cycle
              k2=ms/abs(ms)
              ms=abs(ms)
              s=ms-si+1
              if(r>=s) then
                x=complex_matrix_show(M,r,s)+k1*k2*t(alpha+(spin-1)*site,beta+(spin2-1)*site)
                call complex_matrix_set(M,r,s,x)
                call complex_matrix_set(M,s,r,conjg(x))
              end if
            end do
          end do
        end do
      end do
    end do
    deallocate(t)
  end function

  ! extract the full Hamiltonian matrix
  !
  ! input:
  ! cl - type(hop): cluster coordinations and the hopping 
  !                 matrix
  ! bas - type(basis_type): basis of the cluster
  !
  ! output:
  ! M - type(complex_matrix): complex matrix of the 
  !                           Hamiltonian matrix
  function cluster_Hamiltonian(cl,bas) result(M)
    type(hop),pointer::cl
    type(basis_type),pointer::bas
    type(complex_matrix),pointer::M
    integer::i,j,r,s,n,spin,spin2,site,k1,k2,alpha,beta
    complex(8),pointer::t(:,:)
    complex(8)::x
    n=basis_get_n_basis(bas)
    M=>complex_matrix_init(n)
    site=cl%site
    t=>cluster_update(cl)
    do r=1,n
      x=complex_matrix_show(M,r,r)+basis_double_occupy(r,bas)*cl%U*(-1)**(cl%nambu)
      call complex_matrix_set(M,r,r,x)
      do spin=1,2
        do alpha=1,site
          i=basis_c(bas,alpha,spin,r)
          if(i==0) cycle
          k1=i/abs(i)
          i=abs(i)
          do spin2=1,2
            do beta=1,site
              if(t(alpha+(spin-1)*site,beta+(spin2-1)*site)==0) cycle
              s=basis_cplus(bas,beta,spin2,i)
              if(s==0) cycle
              k2=s/abs(s)
              s=abs(s)
              if(r>=s) then
                x=complex_matrix_show(M,r,s)+k1*k2*t(alpha+(spin-1)*site,beta+(spin2-1)*site)
                call complex_matrix_set(M,r,s,x)
                call complex_matrix_set(M,s,r,conjg(x))
              end if
            end do
          end do
        end do
      end do
    end do
    deallocate(t)
  end function

  ! connect the input matrix to the matrix solver
  ! This is a must before solve the matrix using Lanczos methed
  subroutine cluster_connect_matrix(M)
    type(complex_matrix),pointer::M
    comm_matrix=>M
  end subroutine

  ! disconnect the input matrix from the matrix solver
  ! For most of time it is ok if you don't disconnect it. For
  ! safety, it is better to disconnect it.
  subroutine cluster_disconnect_matrix()
    comm_matrix=>NULL()
  end subroutine

  ! A cluster solver for the whole Hamiltonian using Lanscoz
  ! method.
  !
  ! input:
  ! cluster - type(hubbard_cluster_type): the Hamiltonian of the 
  !                                       cluster
  ! m - integer: the maximum number of eigenstates you want to 
  !              calculate. If the energy range is needed, this
  !              number is not important
  ! D - real(8): a vector pointer (not initialized)
  ! X - real(8): a 2D matrix (not initialized)
  !
  ! output:
  ! m - integer: number of eigenstates for output. If the energy
  !              range is needed, this number may be different 
  !              from input
  ! D - real(8): a vector with length m, containing all the 
  !              eigenvalues for output
  ! X - real(8): an n*m matrix, containing all the eigenvectors 
  !              for output corresponding to these eigenvalues
  subroutine cluster_solver(cluster,m,D,X)
    implicit none
    type(hubbard_cluster_type),pointer,intent(in)::cluster
    real(8),pointer::D(:)
    complex(8),pointer::X(:,:)
    type(complex_matrix),pointer::H
    integer::n,m,orbit,ne,IECODE,q,pinit,iter,i
    
    n=cluster%nbas
    orbit=cluster%nsite
    H=>cluster_Hamiltonian(cluster%hop,cluster%bas)
    call cluster_connect_matrix(H)
    do
      ne=0
      q=m*5
      if(q>n) q=n
      pinit=m

      allocate(D(q))
      allocate(X(n,q))
      iter=n
      do while(ne<m)
        CALL MINVAL(n,q,pinit,m,iter,epsilon,OP,ne,D,X,IECODE)
        !print *,ne
      end do
      if(.not.checkEnergy) exit
      if(D(ne)-D(1)>ediff) exit
      deallocate(D)
      deallocate(X)
      m=m+20
    end do
    call cluster_disconnect_matrix()
    call complex_matrix_del(H)
  end subroutine

  ! A cluster solver for the Hamiltonian with fixed number of 
  ! electrons using Lanscoz method.
  !
  ! input:
  ! cluster - type(hubbard_cluster_type): the Hamiltonian of the 
  !                                       cluster
  ! elec -integer: number of electrons for the calculation
  ! m - integer: the maximum number of eigenstates you want to 
  !              calculate. If the energy range is needed, this
  !              number is not important
  ! D - real(8): a vector pointer (not initialized)
  ! X - real(8): a 2D matrix (not initialized)
  ! dim - integer: none
  !
  ! output:
  ! m - integer: number of eigenstates for output. If the energy
  !              range is needed, this number may be different 
  !              from input
  ! D - real(8): a vector with length m, containing all the 
  !              eigenvalues for output
  ! X - real(8): an n*m matrix, containing all the eigenvectors 
  !              for output corresponding to these eigenvalues
  ! dim - integer: number of column and row dimensions of the 
  !                Hamiltonian matrix
  !
  ! Before using this subroutine, make sure the electron number
  ! is a good quantum number.
  subroutine cluster_solver_elec(cluster,elec,m,D,X,dim)
    implicit none
    type(hubbard_cluster_type),pointer,intent(in)::cluster
    real(8),pointer::D(:)
    complex(8),pointer::X(:,:)
    type(basis_type),pointer::bas
    type(complex_matrix),pointer::H
    integer::dim,elec,n,m,orbit,ne,IECODE,q,pinit,iter,i,si,sf
    n=cluster%nbas
    orbit=cluster%nsite
    bas=>cluster%bas
    si=n
    sf=0
    do i=1,bas%n_block
      if(bas%block(i)%elec==elec) then
        if(si>bas%block(i)%init) si=bas%block(i)%init
        if(sf<bas%block(i)%fin) sf=bas%block(i)%fin
      end if
    end do
    H=>cluster_Hamiltonian_Partition(cluster%hop,cluster%bas,si,sf)
    !call complex_matrix_print(H)
    call cluster_connect_matrix(H)
    n=H.dim
    dim=n
    print *,"dimension:",n
    do
      ne=0
      q=m*5
      if(q>n) q=n
      pinit=m
      allocate(D(q))
      allocate(X(n,q))
      iter=n
      do while(ne<m)
        CALL MINVAL(n,q,pinit,m,iter,epsilon,OP,ne,D,X,IECODE)
        !print *,ne
      end do
      if(.not.checkEnergy) exit
      if(D(ne)-D(1)>ediff) exit
      deallocate(D)
      deallocate(X)
      m=m+20
    end do
    call cluster_disconnect_matrix()
    call complex_matrix_del(H)
  end subroutine

  ! A cluster solver for the Hamiltonian with fixed number of 
  ! electrons and spins using Lanscoz method.
  !
  ! input:
  ! cluster - type(hubbard_cluster_type): the Hamiltonian of the 
  !                                       cluster
  ! elec - integer: number of electrons for the calculation
  ! spin - integer: number of spins for the calculation 
  !                 spin = spin up - spin down
  ! m - integer: the maximum number of eigenstates you want to 
  !              calculate. If the energy range is needed, this
  !              number is not important
  ! D - real(8): a vector pointer (not initialized)
  ! X - real(8): a 2D matrix (not initialized)
  ! dim - integer: none
  !
  ! output:
  ! m - integer: number of eigenstates for output. If the energy
  !              range is needed, this number may be different 
  !              from input
  ! D - real(8): a vector with length m, containing all the 
  !              eigenvalues for output
  ! X - real(8): an n*m matrix, containing all the eigenvectors 
  !              for output corresponding to these eigenvalues
  ! dim - integer: number of column and row dimensions of the 
  !                Hamiltonian matrix
  !
  ! Before using this subroutine, make sure both the electron 
  ! number and spin number are good quantum numbers.
  subroutine cluster_solver_elec_spin(cluster,elec,spin,m,D,X,dim)
    implicit none
    type(hubbard_cluster_type),pointer,intent(in)::cluster
    real(8),pointer::D(:)
    complex(8),pointer::X(:,:)
    type(basis_type),pointer::bas
    type(complex_matrix),pointer::H
    integer::dim,elec,spin,n,m,orbit,ne,IECODE,q,pinit,iter,i,si,sf
    n=cluster%nbas
    orbit=cluster%nsite
    bas=>cluster%bas
    do i=1,bas%n_block
      if(bas%block(i)%elec==elec.and.bas%block(i)%spin==spin) then
        si=bas%block(i)%init
        sf=bas%block(i)%fin
        exit
      end if
    end do
    H=>cluster_Hamiltonian_Partition(cluster%hop,cluster%bas,si,sf)
    call cluster_connect_matrix(H)
    n=H.dim
    dim=n
    do
      ne=0
      q=m*5
      if(q>n) q=n
      pinit=m
      allocate(D(q))
      allocate(X(n,q))
      iter=n
      do while(ne<m)
        CALL MINVAL(n,q,pinit,m,iter,epsilon,OP,ne,D,X,IECODE)
        !print *,ne
      end do
      if(.not.checkEnergy) exit
      if(D(ne)-D(1)>ediff) exit
      deallocate(D)
      deallocate(X)
      m=m+20
    end do
    call cluster_disconnect_matrix()
    call complex_matrix_del(H)
  end subroutine

  ! calculation matrix product of a vector
  subroutine OP(n,x,y)
    integer::n
    complex(8)::x(n),y(n)
    call complex_matrix_product(comm_matrix,x,y)
  end subroutine

  ! return number of sites in the cluster1
  function cluster_get_n_orbit(cluster) result(n)
    implicit none
    type(hubbard_cluster_type),pointer,intent(in)::cluster
    integer::n
    n=cluster%nsite
  end function
  
  
end module hubbard_mod
