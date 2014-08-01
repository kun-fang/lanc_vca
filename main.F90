program hubbard_cluster
  use hubbard_mod
  use rand_mod
  use timer_mod
  use VCA

  implicit none

  type(basis_type),pointer::bas
  type(hop),pointer::cl
  type(hubbard_cluster_type),pointer::cluster
  type(Q_type),pointer::Q

  integer::n,nsite,i,ne,elec,spin
  real(8)::omega,density
  real(8),pointer,dimension(:)::e
  complex(8),pointer,dimension(:,:)::X

  call timer_init()
  cl=>cluster_init()
  nsite=cl%site
  bas=>basis_init(nsite)
  cluster=>hubbard_init(cl,bas)

  write(*,*)"----------------------------------------------------------------"
  write(*,*)"VCA calculation: initial parameters are"
  write(*,*)"U=",cluster%hop%U,"mu=",cluster%hop%lmu
  write(*,*)"t=",cluster%hop%lt
  write(*,*)"Cluster t=",cluster%hop%ct
  write(*,*)"Cluster mu=",cluster%hop%cmu
  write(*,*)"----------------------------------------------------------------"

!------------------ Cluster Eigenvalues and Eigenstates --------------------

  ne=10        ! # of eigenstates you want to calculate
  elec=4      ! # of electrons used in the calculation
  spin=0      ! # of spins used in the calcuation (# of spin up - # of spin down)
  n=cluster%nbas    ! don't change

  ! calculate eigens for all the possible # of electrons and spins
  !call cluster_solver(cluster,ne,e,X)

  ! calculate eigens for specific # of electrons
  call cluster_solver_elec(cluster,elec,ne,e,X,n)

  ! calculate eigens for specific # of electrons and spins
  !call cluster_solver_elec_spin(cluster,elec,spin,ne,e,X,n)

  ! print
  write(*,*) 1,e(1)
  do i=2,ne
    write(*,*) i,e(i),e(i)-e(i-1)
  end do

!----------------------------- VCA calculations ------------------------------
!
! Several setups may be needed to make the calculation smooth
!
!-----------------------------------------------------------------------------

! if(VCA_connect_cluster(cluster)) then       ! connect cluster to VCA solver

!   print *,VCA_potthoff_functional()         ! VCA energy at given parameters
!   print *,VCA_particle_density()            ! electron density at given parameters
!   call VCA_spectral_function()              ! spectral function at given parameters
!   call VCA_DOS(-3.d0,3.d0,40)               ! DOS at given parameters
!   call VCA_fermi_surface()                  ! fermi surface at given parameters
!   if(VCA_optimal(omega)) print *,cluster%hop%cmu,cluster%hop%ct,omega
                                              ! VCA optimization for the initial parameters

!   call VCA_disconnect_cluster()             ! disconnect cluster from VCA solver
! end if


  call hubbard_clean(cluster)
  print *,cputime(),"seconds"


  PAUSE
end program hubbard_cluster


