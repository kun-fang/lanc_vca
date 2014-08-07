!--------------------------------------------------
!
! vca.F90
! module: VCA
! requirement: hubbard_mod,lehmann,qnewton
!
! created by Kun Fang
!
! This module provides several of user interfaces for the VCA 
! calculations, including VCA variational calculation, electron
! density (concentration), spectral function, DOS, fermi surface.
! These interfaces are tested and relatively reliable. I have 
! deleted some old subroutines, because they are written for a 
! single task and not quite reliable.
!
!------------------------------------------------------

module VCA
  use hubbard_mod
  use lehmann
  use qnewton
  complex(8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
  real(8),private,parameter::eps=0.001
  integer,private,parameter::nk=20    ! the number of k-space points in one direction
  integer,private::initm=150    ! the number of eigenstates that will be calculated
  
  type(hubbard_cluster_type),pointer,private::cluster=>NULL()

  contains

  ! connect a cluster object to this module
  ! this is a must before you do any VCA calculation
  ! 
  ! input:
  ! p - type(hubbard_cluster_type): the cluster object
  !                    you want to make connections
  !
  ! output:
  ! OK - logical: whether the connection is successful
  function VCA_connect_cluster(p) result(OK)
    implicit none
    type(hubbard_cluster_type),pointer::p
    logical::OK
    OK=.false.
    if(.not.associated(p)) return
    cluster=>p
    OK=.true.
  end function
  
  ! disconnect the cluster
  subroutine VCA_disconnect_cluster()
    implicit none
    cluster=>NULL()
  end subroutine

  ! variational process of VCA calculation
  ! The variational calculation is performed using quasi-newton
  ! method to find out stationary points in a self-energy space.
  ! The way to connect this subroutine with quasi-newton 
  ! subroutine is really ugly because of the restrictions of 
  ! FORTRAN language. The connection requires a user-defined
  ! type PointArray which is used to create an array of pointers
  ! so that the quasi-newton subroutine can directly change 
  ! values of variational parameters. in subroutine varinit, you
  ! have to assign variational parameters to the pointer array
  ! before beginning of quasi-newton calculation.
  ! 
  ! input:
  ! omega - real(8): the value of this input is not important
  !
  ! output:
  ! omega - real(8): the grand potential of the system obtained
  !                  from variational calculation
  ! OK - logical: indication whether the variational calculation
  !               is successful
  function VCA_optimal(omega) result(OK)
    implicit none
    integer::n_var,i,ie
    real(8)::omega
    real(8),allocatable::realv(:)
    type(PointArray),pointer::v(:)
    logical::OK
    OK=.false.
    call varinit(n_var,v)
    open(unit=9,file='optim.data',status='old',iostat=ie)
    if(ie==0) then
      allocate(realv(n_var))
      do while(ie==0)
        read(9,*,iostat=ie) realv(1:n_var),omega
      end do
      close(9)
      do i=1,n_var
        v(i)%pp=realv(i)
      end do
    end if
    if(qnewton_connect(n_var,v)) then
      call qnewton_run(VCA_potthoff_functional,omega)
      OK=.true.
      call qnewton_disconnect()
    end if
    deallocate(v)
  end function
  
  ! This subroutine is used to assign variational parameters
  ! to an array of pointer
  ! 
  ! input:
  ! n - integer: number of variational parameters
  ! v - type(PointArray),pointer: an initialized array of 
  !               pointers
  ! 
  ! output:
  ! v - type(PointArray),pointer: an array of pointers that holds
  !               memory address of all the variational parameters
  subroutine varinit(n,v)
    implicit none
    type(PointArray),pointer::v(:)
    integer::n
    n=2
    allocate(v(n))
    v(2)%pp=>cluster%hop%cmu
    v(1)%pp=>cluster%hop%ct(1)
    !v(2)=cluster%ty
    !write(7,*) "n=",n
  end subroutine

  ! calculate the potthoff functional of a a system
  ! there is no variational calculation in this subroutine
  !
  ! output:
  ! x - real(8): potthoff functional or grandpotential if the 
  !           system is optimized
  function VCA_potthoff_functional() result(x)
    implicit none
    type(Q_type),pointer::Q
    real(8)::omega,xup,xdn,x
    integer::orbit
    if(.not.associated(cluster)) return
    orbit=cluster%nsite
    Q=>lehmann_init(cluster,omega,initm)
    if(.not.associated(Q)) return
    x=(omega+lehmann_potthoff_functional(cluster,Q)+cluster%hop%shift)/orbit
    !write(*,'(4F20.15)') cluster%hop%ct,cluster%hop%cmu,x
    call lehmann_Q_clean(Q)
  end function

  ! calculate the electron density (concentration) of a 
  ! system. The result is valid only if the system is already
  ! optimized.
  ! 
  ! output:
  ! density - real(8): electron density
  function VCA_particle_density() result(density)
    implicit none
    type(Q_type),pointer::Q
    real(8)::density,omega,dup,ddn
    integer::orbit
    if(.not.associated(cluster)) return
    orbit=cluster_get_n_orbit(cluster)
    Q=>lehmann_init(cluster,omega,initm)
    if(.not.associated(Q)) return
    density=lehmann_particle_density(cluster,Q)/orbit/2
    call lehmann_Q_clean(Q)
  end function

  ! calculate the VCA based spectral function of a system.
  ! The subroutine is set to scan spectral function between
  ! several k-points. You have to define these points first.
  ! change these in the following program.
  !
  ! the output will be a file named "spectr.data"
  subroutine VCA_spectral_function()
    implicit none
    type(Q_type),pointer::Q,p
    real(8)::omega,x,kx,ky,w,Pi,dx,dy,leng,kp(2,4)
    integer::orbit,i,j,k,npoint
    complex(8),pointer,dimension(:,:)::G,gl
    if(.not.associated(cluster)) return
    orbit=cluster_get_n_orbit(cluster)
    allocate(G(orbit*2,orbit*2))
    Q=>lehmann_init(cluster,omega,initm)
    if(.not.associated(Q)) return
    open(unit=7,file='spectr.data',status='replace')
    k=0
    kx=0.d0
    ky=0.d0
    Pi=asin(1.d0)*2

    !------------------------------------
    ! change the following for different lattices
    npoint=4
    kp(1,1)=0.d0
    kp(2,1)=0.d0
    kp(1,2)=0.d0
    kp(2,2)=4.d0*PI/3/sqrt(3.d0)
    kp(1,3)=1.0*PI/3
    kp(2,3)=sqrt(3.d0)*PI/3
    kp(1,4)=0.d0
    kp(2,4)=0.d0
    !-----------------------------------

    do i=1,npoint-1
      kx=kp(1,i)
      ky=kp(2,i)
      dy=kp(2,i+1)-kp(2,i)
      dx=kp(1,i+1)-kp(1,i)
      leng=sqrt(dx*dx+dy*dy)
      dy=dy/leng/10
      dx=dx/leng/10
      do
        if(sqrt((kx-kp(1,i))*(kx-kp(1,i))+(ky-kp(2,i))*(ky-kp(2,i)))>leng) exit
        print *,k,kx,ky
        w=-5
        k=k+1
        do
          if(w>5) exit
          p=>lehmann_clone_Q(Q)
          call lehmann_transform_Q(cluster,p,kx,ky)
          call lehmann_green_function(w+Xi*eps,orbit,p,G)
          gl=>lehmann_periodization(cluster,G,kx,ky)
          x=-imag(gl(1,1)+gl(2,2))/Pi
          write(7,*) k,w,x
          call lehmann_Q_clean(p)
          w=w+0.1
        end do
        kx=kx+dx
        ky=ky+dy
      end do
    end do
    call lehmann_Q_clean(Q)
    deallocate(G)
    close(7)
  end subroutine

  ! Scan the DOS of the system.
  ! 
  ! input:
  ! Emin - real(8): the minimum energy of this calculation
  ! Emax - real(8): the maximum energy of this calculation
  ! bin - integer: number of energy slot that will be calculated
  !
  ! output will a file named "DOS.data"
  subroutine VCA_DOS(Emin,Emax,bin)
    implicit none
    type(Q_type),pointer::Q,p,a
    real(8)::omega,kx,ky,Pi,w,Emin,Emax,dE
    integer,allocatable::S(:)
    real(8),allocatable::E(:),sf(:)
    complex(8),pointer,dimension(:,:)::G,pg
    integer::orbit,i,j,k,bin,l
    logical::OK
    if(.not.associated(cluster)) return
    Pi=asin(1.0)*2
    orbit=cluster_get_n_orbit(cluster)
    Q=>lehmann_init(cluster,omega,initm)
    if(.not.associated(Q)) return
    allocate(S(bin))
    allocate(E(0:bin+1))
    allocate(sf(0:bin+1))
    allocate(G(orbit*2,orbit*2))
    dE=(Emax-Emin)/bin
    E(0)=Emin-dE/2
    do i=1,bin
      E(i)=E(0)+dE*i
      S(i)=0
    end do
    E(bin+1)=Emax+dE/2
    l=0
    do i=0,nk-1
      do j=0,nk-1
        kx=i*2*Pi/nk-Pi
        ky=j*2*Pi/nk-Pi
        p=>lehmann_clone_Q(Q)
        call lehmann_transform_Q(cluster,p,kx,ky)
        do k=0,bin+1
          call lehmann_green_function(E(k)+eps*Xi,orbit,p,G)
          !call matrix_inverse(orbit*2,G,OK)
          pg=>lehmann_periodization(cluster,G,kx,ky)
          sf(k)=-imag(pg(1,1)+pg(2,2))/Pi
          deallocate(pg)
        end do
        do k=1,bin
          if(sf(k)-sf(k-1)>0.and.sf(k)-sf(k+1)>0) then
            S(k)=S(k)+1
            l=l+1
          end if
        end do
        call lehmann_Q_clean(p)
      end do
    end do
    open(unit=10,file='DOS.data',status='replace')
    do i=1,bin
      write(10,*) E(i),S(i)*2.d0/l
    end do
    close(10)
    deallocate(G)
    call lehmann_Q_clean(Q)
  end subroutine

  ! calculate the fermi surface of the system
  ! the calculation is setup to calculation within a square
  ! in the 2D k-space: from (-pi,-pi) to (pi,pi)
  !
  ! output will be a file named "fermi.data"
  subroutine VCA_fermi_surface()
    implicit none
    type(Q_type),pointer::Q,p,a
    real(8)::omega,x,kx,ky,Pi
    complex(8)::gl
    integer::orbit,i,j
    logical::OK
    complex(8),pointer,dimension(:,:)::G,pg
    if(.not.associated(cluster)) return
    orbit=cluster_get_n_orbit(cluster)
    allocate(G(orbit*2,orbit*2))
    Pi=asin(1.d0)*2
    Q=>lehmann_init(cluster,omega,initm)
    if(.not.associated(Q)) return
    open(unit=8,file='fermi.data',status='replace')
    do i=0,nk
      do j=0,nk
        kx=i*2*Pi/nk-Pi
        ky=j*2*Pi/nk-Pi
        p=>lehmann_clone_Q(Q)
        call lehmann_transform_Q(cluster,p,kx,ky)
        call lehmann_green_function(eps*Xi,orbit,p,G)
        !call matrix_inverse(orbit*2,G,OK)
        pg=>lehmann_periodization(cluster,G,kx,ky)
        x=-imag(pg(1,1)+pg(2,2))/Pi
        write(8,*) kx,ky,x
        deallocate(pg)
        call lehmann_Q_clean(p)
      end do
    end do
    deallocate(G)
    call lehmann_Q_clean(Q)
    close(8)
  end subroutine



end module VCA
