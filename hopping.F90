!--------------------------------------------------
!
! hopping.F90
! module: hop_mod
! requirements: none
!
! created by Kun Fang
!
! This module implements a cluster. It provides some interfaces
! for a cluster object including the coordinates and hoppings
! of a cluster.
!
! A cluster is implemented as a set of sites in real space and 
! hopping relation between these sites. The object also includes
! translation vectors that represent the vectors needed to tile
! the clusters to form an infinite lattice.
!
! The information of a cluster is imported from a cluster 
! configuration file and an input file. The cluster configuration
! file includes its name, the coordinates of all the sites and
! the translation vectors. The input file include the hopping 
! relation between cluster sites for both the lattice and the 
! cluster. 
! 
!
!
! types:
! hop
!
!------------------------------------------------------

module hop_mod
  implicit none
  
  complex(8),private,parameter::Zero=(0.0,0.0),Xi=(0.0,1.0)
  
  ! the type represents a cluster
  type,public::hop
    integer::site,dim,nlt,nct,nambu=0,scty=0  ! some parameters: site - number of sites, dim - number of cluster dimensions (usually 2)
                                              ! nlt - number of hopping relation for lattice, nct - number of hopping relation for cluster
                                              ! nambu - flag for nambu representation (work for calculation of SC)
                                              ! scty - type of SC (work with nambu=1)
    real,pointer,dimension(:,:)::coordinate ! coordinates of the cluster sites
    real,pointer,dimension(:,:)::vector     ! translation vector of the cluster
    integer,pointer,dimension(:,:)::lattice ! hopping relation for lattices
    integer,pointer,dimension(:,:)::cluster ! hopping relation for clusters
    real(8),pointer,dimension(:)::lt,ct     ! the final hopping matrix for lattice and cluster
    real(8)::lmu=0.d0,U=0.d0,Tep=0.d0   ! parameters for lattices
    real(8)::cmu=0.d0,shift=0.d0,M=0.d0,delta=0.d0,D=0.d0,h=0.d0,VSO=0.00   ! parameters for clusters
    character(len=20)::SY=""    ! a string to store the types of Weiss field included
  end type hop

  contains

  ! delete the cluster object from the memory
  subroutine cluster_delete(cl)
    type(hop),pointer::cl
    deallocate(cl%coordinate)
    deallocate(cl%vector)
    deallocate(cl%lattice)
    deallocate(cl%cluster)
    deallocate(cl%lt)
    deallocate(cl%ct)
    deallocate(cl)
  end subroutine

  ! initialize a cluster in the memory
  ! two external files are required:
  !   cluster.input
  ! This file include the name of cluster and hopping information 
  ! of the cluster and the lattice
  !   XXX.conf
  ! XXX is the name of the cluster. It must be the same as the 
  ! name that shows up in cluster.input
  function cluster_init() result(cl)
    type(hop),pointer::cl
    character(len=20)::filename,sy,rc,a
    real(8),allocatable::tr(:,:)
    complex(8),pointer::test(:,:)
    integer::site,dim,i,j,k,ierr
    allocate(cl)

    ! read cluster.input file
    open(unit=8,file='cluster.input',status='old',action='read')
    read(8,*) rc
    ! find out the name of cluster XXX
    filename=trim(adjustl(rc))//'.conf'
    ! read XXX.conf file
    open(unit=7,file=trim(adjustl(filename)),status='old',action='read')
    read(7,*)
    read(7,*) dim
    read(7,*)
    read(7,*) site
    cl%dim=dim
    cl%site=site
    allocate(cl%coordinate(dim,site))
    allocate(cl%vector(dim,dim))
    read(7,*)
    ! read coordinates of the cluster sites
    do i=1,site
      read(7,*) cl%coordinate(1:dim,i)
    end do
    read(7,*)
    ! read translation vectors of the cluster
    do i=1,dim
      read(7,*) cl%vector(1:dim,i)
    end do
    close(7)
    read(8,*)
    read(8,*)
    read(8,*)
    i=0
    do
      i=i+1
      read(8,*,iostat=ierr) j
      if(ierr/=0) exit
    end do
    do j=1,i+1
      backspace(8)
    end do
    ! read hopping relations of the lattice
    allocate(cl%lattice(5,i-1))
    k=0
    cl%nlt=i-1
    do j=1,cl%nlt
      read(8,*) cl%lattice(1:5,j)
      if(cl%lattice(3,j)>k) k=cl%lattice(3,j)
    end do
    read(8,*)
    allocate(cl%lt(k))
    do i=1,k
      read(8,*) rc,cl%lt(i)
    end do
    read(8,*)
    read(8,*) rc,cl%lmu
    read(8,*) rc,cl%U
    read(8,*) rc,cl%Tep
    read(8,*)
    read(8,*)
    read(8,*)
    i=0
    do
      i=i+1
      read(8,*,iostat=ierr) j
      if(ierr/=0) exit
    end do
    do j=1,i+1
      backspace(8)
    end do
    ! read hopping relations for the cluster
    allocate(cl%cluster(3,i-1))
    k=0
    cl%nct=i-1
    do j=1,cl%nct
      read(8,*) cl%cluster(1:3,j)
      if(cl%cluster(3,j)>k) k=cl%cluster(3,j)
    end do
    read(8,*)
    allocate(cl%ct(k))
    do i=1,k
      read(8,*) rc,cl%ct(i)
    end do
    read(8,*)
    read(8,*) rc,cl%cmu
    ! deal with any inputs of Weiss fields
    do
      read(8,*,iostat=ierr) sy
      if(ierr/=0) exit
      select case(sy(1:2))
        case("AF")
          read(8,*) rc,cl%M
        case("SC")
          cl%scty=ichar(sy(3:3))-48
          cl%nambu=1
          read(8,*) rc,cl%D
        case("NM")
          read(8,*) rc,cl%delta
        case("SP")
          read(8,*) rc,cl%h
        case("RS")
          read(8,*) rc,cl%VSO
        case default
          cycle
      end select
      if(sy(1:2)=="SC") then
        cl%SY=sy(1:2)//trim(cl%SY)
      else
        cl%SY=trim(cl%SY)//sy(1:2)
      end if
    end do
    close(8)
  end function
  
  ! update the hopping matrix
  !
  ! input:
  ! cl - type(hop): a pointer to a cluster object
  !
  ! output:
  ! tprime - complex(8): a poniter to a hopping matrix
  !                      with dimension = site*2 
  function cluster_update(cl) result(tprime)
    implicit none
    type(hop),pointer::cl
    character(len=2)::rc
    complex(8),pointer::tprime(:,:)
    integer::i,j,k,site
    site=cl%site
    allocate(tprime(site*2,site*2))
    tprime(1:site*2,1:site*2)=Zero
    do i=1,site
      tprime(i,i)=-cl%cmu
      tprime(i+site,i+site)=-cl%cmu
    end do

    do k=1,cl%nct
      i=cl%cluster(1,k)
      j=cl%cluster(2,k)
      tprime(i,j)=-cl%ct(cl%cluster(3,k))
      tprime(i+site,j+site)=tprime(i,j)
      tprime(j,i)=conjg(tprime(i,j))
      tprime(j+site,i+site)=tprime(j,i)
    end do

    do i=1,19,2
      rc=cl%SY(i:i+1)
      if(rc=="") exit
      select case(rc)
        case("AF")
          call cluster_AF(cl,tprime)
        case("SC")
          call cluster_SC(cl,tprime)
        case("NM")
          call cluster_NM(cl,1,2,tprime)
          call cluster_NM(cl,3,4,tprime)
        case("SP")
          call cluster_SP(cl,tprime)
        case("RS")
          call cluster_RS(cl,tprime)
        case default
          cycle
      end select
    end do
  end function

  ! AF Weiss field
  subroutine cluster_AF(cl,tprime)
    implicit none
    type(hop),pointer::cl
    complex(8),pointer::tprime(:,:)
    integer::i,j,spin,site
    site=cl%site
    do i=1,site
      do spin=1,2
        tprime(i+site*(spin-1),i+site*(spin-1))=tprime(i+site*(spin-1),i+site*(spin-1))+cl%M*(-1)**(spin+i)
      end do
    end do
  end subroutine

  subroutine print_matrix(m,n,s1,s2)
    implicit none
    complex(8),pointer::m(:,:)
    integer::n,s1,s2,i,j
    do i=1+s1*n,n+s1*n
      do j=1+s2*n,i-1-s1*n+s2*n
        if(abs(m(i,j))>1.d-5) print *,i-s1*n,j-s2*n,real(m(i,j)),imag(m(i,j))
      end do
    end do
    print *,""
  end subroutine

  ! Rashiba Weiss field
  subroutine cluster_RS(cl,tprime)
    implicit none
    type(hop),pointer::cl
    complex(8),pointer::tprime(:,:)
    real(8)::a,PI,v(2)
    integer::i,j,k,l,site
    site=cl%site
    call print_matrix(tprime,site,0,0)
    do k=1,cl%nct
      i=cl%cluster(1,k)
      j=cl%cluster(2,k)
      if(i>j) then
        i=i+j
        j=i-j
        i=i-j
      end if
      v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
      if(abs(v(1)-3)<1.d-5) v(1)=1.d0
      if(abs(v(1))<1.d-5.or.abs(v(2)-1)<1.d-5) then
        tprime(i,j+site)=tprime(i,j+site)+cl%VSO*Xi
        tprime(i+site,j)=tprime(i+site,j)+cl%VSO*Xi
        tprime(j+site,i)=conjg(tprime(i,j+site))
        tprime(j,i+site)=conjg(tprime(i+site,j))
      else if(abs(v(2))<1.d-5.or.abs(v(1)-1)<1.d-5) then
        tprime(i,j+site)=tprime(i,j+site)+cl%VSO
        tprime(i+site,j)=tprime(i+site,j)-cl%VSO
        tprime(j+site,i)=conjg(tprime(i,j+site))
        tprime(j,i+site)=conjg(tprime(i+site,j))
      end if
    end do
    call print_matrix(tprime,site,0,1)
    call print_matrix(tprime,site,1,0)
  end subroutine

  ! SC Weiss field
  subroutine cluster_SC(cl,tprime)
    implicit none
    type(hop),pointer::cl
    complex(8),pointer::tprime(:,:)
    real(8)::v(2),PI,leng
    integer::i,j,site,ty
    site=cl%site
    PI=asin(1.d0)*2
    do i=1+site,2*site
      do j=1+site,2*site
        tprime(i,j)=-tprime(i,j)
      end do
    end do
    cl%shift=-cl%cmu*site
    do i=1,site
      tprime(i,i)=tprime(i,i)+cl%U
    end do
    select case(cl%scty)
    case(1)
      !-------s-wave SC--------
      do i=1,site
        tprime(i,site+i)=tprime(i,site+i)-cl%D
        tprime(site+i,i)=conjg(tprime(i,site+i))
      end do
      !------------------------
    case(2)
      !---d_(x^2-y^2)-wave SC--
      do i=1,site
        do j=i+1,site
          v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
          tprime(i,j+site)=tprime(i,j+site)-cl%D*(cos(v(1)*PI/2)-cos(v(2)*PI/2))
          tprime(j+site,i)=conjg(tprime(i,j+site))
        end do
      end do
      !------------------------
    case(3)
      !---d_(xy)-wave SC--
      do i=1,site
        do j=i+1,site
          v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
          tprime(i,j+site)=tprime(i,j+site)-cl%D*(sin(v(1)*PI/2)*sin(v(2)*PI/2))
          tprime(j+site,i)=conjg(tprime(i,j+site))
        end do
      end do
      !------------------------
    case(4)
      !---d_(s extend)-wave SC--
      do i=1,site
        do j=i+1,site
          v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
          tprime(i,j+site)=tprime(i,j+site)-cl%D*(cos(v(1)*PI/2)+cos(v(2)*PI/2))
          tprime(j+site,i)=conjg(tprime(i,j+site))
        end do
      end do
      !------------------------
    case(5)
      !---s_(+/-)-wave SC--
      do i=1,site
        do j=i+1,site
          v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
          tprime(i,j+site)=tprime(i,j+site)-cl%D*(cos(v(1)*PI/2)*cos(v(2)*PI/2))
          tprime(j+site,i)=conjg(tprime(i,j+site))
        end do
      end do
      !------------------------
    case(6)
      !---d+id-wave SC--
      do i=1,site
        do j=i+1,site
          v(1:2)=cl%coordinate(1:2,i)-cl%coordinate(1:2,j)
          leng=sqrt(dot_product(v,v))
          tprime(i,j+site)=tprime(i,j+site)-cl%D*(((v(1)/leng)*(v(1)/leng)-(v(2)/leng)*(v(2)/leng))+Xi*((v(1)/leng)*(v(2)/leng)))
          tprime(j+site,i)=conjg(tprime(i,j+site))
        end do
      end do
      !------------------------
    case default
    end select
  end subroutine
  
  ! nematic Weiss field
  subroutine cluster_NM(cl,i,j,tprime)
    implicit none
    type(hop),pointer::cl
    complex(8),pointer::tprime(:,:)
    integer::i,j,spin,site
    site=cl%site
    do spin=1,2
      tprime(i+site*(spin-1),j+site*(spin-1))=tprime(i+site*(spin-1),j+site*(spin-1))*(1+cl%delta)
      tprime(j+site*(spin-1),i+site*(spin-1))=conjg(tprime(i+site*(spin-1),j+site*(spin-1)))
    end do
  end subroutine
  
  ! spiral Weiss field
  subroutine cluster_SP(cl,tprime)
    implicit none
    type(hop),pointer::cl
    complex(8),pointer::tprime(:,:)
    real(8),pointer::e(:,:)
    real(8)::a,PI
    integer::i,j,site
    site=cl%site
    PI=asin(1.d0)*2
    allocate(e(2,site))
    e(1,1)=-sin(PI/3)
    e(2,1)=-cos(PI/3)
    a=-2*PI/3
    do i=2,site
      e(1,i)=cos(a)*e(1,i-1)-sin(a)*e(2,i-1)
      e(2,i)=sin(a)*e(1,i-1)+cos(a)*e(2,i-1)
    end do
    do i=1,site
      tprime(i,i+site)=tprime(i,i+site)+cl%h*(e(1,i)-Xi*e(2,i))
      tprime(i+site,i)=conjg(tprime(i,i+site))
    end do
    deallocate(e)
  end subroutine

  ! return the coordinates of all the sites of the cluster
  subroutine cluster_site(cl,n,a)
    type(hop)::cl
    integer::n,dim
    real,dimension(*)::a
    dim=cl%dim
    a(1:dim)=cl%coordinate(1:dim,n)
  end subroutine

  ! calculate distance between two sites a and b
  ! a and b are 2D coordinate of the two sites
  function cluster_site_distance(dim,a,b) result(r)
    integer::dim,i
    real,dimension(dim)::a,b
    real::r
    r=0
    do i=1,dim
      r=r+(a(i)-b(i))*(a(i)-b(i))
    end do
    r=sqrt(r)
  end function

  ! calculate distance between two sites with
  ! indices i and j
  function cluster_distance(cl,i,j) result(r)
    type(hop)::cl
    integer::dim,i,j
    real(8)::r
    r=-1.0
    if(i>cl%site.or.j>cl%site) return
    dim=cl%dim
    r=cluster_site_distance(dim,cl%coordinate(1:dim,i),cl%coordinate(1:dim,j))
  end function

end module hop_mod

