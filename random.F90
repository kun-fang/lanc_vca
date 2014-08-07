!--------------------------------------------------
!
! random.F90
! module: rand_mod
! requirement: none
!
! created by Kun Fang
!
! This module provides interfaces to create random numbers
! and random vectors. This program is required by underwood 
! module so that the Lanczos method can run correctly.
!
!
!------------------------------------------------------

module rand_mod
  implicit none
  
  ! some parameters
  integer,save,private::seed=0        ! random seed
  real(8),save,private::rand=0.0      ! last random number obtaind
  logical,save,private::init=.false.  ! flag for whether the module is initialized
  
  public::random_n,random_vector
  private::init_seed,new_seed
  
  contains
  
!-----------------------------public-----------------------------

  ! extract a random number within a range
  ! 
  ! input:
  ! a - real(8),optional: lower bound of the random number
  ! b - real(8),optional: upper bound of the random number
  !
  ! output:
  ! x - real(8): the required random number
  function random_n(a,b) result(x)
    real(8),intent(in),optional::a,b
    real(8)::x
    if(.not.init) call init_seed()
    call new_seed()
    if(present(a).and.present(b)) then
      x=a+rand*(b-a)
    else
      x=rand
    end if
  end function

  ! create a random vector
  ! 
  ! input:
  ! n - integer: dimension of the vector
  ! a - real(8): length of the vector
  ! x - real(8): a real vector
  ! 
  ! output:
  ! x - real(8): a real random vector with length
  !              equal to a
  subroutine random_vector(n,a,x)
    integer::n,i
    real(8)::a,sum,x(n)
    if(a<0.or.n<=0) return
    if(.not.init) call init_seed()
    sum=0.0
    do i=1,n
      call new_seed()
      x(i)=rand
      sum=sum+x(i)*x(i)
    end do
    sum=sqrt(sum)
    do i=1,n
      x(i)=x(i)*a/sum
    end do
  end subroutine

!----------------------------private-----------------------------

  !Get initial seed for random routine
  subroutine init_seed()
    integer*4,dimension(3)::today,now
    call itime(now)
    call idate(today)
    seed=today(3)+70*(today(2)+12*(today(1)+31*(now(1)+23*(now(2)+59*now(3)))))
    init=.true.
  end subroutine
  
  !Get a new seed
  subroutine new_seed()
    integer::a,b,m,q,r,Ixx
    a = 16807
    b = 0
    m = 2147483647
    q=int(m/a)
    r=mod(m,a)
    Ixx = a*(mod(seed,q))-r*int(seed/q)
    if (Ixx.lt.0) Ixx = Ixx + m
    seed=Ixx
    rand=1.0*seed/m
  end subroutine

end module rand_mod