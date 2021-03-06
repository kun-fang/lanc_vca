!--------------------------------------------------
!
! timer.F90
! module: timer_mod
! requirement: none
!
! created by Kun Fang
!
! This module create an object to count the runtime.
!
!
!------------------------------------------------------
module timer_mod
  implicit none
  
  real,public,save::itime=0.0
  
  contains

!----------------------public-----------------------

  ! initialize the timer
  subroutine timer_init()
    implicit none
    real,dimension(2)::t
    real::a
    call etime(t,a)
    itime=t(1)
    !print *,p%init
    !print *,t(1)
    !print *,t(2)
  end subroutine

  ! return the cpu runtime from initialization 
  ! of the timer.
  function cputime() result(x)
    implicit none
    real,dimension(2)::t
    real::a,x
    call etime(t,a)
    x=t(1)-itime
    !print *,"------------------------------------------"
    !print *,p%now
    !print *,t(1)
    !print *,t(2)
  end function

end module timer_mod
