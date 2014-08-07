!--------------------------------------------------
!
! qnewton.F90
! module: qnewton
! requirement: rand_mod
!
! created by Kun Fang
!
! This module is used to perform variational calculations. This
! module uses the quasi-newton method to find out the stationary
! points of a function. Please refer to the wikipedia page for 
! details of this method.
!
! This module is setup to find out saddle points in the function
! space.
!
!
! type
! PointArray
!
!------------------------------------------------------

module qnewton
  use rand_mod
  implicit none
  
  real(8),private::tol=3.d-5
  real(8),allocatable,private::dev(:)
  
  type PointArray
    real(8),pointer::pp
  end type PointArray
  
  type(PointArray),pointer,private::vari_para(:)=>NULL()
  integer,private::n=0
  real(8),private,allocatable::var(:)

contains

  ! Connect your functions to this module for variational
  ! calculations. This is required before using this 
  ! module.
  !
  ! input:
  ! n_para - integer: number of parameters in the variational 
  !                   calculation
  ! para - type(PointArray),pointer: an array of pointers to 
  !              the parameters used in the calculation
  ! 
  ! output:
  ! para - type(PointArray),pointer: an array of pointers to
  !         the final results of the variational calculations
  function qnewton_connect(n_para,para) result(OK)
    implicit none
    integer::n_para,i
    type(PointArray),pointer::para(:)
    logical::OK
    OK=.false.
    if(n_para<=0.or..not.associated(para).or.size(para)==0) return
    vari_para=>para
    n=n_para
    allocate(var(n))
    do i=1,n
      var(i)=vari_para(i)%pp
    end do
    OK=.true.
  end function
  
  ! disconnect the function from this module
  ! it is recommended after you finish the variational 
  ! calculation
  subroutine qnewton_disconnect()
    implicit none
    vari_para=>NULL()
    n=0
    deallocate(var)
  end subroutine
  
  ! the function that will be applied to the calculation
  ! 
  ! input:
  ! n - integer: number of parameters
  ! var - real(8): an array to store the parameters
  ! f - real(8),external: the function that will be evaluated
  !
  ! output:
  ! x - real(8): result of the function f.
  function FucVal(n,var,f) result(x)
    integer::n,i
    real(8)::var(n),x
    real(8),external::f
    do i=1,n
      vari_para(i)%pp=var(i)
    end do
    x=f()
  end function
  

  subroutine qnewton_run(f,fx,conv)
  !---------------------------------------------------
  !
  ! quasi-newton method
  !
  ! on input ---
  !
  ! n:   number of parameters to optimize
  ! var: an array of pointers containing pointers to 
  !      n variational parameters, it can be updated
  !      during the optimization
  ! f:   an external function having the function 
  !      form f() gives the equation that needs to 
  !      be optimized
  ! fx:  a real(8) number
  !
  ! on output ---
  ! 
  ! var: a vector containing the optimized parameters
  ! fx:  the value of opitimized function
  !
  !--------------------------------------------------
    integer::i,j,k
    real(8)::x0(n),x1(n),g0(n),g1(n),B(n,n),H(n,n),S(n),q(n),p(n),fx,Z,alpha,d(n),c(n,n)
    real(8),external::f
    real(8),optional::conv
    allocate(dev(n))
    if(present(conv)) tol=conv
    dev(1:n)=1.d-2
    k=0
    B(1:n,1:n)=0.0
    H(1:n,1:n)=0.0
    g0(1:n)=0.0
    do i=1,n
      B(i,i)=1.0
      H(i,i)=1.0
    end do
    x0=var
    open(unit=433,file='optim.data',status='replace')
    call gradient(f,n,x0,g0)
    do
      !-----------the kth iteration-----------
      k=k+1
      !-----------find the next trial----------
      S=-matmul(H,g0)
      Z=sqrt(dot_product(S,S))
      !----------if the criteria is achieved, then exit-----------
      if(Z<tol) exit
      !----------if not, find the next trial------------
      alpha=random_n(5.d-1,1.d0)
      S=alpha*S
      x1=x0+S
      do i=1,n
        if(abs(S(i)/10)>1.d-2) then
          dev(i)=1.d-2
        else
          dev(i)=abs(S(i)/10)
        end if
      end do
      call gradient(f,n,x1,g1)
      q=g1-g0
      p=x1-x0
      call update(n,B,H,p,q)
      x0=x1
      g0=g1
    end do
    var=x0+S
    fx=FucVal(n,var,f)
    deallocate(dev)
    close(433)
  end subroutine

  subroutine update(n,B,H,p,q)
    integer::n
    real(8)::p(n),q(n),B(n,n),H(n,n),d(n),c(n,n)
    !----------update B--------------
    d=q-matmul(B,p)
    c=tensor_product(d,d)/dot_product(d,p)
    B=B+c
    !----------update H--------------
    d=p-matmul(H,q)
    c=tensor_product(d,d)/dot_product(d,q)
    H=H+c
  end subroutine

  subroutine gradient(f,n,x,df)
    integer::n,i
    real(8)::x(n),df(n),f0,dx(n)
    real(8),external::f
    f0=FucVal(n,x,f)
    write(433,*) x,f0
    do i=1,n
      dx=x
      dx(i)=dx(i)-dev(i)
      df(i)=(f0-FucVal(n,dx,f))/dev(i)
    end do
  end subroutine

  function tensor_product(a,b) result(c)
    real(8)::a(:),b(:),c(size(a),size(b))
    integer::n,m,i,j
    n=size(a)
    m=size(b)
    do i=1,n
      do j=1,m
        c(i,j)=a(i)*b(j)
      end do
    end do
  end function

end module qnewton
