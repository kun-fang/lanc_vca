!--------------------------------------------------
!
! matrix.F90
! module: matrix_mod
! requirement: none
!
! created by Kun Fang
!
! This module contains routines for sparse matrix
! creation and calculation. There are two versions
! of sparse matrix: real version and complex version.
! The other part of this program is setup for complex
! version, so it is recommended to always use the
! complex version.
!
! The sparse matrix is stored by rows and only non-zero
! terms are stored. Each non-zero term is stored in a 
! node including its value and column index. Each row is 
! stored as a linked-list of these nodes. The matrix is
! represented as an array of these linked-lists.
!
! types:
! real_matrix
!  |- real_spot
!
! complex_matrix
!  |- complex_spot
!
!------------------------------------------------------

module matrix_mod
	implicit none
	complex(8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)

	! real matrix element node
	type real_spot
		real(8)::val=0.00	! element value
		integer::index=0	! column index
		type(real_spot),pointer::next=>NULL()
	end type

	! complex matrix element node
	type complex_spot
		complex(8)::val=0.00	! element value
		integer::index=0	! column index
		type(complex_spot),pointer::next=>NULL()
	end type

	! real matrix
	type real_matrix
		integer::dim=0
		type(real_spot),pointer::head(:)=>NULL()
	end type real_matrix

	! comple matrix
	type complex_matrix
		integer::dim=0
		type(complex_spot),pointer::head(:)=>NULL()
	end type complex_matrix

	public::real_matrix_init,real_matrix_del,real_matrix_show,real_matrix_set,real_matrix_product
	public::complex_matrix_init,complex_matrix_del,complex_matrix_show,complex_matrix_set,complex_matrix_product

contains

!--------------- Sparse Matrix real version ---------------
!
! It is stored by rows
!
!----------------------------------------------------------

	! initialize the matrix
	function real_matrix_init(n) result(M)
		type(real_matrix),pointer::M
		integer::n,i
		allocate(M)
		M%dim=n
		allocate(M%head(n))
		do i=1,n
			M%head(i)%val=0.0
			M%head(i)%index=i
			M%head(i)%next=>NULL()
		end do
	end function

	! delete the matrix from the memory
	subroutine real_matrix_del(M)
		type(real_matrix),pointer::M
		type(real_spot),pointer::a,b
		integer::i
		do i=1,M%dim
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				b=>a
				a=>a%next
				deallocate(b)
			end do
		end do
		deallocate(M%head)
		M%head=>NULL()
		M%dim=0
		deallocate(M);
	end subroutine

	! print all the non-zero terms of the matrix
	subroutine real_matrix_print(M)
		type(real_matrix),pointer::M
		integer::n,s,i,j
		type(real_spot),pointer::a
		n=M%dim
		do i=1,n
			print *,M%head(i)%index,M%head(i)%index,M%head(i)%val
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				print *,i,a%index,a%val
				a=>a%next
			end do
		end do
	end subroutine

	! return the value of matrix element M(i,j)
	function real_matrix_show(M,i,j) result(x)
		type(real_matrix),pointer::M
		type(real_spot),pointer::a
		integer::i,j,k
		real(8)::x
		x=0.0d0
		if(i>M%dim.or.j>M%dim) return
		if(i==j) then
			x=M%head(i)%val
			return
		end if
		a=>M%head(i)%next
		do
			if(.not.associated(a)) exit
			if(a%index==j) then
				x=a%val
				return
			end if
			if(a%index>j) return
			a=>a%next
		end do
	end function
	
	! calculate y=M*x, where M is a matrix, x and y both are 
	! column vector
	!
	! Note: the subroutine assume that the row dimension of the
	!       matrix is the same as the length of the two vectors.
	!       The subroutine doesn't check the validity of the 
	!       matrix production
	subroutine real_matrix_product(M,x,y)
		type(real_matrix),pointer::M
		type(real_spot),pointer::a
		real(8)::x(:),y(:)
		integer::i,j,k
		do i=1,M%dim
			y(i)=M%head(i)%val*x(i)
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				y(i)=y(i)+a%val*x(a%index)
				a=>a%next
			end do
		end do
	end subroutine

	! assign a value to matrix element M(i,j)
	subroutine real_matrix_set(M,i,j,x)
		type(real_matrix),pointer::M
		integer::i,j
		real(8)::x
		type(real_spot),pointer::a,b,new
		if(i==j) then
			M%head(i)%val=x
			return
		end if
		b=>M%head(i)
		a=>M%head(i)%next
		do
			if(.not.associated(a)) exit
			if(a%index==j) then
				if(abs(x)<1.d-8) then
					b%next=>a%next
					deallocate(a)
				else
					a%val=x
				end if
				return
			end if
			if(a%index>j) exit
			b=>a
			a=>a%next
		end do
		if(abs(x)<1.d-8) return
		allocate(new)
		new%index=j
		new%val=x
		if(associated(a)) then
			b%next=>new
			new%next=>a
		else
			b%next=>new
		end if
	end subroutine

	! convert a sparse matrix to an 2D array
	subroutine real_matrix_convert(M,A)
		type(real_matrix),pointer::M
		type(real_spot),pointer::p
		real(8)::A(:,:)
		integer::i,j,k,n
		n=size(A(1,:))
		A(1:n,1:n)=0.d0
		do i=1,n
			p=>M%head(i)
			A(i,i)=p%val
			do
				p=>p%next
				if(.not.associated(p)) exit
				j=p%index
				A(i,j)=p%val
			end do
		end do
	end subroutine

	! convert an 2D array to a sparse matrix
	subroutine real_matrix_reverse(A,M)
		type(real_matrix),pointer::M
		type(real_spot),pointer::p,new
		real(8)::A(:,:)
		integer::i,j,k,n
		n=size(A(1,:))
		M=>real_matrix_init(n)
		do i=1,n
			p=>M%head(i)
			p%val=A(i,i)
			do j=1,n
				if(abs(A(i,j))<1.d-6) cycle
				allocate(new)
				new%index=j
				new%val=A(i,j)
			end do
		end do
	end subroutine


!--------------- Sparse Matrix complex version ---------------
!
! It is stored in rows
!
!-------------------------------------------------------------

	! initialize the matrix
	function complex_matrix_init(n) result(M)
		type(complex_matrix),pointer::M
		integer::n,i
		allocate(M)
		M%dim=n
		allocate(M%head(n))
		do i=1,n
			M%head(i)%val=Zero
			M%head(i)%index=i
			M%head(i)%next=>NULL()
		end do
	end function

	! delete the matrix from the memory
	subroutine complex_matrix_del(M)
		type(complex_matrix),pointer::M
		type(complex_spot),pointer::a,b
		integer::i
		do i=1,M%dim
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				b=>a
				a=>a%next
				deallocate(b)
			end do
		end do
		deallocate(M%head)
		M%head=>NULL()
		M%dim=0
		deallocate(M);
	end subroutine

	! print all the non-zero terms of the matrix
	subroutine complex_matrix_print(M)
		type(complex_matrix),pointer::M
		integer::n,s,i,j
		type(complex_spot),pointer::a
		n=M%dim
		do i=1,n
			print *,M%head(i)%index,M%head(i)%index,M%head(i)%val
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				print *,i,a%index,a%val
				a=>a%next
			end do
		end do
	end subroutine

	! return the value of matrix element M(i,j)
	function complex_matrix_show(M,i,j) result(x)
		type(complex_matrix),pointer::M
		type(complex_spot),pointer::a
		integer::i,j,k
		complex(8)::x
		x=Zero
		if(i>M%dim.or.j>M%dim) return
		if(i==j) then
			x=M%head(i)%val
			return
		end if
		a=>M%head(i)%next
		do
			if(.not.associated(a)) exit
			if(a%index==j) then
				x=a%val
				return
			end if
			if(a%index>j) return
			a=>a%next
		end do
	end function
	
	! calculate y=M*x, where M is a matrix, x and y both are 
	! column vector
	!
	! Note: the subroutine assume that the row dimension of the
	!       matrix is the same as the length of the two vectors.
	!       The subroutine doesn't check the validity of the 
	!       matrix production
	subroutine complex_matrix_product(M,x,y)
		type(complex_matrix),pointer::M
		type(complex_spot),pointer::a
		complex(8)::x(:),y(:)
		integer::i,j,k
		do i=1,M%dim
			y(i)=M%head(i)%val*x(i)
			a=>M%head(i)%next
			do
				if(.not.associated(a)) exit
				y(i)=y(i)+a%val*x(a%index)
				a=>a%next
			end do
		end do
	end subroutine

	! assign a value to matrix element M(i,j)
	subroutine complex_matrix_set(M,i,j,x)
		type(complex_matrix),pointer::M
		integer::i,j
		complex(8)::x
		type(complex_spot),pointer::a,b,new
		if(i==j) then
			M%head(i)%val=x
			return
		end if
		b=>M%head(i)
		a=>M%head(i)%next
		do
			if(.not.associated(a)) exit
			if(a%index==j) then
				if(abs(x)<1.d-8) then
					b%next=>a%next
					deallocate(a)
				else
					a%val=x
				end if
				return
			end if
			if(a%index>j) exit
			b=>a
			a=>a%next
		end do
		if(abs(x)<1.d-8) return
		allocate(new)
		new%index=j
		new%val=x
		if(associated(a)) then
			b%next=>new
			new%next=>a
		else
			b%next=>new
		end if
	end subroutine

	! convert a sparse matrix to an 2D array
	subroutine complex_matrix_convert(M,A)
		type(complex_matrix),pointer::M
		type(complex_spot),pointer::p
		complex(8)::A(:,:)
		integer::i,j,k,n
		n=size(A(1,:))
		A(1:n,1:n)=Zero
		do i=1,n
			p=>M%head(i)
			A(i,i)=p%val
			do
				p=>p%next
				if(.not.associated(p)) exit
				j=p%index
				A(i,j)=p%val
			end do
		end do
	end subroutine

	! convert an 2D array to a sparse matrix
	subroutine complex_matrix_reverse(A,M)
		type(complex_matrix),pointer::M
		type(complex_spot),pointer::p,new
		complex(8)::A(:,:)
		integer::i,j,k,n
		n=size(A(1,:))
		M=>complex_matrix_init(n)
		do i=1,n
			p=>M%head(i)
			p%val=A(i,i)
			do j=1,n
				if(abs(A(i,j))<1.d-6) cycle
				allocate(new)
				new%index=j
				new%val=A(i,j)
			end do
		end do
	end subroutine

end module matrix_mod