module matrix_mod
	implicit none

	type real_spot
		real(8)::val=0.00
		integer::index=0
		type(real_spot),pointer::next=>NULL()
	end type

	type real_matrix
		integer::dim=0
		type(real_spot),pointer::head(:)=>NULL()
	end type real_matrix

	public::real_matrix_init,real_matrix_del,real_matrix_show,real_matrix_set,real_matrix_product

contains

!-----------------matrix-----------------------

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

end module matrix_mod