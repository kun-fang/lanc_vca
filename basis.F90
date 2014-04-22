module bas_mod
	implicit none

	!define block_type
	type block_type
		integer::init,fin,elec,spin
	end type block_type
	
	type basis_type
		integer::n_basis,n_orbit,n_block
		integer,pointer,dimension(:)::orbit
		type(block_type),pointer,dimension(:)::block
		integer,pointer,dimension(:,:,:)::c_bas
	end type basis_type
	
	public::basis_init,basis_clean,basis_get_n_block,basis_get_n_orbit,basis_get_n_basis,basis_clone_basis,basis_clone_block
	public::basis_clone_cbas,basis_get_basis,basis_get_block,basis_cplus,basis_c,basis_elec,basis_spin,basis_double_occupy

	private::creation,annihilation,occupation,dou_occupy,num_elec,num_spin
	
	contains
	
!---------------------public------------------------

	!initial the basis
	function basis_init(i_site) result(b)
		implicit none
		integer,intent(in)::i_site
		type(basis_type),pointer::b
		integer::n_basis,n_orbit,n_block
		integer,pointer,dimension(:)::orbit
		type(block_type),pointer,dimension(:)::block
		integer,pointer,dimension(:,:,:)::cmat
		
		integer::i,j,k,l,m,n,q,p,r,s,alpha,t,sgn,u
		if(i_site<=0) then
			nullify(b)
			return
		end if
		n_orbit=i_site
		n_basis=4**n_orbit
		n_block=(n_orbit+1)*(n_orbit+1)
		allocate(orbit(n_basis))
		allocate(block(n_block))
		allocate(cmat(n_basis,n_orbit,2))
		cmat(1:n_basis,1:n_orbit,1:2)=0
		k=1
		q=1
		orbit(k)=0
		block(1)%init=1
		block(1)%fin=1
		block(1)%elec=0
		block(1)%spin=0
		m=1
		l=1
		u=0
		do
			if(u==n_orbit*2) exit
			u=u+1
			do r=m,l
				do s=2,1,-1
					t=0
					do i=block(r)%init,block(r)%fin
						do alpha=1,n_orbit
							j=creation(alpha,s,orbit(i),n_orbit)
							if(j==-1) cycle
							p=0
							do n=block(q)%init,k
								if(orbit(n)==j) then
									p=n
									exit
								end if
							end do
							if(p==0) then
								k=k+1
								orbit(k)=j
								t=t+1
								p=k
							end if
							if(s==1) then
								sgn=num_elec(j,n_orbit,alpha+n_orbit,n_orbit*2+1)
							else
								sgn=num_elec(j,n_orbit,alpha,n_orbit*2+1)
							end if
							cmat(p,alpha,s)=i*(-1)**sgn
							cmat(i,alpha,s)=p*(-1)**sgn
						end do
					end do
					if(t==0) cycle
					q=q+1
					block(q)%init=block(q-1)%fin+1
					block(q)%fin=k
					block(q)%elec=u
					block(q)%spin=num_spin(orbit(k),n_orbit)
				end do
			end do
			m=l+1
			l=q
		end do
		allocate(b)
		b%n_orbit=n_orbit
		b%n_basis=n_basis
		b%n_block=n_block
		b%orbit=>orbit
		b%block=>block
		b%c_bas=>cmat
		
	end function
	
	subroutine basis_clean(b)
		implicit none
		type(basis_type),pointer,intent(inout)::b
		if(.not.associated(b)) return
		deallocate(b%block)
		deallocate(b%orbit)
		deallocate(b%c_bas)
		deallocate(b)
		nullify(b)
	end subroutine

	!get number of basis blocks
	!result: n = number of basis blocks
	function basis_get_n_block(b) result(n)
		implicit none
		type(basis_type),pointer,intent(in)::b
		integer::n
		if(associated(b)) then
			n=b%n_block
		else
			n=0
		end if
	end function
	
	!get number of orbits
	!result: n = number of orbits
	function basis_get_n_orbit(b) result(n)
		implicit none
		type(basis_type),pointer,intent(in)::b
		integer::n
		if(associated(b)) then
			n=b%n_orbit
		else
			n=0
		end if
	end function
	
	!get number of basis
	!result: n = number of basis
	function basis_get_n_basis(b) result(n)
		implicit none
		type(basis_type),pointer,intent(in)::b
		integer::n
		if(associated(b)) then
			n=b%n_basis
		else
			n=0
		end if
	end function

	!clone basis
	!result: p = copy of basis
	function basis_clone_basis(b) result(p)
		implicit none
		type(basis_type),pointer,intent(in)::b
		integer,pointer,dimension(:)::p
		integer::i
		if(associated(b)) then
			allocate(p(b%n_basis))
			do i=1,b%n_basis
				p(i)=b%orbit(i)
			end do
		else
			nullify(p)
		end if
	end function

	!clone block
	!result: p = copy of block
	function basis_clone_block(b) result(p)
		implicit none
		type(basis_type),pointer,intent(in)::b
		type(block_type),pointer,dimension(:)::p
		integer::i
		if(associated(b)) then
			allocate(p(b%n_block))
			do i=1,b%n_block
				p(i)=b%block(i)
			end do
		else
			nullify(p)
		end if
	end function

	!clone c_bas
	!result: p = copy of c_bas
	function basis_clone_cbas(b) result(p)
		implicit none
		type(basis_type),pointer,intent(in)::b
		integer,pointer,dimension(:,:,:)::p
		integer::i,j,k
		if(associated(b)) then
			allocate(p(b%n_basis,b%n_orbit,2))
			do k=1,2
				do j=1,b%n_orbit
					do i=1,b%n_basis
						p(i,j,k)=b%c_bas(i,j,k)
					end do
				end do
			end do
		else
			nullify(p)
		end if
	end function

	!get matrix of nth basis
	!result: b = nth basis
	function basis_get_basis(n,p) result(b)
		implicit none
		integer,intent(in)::n
		type(basis_type),pointer,intent(in)::p
		integer,pointer,dimension(:,:)::b
		integer::i,j
		if(associated(p)) then
			allocate(b(p%n_orbit,2))
			do j=1,2
				do i=1,p%n_orbit
					b(i,j)=occupation(i,j,p%orbit(n),p%n_orbit)
				end do
			end do
		else
			nullify(b)
		end if
	end function

	!get nth element of block
	!result: b = nth element of block
	function basis_get_block(n,p) result(b)
		implicit none
		type(basis_type),pointer,intent(in)::p
		integer,intent(in)::n
		integer::b
		if(associated(p)) then
			b=p%block(n)%init
		else
			b=-9999
		end if
	end function

	!calculate |m>=c^{+}(alpha,s)|n>
	function basis_cplus(p,alpha,s,n) result(m)
		implicit none
		integer,intent(in)::alpha,s,n
		type(basis_type),pointer,intent(in)::p
		integer::m
		
		if(associated(p)) then
			m=p%c_bas(n,alpha,s)
			if(abs(m)<n) m=0
		else
			m=0
		end if
	end function
	
	!calculate |m>=c(alpha,s)|n>
	function basis_c(p,alpha,s,n) result(m)
		implicit none
		integer,intent(in)::alpha,s,n
		type(basis_type),pointer,intent(in)::p
		integer::m
		
		if(associated(p)) then
			m=p%c_bas(n,alpha,s)
			if(abs(m)>n) m=0
		else
			m=0
		end if
	end function
	
	!find # of electrons in n state
	!result: elec = # of electron
	function basis_elec(n,p,a,b,s) result(elec)
		implicit none
		integer,intent(in)::n
		type(basis_type),pointer,intent(in)::p
		integer,intent(in),optional::a,b,s
		integer::elec
		if(associated(p)) then
			if(present(a).and.present(b).and.present(s)) then
				if(s==1) then
					elec=num_elec(p%orbit(n),p%n_orbit,a+p%n_orbit,b+p%n_orbit)
				else
					elec=num_elec(p%orbit(n),p%n_orbit,a,b)
				end if
			else
				elec=num_elec(p%orbit(n),p%n_orbit)
			end if
		else
			elec=-9999
		end if
	end function

	!find spin of n state
	!result: elec = spin*2
	function basis_spin(n,p) result(spin)
		implicit none
		integer,intent(in)::n
		type(basis_type),pointer,intent(in)::p
		integer::spin
		
		if(associated(p)) then
			spin=num_spin(p%orbit(n),p%n_orbit)
		else
			spin=-9999
		end if
	end function

	function basis_double_occupy(n,p) result(d)
		implicit none
		integer,intent(in)::n
		type(basis_type),pointer,intent(in)::p
		integer::d
		if(associated(p)) then
			d=dou_occupy(p%orbit(n),p%n_orbit)
		else
			d=-9999
		end if
	end function
	
	
!---------------------private-----------------------


	!create a fermion at alpha-orbital with s-spin of n-state
	!result: new = -1 (error) or new state
	function creation(alpha,s,bas,n_orbit) result(new)
		implicit none
		integer,intent(in)::alpha,s,bas,n_orbit
		integer::i,j,k,new
		if(occupation(alpha,s,bas,n_orbit).eq.1) then
			new=-1
		else
			if(s==1) then
				new=2**(n_orbit+alpha-1)
			else
				new=2**(alpha-1)
			end if
			new=bas+new
		end if
	end function
	
	!annihilate a fermion at alpha-orbital with s-spin of n-state
	!result: new = -1 (error) or new state
	function annihilation(alpha,s,bas,n_orbit) result(new)
		implicit none
		integer,intent(in)::alpha,s,bas,n_orbit
		integer::i,j,k,new
		if(occupation(alpha,s,bas,n_orbit).eq.0) then
			new=-1
		else
			if(s==1) then
				new=2**(n_orbit+alpha-1)
			else
				new=2**(alpha-1)
			end if
			new=bas-new
		end if
	end function

	
	!find out occupation number in alpha-orbital with s-spin of n-state
	!result: r = 1 (occupied) or 0 (empty)
	function occupation(alpha,s,bas,n_orbit) result(r)
		implicit none
		integer,intent(in)::alpha,bas,s,n_orbit
		integer::i,j,k,r
		k=bas
		if(s==1) k=k/(2**n_orbit)
		do i=1,alpha-1
			k=ishft(k,-1)
		end do
		r=iand(k,1)
	end function


	function dou_occupy(bas,n_orbit) result(n)
		implicit none
		integer,intent(in)::bas,n_orbit
		integer::n,a,b
		
		a=ishft(bas,-n_orbit)
		b=mod(bas,2**n_orbit)
		n=num_elec(iand(a,b),n_orbit)
	end function

	!find out number of fermions in n-state
	!result: n_elec = number of fermions
	function num_elec(bas,n_orbit,a,b) result(n_elec)
		implicit none
		integer,intent(in)::bas
		integer,intent(in),optional::a,b
		integer::i,j,k,l,n_orbit,n_elec,ci,cf
		if(present(a).and.present(b)) then
			if(a<b) then
				ci=a
				cf=b
			else
				ci=b
				cf=a
			end if
		else
			ci=0
			cf=n_orbit*2+1
		end if
		j=bas
		l=0
		do i=1,2*n_orbit
			if(j==0) exit
			if((i>ci).and.(i<cf)) l=l+iand(j,1)
			j=ishft(j,-1)
		end do
		n_elec=l
	end function
	
	!find out spin of n-state
	!result: n_spin = spin * 2
	function num_spin(bas,n_orbit) result(n_spin)
		implicit none
		integer,intent(in)::bas,n_orbit
		integer::i,j,k,l,n_spin
		j=bas
		l=0
		do i=1,n_orbit
			k=2**i
			l=l-mod(j,2)
			j=j/2
		end do
		do i=1,n_orbit
			k=2**i
			l=l+mod(j,2)
			j=j/2
		end do
		n_spin=l
	end function
	

end module bas_mod
