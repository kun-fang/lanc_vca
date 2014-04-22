module hubbard_mod
	use hop_mod
	use bas_mod
	use matrix_mod
	use rand_mod
	use underwood_mod
	implicit none

	type hubbard_cluster_type
		integer::nsite=0,nbas=0
		type(hop),pointer::hop=>NULL()
		type(basis_type),pointer::bas=>NULL()
	end type hubbard_cluster_type

	type(real_matrix),private,pointer::comm_matrix=>NULL()
  logical,private,parameter::checkEnergy=.false.
	real(8),private,parameter::ediff=3.5d0
  real(8),private,parameter::epsilon=1.d-5
	

contains

!-----------------hubbard cluster------------------------
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
	
	function cluster_Hamiltonian_Partition(cl,bas,si,sf) result(M)
		type(hop),pointer::cl
		type(basis_type),pointer::bas
    integer::si,sf
		type(real_matrix),pointer::M
		integer::i,j,r,s,n,spin,spin2,site,k1,k2,alpha,beta,mr,ms
		complex(8),pointer::t(:,:)
		real(8)::x
		n=sf-si+1
		M=>real_matrix_init(n)
		site=cl%site
		t=>cluster_update(cl)
		do r=1,n
      mr=r+si-1
			x=real_matrix_show(M,r,r)+basis_double_occupy(mr,bas)*cl%U*(-1)**(cl%nambu)
			call real_matrix_set(M,r,r,x)
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
								x=real_matrix_show(M,r,s)+k1*k2*t(alpha+(spin-1)*site,beta+(spin2-1)*site)
								call real_matrix_set(M,r,s,x)
								call real_matrix_set(M,s,r,x)
							end if
						end do
					end do
				end do
			end do
		end do
		deallocate(t)
	end function

	function cluster_Hamiltonian(cl,bas) result(M)
		type(hop),pointer::cl
		type(basis_type),pointer::bas
		type(real_matrix),pointer::M
		integer::i,j,r,s,n,spin,spin2,site,k1,k2,alpha,beta
		complex(8),pointer::t(:,:)
		real(8)::x
		n=basis_get_n_basis(bas)
		M=>real_matrix_init(n)
		site=cl%site
		t=>cluster_update(cl)
		do r=1,n
			x=real_matrix_show(M,r,r)+basis_double_occupy(r,bas)*cl%U*(-1)**(cl%nambu)
			call real_matrix_set(M,r,r,x)
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
								x=real_matrix_show(M,r,s)+k1*k2*t(alpha+(spin-1)*site,beta+(spin2-1)*site)
								call real_matrix_set(M,r,s,x)
								call real_matrix_set(M,s,r,x)
							end if
						end do
					end do
				end do
			end do
		end do
		deallocate(t)
	end function

	subroutine cluster_connect_matrix(M)
		type(real_matrix),pointer::M
		comm_matrix=>M
	end subroutine

	subroutine cluster_disconnect_matrix()
		comm_matrix=>NULL()
	end subroutine

	subroutine cluster_solver(cluster,m,D,X)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(8),pointer::D(:),X(:,:)
		type(real_matrix),pointer::H
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
		call real_matrix_del(H)
	end subroutine

	subroutine cluster_solver_elec(cluster,elec,m,D,X,dim)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(8),pointer::D(:),X(:,:)
    type(basis_type),pointer::bas
		type(real_matrix),pointer::H
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
		call real_matrix_del(H)
	end subroutine

	subroutine cluster_solver_elec_spin(cluster,elec,spin,m,D,X,dim)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(8),pointer::D(:),X(:,:)
    type(basis_type),pointer::bas
		type(real_matrix),pointer::H
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
		call real_matrix_del(H)
	end subroutine

	subroutine OP(n,x,y)
		integer::n
		real(8)::x(n),y(n)
		call real_matrix_product(comm_matrix,x,y)
	end subroutine

	function cluster_get_n_orbit(cluster) result(n)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		integer::n
		n=cluster%nsite
	end function
	
	
end module hubbard_mod
