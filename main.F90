program hubbard_cluster
	use hubbard_mod
	use rand_mod
	use timer_mod
	use VCA
	implicit none
	type(basis_type),pointer::bas
	type(hop),pointer::cl
	type(hubbard_cluster_type),pointer::cluster
	type(real_matrix),pointer::M
	type(Q_type),pointer::Q
	integer::n,nsite,norb,i
	real(8)::omega,density
	complex(8),parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	character(len=30)::filename
	call timer_init()
	filename="betts8.conf"
	cl=>cluster_init()
	nsite=cl%site
	norb=nsite
	bas=>basis_init(norb)
	cluster=>hubbard_init(cl,bas)
	write(*,*)"----------------------------------------------------------------"
	write(*,*)"VCA calculation: initial parameters are"
	write(*,*)"U=",cluster%hop%U,"mu=",cluster%hop%lmu
	write(*,*)"t=",cluster%hop%lt
	write(*,*)"Cluster t=",cluster%hop%ct
	write(*,*)"Cluster mu=",cluster%hop%cmu
	write(*,*)"----------------------------------------------------------------"
	!call VCA_fermi_surface(cluster)
	if(VCA_connect_cluster(cluster)) then
		call VCA_optimal(omega,density)
		!do i=0,10
		!	print *,cluster%hop%M,VCA_potthoff_functional()
		!	cluster%hop%M=cluster%hop%M+0.03
		!end do
		!call VCA_spectral_function()
		!call VCA_DOS(-3.d0,3.d0,40)
		!call VCA_fermi_surface()
		call VCA_disconnect_cluster()
	end if
	call hubbard_clean(cluster)
	print *,cputime(),"seconds"
end program hubbard_cluster


