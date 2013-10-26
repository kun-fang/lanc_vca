program hubbard_cluster
	use hubbard_mod
	use rand_mod
	use timer_mod
	use lehmann
	implicit none
	type(basis_type),pointer::bas
	type(hop),pointer::cl
	type(hubbard_cluster_type),pointer::cluster
	type(real_matrix),pointer::M
	type(Q_type),pointer::Q
	integer::n,nsite,norb,r
	real(8)::omega
	real(8),pointer::D(:),X(:,:)
	complex(8),parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	character(len=30)::filename
	call timer_init()
	filename="betts8.conf"
	cl=>cluster_init()
	nsite=cl%site
	norb=nsite
	bas=>basis_init(norb)
	cluster=>hubbard_init(cl,bas)
	r=20
	Q=>lehmann_init(cluster,omega)
	print *,lehmann_potthoff_functional(cluster,Q)
	call hubbard_clean(cluster)
	print *,cputime(),"seconds"
end program hubbard_cluster


