module VCA
	use hubbard_mod
	use lehmann
	!use optimal
	!use cmatrix
	!use MCDOS
	!use cmatrix
	complex(kind=8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	real(kind=8),private,parameter::eps=0.001
	integer,private,parameter::nk=20


	contains



	function VCA_potthoff_functional(cluster) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q
		real(kind=8)::omega,xup,xdn,x
		integer::orbit
		if(.not.associated(cluster)) return
		orbit=cluster%nsite
		!---> spin up
		Q=>lehmann_init(cluster,omega)
		if(.not.associated(Q)) return
		x=(omega+lehmann_potthoff_functional(cluster,Q)+cluster%hop%shift)/orbit
		!write(*,'(4F20.15)') cluster%hop%ct,cluster%hop%cmu,x
		call lehmann_Q_clean(Q)
	end function



end module VCA
