program Dentate
integer n_ec,n_gc,n_gc1,n_gc2,n_hil,n_bc,n_mc,n0,n1,n2
integer Nf
integer inn,ijj
integer, parameter :: nn=1200								! Total number of neurons 
double precision :: dt
double precision, dimension (nn) :: a,b,c,d
double precision, dimension (nn) :: v,u,vn,un
double precision :: rate,Pp,E_po,tau_po
double precision,dimension (nn) :: p_po,g_po,rP,IP
integer, dimension (nn) :: 	    jP
double precision :: kit
double precision, dimension (5,5) :: pr
integer, dimension (nn,nn):: CM
integer, dimension (nn) :: nc_a,nc_n,nc_g
double precision :: tau_n,tau_a,tau_g,g_n,g_a,g_g,E_n,E_a,E_g
double precision, dimension (nn) :: r_a,r_n,r_g
double precision, dimension (nn) :: r_an,r_nn,r_gn
double precision, dimension (nn) :: rsum_a,rsum_n,rsum_g
double precision, dimension (nn,nn) :: p_a,p_n,p_g
double precision, dimension (nn) :: I_a,I_g,I_n,I_tot
double precision :: v_ec,v_gc,v_hil,v_bc,v_mc
double precision :: pcmi
double precision :: rec_a,Iec_a	
double precision :: p_ec_g
double precision :: Inh,Exc,Nmda

!Generating seed for random number generators
ns=1 
call dran_ini(58115+25*ns)
!Time Series parameters
dt=0.01d0 ! time step (ms)	
Nf=700000 ! Length time series

!Number of neurons for each population
n_ec=200							!Total neurons EC (excitatory+inhibitory neurons)
n_ec1=160							!Excitatory neurons EC
n_gc=500							!
n_gc1=n_ec+int(0.65d0*dble(n_gc))	!Burst group
n_gc2=n_gc1+int(0.35d0*dble(n_gc))	!Regular Spiking group
n_hil=100
n_bc=100
n_mc=300
n0=n_ec+n_gc
n1=n0+n_hil
n2=n1+n_bc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!IZHIKEVICH PARAMETERS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Excitatory EC neurons
do inn=1,n_ec1
	a(inn)=0.02d0
	b(inn)=0.2d0
	c(inn)=-60.d0
	d(inn)=8.d0
enddo
!Inhibitory EC neurons
do inn=n_ec1+1,n_ec
	a(inn)=0.1d0
	b(inn)=0.26d0
	c(inn)=-65.d0
	d(inn)=1.d0
enddo
!Bursting GC neurons
do inn=n_ec+1,n_gc1
	a(inn)=0.02d0
	b(inn)=0.2d0
	c(inn)=-50.d0
	d(inn)=2.d0
enddo
!Regular Spiking GC neurons
do inn=n_gc1+1,n_gc2
	a(inn)=0.02d0
	b(inn)=0.2d0
	c(inn)=-60.d0
	d(inn)=8.d0
enddo
!HIL interneurons
do inn=n0+1,n1
	a(inn)=0.1d0
	b(inn)=0.26d0
	c(inn)=-65.d0
	d(inn)=1.d0
enddo
!BC interneurons
do inn=n1+1,n2
	a(inn)=0.1d0
	b(inn)=0.26d0
	c(inn)=-65.d0
	d(inn)=1.d0
enddo
!Mossy Cells
do inn=n2+1,nn
	a(inn)=0.1d0
	b(inn)=0.26d0
	c(inn)=-65.d0
	d(inn)=1.d0
enddo

!Initialization Izhikevih's model
do inn=1,nn
	v(inn)=-70.d0+5.d0*dran_u()
	u(inn)=v(inn)*b(inn)
enddo
vn=0.d0;un=0.d0

! Poisson Noise Parameters
rate=0.08d0
Pp=exp(-rate*dt)*rate*dt

E_po=0d0
tau_po=5.26d0	
p_po=0.d0
g_po=0.d0
jP=0.d0
rP=0.d0

!Maximum neurotransmitter probabilities
do inn=1,n0		!slow synapses
	p_po(inn)=0.1d0	
	g_po(inn)=0.3d0
enddo
do inn=n0+1,nn	!fast synapses
	p_po(inn)=0.1d0		 
	g_po(inn)=0.2d0
enddo	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CONNECTIVITY MATRIX
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!matrix of probabilities: ec-gc-hipp-hil-bc-mc
!matrix reference
!EC:1 - GC:2 - HIL:3 - BC:4 - MC:5
pr=0.d0
!From EC to
pr(1,1)=0.5d0		!EC
pr(1,2)=0.05d0 		!GC
pr(1,4)=0.02d0  	!BC
!From GC to
pr(2,3)=0.02d0		!HIL
pr(2,4)=0.02d0		!BC		
pr(2,5)=0.02d0		!MC
!From HIL to
pr(3,3)=0.3d0		!HIL
pr(3,4)=0.6d0		!BC
pr(3,5)=0.2d0		!MC
!From BC to
pr(4,4)=0.1d0		!BC
pr(4,2)=0.12d0		!GC
pr(4,5)=0.15d0		!MC	
!From MC to
pr(5,5)=0.15d0		!MC
pr(5,2)=0.1d0		!GC
pr(5,3)=0.2d0		!HIL
pr(5,4)=0.2d0		!BC


CM=0
p_ec_g=0.3d0
!Connectivity EC
do inn=1,n_ec1
	!EC
	do ijj=1,n_ec
		if(dran_u().le.pr(1,1)) CM(inn,ijj)=1
	enddo
	!GC
	do ijj=n_ec+1,n0
		if(dran_u().le.pr(1,2)) CM(inn,ijj)=1
	enddo
	!BC
	do ijj=n1+1,n2
		if(dran_u().le.pr(1,4)) CM(inn,ijj)=1
	enddo
enddo
do inn=n_ec1+1,n_ec
	!EC GABA
	do ijj=1,n_ec
		if(dran_u().le.p_ec_g) CM(inn,ijj)=1
	enddo
enddo

!Connectivity GC
do inn=n_ec+1,n0
	!HIL
	do ijj=n0+1,n1
		if(dran_u().le.pr(2,3)) CM(inn,ijj)=1
	enddo
	!BC		
	do ijj=n1+1,n2
		if(dran_u().le.pr(2,4)) CM(inn,ijj)=1
	enddo
	!MC
	do ijj=n2+1,nn
		if(dran_u().le.pr(2,5)) CM(inn,ijj)=1
	enddo
enddo

!Connectivity Hil
do inn=n0+1,n1
	!HIL
	do ijj=n0+1,n1
		if(dran_u().le.pr(3,3)) CM(inn,ijj)=1	
	enddo
	!BC
	do ijj=n1+1,n2
		if(dran_u().le.pr(3,4)) CM(inn,ijj)=1	
	enddo			
	!MC
	do ijj=n2+1,nn
		if(dran_u().le.pr(3,5)) CM(inn,ijj)=1	
	enddo
enddo

!Connectivity BC
do inn=n1+1,n2
	!GC
	do ijj=n_ec+1,n0
		if(dran_u().le.pr(4,2)) CM(inn,ijj)=1	
	enddo
	!BC
	do ijj=n1+1,n2
		if(dran_u().le.pr(4,4)) CM(inn,ijj)=1	
	enddo
	!MC
	do ijj=n2+1,nn
		if(dran_u().le.pr(4,5)) CM(inn,ijj)=1	
	enddo
enddo	

!Connectivity MC
do inn=n2+1,nn
	!GC
	do ijj=n_ec+1,n0
		if(dran_u().le.pr(5,2)) CM(inn,ijj)=1
	enddo
	!HIL
	do ijj=n0+1,n1
		if(dran_u().le.pr(5,3)) CM(inn,ijj)=1
	enddo
	!BC
	do ijj=n1+1,n2
		if(dran_u().le.pr(5,4)) CM(inn,ijj)=1	
	enddo
	!MC
	do ijj=n2+1,nn
		if(dran_u().le.pr(5,5)) CM(inn,ijj)=1
	enddo
enddo

!Removing recurrent connections
do inn=1,nn
	CM(inn,inn)=0
enddo


!Connections count
nc_a=0;nc_n=0;nc_g=0	
do inn=1,nn
	!AMPA
	!Pre-EC
	do ijj=1,n_ec1			
		if(CM(ijj,inn).eq.1) then
			nc_a(inn)=nc_a(inn)+1
		endif		
	enddo
	!Pre-GC
	do ijj=n_ec+1,n0
		if(CM(ijj,inn).eq.1) then
			nc_a(inn)=nc_a(inn)+1
		endif		
	enddo
	!Pre-MC
	do ijj=n2+1,nn
		if(CM(ijj,inn).eq.1) then
			nc_a(inn)=nc_a(inn)+1
		endif		
	enddo			
	!NMDA
	!Pre-EC
	if(inn.gt.n_ec) then
		do ijj=1,n_ec1			
			if(CM(ijj,inn).eq.1) then
				nc_n(inn)=nc_n(inn)+1
			endif		
		enddo
	endif
	!Pre-GC
	do ijj=n_ec+1,n0
		if(CM(ijj,inn).eq.1) then
			nc_n(inn)=nc_n(inn)+1
		endif		
	enddo
	!GABA
	!Pre-HIL,BC
	do ijj=n0+1,n2 
		if(CM(ijj,inn).eq.1) then
			nc_g(inn)=nc_g(inn)+1
		endif	
	enddo
enddo

!Setting to 1 the null connection count
do inn=1,nn
	if(nc_a(inn).eq.0) nc_a(inn)=1
	if(nc_n(inn).eq.0) nc_n(inn)=1
	if(nc_g(inn).eq.0) nc_g(inn)=1
enddo

!Synapses
jt=0
tau_n=200.d0
tau_a=5.6d0
tau_g=5.6d0	
g_n=0.5d0
g_a=0.5d0
g_g=2.d0
E_n=0.d0
E_a=0.d0
E_g=-85.d0
r_a=0.d0
r_n=0.d0
r_g=0.d0
!probability neurotransmitter matrix
p_n=0.d0;p_a=0.d0;p_g=0.d0
!EC	
do inn=1,n_ec1
	!EC	
	do ijj=1,n_ec
		p_a(inn,ijj)=0.5d0
	enddo			
	!GC	
	do ijj=n_ec+1,n0
		!Potenciación LTP->PP
		p_a(inn,ijj)=0.05d0!+0.2d0*0.05d0!0.03d0!0.02d0
		p_n(inn,ijj)=0.05d0!+0.2d0*0.05d0!0.035d0!0.03d0
	
	enddo
	!BC
	do ijj=n1+1,n2
		p_a(inn,ijj)=0.01d0!-0.1d0*0.01d0!0.005d0!0.01d0
		p_n(inn,ijj)=0.01d0!-0.1d0*0.01d0!0.005d0!0.01d0
	enddo
enddo	
!EC-GABA
do inn=n_ec1+1,n_ec
	!EC
	do ijj=1,n_ec
		p_g(inn,ijj)=0.025d0
	enddo
enddo
!GC	
do inn=n_ec+1,n0
	!HIL
	do ijj=n0+1,n1
		p_a(inn,ijj)=0.1d0
		p_n(inn,ijj)=0.1d0
	enddo
	!BC
	do ijj=n1+1,n2
		p_a(inn,ijj)=0.02d0
		p_n(inn,ijj)=0.02d0
	enddo
	!MC
	do ijj=n2+1,nn
		p_a(inn,ijj)=0.01d0
		p_n(inn,ijj)=0.03d0
	enddo
enddo
!HIL	
do inn=n0+1,n1
	!HIL
	do ijj=n0+1,n1
		p_g(inn,ijj)=0.3d0
	enddo
	!BC
	do ijj=n1+1,n2
		p_g(inn,ijj)=0.11d0
	enddo
	!MC
	do ijj=n2+1,nn
		p_g(inn,ijj)=0.4d0
	enddo
enddo	
!BC
do inn=n1+1,n2
	!GC
	do ijj=n_ec+1,n0
		p_g(inn,ijj)=0.4d0+0.1d0*0.4d0
	enddo
	!BC	
	do ijj=n1+1,n2
		p_g(inn,ijj)=0.2d0
	enddo
	!MC
	do ijj=n2+1,nn
		p_g(inn,ijj)=0.3d0
	enddo
enddo
!MC
do inn=n2+1,nn
	!GC
	do ijj=n_ec+1,n0
		p_a(inn,ijj)=0.3d0
	enddo
	!BC	
	do ijj=n1+1,n2
		p_a(inn,ijj)=0.2d0
	enddo
	!HIL
	do ijj=n0+1,n1
		p_a(inn,ijj)=0.4d0
	enddo
	!MC
	do ijj=n2+1,nn
		p_a(inn,ijj)=0.4d0
	enddo
enddo
rec_a=0.d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!DENTATE GYRUS SIMULATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do it=1,Nf
!Currents current
rsum_a=0.d0;rsum_g=0.d0;rsum_n=0.d0

	do inn=1,nn
		if(dran_u().le.Pp) then
			jP(inn)=0
			kit=1.d0
		endif
		rP(inn)=rP(inn)*kit+p_po(inn)*exp(-jP(inn)*dt/tau_po)
		IP(inn)=-g_po(inn)*rP(inn)*(v(inn)-E_po)
		jP(inn)=jP(inn)+1
		kit=0.d0
		do ijj=1,nn
			rsum_a(inn)=rsum_a(inn)+p_a(ijj,inn)*r_a(ijj)*dble(CM(ijj,inn))
			rsum_n(inn)=rsum_n(inn)+p_n(ijj,inn)*r_n(ijj)*dble(CM(ijj,inn))		
			rsum_g(inn)=rsum_g(inn)+p_g(ijj,inn)*r_g(ijj)*dble(CM(ijj,inn))
		enddo
		if(inn.le.n_ec1) then
			I_a(inn)=-g_a*rsum_a(inn)*(v(inn)-E_a)/dble(nc_a(inn))+4.d0
		else
			I_a(inn)=-g_a*rsum_a(inn)*(v(inn)-E_a)/dble(nc_a(inn))
		endif
		I_n(inn)=-g_n*rsum_n(inn)*(v(inn)-E_n)/dble(nc_n(inn))
		I_g(inn)=-g_g*rsum_g(inn)*(v(inn)-E_g)/dble(nc_g(inn))			
	enddo
	
	!print*, IP(1)
	I_tot=IP+I_a+I_g+I_n
!Dynamic Neuron


	call Runge_AMPA(r_a,tau_a,dt,r_an,nn)
	call Runge_NMDA(r_n,tau_n,dt,r_nn,nn)
	call Runge_GABA(r_g,tau_g,dt,r_gn,nn)
	r_a=r_an;r_g=r_gn;r_n=r_nn


	call Runge_Izhikevich(v,u,a,b,I_tot,dt,vn,un,nn)

	do inn=1,nn
		if(vn(inn).ge.30.d0) then
					write(100,*) dble(it)*dt,inn
					v(inn)=30.d0
			vn(inn)=c(inn)
			un(inn)=un(inn)+d(inn)
			r_a(inn)=r_a(inn)+1.d0
			r_n(inn)=r_n(inn)+1.d0
			r_g(inn)=r_g(inn)+1.d0
		endif
	enddo
!LFP-Write
	if(it.gt.200000) then

		v_gc=0.d0;v_hil=0.d0;v_bc=0.d0;v_mc=0.d0;v_ec=0.d0
!LFP		
		!EC
		do inn=1,n_ec1
			v_ec=v_ec+v(inn)/dble(n_ec1)

		enddo
		write(200,*) v_ec
		!GC
		do inn=n_ec+1,n0
			v_gc=v_gc+v(inn)/dble(n_gc)
		enddo
		write(300,*) v_gc
		!HIL
		do inn=n0+1,n1
			v_hil=v_hil+v(inn)/dble(n_hil)
		enddo
		write(400,*) v_hil
		!BC
		do inn=n1+1,n2
			v_bc=v_bc+v(inn)/dble(n_bc)
		enddo
		write(500,*) v_bc
		!MC
		do inn=n2+1,nn
			v_mc=v_mc+v(inn)/dble(n_mc)
		enddo
		write(600,*) v_mc
!Currents
		!IPSC EPSC
		Inh=0.d0
		Exc=0.d0
		Nmda=0.d0
		do inn=n_ec+1,n0
			Inh=Inh+I_g(inn)/dble(n_gc)
			Exc=Exc+I_a(inn)/dble(n_gc)
			Nmda=Nmda+I_n(inn)/dble(n_gc)
		enddo
		write(700,*) Inh,Exc,Nmda

	endif
	
!Update Neuron
	v=vn
	u=un
	
enddo


end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!RUNGE-KUTTA SUBROUTINES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
subroutine Runge_Izhikevich(v,u,a,b,I,dt,vn,un,nn)
integer :: nn,inn
double precision :: v(nn),u(nn),vn(nn),un(nn)
double precision :: a(nn),b(nn),I(nn)
double precision :: dt
double precision :: k1(nn),k2(nn),k3(nn),k4(nn)
double precision :: l1(nn),l2(nn),l3(nn),l4(nn)

k1=0.d0;k2=0.d0;k3=0.d0;k4=0.d0	
l1=0.d0;l2=0.d0;l3=0.d0;l4=0.d0
vn=0.d0;un=0.d0

do inn=1,nn
	k1(inn)=0.04d0*v(inn)**2+5.d0*v(inn)-u(inn)+140.d0+I(inn)
	l1(inn)=a(inn)*(b(inn)*v(inn)-u(inn))
	k2(inn)=0.04d0*(v(inn)+k1(inn)/2.d0)**2+5.d0*(v(inn)+k1(inn)/2.d0)-u(inn)-l1(inn)/2.d0+140.d0+I(inn)
	l2(inn)=a(inn)*(b(inn)*(v(inn)+k1(inn)/2.d0)-u(inn)-l1(inn)/2.d0)
	k3(inn)=0.04d0*(v(inn)+k2(inn)/2.d0)**2+5.d0*(v(inn)+k2(inn)/2.d0)-u(inn)-l2(inn)/2.d0+140.d0+I(inn)
	l3(inn)=a(inn)*(b(inn)*(v(inn)+k2(inn)/2.d0)-u(inn)-l2(inn)/2.d0)
	k4(inn)=0.04d0*(v(inn)+k3(inn))**2+5.d0*(v(inn)+k3(inn))-u(inn)-l3(inn)+140.d0+I(inn)
	l4(inn)=a(inn)*(b(inn)*(v(inn)+k3(inn))-u(inn)-l3(inn))
	vn(inn)=v(inn)+dt/6.d0*(k1(inn)+2.d0*k2(inn)+2.d0*k3(inn)+k4(inn))
	un(inn)=u(inn)+dt/6.d0*(l1(inn)+2.d0*l2(inn)+2.d0*l3(inn)+l4(inn))
enddo

end subroutine

subroutine Runge_AMPA(r_a,tau_a,dt,r_an,nn)
integer :: nn,inn
double precision, dimension (nn) :: r_a,r_an
double precision, dimension (nn) :: k1,k2,k3,k4
double precision :: dt,tau_a
k1=0.d0;k2=0.d0;k3=0.d0;k4=0.d0	

r_an=0.d0;

do inn=1,nn
	k1(inn)=-(r_a(inn))/tau_a
	k2(inn)=-(r_a(inn)+k1(inn)/2.d0)/tau_a
	k3(inn)=-(r_a(inn)+k2(inn)/2.d0)/tau_a
	k4(inn)=-(r_a(inn)+k3(inn))/tau_a
	r_an(inn)=r_a(inn)+dt/6.d0*(k1(inn)+2.d0*k2(inn)+2.d0*k3(inn)+k4(inn))
enddo
end subroutine


subroutine Runge_NMDA(r_n,tau_n,dt,r_nn,nn)
integer :: nn,inn
double precision, dimension (nn) :: r_n,r_nn
double precision, dimension (nn) :: k1,k2,k3,k4
double precision :: dt,tau_n
k1=0.d0;k2=0.d0;k3=0.d0;k4=0.d0	

r_nn=0.d0;

do inn=1,nn
	k1(inn)=-(r_n(inn))/tau_n
	k2(inn)=-(r_n(inn)+k1(inn)/2.d0)/tau_n
	k3(inn)=-(r_n(inn)+k2(inn)/2.d0)/tau_n
	k4(inn)=-(r_n(inn)+k3(inn))/tau_n
	r_nn(inn)=r_n(inn)+dt/6.d0*(k1(inn)+2.d0*k2(inn)+2.d0*k3(inn)+k4(inn))
enddo
end subroutine

	
subroutine Runge_GABA(r_g,tau_g,dt,r_gn,nn)
integer :: nn,inn
double precision, dimension (nn) :: r_g,r_gn
double precision, dimension (nn) :: k1,k2,k3,k4
double precision :: dt,tau_g
k1=0.d0;k2=0.d0;k3=0.d0;k4=0.d0	

r_gn=0.d0;

do inn=1,nn
	k1(inn)=-(r_g(inn))/tau_g
	k2(inn)=-(r_g(inn)+k1(inn)/2.d0)/tau_g
	k3(inn)=-(r_g(inn)+k2(inn)/2.d0)/tau_g
	k4(inn)=-(r_g(inn)+k3(inn))/tau_g
	r_gn(inn)=r_g(inn)+dt/6.d0*(k1(inn)+2.d0*k2(inn)+2.d0*k3(inn)+k4(inn))
enddo
end subroutine

subroutine Runge_input_gc(r_in_gc,tau_in_gc,dt,r_in_gcn)
double precision :: r_in_gc,r_in_gcn
double precision :: k1,k2,k3,k4
double precision :: dt,tau_in_gc
k1=0.d0;k2=0.d0;k3=0.d0;k4=0.d0	

r_in_gcn=0.d0;
	k1=-(r_in_gc)/tau_in_gc
	k2=-(r_in_gc+k1/2.d0)/tau_in_gc
	k3=-(r_in_gc+k2/2.d0)/tau_in_gc
	k4=-(r_in_gc+k3)/tau_in_gc
	r_in_gcn=r_in_gc+dt/6.d0*(k1+2.d0*k2+2.d0*k3+k4)
end subroutine





