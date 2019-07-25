module eqcalc
  implicit none

  real, parameter   :: pi = 3.141592653589793238462643383279502884197

contains


!! Define kappa-shape profile
real function k(rho)  result(y)
  real, intent(in)  :: rho  !input
  y   = 1. !+ (rho**2.0) / 100.0
end function k

!! Define derivative of kappa-shape profile
real function dk(rho)  result(y)
  real, intent(in)  :: rho   !input
  y   = 0. !+ rho / 50.0
end function dk

!! Define heat deposition shape profile
real function exheat(rho)     result(y)
  real, intent(in)  :: rho !input
  y   = 1 - rho**2
end function exheat

!! int( rho*p(rho)*drho,rho)/rho
real function intexheat(rho)     result(y)
  real, intent(in)  :: rho !input
  y   = 0.5*(rho - 0.5*rho**3.0)
end function intexheat

!! Define the q-profile
real function qprof(rho,epsilon)     result(y)
  real, intent(in)  :: rho, epsilon !input
  y   = 10 * 2.1 * epsilon / (2 - 1.0*rho**2.0)
end function qprof

!! Define the viscosity profile
real function m(rho)     result(y)
  real, intent(in)  :: rho !input
  y   = 1
end function m

!! Define the derivative of the viscosity profile
real function dm(rho)    result(y)
  real, intent(in)  :: rho !input
  y   = 0
end function dm

!! Define the momentum input profile
real function f_m(rho)     result(y)
  real, intent(in)  :: rho !input
  real              :: a_div_sigma
  a_div_sigma   = 5.
  y   = exp( -0.5 * a_div_sigma**2.0 * rho**2.0)
end function f_m









!! Routine that calculates rho0 as a function of gamma
real function rho0(gamma,beta,drho,t_end) result(y)
  real, intent(in)  :: gamma, beta, drho, t_end
  real              :: rho
  real, dimension(2):: t_bucket, g_bucket
  integer           :: rho_idx

  !! Initialize from boundary conditions
  g_bucket(1)     =   0
  g_bucket(2)     =   -(gamma + gamma*beta*exheat(0.0))/(2*k(0.0))*drho
  t_bucket(1)     =   1
  t_bucket(2)     =   1
  rho_idx         =   3
  !! Loop through all nodes untill the temperature is negative
  do while (t_bucket(2) .gt. t_end)
    rho           = (rho_idx-1)*drho
    t_bucket(1)   = t_bucket(2)
    g_bucket(1)   = g_bucket(2)
    g_bucket(2)   = g_bucket(1) - 1/k(rho-drho) * ( gamma*t_bucket(1)**1.5 + &
                    gamma*beta*exheat(rho-drho) + dk(rho-drho)*g_bucket(1) + &
                    g_bucket(1)/(rho-drho))*drho
    t_bucket(2)   = t_bucket(1) + g_bucket(1)*drho
    rho_idx       = rho_idx + 1
  end do
  !! Return the final coordinate
  y             = rho
end function rho0

!! Routine that finds rho0(gamma_final)=1
real function gamma_final(beta,drho,t_end) result(y)
  real, intent(in)  :: beta, drho,t_end
  real, dimension(2):: gamma0, rhoarr

  !! Initialize two gamma guesses and calculate corresponding rho0
  gamma0        = (/ 1, 10 /)
  rhoarr(1)     = rho0(gamma0(1),beta,drho,t_end)
  rhoarr(2)     = rho0(gamma0(2),beta,drho,t_end)

  !! Go through Newton-Rhpason method to approximate root.
  do while (abs(rhoarr(2)-1) .gt. 1E-6)
    !! Abs as gamma > 0
    gamma0(1) = abs(gamma0(1) - (rhoarr(1) - 1) * (gamma0(2) - gamma0(1)) / (rhoarr(2) - rhoarr(1)))
    gamma0    = gamma0(2:1:-1)

    rhoarr(1)   = rho0(gamma0(2),beta,drho,t_end)
    rhoarr      = rhoarr(2:1:-1)
  end do
  ! Return final gamma
  y             = gamma0(2)
end function gamma_final

!! Returns t'(1)
real function dtdrhoa(gamma,beta,drho,t_end) result(y)
  real, intent(in)  :: gamma, beta, drho,t_end
  real              :: rho
  real, dimension(2):: t_bucket, g_bucket
  integer           :: rho_idx

  ! Initialize
  g_bucket(1)     =   0
  g_bucket(2)     =   -(gamma + beta*exheat(0.0))/(2*k(0.0))*drho
  t_bucket(1)     =   1
  t_bucket(2)     =   1
  rho_idx         =   3

  ! Loop over all nodes
  do while (t_bucket(2) .gt. t_end)
    rho           = (rho_idx-1)*drho
    t_bucket(1)   = t_bucket(2)
    g_bucket(1)   = g_bucket(2)
    g_bucket(2)   = g_bucket(1) - ( gamma*(t_bucket(1)**1.5 + &
                    beta*exheat(rho))/k(rho) + &
                    (1/rho + dk(rho)/k(rho))*g_bucket(1) )*drho
    t_bucket(2)   = t_bucket(1) + g_bucket(1)*drho
    rho_idx       = rho_idx + 1
  end do
  ! return final slope
  y = (t_bucket(2)-t_bucket(1))/drho
end function dtdrhoa

!! Routine that calculates p_edge as a function of beta_p
real function p_edge(betap,epsilon,drho,btheta_prof,bz_prof,n_nodes) result(y)
  real, intent(in)                    :: betap, epsilon, drho
  integer, intent(in)                 :: n_nodes
  real, dimension(n_nodes), intent(in):: btheta_prof, bz_prof
  real                                :: rho, dbtheta2drho, dbz2drho, dpdrho
  real, dimension(n_nodes)            :: p_prof
  integer                             :: rho_idx

  !! First two by hand
  p_prof(1)       = 1
  p_prof(2)       = 1

  !! Calculate the profile
  do rho_idx = 3, n_nodes
    rho             = (rho_idx-1)*drho
    dbtheta2drho    = (btheta_prof(rho_idx)**2.0 - btheta_prof(rho_idx-1)**2.0)/drho
    dbz2drho        = ( bz_prof(rho_idx)**2.0 - bz_prof(rho_idx-1)**2.0 ) / drho

    dpdrho          = -betap*(dbtheta2drho + dbz2drho/epsilon**2.0 + &
                      2.0*btheta_prof(rho_idx)**2.0/rho)
    p_prof(rho_idx) = p_prof(rho_idx-1) + dpdrho*drho
  end do
  !! return final pressure
  y             = p_prof(n_nodes)
end function p_edge

!! Routine that solves n_edge(betap) = 0 for betap
real function beta_p0(epsilon,drho,btheta_prof,bz_prof,n_nodes,p_inf) result(y)
  real, intent(in)        :: epsilon, drho, p_inf
  integer, intent(in)     :: n_nodes
  real, dimension(n_nodes), intent(in):: btheta_prof, bz_prof
  real, dimension(2)      :: betap_arr, p_edge_arr

  ! Make two guesses
  betap_arr     = (/ 1.0E-15, 1.0E-12  /)
  p_edge_arr(1) = p_edge(betap_arr(1),epsilon,drho,btheta_prof,bz_prof,n_nodes)
  p_edge_arr(2) = p_edge(betap_arr(2),epsilon,drho,btheta_prof,bz_prof,n_nodes)

  ! Newton Rhapson
  do while (abs(p_edge_arr(2)-p_inf) .gt. 1E-6)
    betap_arr(1)  = betap_arr(1) - ( p_edge_arr(1) - p_inf) * ( betap_arr(2) -&
                    betap_arr(1) ) / ( p_edge_arr(2) - p_edge_arr(1) + 1E-8 )
    betap_arr     = betap_arr(2:1:-1)

    p_edge_arr(1) = p_edge(betap_arr(2),epsilon,drho,btheta_prof,bz_prof,n_nodes)
    p_edge_arr    = p_edge_arr(2:1:-1)
  end do

  ! Return final beta
  y             =   betap_arr(2)
end function beta_p0

!! Routine that calculates v_end as a function of alpha
real function v_end(alpha,drho) result(y)
  real, intent(in)  :: alpha, drho
  real              :: rho, dhdrho
  real, dimension(:), allocatable :: v_arr, h_arr
  integer           :: rho_idx, n_nodes

  n_nodes = ceiling(1/drho)+1
  allocate ( v_arr(n_nodes) )
  allocate ( h_arr(n_nodes) )

  !! First two by hand
  h_arr(1)        = 0
  v_arr(1)        = 1
  h_arr(2)        = h_arr(1) - alpha*f_m(0.0)/2.0 * drho
  v_arr(2)        = v_arr(1) + h_arr(1) * drho


  !! Calculate the profile
  do rho_idx = 3, n_nodes
    rho             = (rho_idx-1)*drho
    dhdrho          = -( alpha*f_m(rho)/m(rho) + h_arr(rho_idx-1)/rho + dm(rho)/m(rho)*h_arr(rho_idx-1) )
    h_arr(rho_idx)  = h_arr(rho_idx-1) + dhdrho * drho
    v_arr(rho_idx)  = v_arr(rho_idx-1) + h_arr(rho_idx-1) * drho
  end do

  !! Return v_end
  y             = v_arr(n_nodes)
end function v_end

!! Routine that solves v_end(alpha) = val. for alpha
real function alpha_fin(drho,v_fin) result(y)
  real, intent(in)  :: drho, v_fin
  real, dimension(2):: alpha_arr, v_end_arr


  ! Make two guesses
  alpha_arr     = (/ 1.0E-1, 1.0E0  /)
  v_end_arr(1)  = v_end(alpha_arr(1),drho)
  v_end_arr(2)  = v_end(alpha_arr(2),drho)

  ! Newton Rhapson
  do while (abs(v_end_arr(2)-v_fin) .gt. 1E-6)
    alpha_arr(1)  = alpha_arr(1) - ( v_end_arr(1) - v_fin) * ( alpha_arr(2) -&
                    alpha_arr(1) ) / ( v_end_arr(2) - v_end_arr(1) )
    alpha_arr     = alpha_arr(2:1:-1)

    v_end_arr(1) = v_end(alpha_arr(2),drho)
    v_end_arr    = v_end_arr(2:1:-1)
  end do

  ! Return final beta
  y             =   alpha_arr(2)
end function alpha_fin




!! Calculates all profiles
subroutine all_profs(drho,n_nodes,t_prof,g_prof,beta,gamma,btheta_prof,p_prof,&
     betap,epsilon,bz_arr,rho_arr,b_flat,t_end,p_inf,alpha,v_arr,v_fin)
  real, intent(in)        :: drho, beta, gamma, epsilon, t_end, p_inf
  integer, intent(in)     :: n_nodes, b_flat
  real                    :: rho, betap, dbtheta2drho, dbz2drho, dpdrho, alpha, &
                             dhdrho, v_fin
  real, dimension(n_nodes):: t_prof, g_prof, btheta_prof, k_arr, &
                             p_prof, intexheat_arr, q_arr, rho_arr, bz_arr, &
                             v_arr, h_arr
  integer                 :: rho_idx

  !! Do first two steps by hand due to limits
  rho_idx         = 1
  t_prof(rho_idx) = 1
  g_prof(rho_idx) = 0
  rho_idx         = rho_idx + 1
  t_prof(rho_idx) = 1
  g_prof(rho_idx) = -(gamma + gamma*beta*exheat(0.0))/(2*k(0.0))*drho

  !! Loop over all nodes for t and t' profiles
  do rho_idx = 3, n_nodes
    rho             = (rho_idx-1)*drho
    g_prof(rho_idx) = g_prof(rho_idx-1) - 1/k(rho-drho) * ( gamma*t_prof(rho_idx-1)**1.5 + &
                      gamma*beta*exheat(rho-drho) + dk(rho-drho)*g_prof(rho_idx-1) + &
                      g_prof(rho_idx-1)/(rho-drho))*drho
    t_prof(rho_idx) = t_prof(rho_idx-1) + g_prof(rho_idx-1)*drho
  end do

  !! Calculate btheta profile
  !! Convert kappa and intexheat in array
  rho_idx         = 1
  do rho_idx = 1, n_nodes
    rho                   = (rho_idx-1)*drho
    rho_arr(rho_idx)      = rho
    k_arr(rho_idx)        = k(rho)
    intexheat_arr(rho_idx)= intexheat(rho)
    q_arr(rho_idx)        = qprof(rho, epsilon)
  end do
  btheta_prof   = (beta * gamma * intexheat_arr    + k_arr * g_prof          )/ &
                  (beta * gamma * intexheat(1.0)   + k(1.0)* g_prof(n_nodes) )


  !! Calculate the number density profile

  !! First two by hand
  p_prof(1)       = 1
  p_prof(2)       = 1

  !! make bz profile and adjust first few indices
  bz_arr          = q_arr * btheta_prof / rho_arr
  bz_arr(1)       = qprof(0.0,epsilon)    * ( btheta_prof(2) - btheta_prof(1) ) / drho
  bz_arr          = bz_arr*b_flat
  betap = beta_p0(epsilon,drho,btheta_prof,bz_arr,n_nodes,p_inf)

  !! The beginning is quite sensitive to noise, so we rewrite this as a linear slope


  !! Calculate the profile
  do rho_idx = 3, n_nodes
    rho             = (rho_idx-1)*drho
    dbtheta2drho    = (btheta_prof(rho_idx)**2.0 - btheta_prof(rho_idx-1)**2.0)/drho
    dbz2drho        = (bz_arr(rho_idx)**2.0  - bz_arr(rho_idx-1)**2.0)/drho
    dpdrho          = -betap*(dbtheta2drho  + dbz2drho/(epsilon**2.0) + &
                      2.0*btheta_prof(rho_idx)**2.0/rho)
    p_prof(rho_idx) = p_prof(rho_idx-1) + dpdrho*drho
  end do

  !! Calculate the velocity profile in the axial direction
  !! Find dimensionless constant alpha for given profile and v_end
  alpha = alpha_fin(drho,v_fin)

  !! Calculate profiles using forward Euler.
  !! First two by hand
  h_arr(1)        = 0
  v_arr(1)        = 1
  h_arr(2)        = h_arr(1) - alpha*f_m(0.0)/2.0 * drho
  v_arr(2)        = v_arr(1) + h_arr(1) * drho


  !! Calculate the profile
  do rho_idx = 3, n_nodes
    rho             = (rho_idx-1)*drho
    dhdrho          = -( alpha*f_m(rho)/m(rho) + h_arr(rho_idx-1)/rho + dm(rho)/m(rho)*h_arr(rho_idx-1) )
    h_arr(rho_idx)  = h_arr(rho_idx-1) + dhdrho * drho
    v_arr(rho_idx)  = v_arr(rho_idx-1) + h_arr(rho_idx-1) * drho
  end do


end subroutine all_profs


end module eqcalc
















program eqcalc_2

  use eqcalc
  implicit none
  real              :: drho, beta, gamma,t1, betap, b_zval, t_end, p_inf
  real              :: eta0, a, Ip, kappa0, R0, T0, VL, mu0, epsilon, &
                       B1, p_tot, p0, eta_fac, ln_lambda, z_ion, mu_fac
  integer           :: b_flat
  real, dimension(:), allocatable  :: t_prof, g_prof, btheta_prof,rho_array, p_prof, bz_arr, rho_arr, v_arr
  real,dimension(:,:),allocatable  :: profile_j, profile_p, profile_ne, profile_f, profile_q, profile_t, profile_omega
  integer           :: n_nodes, write_idx
  real              :: v_fin, alpha, F0, mu_v, q0, F0_real, mu_v_real
  namelist /mackenbach_input/ ln_lambda, z_ion, drho, beta, t_end, p_inf, eta_fac, mu_fac, a, Ip, kappa0, R0, b_flat, v_fin, F0_real, mu_v_real, q0

  open(12,file='mackenbach_input')
  read(12,nml=mackenbach_input)
  write(*,nml=mackenbach_input)
  
  !! Declare main varibales
  !ln_lambda = 17.
  !z_ion   = 1.
  !drho    = 1E-5    ! step size in rho space
  !beta    = 0       ! fraction of non-ohmic heating
  !t_end   = 5e-2    ! Fraction of T0 at edge
  !p_inf   = 1e-3    ! Fraction of p0 at the edge
  gamma   = gamma_final(beta, drho, t_end) ! the value of dimless const. gamma
  t1      = abs(dtdrhoa(gamma,beta,drho,t_end)) ! t profile derivative at t1
  p_tot   = intexheat(1.0)

  !! Statement to check if t_profile gets ohmic heating (necessary for current constraint)
  IF ( beta .gt. abs( k(1.0)*t1/( gamma*p_tot ) ))  THEN
   write(*,*) 'Beta too large. Decrease beta.'
   stop
  end if


  !eta_fac = 100.
  !mu_fac  = 10.
!  eta0    = 28.05E-9 * (1.6022e-19)**1.5 * eta_fac ! Prefactor for spitzer conductivity
  eta0    = 1.03e-4 * z_ion * ln_lambda * (1.6022e-19)**1.5 * eta_fac
  !a       = 1.0       ! Minor radius
  !Ip      = 3E7       ! Total plasma current
  !kappa0  = 1E21      ! Heat diffusivity
  !R0      = 10.0      ! Major radius
  mu0     = 4*pi*1E-7 ! Vacuum permeability
  B1      = mu0*Ip/(2*pi*a) ! Strength of poloidal magnetic field at edge
  epsilon = a/R0      ! Minor over Major radius
  !b_flat  = 0         ! Parameter indiciting if the pressure will be calculated
                      ! consistently via the q profile (b_flat = 1),
                      ! or if one will assume a flat b_z profile (b_flat = 0).
  !v_fin   = 0.0       ! Edge value of velocity profile.
  F0      = mu_fac*F0_real    ! Total toroidal force (in Newtons/m3)
  mu_v    = mu_fac*mu_v_real      ! viscosity value at centre

  !! Allocate memory for profiles
  n_nodes = ceiling(1/drho)+1
  allocate ( t_prof(n_nodes) )
  allocate ( g_prof(n_nodes) )
  allocate ( btheta_prof(n_nodes) )
  allocate ( p_prof(n_nodes) )
  allocate ( rho_array(n_nodes) )
  allocate ( bz_arr(n_nodes) )
  allocate ( rho_arr(n_nodes) )
  allocate ( v_arr(n_nodes) )
  call all_profs(drho,n_nodes,t_prof,g_prof,beta,gamma,btheta_prof,p_prof,betap, &
                  epsilon,bz_arr,rho_arr,b_flat,t_end,p_inf,alpha,v_arr,v_fin)

  !! Calculate dimensionfull quantities
  T0 =  (eta0**0.4*gamma**0.4*Ip**0.8)/(a**0.8*kappa0**0.4*(2*pi)**0.8*(abs(beta*gamma*p_tot - k(1.0)*t1))**0.8)
  VL =  (2*2**0.2*eta0**0.4*gamma**0.4*kappa0**0.6*pi**1.2*R0*(abs(beta*gamma*p_tot - k(1.0)*t1))**0.2)/(a**0.8*Ip**0.2)
  p0 =  mu0*Ip**2.0/(8.0 * pi**2.0 * a**2.0 * betap )
  !q0      = 2.0     !! Centre value for q
  b_zval  = q0*B1*(btheta_prof(2) - btheta_prof(1))/epsilon/drho
            ! Value for const. magnetic field such that q(0) = q0


  !! Make to be written matrices
  allocate ( profile_j(2,n_nodes) )
  profile_j(1,:) = a*rho_arr
  profile_j(2,:) = VL*T0**(1.5)/(2.0*pi*R0*eta0)*t_prof**(1.5)

  allocate ( profile_p(2,n_nodes) )
  profile_p(1,:) = a*rho_arr
  profile_p(2,:) = p_prof*p0 + p_inf

  allocate ( profile_ne(2,n_nodes) )
  profile_ne(1,:) = a*rho_arr
  profile_ne(2,:) = p_prof/(2.0*t_prof) * p0/T0

  allocate ( profile_t(2,n_nodes) )
  profile_t(1,:) = a*rho_arr
  profile_t(2,:) = t_prof * T0 / (1.6022e-19) / 1e3

  allocate ( profile_f(2,n_nodes) )
  profile_f(1,:)  = a*rho_arr
  profile_f(2,:)  = B1*bz_arr/epsilon
  if ( b_flat  == 0 ) then
    profile_f(2,:) = (rho_arr+1)/(rho_arr+1) * b_zval
  end if

  allocate ( profile_q(2,n_nodes) )
  profile_q(1,:)  = a*rho_arr
  profile_q(2,:)  = rho_arr*bz_arr/(epsilon*R0*btheta_prof)
  if ( b_flat  == 0 ) then
    profile_q(2,:) = a*rho_arr*b_zval/(R0*btheta_prof*B1)
  end if

  allocate ( profile_omega(2,n_nodes) )
  profile_omega(1,:)  = a*rho_arr
  profile_omega(2,:)  = v_arr * ( F0 * a**2.0 /(mu_v * alpha) )
  print *, 'Core rotation (m/s) :', profile_omega(2,1)
  ! convert from m/s to krad/s
  profile_omega(2,:) = profile_omega(2,:) / R0 / 1e3


  !! Write output
  open (unit=1,file="profile_j",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(1,*) profile_j(:,write_idx)
  end do
  close(unit=1)

  open (unit=2,file="profile_p",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(2,*) profile_p(:,write_idx)
  end do
  close(unit=2)

  open (unit=3,file="profile_ne",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(3,*) profile_ne(:,write_idx)
  end do
  close(unit=3)

  open (unit=4,file="profile_f",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(4,*) profile_f(:,write_idx)
  end do
  close(unit=4)

  open (unit=5,file="profile_kappa",action="write",status="replace")
  do write_idx = 2, n_nodes-1
     write(5,*) a*rho_arr(write_idx), kappa0*k(rho_arr(write_idx))
  end do
  close(unit=5)

  open (unit=6,file="vloop",action="write",status="replace")
  write(6,*) VL
  close(unit=6)

  open (unit=7,file="profile_te",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(7,*) profile_t(:,write_idx)
  end do
  close(unit=7)

  open (unit=8,file="profile_omega",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(8,*) profile_omega(:,write_idx)
  end do
  close(unit=8)

  open (unit=9,file="profile_q",action="write",status="replace")
  do write_idx = 2, n_nodes-1
    write(9,*) profile_q(:,write_idx)
  end do
  close(unit=9)

  open (unit=10,file="profile_mu",action="write",status="replace")
  do write_idx = 2, n_nodes-1
     write(10,*) a*rho_arr(write_idx), mu_v*m(rho_arr(write_idx))
  end do
  close(unit=10)

  open (unit=11,file="profile_force",action="write",status="replace")
  do write_idx = 2, n_nodes-1
     write(11,*) a*rho_arr(write_idx), F0*f_m(rho_arr(write_idx))
  end do
  close(unit=11)


  deallocate(profile_j, profile_p, profile_ne, profile_f, profile_t,profile_omega)


  open (unit=13,file="Physical_params",action="write",status="replace") 
  write(13,*) a
  write(13,*) R0
  write(13,*) T0
  write(13,*) mu_v
  write(13,*) mu0
  write(13,*) eta0
  write(13,*) p0
  write(13,*) VL
  write(13,*) ( F0 * a**2.0 /(mu_v * alpha) )
  write(13,*) b_zval

  close(unit=13)


  print *, "eta_fac = ", eta_fac
  print *, "kappat = ", kappa0/ 2.18122e+26
  print *, "amu = ", mu_v / 0.364831
  print *, "vloop = ", VL/2.18122e+06 
  print *, "tcur = ", IP/795217.
  print *, "beam_v =", 1e5
  print *, "beam_rate = ", F0*8.595e21
  print *, "beam_dr = ", 0.2


  close(12)

end program eqcalc_2
