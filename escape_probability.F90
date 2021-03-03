!T.Bisbas, T.Bell

subroutine escape_probability(transition, dust_temperature, nrays, nlev,nfreq, &
                   &A_COEFFS, B_COEFFS, C_COEFFS, &
                   &frequencies,s_evalpop, maxpoints,&
                   & T_evalpoint, vel_evalpoint, Tguess, v_turb, v_gas, &
                   &s_jjr, s_pop, s_evalpoint, weights,cooling_rate,line,tau,&
                   &coolant,density,metallicity,bbeta)


use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d

implicit none

integer(kind=i4b), intent(in) :: nrays
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: nfreq
integer(kind=i4b), intent(in) :: maxpoints
integer(kind=i4b), intent(in) :: s_jjr(0:nrays-1)
integer(kind=i4b), intent(in) :: coolant
real(kind=dp), intent(in) :: A_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: B_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: C_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: frequencies(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: s_evalpop(0:nrays-1,0:maxpoints,1:nlev)
real(kind=dp), intent(in) :: s_evalpoint(1:3,0:nrays-1,0:maxpoints)
real(kind=dp), intent(in) :: T_evalpoint(0:nrays-1,0:maxpoints)
real(kind=dp), intent(in) :: vel_evalpoint(0:nrays-1,0:maxpoints)
real(kind=dp) :: Tguess, v_turb, v_gas !, intent(in)
real(kind=dp), intent(in) :: weights(1:nlev)
real(kind=dp), intent(in) :: s_pop(1:nlev)
real(kind=dp), intent(in) :: dust_temperature,density,metallicity

real(kind=dp), intent(out) :: line(1:nlev,1:nlev)
real(kind=dp), intent(out) :: cooling_rate
real(kind=dp), intent(out) :: tau(1:nlev,1:nlev,0:nrays-1)
real(kind=dp),intent(out) :: bbeta(1:nlev,1:nlev,0:nrays-1)
real(kind=dp), intent(inout) :: transition(1:nlev,1:nlev)

integer(kind=i4b) :: i, j
integer(kind=i4b) :: ifreq
integer(kind=i4b) :: ilevel, jlevel
real(kind=dp) :: beta_ij, beta_ij_sum
real(kind=dp) :: frac1, frac2, frac3, rhs2 
real(kind=dp) :: thermal_velocity
real(kind=dp) :: tpop, tmp2
real(kind=dp) :: S_ij, BB_ij
real(kind=dp) :: tau_increment
real(kind=dp) :: frequency(0:nfreq-1),velocities(0:nfreq-1)
real(kind=dp) :: doppler_profile(0:nfreq-1)
real(kind=dp) :: tau_ij(0:nrays-1)
real(kind=dp) :: tau_ij_profile(0:nfreq-1,0:nrays-1)
real(kind=dp) :: tau_increment_profile(0:nfreq-1)
real(kind=dp) :: tau_profile(1:nlev,1:nlev,0:nfreq-1,0:nrays-1)
real(kind=dp) :: beta_ij_ray(0:nrays-1)
real(kind=dp) :: beta_ij_ray_profile(0:nfreq-1,0:nrays-1)
real(kind=dp) :: field(1:nlev,1:nlev)
real(kind=dp) :: emissivity, bb_ij_dust, ngrain, rho_grain

v_gas = v_gas*1.
line=0.0D0
cooling_rate = 0.0D0
field=0.0D0
frac2=1.0D0/sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)
thermal_velocity=sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)

    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit	 
	 !init frequency array
	 do ifreq=0,nfreq-1
	   frequency(ifreq)=frequencies(ilevel,jlevel)*C/(C+v_gas*1d5)-3*thermal_velocity+ifreq*3*thermal_velocity*2/(nfreq-1)
	 enddo
	 doppler_profile=exp(-((1+v_gas*1d5/C)*frequency(:)-frequencies(ilevel,jlevel))**2/(2.0*thermal_velocity**2))/&
                                    &(thermal_velocity*sqrt(2.*pi))
         tau_ij=0.0D0
	 tau_ij_profile(0:nfreq-1,0:nrays-1)=0.0D0
         beta_ij=0.0D0; beta_ij_ray=0.0D0
	 beta_ij_ray_profile=0.0D0
	 beta_ij_sum=0.0D0

         frac1=(A_COEFFS(ilevel,jlevel)*(C**3))/(8.0*pi*(frequencies(ilevel,jlevel)**3))
         TMP2=2.0D0*HP*(FREQUENCIES(ilevel,jlevel)**3)/(C**2)
         BB_ij = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/2.7D0)-1.0D0)) !Planck function !2.7D0 is the CMBR temperature
         NGRAIN=2.0D-12*density*metallicity*100./g2d
         rho_grain=2.0D0
         EMISSIVITY=(RHO_GRAIN*NGRAIN)*(0.01*(1.3*FREQUENCIES(ilevel,jlevel)/3.0D11))
         BB_ij_dust = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/DUST_TEMPERATURE)-1.D0)*EMISSIVITY)
         BB_ij = BB_ij + BB_ij_dust
         if (s_pop(ilevel).eq.0) then
            S_ij=0.0D0
            beta_ij=1.0D0
            goto 2
         endif
         TPOP=(s_pop(jlevel)*WEIGHTS(ilevel))/(s_pop(ilevel)*WEIGHTS(jlevel))-1.0D0
         if(abs(TPOP).lt.1.0D-50) then
              S_ij=HP*FREQUENCIES(ilevel,jlevel)*s_pop(ilevel)*A_COEFFS(ilevel,jlevel)/4./pi
              beta_ij=1.0D0
              goto 1
         else
         !calculation of source function (taken from UCL_PDR)
              S_ij=TMP2/TPOP
         endif
         do j=0,nrays-1
#ifdef PSEUDO_1D
         if (j.ne.6) then
           tau_ij(j) = 1.0D50
	   tau_ij_profile(0:nfreq-1,j) = 1.0D50
         else
#endif
#ifdef PSEUDO_2D
         if (abs(vectors(3,j).gt.1d-10) then
	     tau_ij(j) = 1.0D50 !Not in Equator
	     tau_ij_profile(0:nfreq-1,j) = 1.0D50
#endif

 do i=1,s_jjr(j)
             !calculations of tau_ij
     frac3=((s_evalpop(j,i-1,jlevel)*weights(ilevel)-s_evalpop(j,i-1,ilevel)*weights(jlevel))+&
      &(s_evalpop(j,i,jlevel)*weights(ilevel)-s_evalpop(j,i,ilevel)*weights(jlevel)))/2./weights(jlevel)
     rhs2=sqrt((s_evalpoint(1,j,i-1)-s_evalpoint(1,j,i))**2+&
              &(s_evalpoint(2,j,i-1)-s_evalpoint(2,j,i))**2+&
              &(s_evalpoint(3,j,i-1)-s_evalpoint(3,j,i))**2) !adaptive step
     tau_increment=frac1*frac2*frac3*rhs2*PC
     tau_increment_profile(0:nfreq-1)=frac1*doppler_profile(0:nfreq-1)*frac3*rhs2*PC
     tau_ij(j)=tau_ij(j)+tau_increment !optical depth
     tau_ij_profile(0:nfreq-1,j)=tau_ij_profile(0:nfreq-1,j)+tau_increment_profile(0:nfreq-1) !optical depth depending on frequency
 enddo !i=1,jr(j)
#ifdef PSEUDO_1D
         endif
#endif
#ifdef PSEUDO_2D
         endif
#endif

           ! Prevent exploding beta values caused by strong masing (tau < -5)
           ! Assume tau = -5 and calculate the escape probability accordingly
           if (tau_ij(j).lt.-5.0D0) then
              beta_ij_ray(j)=(1.0D0-EXP(5.0D0))/(-5.0D0)
!           ! Treat weak masing using the standard escape probability formalism
!           else if (tau_ij(j).lt.0.0D0) then
!              beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
           ! Prevent floating point overflow caused by very low opacity (tau < 1e-8)
           else if (abs(tau_ij(j)).lt.1.0D-8) then
              beta_ij_ray(j)=1.0D0
           ! For all other cases use the standard escape probability formalism
           else
              beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
           endif
        !if ((j.eq.6).and.(coolant.eq.1)) then
        !write(*,*)sum(beta_ij_ray),'before'
        !endif

           beta_ij_ray_profile(:,j)=doppler_profile*(1.0D0-EXP(-tau_ij_profile(:,j)))/tau_ij_profile(:,j)
           where (tau_ij_profile(:,j).lt.-5.0D0)
              beta_ij_ray_profile(:,j)=doppler_profile*(1.0D0-EXP(5.0D0))/(-5.0D0)
           elsewhere (abs(tau_ij_profile(:,j)).lt.1.0D-8)
              beta_ij_ray_profile(:,j)=doppler_profile
           end where
           beta_ij_ray(j) = 0.0d0
           do ifreq=0,nfreq-2
             beta_ij_ray(j) = beta_ij_ray(j)+&
              &(beta_ij_ray_profile(ifreq+1,j)+beta_ij_ray_profile(ifreq,j))*&
              &abs(frequency(ifreq+1)-frequency(ifreq))/2.
           enddo
        !if ((j.eq.6).and.(coolant.eq.1)) then
        !write(*,*)sum(beta_ij_ray),'after'
        !endif
	!=============
	tau(ilevel,jlevel,j)=tau_ij(j)
	bbeta(ilevel,jlevel,j)=beta_ij_ray(j)
	!=============
         enddo !j=0,nrays-1

         beta_ij_sum=sum(beta_ij_ray)
         !calculation of average beta_ij in the origin grid point
#ifdef PSEUDO_1D
         beta_ij = beta_ij_sum
#elif PSEUDO_2D
         beta_ij = beta_ij_sum / 4.
#else
         beta_ij = beta_ij_sum / real(nrays,kind=DP) 
#endif


1 continue
         line(ilevel,jlevel) = A_COEFFS(ilevel,jlevel)*HP*frequencies(ilevel,jlevel) * &
                             & s_pop(ilevel)*beta_ij*(S_ij-BB_ij)/S_ij
         cooling_rate = cooling_rate + line(ilevel,jlevel)
2 continue
         !<J_ij>
         field(ilevel,jlevel) = (1.0D0-beta_ij)*S_ij + beta_ij*BB_ij
         field(jlevel,ilevel) = field(ilevel,jlevel)
         !J_ij(p)
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev
 

    !R_IJ CALCULATIONS
    !Update the transition matrix: Rij = Aij + Bij.<J> + Cij				    	 
    DO ilevel=1,NLEV
      DO jlevel=1,NLEV
        TRANSITION(ilevel,jlevel)=A_COEFFS(ilevel,jlevel)&
        & +B_COEFFS(ilevel,jlevel)*FIELD(ilevel,jlevel)&
        & +C_COEFFS(ilevel,jlevel)
        IF(ABS(TRANSITION(ilevel,jlevel)).LT.1.0D-50) TRANSITION(ilevel,jlevel)=0.0D0

      ENDDO !jlevel=1,nlev
    ENDDO !ilevel=1,nlev

  return

end subroutine escape_probability
