subroutine radiation_transfer(pdr_ptot,nlev,nfreq,&
                     &frequencies,density,metallicity,dust_temperature,s_pop_array,weights,A_COEFFS,&
                     &tau_profile_array,bright_temperature)

use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d
implicit none

integer(kind=i4b), intent(in) :: pdr_ptot
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: nfreq
real(kind=dp), intent(in) :: frequencies(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: metallicity
real(kind=dp), intent(in) :: density(1:pdr_ptot),dust_temperature(1:pdr_ptot)
real(kind=dp), intent(in) :: s_pop_array(1:nlev,1:pdr_ptot)
real(kind=dp), intent(in) :: weights(1:nlev)
real(kind=dp), intent(in) :: A_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: tau_profile_array(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
integer(kind=i4b) :: ilevel, jlevel
integer(kind=i4b) :: pp
real(kind=dp) :: TMP
real(kind=dp) :: BB_ij(1:pdr_ptot), TPOP(1:pdr_ptot), S_ij(1:pdr_ptot)
real(kind=dp) :: rho_grain
real(kind=dp) :: ngrain(1:pdr_ptot), emissivity(1:pdr_ptot), BB_ij_dust(1:pdr_ptot)
real(kind=dp) :: dtau(0:nfreq-1)
real(kind=dp) :: beta(0:nfreq-1,1:pdr_ptot)
real(kind=dp) :: intensity_profile(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
real(kind=dp), intent(out) ::bright_temperature(1:nlev,1:nlev,0:nfreq-1)

    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit

         TMP=2.0D0*hp*(frequencies(ilevel,jlevel)**3)/(c**2)
         BB_ij = TMP*(1.0D0/(EXP(hp*frequencies(ilevel,jlevel)/kb/2.7D0)-1.0D0))!Planck function !2.7D0 is the CMBR temperature
         ngrain(1:pdr_ptot)=2.0D-12*density(:)*metallicity*100./g2d
         rho_grain=2.0D0
         emissivity(1:pdr_ptot)=(rho_grain*ngrain(:))*(0.01*(1.3*frequencies(ilevel,jlevel)/3.0D11))
         BB_ij_dust(1:pdr_ptot) = TMP*(1.0D0/(EXP(hp*frequencies(ilevel,jlevel)/kb/dust_temperature(:))-1.D0)*emissivity(:))
         BB_ij(1:pdr_ptot) = BB_ij + BB_ij_dust
         where (s_pop_array(ilevel,:).eq.0)
           S_ij(1:pdr_ptot)=0.0D0
         elsewhere
           TPOP(1:pdr_ptot)=(s_pop_array(jlevel,:)*weights(ilevel))/(s_pop_array(ilevel,:)*weights(jlevel))-1.0D0
           where(abs(TPOP).lt.1.0D-50)
             S_ij(1:pdr_ptot)=hp*frequencies(ilevel,jlevel)*s_pop_array(ilevel,:)*A_COEFFS(ilevel,jlevel)/4./pi
           elsewhere
	     S_ij(1:pdr_ptot)=TMP/TPOP
           end where
         end where
         

         intensity_profile(ilevel,jlevel,:,1) = 0.0D0
	 do pp=1,pdr_ptot-2
           dtau(0:nfreq-1) = abs(tau_profile_array(ilevel,jlevel,:,pp+1)-tau_profile_array(ilevel,jlevel,:,pp))
           beta(:,pp) = (1-EXP(-tau_profile_array(ilevel,jlevel,:,pp)))/tau_profile_array(ilevel,jlevel,:,pp)
           intensity_profile(ilevel,jlevel,:,pp) = intensity_profile(ilevel,jlevel,:,pp)*beta(:,pp)+&
                                  &S_ij(pp+1)*(1-beta(:,pp))+BB_ij(pp)*(beta(:,pp)-EXP(-dtau))
	   where(dtau(:).gt.1d10)
             intensity_profile(ilevel,jlevel,:,pp) = BB_ij(pp+1)
           elsewhere(dtau(:).lt.1d-6)
             intensity_profile(ilevel,jlevel,:,pp) = (1-dtau)*intensity_profile(ilevel,jlevel,:,pp)+&
                                                      &dtau*(BB_ij(pp)+BB_ij(pp+1))/2
           end where
         enddo
bright_temperature(ilevel,jlevel,:) = intensity_profile(ilevel,jlevel,:,pdr_ptot-2)
!bright_temperature(ilevel,jlevel,:) = S_ij(2)*(1-EXP(-(tau_profile_array(ilevel,jlevel,:,2)-tau_profile_array(ilevel,jlevel,:,1))))
       !bright_temperature(ilevel,jlevel,:) = intensity_profile(ilevel,jlevel,:,2)                                     
       !bright_temperature(ilevel,jlevel,:) = intensity_profile(ilevel,jlevel,:,pdr_ptot)*c**2/(2*kb*frequencies(ilevel,jlevel)**2)
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev

  return

end subroutine radiation_transfer
