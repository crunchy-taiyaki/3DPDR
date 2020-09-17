subroutine radiation_transfer(intensity_profile_array,pdr_ptot,nlev,nfreq,tau_profile_array,intensity_profile)

use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d
implicit none

integer(kind=i4b), intent(in) :: pdr_ptot
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: nfreq
real(kind=dp), intent(in) :: intensity_profile_array(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
real(kind=dp), intent(in) :: tau_profile_array(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
integer(kind=i4b) :: ilevel, jlevel
integer(kind=i4b) :: pp
real(kind=dp), intent(out) :: intensity_profile(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)

    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit
         intensity_profile(ilevel,jlevel,:,1) = intensity_profile_array(ilevel,jlevel,:,1)
	 do pp=2,pdr_ptot
           intensity_profile(ilevel,jlevel,:,pp) = EXP(-(tau_profile_array(ilevel,jlevel,:,pp)-&
                &tau_profile_array(ilevel,jlevel,:,pp-1)))*intensity_profile(ilevel,jlevel,:,pp-1)+&
		&intensity_profile_array(ilevel,jlevel,:,pp)
         enddo
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev
  return

end subroutine radiation_transfer
