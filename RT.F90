program RT
implicit none
character(len=20) :: directory, output
character(len=50)::filepdr,fileline,filetau,filepop,prefix,spec(1:17),fileout
character(len=50) :: fileout_CII,fileout_CI,fileout_OI,fileout_CO,&
                     fileout_CII_all,fileout_CI_all,fileout_OI_all,fileout_CO_all
integer :: nfreq
integer::ptot,p,id,vv,j,ifreq
real :: dummy
double precision :: sigma, sigma_p, frac
double precision :: TMP, rho_grain
double precision,allocatable :: BB(:,:), ngrain(:), BB_dust(:,:),TPOP(:),emissivity(:),S(:,:)
double precision,allocatable::freq(:,:),velocities(:,:),phi(:),tau_test(:,:),dtau(:)
double precision,allocatable::tr_incr(:,:,:)
double precision,allocatable::x(:),av(:),tau(:,:,:)
double precision,allocatable::N(:),Tex(:,:),Tgas(:),Tdust(:),rho(:)
double precision,allocatable::pop(:,:,:),abun(:,:)
double precision,allocatable::rx(:),rav(:),rtau(:,:)
double precision,allocatable::rTgas(:),rTdust(:),rrho(:)
double precision,allocatable::rpop(:,:,:),rabun(:,:)
double precision::c,kb,pc,hp,pi,T,mhp,freq0(1:17),g(1:17,1:2)
double precision:: mh(1:17),Ntot,Ntgas
double precision::uv,Nwrite(1:17)
double precision::A(1:17)
double precision :: v_turb, v_gas, metallicity, gas_to_dust
v_gas=0     !cm/s
v_turb=1*1d5
gas_to_dust = 100.0
write(6,*) 'prefix?'
read(5,*) prefix

call constants
call readfile

write(6,*) 'Z='
read(5,*) metallicity

allocate(freq(1:17,0:nfreq-1))
allocate(TPOP(1:ptot),ngrain(1:ptot),emissivity(1:ptot))
allocate(BB(1:17,1:ptot))
allocate(BB_dust(1:17,1:ptot))
allocate(S(1:17,1:ptot))
allocate(velocities(1:17,0:nfreq-1))
allocate(phi(0:nfreq-1))
allocate(tau_test(1:17,0:nfreq-1))
allocate(tau(1:17,1:ptot,0:nfreq-1))
allocate(dtau(0:nfreq-1))
allocate(tr_incr(1:17,1:ptot,0:nfreq-1))
allocate(Tex(1:17,1:ptot))
allocate(N(1:33))

!source function calculation
do j=1,17         
  TMP=2.0D0*hp*(freq0(j)**3)/(c**2)
  BB(j,:) = TMP*(1.0D0/(exp(hp*freq0(j)/kb/2.7D0)-1.0D0))!Planck function !2.7D0 is the CMBR temperature
  ngrain=2.0D-12*rho(:)*metallicity*100./gas_to_dust
  rho_grain=2.0D0
  emissivity=(rho_grain*ngrain(:))*(0.01*(1.3*freq0(j)/3.0D11))
  BB_dust(j,:) = TMP*(1.0D0/(exp(hp*freq0(j)/kb/Tdust(:))-1.D0)*emissivity(:))
  BB = BB + BB_dust
  where (pop(j,2,:).eq.0)
    S(j,:)=0.0D0
  elsewhere
    TPOP=(pop(j,1,:)*g(j,2))/(pop(j,2,:)*g(j,1))-1.0D0
    where (abs(TPOP).lt.1.0D-50)
      S(j,:)=hp*freq0(j)*pop(j,2,:)*A(j)/4./pi
    elsewhere
      S(j,:)=TMP/TPOP
    end where
  end where
enddo

!optical depth calculation
tau_test=0
do j=1,17
  do p=1,ptot-1
  sigma=(freq0(j)/c)*sqrt(kb*Tgas(p)/mh(j)+v_turb**2/2.)
  sigma_p=(freq0(j)/c)*sqrt(kb*Tgas(ptot-1)/mh(j)+v_turb**2/2.)  
  do ifreq=0,nfreq-1 !frequencies calcuation
   freq(j,ifreq)=freq0(j)-3*sigma_p+ifreq*2*3*sigma_p/(nfreq-1)
  enddo
  phi=exp(-((1+v_gas/c)*freq(j,:)-freq0(j))**2/sigma**2/2.)/sigma/sqrt(2.*pi)
  frac=0.5*((pop(j,1,p)+pop(j,1,p+1))*g(j,2)/g(j,1)-(pop(j,2,p)+pop(j,2,p+1)))
  tau_test(j,:)=tau_test(j,:)+phi(:)*(A(j)*c**2/8./pi/freq0(j)**2)*frac*abs(x(p+1)-x(p))*pc
  tau(j,p,:)=tau_test(j,:)
  enddo

velocities(j,:) = C*(freq(j,:)/freq0(j) - 1.0D0)*1d-5
Tex(j,:)=(hp*freq0(j)/kb)/log(g(j,2)*pop(j,1,:)/pop(j,2,:)/g(j,1))
enddo

!radiation transfer solving
tr_incr(j,0,:) = 0
do p=1,ptot-2
  do j=1,17
  !dtau=abs(tau(j,p+1,:)-tau(j,p,:))
  dtau=tau(j,p,:)
    where (dtau(:).gt.1d10)
      tr_incr(j,p+1,:)=BB(j,p+1)
    elsewhere (dtau(:).gt.1d-6)
      tr_incr(j,p+1,:)=tr_incr(j,p,:)*exp(-dtau)+&
             &S(j,p)*((1-exp(-dtau))/dtau-exp(-dtau))+&
             &S(j,p+1)*(1.-(1.-exp(-dtau))/dtau)
    elsewhere (dtau(:).le.1d-6)
      tr_incr(j,p+1,:)=tr_incr(j,p,:)*(1-dtau)+(BB(j,p)+BB(j,p+1))*dtau/2.
    end where

    if (tr_incr(j,p+1,20).le.tr_incr(j,p,20)) then
      tr_incr(j,p+1,:)=tr_incr(j,p,:)
    end if
  enddo
enddo

!fileout_CII=trim(adjustl(directory))//'/'//trim(adjustl(output))//'.CII_br_temp.dat'
!fileout_CI=trim(adjustl(directory))//'/'//trim(adjustl(output))//'.CI_br_temp.dat'
!fileout_OI=trim(adjustl(directory))//'/'//trim(adjustl(output))//'.OI_br_temp.dat'
!fileout_CO=trim(adjustl(directory))//'/'//trim(adjustl(output))//'.C12O_br_temp.dat'
!fileout_CII_all=trim(adjustl(directory))//'/'//trim(adjustl(output))//'CII.dat'
!fileout_CI_all=trim(adjustl(directory))//'/'//trim(adjustl(output))//'CI.dat'
!fileout_OI_all=trim(adjustl(directory))//'/'//trim(adjustl(output))//'OI.dat'
!fileout_CO_all=trim(adjustl(directory))//'/'//trim(adjustl(output))//'C12O.dat'

fileout_CII='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'.CII_br_temp.dat'
fileout_CI='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'.CI_br_temp.dat'
fileout_OI='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'.OI_br_temp.dat'
fileout_CO='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'.C12O_br_temp.dat'
fileout_CII_all='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'CII.dat'
fileout_CI_all='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'CI.dat'
fileout_OI_all='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'OI.dat'
fileout_CO_all='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//'C12O.dat'

open(unit=100,file=fileout_CII,status='replace')
open(unit=101,file=fileout_CI,status='replace')
open(unit=102,file=fileout_OI,status='replace')
open(unit=103,file=fileout_CO,status='replace')

open(unit=121,file=fileout_CII_all,status='replace')
open(unit=122,file=fileout_CI_all,status='replace')
open(unit=123,file=fileout_OI_all,status='replace')
open(unit=124,file=fileout_CO_all,status='replace')
  write(100,*) 0, 0, velocities(1,:)
  write(101,*) 0, 0, velocities(2,:)
  write(102,*) 0, 0, velocities(5,:)
  write(103,*) 0, 0, velocities(8,:)

N=0;Ntot=0;Ntgas=0
do p=1,ptot-2
  N=N+0.5*(rho(p)*abun(:,p)+rho(p+1)*abun(:,p+1))*abs(x(p+1)-x(p))*pc
  Ntot=Ntot+0.5*(rho(p)+rho(p+1))*abs(x(p+1)-x(p))*pc
  Ntgas=Ntgas+0.5*(rho(p)*tgas(p)+rho(p+1)*tgas(p+1))*abs(x(p+1)-x(p))*pc
  write(100,*) p, 0, tr_incr(1,ptot-2,:)*c**2/2./kb/freq0(1)**2
  write(101,*) p, 0, tr_incr(2,ptot-2,:)*c**2/2./kb/freq0(2)**2
  write(102,*) p, 0, tr_incr(5,ptot-2,:)*c**2/2./kb/freq0(5)**2
  write(103,*) p, 0, tr_incr(8,ptot-2,:)*c**2/2./kb/freq0(8)**2

  write(121,'(100ES11.3)') x(p),N(11),Tex(1,p),tau(1,p,nfreq/2),tr_incr(1,p,nfreq/2)*c**2/2./kb/freq0(1)**2
  write(122,'(100ES11.3)') x(p),N(25),Tex(2,p),tau(2,p,nfreq/2),tr_incr(2,p,nfreq/2)*c**2/2./kb/freq0(2)**2
  write(123,'(100ES11.3)') x(p),N(30),Tex(5,p),tau(5,p,nfreq/2),tr_incr(5,p,nfreq/2)*c**2/2./kb/freq0(5)**2
  write(124,'(100ES11.3)') Ntot,N(28),Tex(8,p),tau(8,p,nfreq/2),tr_incr(8,p,nfreq/2)*c**2/2./kb/freq0(8)**2,rho(p)
enddo


write(6,*) '<Tgas>=',NTgas/Ntot
write(6,*) 'Linewidth (CO) = ',2.*sqrt(2.*log(2.))*sqrt(kb*(NTgas/Ntot)/mh(8)+v_turb**2/2.)/1d5,' [km/s]'
Nwrite(1)=N(11);Nwrite(2:4)=N(25);Nwrite(5:7)=N(30);Nwrite(8:17)=N(28)
write(6,*) 'Species || Column density || Tex || tau (pdr_tot-2) || Tr (pdr_tot-2)'
do j=1,17
  write(6,'(A10,2X,5ES11.3)') spec(j),Nwrite(j),Tex(j,ptot),tau(j,ptot-2,nfreq/2),&
                              &tr_incr(j,ptot-2,nfreq/2)*c**2/2./kb/freq0(j)**2
enddo


contains
subroutine readfile

!open(unit=12,file='params.dat',status='old')
  !read(12,*); read(12,*); read(12,*)
  !read(12,'(A)')
  !read(12,*)
  !read(12,'(A)') directory
  !read(12,*) output
  !read(12,*)
  !read(12,*)
  !read(12,*)
  !read(12,*)
  !read(12,*)
  !read(12,*)
  !read(12,*)
  !read(12,*) v_turb
  !read(12,*)
  !read(12,*)
  !read(12,*) gas_to_dust
  !read(12,*) metallicity
  !close(12)

filepdr='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//".pdr.fin"
filetau='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//".opdp.fin"
filepop='./ALL_TESTS/'//trim(adjustl(prefix))//'p/'//trim(adjustl(prefix))//".spop.fin"
open(unit=2,file=filepdr,status='old')
open(unit=3,file=filepop,status='old')
open(unit=4,file=filetau,status='old')

nfreq=40
ptot=0
do 
  read(2,*,end=100) dummy
  ptot=ptot+1
enddo
100 continue
rewind(2)
allocate(x(1:ptot),av(1:ptot),abun(1:33,1:ptot))
allocate(Tgas(1:ptot),Tdust(1:ptot),rho(1:ptot),pop(1:17,1:2,1:ptot))
allocate(rx(1:ptot),rav(1:ptot),rabun(1:33,1:ptot))
allocate(rTgas(1:ptot),rTdust(1:ptot),rrho(1:ptot),rpop(1:17,1:2,1:ptot),rtau(1:17,1:ptot))
do p=1,ptot
  read(2,*) id,rx(p),rav(p),rTgas(p),rTdust(p),dummy,rrho(p),uv,rabun(1:33,p)
  read(3,*) id,dummy,rpop(1,1,p),rpop(1,2,p),dummy,dummy,dummy,&
       &rpop(2,1,p),rpop(2,2,p),rpop(3,2,p),dummy,dummy,&
       &rpop(5,1,p),rpop(5,2,p),rpop(6,2,p),dummy,dummy,&
       &rpop(8,1,p),rpop(8,2,p),rpop(9,2,p),rpop(10,2,p),rpop(11,2,p),&
       &rpop(12,2,p),rpop(13,2,p),rpop(14,2,p),rpop(15,2,p),rpop(16,2,p),rpop(17,2,p)
  rpop(3,1,p)=rpop(2,1,p);  rpop(4,2,p)=rpop(3,2,p);  rpop(4,1,p)=rpop(2,2,p)
  rpop(6,1,p)=rpop(5,1,p);  rpop(7,2,p)=rpop(6,2,p);  rpop(7,1,p)=rpop(5,2,p)
  rpop(9,1,p)=rpop(8,2,p);  rpop(10,1,p)=rpop(9,2,p);  rpop(11,1,p)=rpop(10,2,p)
  rpop(12,1,p)=rpop(11,2,p);  rpop(13,1,p)=rpop(12,2,p);  rpop(14,1,p)=rpop(13,2,p)
  rpop(15,1,p)=rpop(14,2,p);  rpop(16,1,p)=rpop(15,2,p);  rpop(17,1,p)=rpop(16,2,p)
  read(4,*) id,dummy,rtau(1:17,p)
enddo

x=rx;av=rav;abun=rabun;Tgas=rTgas;Tdust=rTdust;rho=rrho;pop=rpop
!tau(1:17,1:ptot,nfreq/2)=rtau

write(6,*) 'Density=',rho(1)
write(6,*) 'Temperature=',Tgas(1)

return
end subroutine

subroutine constants

 c=2.9979246d10 !cm/s
kb=1.380650d-16 !erg / K  *or*  g cm^2 / K s^2
hp=6.6260696e-27 !erg s  *or*  g cm^2 / s
pc=3.085677d+18
pi=3.1415927
mhp=1.6726218d-24

freq0(1)=1900.5369d9     ;g(1,1)=2.0  ;g(1,2)=4.0   ; spec(1)="CII 158um"
freq0(2)=492.16065d9     ;g(2,1)=1.0  ;g(2,2)=3.0   ; spec(2)="CI (1-0)"
freq0(3)=1301.50262d9    ;g(3,1)=1.0  ;g(3,2)=5.0   ; spec(3)="CI (2-0)"
freq0(4)=809.34197d9     ;g(4,1)=3.0  ;g(4,2)=5.0   ; spec(4)="CI (2-1)"
freq0(5)=4744.77749d9    ;g(5,1)=5.0  ;g(5,2)=3.0   ; spec(5)="OI  1-0 "
freq0(6)=6804.84658d9    ;g(6,1)=5.0  ;g(6,2)=1.0   ; spec(6)="OI  2-0 "
freq0(7)=2060.06909d9    ;g(7,1)=3.0  ;g(7,2)=1.0   ; spec(7)="OI  2-1 "
freq0(8)=115.2712018d9   ;g(8,1)=1.0  ;g(8,2)=3.0   ; spec(8)="CO (1-0)"
freq0(9)=230.538d9       ;g(9,1)=3.0  ;g(9,2)=5.0   ; spec(9)="CO (2-1)"
freq0(10)=345.7959899d9  ;g(10,1)=5.0 ;g(10,2)=7.0  ; spec(10)="CO (3-2)"
freq0(11)=461.040768d9   ;g(11,1)=7.0 ;g(11,2)=9.0  ; spec(11)="CO (4-3)"
freq0(12)=576.2679305d9  ;g(12,1)=9.0 ;g(12,2)=11.0 ; spec(12)="CO (5-4)"
freq0(13)=691.4730763d9  ;g(13,1)=11.0;g(13,2)=13.0 ; spec(13)="CO (6-5)"
freq0(14)=806.6518060d9  ;g(14,1)=13.0;g(14,2)=15.0 ; spec(14)="CO (7-6)"
freq0(15)=921.7997000d9  ;g(15,1)=15.0;g(15,2)=17.0 ; spec(15)="CO (8-7)"
freq0(16)=1036.9123930d9 ;g(16,1)=17.0;g(16,2)=19.0 ; spec(16)="CO (9-8)"
freq0(17)=1151.9854520d9 ;g(17,1)=19.0;g(17,2)=21.0 ; spec(17)="CO (10-9)"

!Einstein A coefficients
A(1)=2.321d-06 ;  mh(1)=12.*mhp  
A(2)=7.880d-08 ;  mh(2)=12.*mhp  
A(3)=1.810d-14 ;  mh(3)=12.*mhp  
A(4)=2.650d-07 ;  mh(4)=12.*mhp  
A(5)=8.910d-05 ;  mh(5)=16.*mhp  
A(6)=1.340d-10 ;  mh(6)=16.*mhp  
A(7)=1.750d-05 ;  mh(7)=16.*mhp  
A(8)=7.203d-08 ;  mh(8)=28.*mhp  
A(9)=6.910d-07 ;  mh(9)=28.*mhp  
A(10)=2.497d-06;  mh(10)=28.*mhp 
A(11)=6.126d-06;  mh(11)=28.*mhp 
A(12)=1.221d-05;  mh(12)=28.*mhp 
A(13)=2.137d-05;  mh(13)=28.*mhp 
A(14)=3.422d-05;  mh(14)=28.*mhp 
A(15)=5.134d-05;  mh(15)=28.*mhp 
A(16)=7.330d-05;  mh(16)=28.*mhp 
A(17)=1.006d-04;  mh(17)=28.*mhp 

return
end subroutine

end program
