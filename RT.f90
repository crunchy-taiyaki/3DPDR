program RT
use line_transfer
implicit none

type(line) :: CII_10, CI_10, CI_21, OI_10, OI_21, CO_10, CO_21
character(len=40) :: prefix, root, velocity_flag
character(len=100) :: directory
character(len=150) :: filevel,filepdr,filepop
integer::ptot,p,id,vv,j,ifreq
real :: dummy
double precision,allocatable::x(:),av(:),init_velocities(:),Tgas(:),Tdust(:)
double precision,allocatable::rho(:),abun(:,:),CIIpop(:,:),CIpop(:,:),OIpop(:,:),COpop(:,:)
double precision :: metallicity, gas_to_dust

root = './ALL_TESTS/vel'
init_velocities = 0.0

write(6,*) 'prefix?'
read(5,*) prefix
write(6,*) 'velocities on ("y" or "n")?'
read(5,*) velocity_flag

call readfile
write(6,*) 'Z='
read(5,*) metallicity
gas_to_dust = 100.

 CII_10 = line("CII_158um",1900.5369d9,2.0,4.0,2.321d-06,12)
 CI_10 = line("CI_(1-0)",492.16065d9,1.0,3.0,7.880d-08,12)
 CI_21 = line("CI_(2-1)",809.34197d9,3.0,5.0,2.650d-07,12)
 OI_10 = line("OI_(1-0)",4744.77749d9,5.0,3.0,8.910d-05,16)
 OI_21 = line("OI_(2-1)",2060.06909d9,3.0,1.0,1.750d-05,16)
 CO_10 = line("CO_(1-0)",115.2712018d9,1.0,3.0,7.203d-08,28)
 CO_21 = line("CO_(2-1)",230.538d9,3.0,5.0,6.910d-07,28)

directory = trim(adjustl(root))//'/'//trim(adjustl(prefix))
if (velocity_flag.eq.'y') then
 directory= trim(adjustl(root))//'/'//trim(adjustl(prefix))//'/'//'vel'
endif
write(6,*) 'Species || Column density || Tex || tau || Tr '

  call CII_10%solve_rt(directory,x,av,init_velocities,abun(11,:),Tgas,Tdust,&
                       &rho,CIIpop(1,:),CIIpop(0,:),metallicity,gas_to_dust)
  call CI_10%solve_rt(directory,x,av,init_velocities,abun(25,:),Tgas,Tdust,&
                      &rho,CIpop(1,:),CIpop(0,:),metallicity,gas_to_dust)
  call CI_21%solve_rt(directory,x,av,init_velocities,abun(25,:),Tgas,Tdust,&
                      &rho,CIpop(2,:),CIpop(1,:),metallicity,gas_to_dust)
  call OI_10%solve_rt(directory,x,av,init_velocities,abun(30,:),Tgas,Tdust,&
                      &rho,OIpop(1,:),OIpop(0,:),metallicity,gas_to_dust)
  call OI_21%solve_rt(directory,x,av,init_velocities,abun(30,:),Tgas,Tdust,&
                      &rho,OIpop(2,:),OIpop(1,:),metallicity,gas_to_dust)
  call CO_10%solve_rt(directory,x,av,init_velocities,abun(28,:),Tgas,Tdust,&
                      &rho,COpop(1,:),COpop(0,:),metallicity,gas_to_dust)
  call CO_21%solve_rt(directory,x,av,init_velocities,abun(28,:),Tgas,Tdust,&
                      &rho,COpop(2,:),COpop(1,:),metallicity,gas_to_dust)

contains

subroutine readfile
filepdr=trim(adjustl(root))//'/'//trim(adjustl(prefix))//'/'//"pdr.fin"
filepop=trim(adjustl(root))//'/'//trim(adjustl(prefix))//'/'//"spop.fin"
open(unit=1,file=filepdr,status='old')
open(unit=2,file=filepop,status='old')

if (velocity_flag.eq.'y') then
  filevel=trim(adjustl(root))//'/'//trim(adjustl(prefix))//'/'//trim(adjustl(prefix))//".dat"
  open(unit=3,file=filevel,status='old')
endif
ptot=0
do 
  read(1,*,end=100) dummy
  ptot=ptot+1
enddo
100 continue
rewind(1)
allocate(x(1:ptot),av(1:ptot),init_velocities(1:ptot),Tgas(1:ptot),Tdust(1:ptot))
allocate(rho(1:ptot),abun(1:33,1:ptot))
allocate(CIIpop(0:4,1:ptot),CIpop(0:4,1:ptot),OIpop(0:4,1:ptot),COpop(0:40,1:ptot))

do p=1,ptot
  read(1,*) id,x(p),av(p),Tgas(p),Tdust(p),dummy,rho(p),dummy,abun(1:33,p)
  read(2,*)id,dummy,CIIpop(0:4,p),CIpop(0:4,p),OIpop(0:4,p),COpop(0:40,p) 
enddo

if (velocity_flag.eq.'y') then
do p=1,ptot
  read(3,*)dummy,dummy,dummy,dummy,init_velocities(p)
  write(6,*)init_velocities(p)
enddo
endif 

write(6,*) 'Density=',rho(1)
write(6,*) 'Temperature=',Tgas(1)
 end subroutine

end program
