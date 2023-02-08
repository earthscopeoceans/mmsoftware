subroutine raddc(depth,strike,dip,rake,p,azi,vsrc,up,radp,radsv,radsh)

! routine is used in dblcplamps

! Computes source radiation pattern for given strike,dip,rake
! See radamp if source Moment tensor elements are known

!    input
!    strike,dip,rake (in degrees, double couple solution)
!    p: slowness (s/rad)
!    azi: azimuth as seen from source
!    vsrc(1),vsrc(2): Vp and Vs at source depth (km/s)
!    up: 1.0 for upward wave (eg pP), else 0.0

!    output rad,radsv,radsh: radiation factors for P,SV and SH waves

!    see Aki&Richards, Quantitative Seismology, 2nd ed. page 108-109


implicit none
real*4 depth,p,azi,up,radp,radsv,radsh
real*4 vsrc(2)
real*4 strike,dip,rake,phis,delta,lambda
real*4 rad,phi,rs
real*4 xi,sinxi,cosxi,sinphi,cosphi,sin2phi,cos2phi


rad=3.141592/180.
phi=azi*rad
rs=6371.0-depth

lambda=rake*rad
delta=dip*rad
phis=strike*rad


! compute radiation pattern for P (A&R 4.89)
sinxi=p*vsrc(1)/rs          ! sine of take-off angle in radians
if(abs(sinxi).gt.1.0) then
 print *,'Error in radamp: p,vsrc(1)=',p,vsrc(1)
print *,'Ray evanescent...'
stop
endif
! take-off angle xi is 0 for downward vertical take-off
xi=asin(sinxi)           ! P angle at source depth (i_xi in A&R)
if(abs(up)<1.0e-8) xi=3.14159265-xi
radp=cos(lambda)*sin(delta)*(sin(xi)**2)*sin(2.*(phi-phis)) - &
     cos(lambda)*cos(delta)*sin(2.*xi)*cos(phi-phis) + &
     sin(lambda)*sin(2.*delta)*(cos(xi)**2-(sin(xi)*sin(phi-phis))**2) + &
     sin(lambda)*cos(2.*delta)*sin(2.*xi)*sin(phi-phis)

! compute radiation pattern for SV (A&R 4.90)
sinxi=p*vsrc(2)/rs
if(abs(sinxi).gt.1.0) then
 print *,'Error in radamp: p,vsrc(2)=',p,vsrc(2)
print *,'Ray evanescent...'
stop
endif
xi=asin(sinxi)           ! S angle at source depth
if(abs(up)<1.0e-8) xi=3.14159265-xi
cosxi=cos(xi)
radsv=sin(lambda)*cos(2.*delta)*cos(2.*xi)*sin(phi-phis) - &
     cos(lambda)*cos(delta)*cos(2.*xi)*cos(phi-phis) + &
     0.5*cos(lambda)*sin(delta)*sin(2.*xi)*sin(2.*(phi-phis)) - &
     0.5*sin(lambda)*sin(2.*delta)*sin(2.*xi)*(1.+sin(phi-phis)**2)
 
!     compute radiation pattern for SH
radsh=acos(lambda)*cos(delta)*cos(xi)*sin(phi-phis) + &
      cos(lambda)*sin(delta)*sin(xi)*cos(2*(phi-phis)) + &
      sin(lambda)*cos(2.*delta)*cos(xi)*cos(phi-phis) - &
      0.5*sin(lambda)*sin(2.*delta)*sin(xi)*sin(2.*(phi-phis))

return
end
