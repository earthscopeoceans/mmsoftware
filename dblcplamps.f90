program dblcplamps

! compile: g7 dblcplamps delaz raddc ttak135 bin

! Construct a GMT file to map P wave amplitudes on a globe

! Input is:
! Earthquake lat,lon,depth (deg, km)
! Lat/lon limits of the map in degrees
! strike,dip,rake in degrees

! Get moment tensor elements from the CMT catalogue, eg for older events:
! http://www.globalcmt.org/CMTfiles.html
! for more recent ones:
! http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/
! or the quick solutions for very recent ones:
! http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/

use ttak135

implicit none

integer, parameter :: NGRID=100                 ! map size in pixels
real, parameter :: FG=0.02                      ! lat/lon subdivide
real*4 :: slat,slon,depth                       ! source coordinates
real*4 :: lat1,lat2,lon1,lon2,dlat,dlon         ! map limits
real*4 :: p,azi,vsrc,up,radp,radsv,radsh,dist  ! radiation
real*4 :: strike,dip,rake
real*4 :: azie,azis,xlat,xlon,rmax,pvel,pmax,r2d
real*4 :: prad(NGRID,NGRID)                     ! map
integer :: i,j,k,nlat,nlon,ios

r2d=180./3.14159265

print *,'Give quake lat,lon,depth (deg, km):'
read(5,*,iostat=ios)  slat,slon,depth
if(ios.ne.0) then
  print *,'ERROR in plotamps first line input'
  stop
endif  

print *,'Give plot min max latitude:'
read *,lat1,lat2

print *,'Give plot min max longitude:'
read *,lon1,lon2

print *,'Give strike,dip,rake (deg):'
read(5,*,iostat=ios) strike,dip,rake
if(ios.ne.0) then
  print *,'ERROR in rake,dip strike input'
  stop
endif  

vsrc=Pvel(depth)
up=0.           ! downward directed P wave

dlat=FG*(lat2-lat1)
dlon=FG*(lon2-lon1)
nlat=nint((lat2-lat1)/dlat)
nlon=nint((lon2-lon1)/dlon)
print *,'Plotting ',nlat,'X',nlon,' values spaced:',dlat,dlon
if(nlat>NGRID .or. nlon>NGRID) stop 'Increase NGRID'

open(1,file='Prad.xy')

rmax=0.
pmax=0.999*(6371.-depth)/vsrc            ! evanescent wave if larger p
print *,'lon1, lon2=',lon1,lon2,'  lat1,lat2=',lat1,lat2
! write(13,'(a)') '     lat     lon    dist   angle     Prad'
do j=1,nlon
  xlon=lon1+(j-1)*dlon
  do i=1,nlat
    xlat=lat1+(i-1)*dlat
    call delaz(xlat,xlon,slat,slon,dist,azis,azie)
    p=min(pmax,slw(dist,depth))
    call raddc(depth,strike,dip,rake,p,azie,vsrc,up,radp,radsv,radsh)
    prad(i,j)=radp
    ! write(13,'(5f8.2)') xlat,xlon,dist,r2d*p*vsrc/(6371.-depth),radp
    rmax=max(rmax,abs(radp))
  enddo
enddo

prad=prad/rmax          ! normalize to max=+1

do j=1,nlon
  xlon=lon1+(j-1)*dlon
  do i=1,nlat
    xlat=lat1+(i-1)*dlat
    write(1,fmt='(2f8.2,f7.3)') xlon,xlat,prad(i,j)
  enddo
enddo

! Now repeat for smallP (pP)
close(1)

! Check for deep quakes if pP is different from direct P

open(1,file='smallPrad.xy')

up=1.0          ! indicates pP
rmax=0.
! write(13,'(a)') '     lat     lon    dist   angle    pPrad'
do j=1,nlon
  xlon=lon1+(j-1)*dlon
  do i=1,nlat
    xlat=lat1+(i-1)*dlat
    call delaz(xlat,xlon,slat,slon,dist,azis,azie)
    p=min(pmax, slwsmallP(dist,depth))
    call raddc(depth,strike,dip,rake,p,azie,vsrc,up,radp,radsv,radsh)
    ! set radp to zero if pP is in reality direct P
    if(dist<14+0.02*depth .and. depth>35.) radp=0.
    if(depth>600. .and. dist<26.+0.12*(depth-600) ) radp=0.
    prad(i,j)=radp
    rmax=max(rmax,abs(radp))
    ! write(13,'(5f8.2)') xlat,xlon,dist,r2d*p*vsrc/(6371.-depth),radp
  enddo
enddo

prad=prad/rmax          ! normalize to max=+1

do j=1,nlon
  xlon=lon1+(j-1)*dlon
  do i=1,nlat
    xlat=lat1+(i-1)*dlat
    write(1,fmt='(2f8.2,f7.3)') xlon,xlat,prad(i,j)
  enddo
enddo

print *,'Radiation pattern is in Prad.xy, plot with gmtPrad'

end

real function Pvel(depth)

! returns P velocity at depth

implicit none

real*4, intent(in) :: depth
real*4 :: vp(13),dvp(13)
integer :: i,j,k

data dvp/0.,15.,15.,34.,34.,220.,220.,400.,400.,510.,600.,670.,770./
data vp/5.80,5.80,6.80,6.80,8.02,7.80,8.56,8.91,9.13,9.70,10.16,10.27, &
  10.27/

i=2
do while(dvp(i)<depth .and. i<13)
  i=i+1
enddo

if(abs(dvp(i)-depth)<0.001) then
  Pvel=vp(i)
else  
  Pvel=vp(i-1)+(vp(i)-vp(i-1))*(depth-dvp(i-1))/(dvp(i)-dvp(i-1))
endif

return
end function Pvel
