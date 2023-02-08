      subroutine radamp(depth,xmom,p,azi,vsrc,up,radp,radsv,radsh)

c     routine is used in pltamps
      
c     input xmom: moment tensor (Mrr, Mss, Mee, Mrs, Mre, Mse)
c                            or (m1 , m2,  m3,  m4,  m5,  m6)
c                            or (Mzz, Mxx, Myy, Mzx,-Mzy,-Mxy)
c           depth: source depth in km
c           azi: azimuth (deg) from source
c	    p: slowness (s/rad)
c           vsrc(1),vsrc(2): Vp and Vs at source depth (km/s)
c           up: 1.0 for upward wave (eg pP), else 0.0
c     output rad,radsv,radsh: radiation factors for P,SV and SH waves

c     see Aki&Richards, Quantitative Seismology, Box 9.10 page 467 [and 
c     page 118, both first edition - if you have the 2nd edition explicit
c     eqs can be derived from (4.88), (4.96) and Box 4.4]

c     Caution: the elements of the HRVD moment tensor are ranked and
c     signed differentl from Mxx,Myy etc!
c     HRVD coordinate system is r-up, s-South, e-East, and the 6 elements 
c         are ordered Mrr, Mss, Mee, Mrs, Mre, Mse.

      implicit none
      real*4 depth,p,azi,up,radp,radsv,radsh
      real*4 vsrc(2),xmom(6)
      real*4 Mxx,Myy,Mzz,Mxy,Mxz,Myz
      real*4 rad,phi,rs
      real*4 xi,sinxi,cosxi,sinphi,cosphi,sin2phi,cos2phi
      

      rad=3.141592/180.
      phi=azi*rad
      rs=6371.0-depth

c     translate spherical moment tensor into Cartesian
      Mxx=xmom(2)
      Myy=xmom(3)
      Mzz=xmom(1)
      Mxy=-xmom(6)
      Mxz=xmom(4)
      Myz=-xmom(5)

      
c     compute radiation pattern for P
      sinxi=p*vsrc(1)/rs          ! sine of take-off angle in radians
      if(abs(sinxi).gt.1.0) then
        print *,'Error in radamp: p,vsrc(1)=',p,vsrc(1)
	print *,'Ray evanescent...'
	stop
      endif	
      ! take-off angle xi is 0 for downward vertical take-off
      xi=asin(sinxi)           ! P angle at source depth (i_xi in A&R)
      if(abs(up)<1.0e-8) xi=3.14159265-xi
      cosxi=cos(xi)
      sinphi=sin(phi)
      cosphi=cos(phi)
      sin2phi=sin(2.0*phi)
      cos2phi=cos(2.0*phi)
      radp=sinxi**2*(cosphi**2*Mxx+sin2phi*Mxy+
     #   sinphi**2*Myy-Mzz)+2.0*sinxi*cosxi*
     #  (cosphi*Mxz+sinphi*Myz)+Mzz

c     compute radiation pattern for SV
      sinxi=p*vsrc(2)/rs
      if(abs(sinxi).gt.1.0) then
        print *,'Error in radamp: p,vsrc(2)=',p,vsrc(2)
	print *,'Ray evanescent...'
	stop
      endif	
      xi=asin(sinxi)           ! S angle at source depth
      if(abs(up)<1.0e-8) xi=3.14159265-xi
      cosxi=cos(xi)
      radsv=sinxi*cosxi*(cosphi**2*Mxx+sin2phi*Mxy+
     #     sinphi**2*Myy-Mzz)+(1.0-2.0*sinxi**2)*
     #     (cosphi*Mxz+sinphi*Myz)
        
c     compute radiation pattern for SH
      radsh=sinxi*(0.5*sin2phi*(Myy-Mxx)+cos2phi*Mxy)+
     #     cosxi*(cosphi*Myz-sinphi*Mxz) 

      return
      end
