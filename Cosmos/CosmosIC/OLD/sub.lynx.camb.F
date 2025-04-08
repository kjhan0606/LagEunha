c---  Reads position or position & velocity of particles
c---  npkid = number of particles in the subset
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.

C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.

      SUBROUTINE qsimp(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-5, JMAX=40)
CU    USES trapzd
      INTEGER j
      REAL os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        if (abs(s-os).lt.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      pause 'too many steps in qsimp'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.



      function power(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  biquadratic cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and 
c---  transfer function
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      if(rk.gt.rknyq.or.rk.eq.0.)then
      power=0.
      return
      else
c      power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if



c***  Transfer function (BBKS 1986, p60) 
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of 
c---  relativistic neutrinos plus photons and negligible baryon
      if(dark.eq.'cdm')then
          power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)

c---  Hot dark matter, adiabatic fluctuations with 2 species of
c---  relativiatic neutrinos
      else if(dark.eq.'hdm')then 
          damp=2.6*qq
          power=power*exp(-0.16*damp-0.5*damp**2)/(1+1.6*qq
     &     +(4.*qq)**1.5+(0.92*qq)**2)
      end if
ccccccccccccccccccccc a cutoff at k=8
c     if(rk/2./3.141592654.gt.8)power = 0
      return
      end


c---  Integrand for calculating (dM/M)**2 in top hat window
      function powint(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum(without knyquist cutoff) times volumn factor
      powint = fopi*xk*xk*(power1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      powint = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      powint = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end

c---  Integrand for calculating (dM/M)**2 in top hat window
      function bpowint(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/


c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum(without knyquist cutoff) times volumn factor
      powint = fopi*xk*xk*(bpower1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      bpowint = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      bpowint = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end
c---  Integrand for calculating integral(P(k)*k^2) in top hat window
      function bpow2int(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum times volumn factor
      powint = fopi*xk**2*rk**2*(bpower1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      bpow2int = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      bpow2int = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end


c---  Integrand for calculating integral(P(k)*k^2) in top hat window
      function pow2int(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum times volumn factor
      powint = fopi*xk**2*rk**2*(power1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      pow2int = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      pow2int = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end

      function power1(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  biquadratic cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and 
c---  transfer function
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      power=rk**powh*exp(-(rk*smooth)**2)

c***  Transfer function (BBKS 1986, p60) 
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of 
c---  relativistic neutrinos plus photons and negligible baryon
      if(dark.eq.'cdm')then
      if(qq.lt.0.002)then
          power=power*(1-1.17*qq+(1.351*qq)**2)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)
      else
          power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)
      end if

c---  Hot dark matter, adiabatic fluctuations with 2 species of
c---  relativiatic neutrinos
      else if(dark.eq.'hdm')then 
          damp=2.6*qq
          power=power*exp(-0.16*damp-0.5*damp**2)/(1+1.6*qq
     &     +(4*qq)**1.5+(0.92*qq)**2)
      end if
      power1 = power
ccccccccccccccccccccccccccc
c     if(rk/2./3.141592654.gt.8)power = 0
      return
      end

#ifdef EH98
      function bpower(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and
c---  transfer function
c---  Ref. Eisenstein and Hu, 1998, ApJ, 496,605
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      save in1,cmbt27,o,ob,oc,zeq,rkeq,b1,b2,zd,Rd,Req,sound,rksk
      save y,gy,alphab,a11,a22,alphac,b11,b22,betac,betab,betan
      data in1/0/

c---  Important cosmological scales and epochs
c---  Tcmb = 2.728+-0.0004K (Fixen et al. 1996) adopted
c---  rkeq = scale of particle horizon at equality epoch in unit of [h Mpc^-1]
c---  zd = drag epoch (baryon released from the drag of photons)
c---  sound = sound horizon at drag epoch in [h^-1 Mpc]
c---  rksk = Silk damping scale in unit of [h Mpc^-1]
      if(in1.eq.0)then
      cmbt27 = 2.728/2.7
      o   = omep*h**2
      ob  = omepb*h**2
      oc  = (omep-omepb)*h**2

      zeq = 2.50e4*o/cmbt27**4
      rkeq = 7.46e-2*omep*h/cmbt27**2
      b1  = 0.313/o**0.419*(1+0.607*o**0.674)
      b2  = 0.238*o**0.223
      zd  = 1291*o**0.251/(1+0.659*o**0.828)*(1+b1*ob**b2)
      Rd  = 31.5*ob/cmbt27**4/(zd/1000.)
      Req = 31.5*ob/cmbt27**4/(zeq/1000.)
      sound=2./3./rkeq*sqrt(6./Req)
     &      *alog((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
      rksk = 1.6*ob**0.52*o**0.73*(1+(10.4*o)**(-0.95))/h

      y = (1+zeq)/(1+zd)
      gy = y*(-6*sqrt(1+y)+(2+3*y)*alog((sqrt(1+y)+1)/(sqrt(1+y)-1)))
      alphab = 2.07*rkeq*sound/(1+Rd)**(3./4.)*gy

      a11  = (46.9*o)**0.670*(1+(32.1*o)**(-0.532))
      a22  = (12.0*o)**0.424*(1+(45.0*o)**(-0.582))
      alphac = a11**(-ob/o)*a22**(-ob/o*ob/o*ob/o)
      b11  = 0.944/(1+(458*o)**(-0.708))
      b22  = 1./(0.395*o)**0.0266
      betac = 1./(1+b11*((oc/o)**b22-1))
      betab = 0.5+ob/o+(3-2*ob/o)*sqrt((17.2*o)**2+1)
      betan = 8.41*o**0.435

      in1 = 1
      end if

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      if(rk.gt.rknyq.or.rk.eq.0.)then
      bpower=0.
      return
      else
c     power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if

c***  Baryonic effects included
      q   = rk*omehs*cmbt27**2
      rks = rk/ssize*sound
      f   = 1./(1+(rks/5.4)**4)
      transc = f*tnot(q,1.0,betac)+(1-f)*tnot(q,alphac,betac)

      soundt = sound/(1+(betan/rks)**3)**(1./3.)
      rkstilde = rk/ssize*soundt
      rkrksk = rk/ssize/rksk
      transb = (tnot(q,1.0,1.0)/(1+(rks/5.2)**2)
     & +alphab/(1.+(betab/rks)**3)*exp(-rkrksk**1.4))
     & *sin(rkstilde)/rkstilde
      
      bpower = power*abs(ob*transb+oc*transc)/o

c***  Transfer function (BBKS 1986, p60)
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of
c---  relativistic neutrinos plus photons and negligible baryon
c         power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
c    &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)

      return
      end

      function tnot(q,alphac,betac)
c     data e/2.718281828/
      e = 2.718281828
      c = 14.2/alphac+386./(1+69.9*q**1.08)
      tnot = alog(e+1.8*betac*q)/(alog(e+1.8*betac*q)+c*q**2)
      return
      end

      function bpower1(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and
c---  transfer function
c---  Ref. Eisenstein and Hu, 1998, ApJ, 496,605
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      save in1,cmbt27,o,ob,oc,zeq,rkeq,b1,b2,zd,Rd,Req,sound,rksk
      save y,gy,alphab,a11,a22,alphac,b11,b22,betac,betab,betan
      data in1/0/

c---  Important cosmological scales and epochs
c---  Tcmb = 2.728+-0.0004K (Fixen et al. 1996) adopted
c---  rkeq = scale of particle horizon at equality epoch in unit of [h Mpc^-1]
c---  zd = drag epoch (baryon released from the drag of photons)
c---  sound = sound horizon at drag epoch in [h^-1 Mpc]
c---  rksk = Silk damping scale in unit of [h Mpc^-1]
      if(in1.eq.0)then
      cmbt27 = 2.728/2.7
      o   = omep*h**2
      ob  = omepb*h**2
      oc  = (omep-omepb)*h**2

      zeq = 2.50e4*o/cmbt27**4
      rkeq = 7.46e-2*omep*h/cmbt27**2
      b1  = 0.313/o**0.419*(1+0.607*o**0.674)
      b2  = 0.238*o**0.223
      zd  = 1291*o**0.251/(1+0.659*o**0.828)*(1+b1*ob**b2)
      Rd  = 31.5*ob/cmbt27**4/(zd/1000.)
      Req = 31.5*ob/cmbt27**4/(zeq/1000.)
      sound=2./3./rkeq*sqrt(6./Req)
     &      *alog((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
      rksk = 1.6*ob**0.52*o**0.73*(1+(10.4*o)**(-0.95))/h

      y = (1+zeq)/(1+zd)
      gy = y*(-6*sqrt(1+y)+(2+3*y)*alog((sqrt(1+y)+1)/(sqrt(1+y)-1)))
      alphab = 2.07*rkeq*sound/(1+Rd)**(3./4.)*gy

      a11  = (46.9*o)**0.670*(1+(32.1*o)**(-0.532))
      a22  = (12.0*o)**0.424*(1+(45.0*o)**(-0.582))
      alphac = a11**(-ob/o)*a22**(-ob/o*ob/o*ob/o)
      b11  = 0.944/(1+(458*o)**(-0.708))
      b22  = 1./(0.395*o)**0.0266
      betac = 1./(1+b11*((oc/o)**b22-1))
      betab = 0.5+ob/o+(3-2*ob/o)*sqrt((17.2*o)**2+1)
      betan = 8.41*o**0.435

      in1 = 1
      end if

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
c     if(rk.gt.rknyq.or.rk.eq.0.)then
      if(rk.eq.0.)then
      bpower1=0.
      return
      else 
c     power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if

c***  Baryonic effects included
      q   = rk*omehs*cmbt27**2
      rks = rk/ssize*sound
      f   = 1./(1+(rks/5.4)**4)
      transc = f*tnot(q,1.0,betac)+(1-f)*tnot(q,alphac,betac)

      soundt = sound/(1+(betan/rks)**3)**(1./3.)
      rkstilde = rk/ssize*soundt
      rkrksk = rk/ssize/rksk
      transb = (tnot(q,1.0,1.0)/(1+(rks/5.2)**2)
     &         +alphab/(1.+(betab/rks)**3)*exp(-rkrksk**1.4))
     &         *sin(rkstilde)/rkstilde
      
      bpower1 = power*abs(ob*transb+oc*transc)/o
      return
      end
#else
c---------------
c     rk is the 2pi * u where u is the wave number
c     rks is the k = 2pi * u /boxsize, real wave vector
      function bpower1(rk)
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      data in1/0/
      real matterpk(8192),rkarray(8192),Trans(8192),minkh,dlnkh
      integer point,nmatterpk,nn,nmin,nmax
      save in1
      common /cambpower/ matterpk,rkarray,nmatterpk,minkh,dlnkh
      real asciirk(8192), asciipk(8192)
      integer nascii,powreadflag
      common /ASCIIPOWER/powreadflag,nascii,asciirk,asciipk
      rks = rk/ssize
      if(in1 .eq.0) then
         nmatterpk = 8192
         call GetPower(minkh,dlnkh,matterpk,nmatterpk,Trans)
         if(powreadflag .eq.0 .or. powreadflag .eq. 1) then
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
               matterpk(i) = sqrt(matterpk(i))
            enddo
         else if(powreadflag .eq. 2) then
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
            enddo
            call getpkfromascii(rkarray,matterpk,nmatterpk,asciirk,
     &               asciipk,nascii,Trans)
            do i = 1, nmatterpk
               matterpk(i) = sqrt(matterpk(i))
            enddo
         endif
         in1 = 1
      endif

      if(rk .eq. 0 ) then
          bpower1 = 0
          return 
      endif 
      if(rks .lt. rkarray(1)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Less than lowest boundary ', rks, rkarray(1)
         print *,'############################################'
         print *,'############################################'
         stop
      else if(rks .gt. rkarray(nmatterpk)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Larger than highest boundary ',rks,rkarray(nmatterpk)
         print *,'############################################'
         print *,'############################################'
         stop
      endif
      point = (alog(rks/minkh)/dlnkh+1)
      nmin = max(1,point-2)
      nmax = min(nmatterpk,point+2)
      nn = nmax-nmin+1
      call polint(rkarray(nmin),matterpk(nmin),nn,rks,y,dy)
      bpower1 =y
      return
      end
c---------------
c     bpower = sqrt(PK)
c---------------
      function bpower(rk)
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      data in1/0/
      real matterpk(8192),rkarray(8192),Trans(8192),minkh,dlnkh
      integer point,nmatterpk,nn,nmin,nmax
      save in1
      common /cambpower/ matterpk,rkarray,nmatterpk,minkh,dlnkh
      real asciirk(8192), asciipk(8192)
      integer nascii,powreadflag
      common /ASCIIPOWER/powreadflag,nascii,asciirk,asciipk
      rks = rk/ssize
      if(in1 .eq.0) then
         nmatterpk = 8192
         call GetPower(minkh,dlnkh,matterpk,nmatterpk,Trans)
         if(powreadflag .eq.0 .or. powreadflag .eq. 1) then
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
               matterpk(i) = sqrt(matterpk(i))
            enddo
         else if(powreadflag .eq. 2) then
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
            enddo
            call getpkfromascii(rkarray,matterpk,nmatterpk,asciirk,
     &               asciipk,nascii,Trans)
            do i = 1, nmatterpk
               matterpk(i) = sqrt(matterpk(i))
            enddo
         endif
         in1 = 1
      endif
      if(rk .eq. 0 .or. rk .gt. rknyq) then
          bpower = 0
          return 
      endif 
      if(rks .lt. rkarray(1)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Less than lowest boundary ', rks, rkarray(1)
         print *,'############################################'
         print *,'############################################'
         stop
      else if(rks .gt. rkarray(nmatterpk)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Larger than highest boundary ',rks,rkarray(nmatterpk)
         print *,'############################################'
         print *,'############################################'
         stop
      endif
      point = (alog(rks/minkh)/dlnkh+1)
      nmin = max(1,point-2)
      nmax = min(nmatterpk,point+2)
      nn = nmax-nmin+1
      call polint(rkarray(nmin),matterpk(nmin),nn,rks,y,dy)
      bpower = y
      return
      end
c###############
      subroutine getpkfromascii(rk,pk,npk,ark,apk,napk,Ts2)
      integer npk,napk,n1,n2,nn
      real rk(npk),pk(npk),Ts2(npk),ark(napk),apk(napk)
      real y
c --- change scaling: Be noted that the number of effective digits might be changed.
      do i = 1, npk
         rk(i) = alog10(rk(i))
      enddo
      do i = 1, napk
         ark(i) = alog10(ark(i))
         apk(i) = alog10(apk(i))
      enddo
      do i = 1, npk
         if(rk(i) .lt. ark(1)) then
            y = (apk(2)-apk(1))/(ark(2)-ark(1))*(rk(i)-ark(1)) + apk(1)
         else if(rk(i) .gt. ark(napk)) then
            y = (apk(napk)-apk(napk-1))/(ark(napk)-ark(napk-1))*
     &             (rk(i)-ark(napk)) + apk(napk)
         else
            do j = 1, napk-1
               if(rk(i).ge. ark(j) .and. rk(i).le.ark(j+1)) goto 10
            enddo
10          continue
            n1 = j-2
            n2 = j+2
            n1 = max(n1,1)
            n2 = min(n2,napk)
            nn = n2 - n1 + 1
            call polint(ark(n1),apk(n1),nn,rk(i),y,dy)
         endif
         pk(i) = 10.**(y+rk(i))*Ts2(i)
c        print *, 10**rk(i), pk(i),y,Ts2(i)
      enddo
c     stop
c --- change scaling: Be noted that the number of effective digits might be changed.
      do i = 1, npk
         rk(i) = 10.**(rk(i))
      enddo
      do i = 1, napk
         ark(i) = 10.**(ark(i))
         apk(i) = 10.**(apk(i))
      enddo
  
      return
      end
#endif
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran30(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran30,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran30=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran31(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran31,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran31=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran32(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran32,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran32=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran33(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran33,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran33=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran34(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran34,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran34=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran35(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran35,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran35=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran36(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran36,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran36=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran37(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran37,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran37=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION ran38(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran38,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran38=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION hypgeo(a,b,c,z)
      COMPLEX hypgeo,a,b,c,z
      REAL EPS
      PARAMETER (EPS=1.e-6)
CU    USES bsstep,hypdrv,hypser,odeint
      INTEGER kmax,nbad,nok,kount
      EXTERNAL bsstep,hypdrv
      COMPLEX z0,dz,aa,bb,cc,y(2)
      COMMON /hypg/ aa,bb,cc,z0,dz

      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      REAL dxsav,xp(KMAXX),yp(NMAX,KMAXX)
      COMMON /path/ kmax,kount,dxsav,xp,yp

      kmax=0
      if (real(z)**2+aimag(z)**2.le.0.25) then
        call hypser(a,b,c,z,hypgeo,y(2))
        return
      else if (real(z).lt.0.) then
        z0=cmplx(-0.5,0.)
      else if (real(z).le.1.0) then
        z0=cmplx(0.5,0.)
      else
        z0=cmplx(0.,sign(0.5,aimag(z)))
      endif
      aa=a
      bb=b
      cc=c
      dz=z-z0
      call hypser(aa,bb,cc,z0,y(1),y(2))
      call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
      hypgeo=y(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE hypser(a,b,c,z,series,deriv)
      INTEGER n
      COMPLEX a,b,c,z,series,deriv,aa,bb,cc,fac,temp
      deriv=cmplx(0.,0.)
      fac=cmplx(1.,0.)
      temp=fac
      aa=a
      bb=b
      cc=c
      do 11 n=1,1000
        fac=fac*aa*bb/cc
        deriv=deriv+fac
        fac=fac*z/n
        series=temp+fac
        if (series.eq.temp) return
        temp=series
        aa=aa+1.
        bb=bb+1.
        cc=cc+1.
11    continue
      pause 'convergence failure in hypser'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE hypdrv(s,y,dyds)
      REAL s
      COMPLEX y(2),dyds(2),aa,bb,cc,z0,dz,z
      COMMON /hypg/ aa,bb,cc,z0,dz
      z=z0+s*dz
      dyds(1)=y(2)*dz
      dyds(2)=(aa*bb*y(1)-(cc-(aa+bb+1.)*z)*y(2))*dz/(z*(1.-z))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
