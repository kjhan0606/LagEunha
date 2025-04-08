      function growthgen(a,amax,omep,omeplam,wlam)
      complex a1,b1,c1,z,a2,b2,c2,hypgeo


      x = (1-omep)/omep * (a/amax)**(-3*wlam)


      a1 = cmplx(1.5,0)
      b1 = cmplx(-5./6./wlam,0)
      c1 = cmplx(1.-5./6./wlam,0)
      z = cmplx(-x,0)
      a2 = cmplx(1.5,0)
      b2 = cmplx(1.-5./6./wlam,0)
      c2 = cmplx(2.-5./6./wlam,0)


      growthgen = (a/amax)*sqrt(1+x)*(real(hypgeo(a1,b1,c1,z))
     &     + x*5*(1+wlam)/(5-6*wlam)*real(hypgeo(a2,b2,c2, z)) )

      return
      end
c---  Growth factor in the universe with zero curvature
c---  omep= present matter density parameter
c---  omeplam = 1-omep is assumed (zero curvature inflation model)
      function growth(a,amax,omep)
      external coscon
      common/pass/x

      x0 = (1./omep-1)**(1./3.)
      red = amax/a-1
      x = (1./omep-1)**(1./3.)/(1+red)
      call qsimp(coscon,0.,1.,gfac)
      growth = (x/x0)*sqrt(x**3+1)*gfac
      return
      end


      function coscon(y)
      common/pass/x
      coscon = 1./(sqrt(1+x**3*y**1.2))**3
      return
      end


c--- growth factor for non-flat and non-zero cosmological constant
c---
      function growth2(a,amax,om,ol)
      external coscon2
      common/cospar/omep,omeplam

      omep = om
      omeplam = ol

      omepk = 1-omep-omeplam
      aa = a/amax
      growth2 = sqrt(omeplam*aa**3+omepk*aa+omep)/aa**1.5

      call qsimp(coscon2,0.,aa,gfac)
      growth2 = growth2*gfac

      return
      end

      function coscon2(y)
      real omep,omepl,omepk
      common/cospar/omep,omeplam
      omepk = 1-omep-omeplam
      coscon2 = y**1.5/(omeplam*y**3+omepk*y+omep)**1.5
      return
      end

