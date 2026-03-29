c234567
c     N(>40000 ?) is the size of array for spline interpolation.
c     range2 is the range of distance between particles in pixel scale,
      subroutine pspline
	  include 'spline.inc'
      integer i
      real e
      double precision a(4)
c     data a/0.96473D-06,0.61076D+00,0.14498D-07,-0.30355D+00/
c	  data a/-1.0D0,0.61076D+00,-0.0150280389D0,-0.30355D+00/
c     data a/0.60491D-04  0.95870D+00  0.39900D-05 -0.60338D+00/
 	  data a/-1.0D0,0.9587D+00,-0.0659602255D0,-0.60338D0/
c     real diff(3,0:Nrange),slope(3,0:Nrange)
      real diff(2,3,0:Nrange)
	  real ran2nran
      double precision a2x,dtanha2x,expdxx
      double precision f_a,f_b,fp_a,fp_b,x
      double precision ar,br,arprime,brprime
      double precision plmf_1,plmf_2,plmf_3,plmf_4
      double precision pmf1_1,pmf1_2,pmf1_3,pmf1_4
      double precision pmf2_1,pmf2_2,pmf2_3,pmf2_4
      double precision pmf3_1,pmf3_2,pmf3_3,pmf3_4
      double precision xstep
	  common /spline/diff,ran2nran

	  ran2nran = range2/Nrange

      diff(1,1,0) = -0.95870E+00
      diff(1,2,0) =  0.22775E+00
      diff(1,3,0) = -0.27274E+00

      do i = 1, Nrange
         x = (dble(i)*dble(range2)/dble(Nrange)) ! here x is distance
c        if(i.eq.0) x = (1.D-10) ! to avoid 1./0.
c----------------------------------------------------------------------------
c---     pmf1_1 --> multipole expansion of the force of dtanh(x)/x potential
c        if(i .eq. 3709028) print *, i,x,a(1),a(2),dtanh(a(2)*x)
c        if(i .eq. 3709000) print *, i,x,dtanh(a(2)*x),dcosh(a(2)*x),
c    &                          dsinh(a(2)*x)
         pmf1_1 = a(1)*dtanh(a(2)*x)/x
         pmf1_2 = pmf1_1*f_a(a,x)/x/x
         pmf1_3 = pmf1_1*4.d0*f_b(a,x)/x/x/x/x
c---     pmf2_1 --> multipole expansion of the force of (x+x^3)*exp(x^2)
c                   potential
         expdxx = dexp(a(4)*x*x)
         pmf2_1 = a(3)*x*x*expdxx
         pmf2_2 = a(3)*(1.d0+a(4)*x*x)*expdxx
         pmf2_3 = 4.d0*a(3)*a(4)*(1.d0+0.5d0*a(4)*x*x)*expdxx
c---     diff is the potential correction term
         diff(1,1,i) = pmf1_1+pmf2_1
         diff(1,2,i) = pmf1_2+pmf2_2
         diff(1,3,i) = pmf1_3+pmf2_3
c        diff(1,1,i) = -1./sqrt(x*x+1.10)
c        diff(1,2,i) = +0.5/(x*x+1.10)/sqrt(x*x+1.10)
c        diff(1,3,i) = -1.5/(x*x+1.10)**2/sqrt(x*x+1.10)
c        diff(1,1,i) = -1./sqrt(x*x+0.10)
c        diff(1,2,i) = +0.5/(x*x+0.10)/sqrt(x*x+0.10)
c        diff(1,3,i) = -1.5/(x*x+0.10)**2/sqrt(x*x+0.10)
c----------------------------------------------------------------------------
      enddo
c	  open(1,file='scale.dat')
c	  do i = 0,Nrange
c	     x = (dble(i)*dble(range2)/dble(Nrange))
c		 write(1,*) x,diff(1,1,i),diff(1,2,i),diff(1,3,i)
c	  enddo
c	  close(1)
c	  call MPI_Finalize()
c	  stop
100   format(6(1x,e12.5))
      do i = 0, Nrange-1
         x = (dble(i)*dble(range2)/dble(Nrange))
c        if(i.eq.0) x = (1.D-10)
         xstep = dble(i+1)*dble(range2)/dble(Nrange)-x
         diff(2,1,i) = (diff(1,1,i+1)-diff(1,1,i))/xstep
         diff(2,2,i) = (diff(1,2,i+1)-diff(1,2,i))/xstep
         diff(2,3,i) = (diff(1,3,i+1)-diff(1,3,i))/xstep
c        slope(4,i) = (diff(4,i+1)-diff(4,i))/xstep
      enddo
      return
      end
ccc
      double precision function f_a(a,x)
      double precision a(4)
      double precision x
	  if(x .gt. 100) then
    	  f_a = -0.5d0
	  else
          f_a = 0.5d0*( a(2)*x/dsinh(a(2)*x)/dcosh(a(2)*x) - 1.d0 )
	  endif
      return
      end
ccc
      double precision function fp_a(a,x)
      double precision a(4)
      double precision x
	  if(x .gt. 100.) then
	      fp_a = 0.D0
	  else 
          fp_a = a(2)/2.d0/dsinh(a(2)*x)/dcosh(a(2)*x)*
     &       ( 1.d0-a(2)*x*(dtanh(a(2)*x)+1.d0/dtanh(a(2)*x)) )
	  endif
      return
      end
ccc
      double precision function f_b(a,x)
      double precision a(4)
      double precision x
	  if(x .gt. 100.) then
	      f_b = 0.125D0* 3.D0
	  else 
          f_b = 0.125d0*( 3.d0-a(2)*x*(2.d0*a(2)*x*dtanh(a(2)*x)+3.d0)/
     &                dsinh(a(2)*x)/dcosh(a(2)*x) )
	  endif
      return
      end
      double precision function fp_b(a,x)
      double precision a(4)
      double precision x
	  if(x .gt. 100.) then
    	  fp_b = 0.D0
	  else
          fp_b = -a(2)*0.125d0/dsinh(a(2)*x)/dcosh(a(2)*x)*
     &      ( a(2)*x*dtanh(a(2)*x)-(2.d0*a(2)*x*dtanh(a(2)*x))**2+
     &        3.d0-3.d0*a(2)*x/dtanh(a(2)*x) )
	  endif
      return
      end
