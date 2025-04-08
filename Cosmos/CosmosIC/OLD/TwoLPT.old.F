c234567
      subroutine SelfDiffPhi(phi,ixyz)
      IMPLICIT none
      include 'mpif.h'
      include 'fftw_f77.i'
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid,ierror
      real phi(2*(nx/2+1),ny,local_nz)
      integer ngy,ngh,nghx,i,j,k,ixyz,ng,ngx
      real pi,twopi,fourpi2,wrk
      double precision rincube
      integer ky,kz
      real xyz,xyy,tmp1,tmp2
      real rngx,rngy,rngz

      rngx = nx
      rngy = ny
      rngz = nz
      xyz = rngx/rngz
      xyy = rngx/rngy

      rincube = 1.d0/dble(nx)/dble(ny)/dble(nz)

c     if(kflag .ne. 1) then
c     call rfftwnd_f77_mpi(plan,1,phi,wrk,USE_WORK,FFTW_NORMAL_ORDER)
c     endif

      pi = 4.*atan(1.);

      twopi = 2.*pi
      fourpi2 = twopi**2

      ng = nx
      ngx = nx
      nghx = nx/2
      ngy = ny
      ngh = nghx

      if(ixyz .eq. 1) then
c$omp parallel do private(i,j,k,ky,kz,tmp1,tmp2)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         tmp1 = phi(2*i+1,j,k)
         tmp2 = phi(2*i+2,j,k)
         phi(2*i+1,j,k) = -twopi*i*tmp2
         phi(2*i+2,j,k) =  twopi*i*tmp1
      enddo
      enddo
      enddo
      else if(ixyz .eq.2) then
c$omp parallel do private(i,j,k,ky,kz,tmp1,tmp2)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         tmp1 = phi(2*i+1,j,k)
         tmp2 = phi(2*i+2,j,k)
         phi(2*i+1,j,k) = -twopi*ky*tmp2
         phi(2*i+2,j,k) =  twopi*ky*tmp1
      enddo
      enddo
      enddo
      else if(ixyz .eq.3) then
c$omp parallel do private(i,j,k,ky,kz,tmp1,tmp2)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         tmp1 = phi(2*i+1,j,k)
         tmp2 = phi(2*i+2,j,k)
         phi(2*i+1,j,k) = -twopi*kz*tmp2
         phi(2*i+2,j,k) =  twopi*kz*tmp1
      enddo
      enddo
      enddo
      else
         print *, 'Please check input ixyz flag 1/2/3 ', ixyz
         call MPI_Finalize(MPI_COMM_WORLD,ierror)
         stop
      endif

c     call rfftwnd_f77_mpi(iplan,1,phi,wrk,USE_WORK,
c    &                  FFTW_NORMAL_ORDER)
      call fftwbackward_(phi)

      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi(i,j,k) = phi(i,j,k)*rincube
      enddo
      enddo
      enddo

      return
      end

c234567
      subroutine DiffPhi(phi,dphidx,ixyz)
      IMPLICIT none
      include 'mpif.h'
      include 'fftw_f77.i'
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid,ierror
      real phi(2*(nx/2+1),ny,local_nz)
      real dphidx(2*(nx/2+1),ny,local_nz)
      integer ngy,ngh,nghx,i,j,k,ixyz,ng,ngx
      real pi,twopi,fourpi2,wrk
      double precision rincube
      integer ky,kz
      real xyz,xyy,tmp1,tmp2
      real rngx,rngy,rngz

      rngx = nx
      rngy = ny
      rngz = nz
      xyz = rngx/rngz
      xyy = rngx/rngy

      rincube = 1.d0/dble(nx)/dble(ny)/dble(nz)

c     if(kflag .ne. 1) then
c     call rfftwnd_f77_mpi(plan,1,phi,wrk,USE_WORK,FFTW_NORMAL_ORDER)
c     endif

      pi = 4.*atan(1.);

      twopi = 2.*pi
      fourpi2 = twopi**2

      ng = nx
      ngx = nx
      nghx = nx/2
      ngy = ny
      ngh = nghx

      if(ixyz .eq. 1) then
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         dphidx(2*i+1,j,k) = -twopi*i*phi(2*i+2,j,k)
         dphidx(2*i+2,j,k) =  twopi*i*phi(2*i+1,j,k)
      enddo
      enddo
      enddo
      else if(ixyz .eq.2) then
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         dphidx(2*i+1,j,k) = -twopi*ky*phi(2*i+2,j,k)
         dphidx(2*i+2,j,k) =  twopi*ky*phi(2*i+1,j,k)
      enddo
      enddo
      enddo
      else if(ixyz .eq.3) then
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         dphidx(2*i+1,j,k) = -twopi*kz*phi(2*i+2,j,k)
         dphidx(2*i+2,j,k) =  twopi*kz*phi(2*i+1,j,k)
      enddo
      enddo
      enddo
      else
         print *, 'Please check input ixyz flag 1/2/3 ', ixyz
         call MPI_Finalize(MPI_COMM_WORLD,ierror)
         stop
      endif

c     call rfftwnd_f77_mpi(iplan,1,dphidx,wrk,USE_WORK,
c    &                  FFTW_NORMAL_ORDER)
      call fftwbackward_(dphidx)

      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         dphidx(i,j,k) = dphidx(i,j,k)*rincube
      enddo
      enddo
      enddo

      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine twoLPT(phi1,work1,work2,phi2,myid,nid,fNL,gNL,
     &      boxsize,amax,omep,Growthfactor)
      IMPLICIT none
      include 'mpif.h'
      include 'fftw_f77.i'
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid
      real phi1(2*(nx/2+1),ny,local_nz),fNL,gNL,boxsize,amax,omep
      real phi2(2*(nx/2+1),ny,local_nz)
      real work1(2*(nx/2+1),ny,local_nz)
      real work2(2*(nx/2+1),ny,local_nz)
      integer nghx,i,j,k
      real pi,twopi,fourpi2,wrk,Growthfactor
      double precision rincube,sqrrincube
      real rngx,rngy,rngz,xyz,xyy
      integer ngx,ngy,ngz,ngh,ky,ng,kz,ierror
      double precision wave2,invwave2,phistd,phimean,tphistd,tphimean
      double precision VelPot2RealPot,Cvel,iVelPot2RealPot
      real fNLzcorr,gNLzcorr



      pi = 4.*atan(1.);
      Cvel = 2.99792458D5


      twopi = 2.*pi
      fourpi2 = twopi**2;

      nghx = nx/2
      ngy = ny
      ng = nx

      rngx = nx
      rngy = ny
      rngz = nz
      xyz = rngx/rngz
      xyy = rngx/rngy
      ngx = nx
      ngy = ny
      ngh = nghx



      rincube = 1.d0/dble(nx)/dble(ny)/dble(nz)
      sqrrincube = rincube**2

c----
      if(fNL .ne.0 .or. gNL .ne.0) then 
      VelPot2RealPot = 3./2.*(100*boxsize/Cvel)**2*omep*amax
      iVelPot2RealPot = 1.d0/VelPot2RealPot
      do k = 1, local_nz
      do j = 1, ny
      do i = 1, nx
         phi1(i,j,k) = phi1(i,j,k)*VelPot2RealPot
      enddo
      enddo
      enddo
      call rfftwnd_f77_mpi(iplan,1,phi1,wrk,USE_WORK, 
     &            FFTW_NORMAL_ORDER)
      do k=1,local_nz
      do j=1,ny
      do i=1,nx
      phi1(i,j,k) = phi1(i,j,k)*rincube
      enddo
      enddo
      enddo
      phistd = 0
      phimean = 0
      do k=1,local_nz
      do j=1,ny
      do i=1,nx
      phistd = phistd + phi1(i,j,k)**2
      phimean = phimean + phi1(i,j,k)
      enddo
      enddo
      enddo
      call MPI_Reduce(phistd,tphistd,1,MPI_DOUBLE_PRECISION,
     &         MPI_SUM,0,MPI_COMM_WORLD,ierror)
      phistd = tphistd
      call MPI_Bcast(phistd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           ierror)

      call MPI_Reduce(phimean,tphimean,1,MPI_DOUBLE_PRECISION,
     &         MPI_SUM,0,MPI_COMM_WORLD,ierror)
      phimean = tphimean
      phimean = phimean*rincube
      call MPI_Bcast(phimean,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &          ierror)
      phistd = dsqrt(phistd*rincube - phimean**2)
      fNLzcorr = 1./growthfactor/amax
      gNLzcorr = (1./growthfactor/amax)**2
      do k=1,local_nz
      do j=1,ny
      do i=1,nx
      phi1(i,j,k) = phi1(i,j,k)
     &        +fNLzcorr*fNL*(phi1(i,j,k)**2-phistd**2)
     &        +gNLzcorr*gNL*(phi1(i,j,k)**3-3*phistd**2*phi1(i,j,k))
      enddo
      enddo
      enddo
      call rfftwnd_f77_mpi(plan,1,phi1,wrk,USE_WORK,
     &            FFTW_NORMAL_ORDER)

      do k = 1, local_nz
      do j = 1, ny
      do i = 1, nx
         phi1(i,j,k) = phi1(i,j,k)*iVelPot2RealPot
      enddo
      enddo
      enddo
      endif
c----



c-----   phi,11 * phi,22 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work1(2*i+1,j,k) = -fourpi2*i*i*phi1(2*i+1,j,k) ! phi,11
         work1(2*i+2,j,k) = -fourpi2*i*i*phi1(2*i+2,j,k) ! phi,11
         work2(2*i+1,j,k) = -fourpi2*ky*ky*phi1(2*i+1,j,k) ! phi,22
         work2(2*i+2,j,k) = -fourpi2*ky*ky*phi1(2*i+2,j,k) ! phi,22
      enddo
      enddo
      enddo
      call fftwbackward_(work1)
      call fftwbackward_(work2)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = work1(i,j,k)*work2(i,j,k)*sqrrincube
      enddo
      enddo
      enddo




c-----   phi,11 * phi,33 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work2(2*i+1,j,k) = -fourpi2*kz*kz*phi1(2*i+1,j,k) ! phi,33
         work2(2*i+2,j,k) = -fourpi2*kz*kz*phi1(2*i+2,j,k) ! phi,33
      enddo
      enddo
      enddo
      call fftwbackward_(work2)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = phi2(i,j,k)+work1(i,j,k)*work2(i,j,k)*sqrrincube
      enddo
      enddo
      enddo


c-----   phi,22 * phi,33 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work1(2*i+1,j,k) = -fourpi2*ky*ky*phi1(2*i+1,j,k) ! phi,22
         work1(2*i+2,j,k) = -fourpi2*ky*ky*phi1(2*i+2,j,k) ! phi,22
      enddo
      enddo
      enddo
      call fftwbackward_(work1)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = phi2(i,j,k)+work1(i,j,k)*work2(i,j,k)*sqrrincube
      enddo
      enddo
      enddo

c-----   phi,21 * phi,21 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work1(2*i+1,j,k) = -fourpi2*ky*i*phi1(2*i+1,j,k) ! phi,21
         work1(2*i+2,j,k) = -fourpi2*ky*i*phi1(2*i+2,j,k) ! phi,21
      enddo
      enddo
      enddo
      call fftwbackward_(work1)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = phi2(i,j,k)-work1(i,j,k)*work1(i,j,k)*sqrrincube
      enddo
      enddo
      enddo

c-----   phi,31 * phi,31 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work1(2*i+1,j,k) = -fourpi2*kz*i*phi1(2*i+1,j,k) ! phi,31
         work1(2*i+2,j,k) = -fourpi2*kz*i*phi1(2*i+2,j,k) ! phi,31
      enddo
      enddo
      enddo
      call fftwbackward_(work1)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = phi2(i,j,k)-work1(i,j,k)*work1(i,j,k)*sqrrincube
      enddo
      enddo
      enddo

c-----   phi,32 * phi,32 compenent
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         work1(2*i+1,j,k) = -fourpi2*kz*ky*phi1(2*i+1,j,k) ! phi,32
         work1(2*i+2,j,k) = -fourpi2*kz*ky*phi1(2*i+2,j,k) ! phi,32
      enddo
      enddo
      enddo
      call fftwbackward_(work1)
      do k=1,local_nz
c$omp parallel do private(i,j)
      do j=1,ngy
      do i=1,ngx
         phi2(i,j,k) = phi2(i,j,k)-work1(i,j,k)*work1(i,j,k)*sqrrincube
      enddo
      enddo
      enddo



c---  go back to the Fourier Space
      call fftwforward_(phi2)
c-----   Final factor -1/k^2
c$omp parallel do private(i,j,k,ky,kz,wave2,invwave2)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         wave2 = -(i*i + ky*ky + kz*kz)*fourpi2
         invwave2 = 1.d0/wave2
         if(wave2 .lt. 0) then ! be aware of the negative value of wave2 
            phi2(2*i+1,j,k) = phi2(2*i+1,j,k)*invwave2 
            phi2(2*i+2,j,k) = phi2(2*i+2,j,k)*invwave2 
         else
            phi2(2*i+1,j,k) = 0
            phi2(2*i+2,j,k) = 0
         endif
      enddo
      enddo
      enddo

c-----
c=     Checking the data output
c     call fftwbackward(phi1)
c     call fftwbackward(phi2)
c     do k=1,local_nz
c     do j=1,ngy
c     do i=1,ngx
c        phi2(i,j,k) = phi2(i,j,k)*rincube
c        phi1(i,j,k) = phi1(i,j,k)*rincube
c     enddo
c     enddo
c     enddo
c     print *, phi1(20,20,4), phi1(20,20,3), phi2(20,20,4),phi2(20,20,3)
c     call MPI_Finalize(ierror)
c     stop
c----

      return
      end

