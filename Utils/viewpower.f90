     !Simple test program to print out sigma_8 as a function of the CDM density
      module MyModule
         use CAMB
         type (CAMBdata) OutData,InData
      end module MyModule

      program viewpower
        use CAMB
        use MyModule
        use InitialPower
        use LambdaGeneral
        implicit none
        integer i,j,k,nmatterpk
        real PI2,rng,maxkhfactor
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
    
        real matterpk(8192),cmatterpk(8192),bmatterpk(8192),rkarray(8192),Ts(8192)
        type(CAMBparams) P !defined in ModelParams in modules.f90
        !type (CAMBdata) OutData
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in
        integer error
        real minkh,dlnkh,pamp,bias,khmax,opamp
        integer npoints
        character(LEN=80) fmt,inputfile,outfile
        logical OK
        double precision cubicbox
        real (dl) rk



        call getarg(1,inputfile)
       open(1,file=inputfile,form='unformatted')
!        read(1) pamp
        read(1) InData%Params
        CP = Indata%Params
        read(1) InData%MTrans%num_q_trans
        call Transfer_Allocate(InData%MTrans)
        read(1) (InData%MTrans%q_trans(i),i=1,InData%MTrans%num_q_trans)
        read(1) (((InData%MTrans%TransferData(i,j,k),&
                i=1,Transfer_max), &
                j=1,InData%MTrans%num_q_trans), &
                k=1,CP%Transfer%num_redshifts)
        read(1) ((InData%MTrans%sigma_8(i,j), &
                i=1,CP%Transfer%num_redshifts), &
                j=1,CP%InitPower%nn)
        close(1)
        pamp = InData%MTrans%sigma_8(1,1)**2

!        print *,CP%Transfer%redshifts(1),Indata%Params%Transfer%redshifts(1)
!        boxsize = 128
!        print *,pamp, InData%MTrans%sigma_8(1,1)**2*boxsize**3

       
        omepb = CP%omegab 
        omep = CP%omegac  + omepb
        omeplam  = CP%omegav  
        hubble = CP%H0  /100.
        write(*,*)
        write(*,*)
        print *, '#####################################################'
        print *, '####      POWER SPECTRUM PARAMETER        ###########'
        print *, '#####################################################'
        print *, 'omep                  = ',omep
        print *, 'omepb                 = ',omepb
        print *, 'omeplam               = ',omeplam
        print *, 'npower                = ',CP%InitPower%an(1)
        print *, 'Hubble (km/sec.)      = ',CP%H0
!       if(CP%OutputNormalization) then
        if(CP%OutputNormalization .ne. 0) then
        print *, 'Normalization applied   '
        else 
        print *, 'Normalization unapplied '
        endif
        print *, 'simulation amplitude/boxsize**3  = ',pamp
        print *, 'power Amplitude       = ',CP%InitPower%ScalarPowerAmp(1)
        print *, 'number of redshift    = ',CP%InitPower%nn
        print *, 'Redshift              = ',CP%Transfer%redshifts(1)
        print *, 'Highprecision         = ',CP%Transfer%high_precision
        print *, 'kmax                  = ',CP%Transfer%kmax
        print *, 'k per logint          = ',CP%Transfer%k_per_logint
        print *, '#####################################################'
        print *, '####      POWER SPECTRUM PARAMETER        ###########'
        print *, '#####################################################'
        write(*,*)
        write(*,*)
        print *, '#    k           P_t(k)        P_cdm(k)     P_b(k)       T^2(k)  Amp of k^(n-1)'

        OK= CAMB_ValidateParams(CP)
!       if(OK .eq. .false.) then
        if(.not. OK) then
           print *,'###################################'
           print *,'###################################'
           print *,'Wrong set of parameters in the CAMB'
           print *,'###################################'
           print *,'###################################'
           print *,'Ok = ',OK
           stop
        endif

        P = CP
        call InitializePowers(CP%InitPower,CP%curv)

         nmatterpk = 8192
         call GetPower(minkh,dlnkh,matterpk,cmatterpk,bmatterpk,nmatterpk,Ts)
         do i = 1, nmatterpk
            rkarray(i) = minkh*exp((i-1)*dlnkh)
            rk = rkarray(i)
!            matterpk(i) = sqrt(matterpk(i))
            write(*,100) rkarray(i), matterpk(i),cmatterpk(i),bmatterpk(i),Ts(i),ScalarPower(rk,1)
         enddo
100      format(6(e15.8,1x))


        !call CAMB_GetResults(P) 
!        call CAMB_GetTransfers(P,OutData)
             
!       This cubic box multiplication is to normalize the amplitude in the unit of 2pi*u not k.
!        pamp = (MT%sigma_8(1,1))**2*cubicbox

!       open(1,file='Check.dat',form='unformatted')
!       read(1) P,opamp
!        read(1) InData%Params
!        read(1) InData%MTrans%num_q_trans
!        call Transfer_Allocate(InData%MTrans)
!        read(1) (InData%MTrans%q_trans(i),i=1,InData%MTrans%num_q_trans)
!        read(1) (((InData%MTrans%TransferData(i,j,k),i=1,Transfer_max),j=1,InData%MTrans%num_q_trans), &
!                k=1,CP%Transfer%num_redshifts)
!        read(1) ((InData%MTrans%sigma_8(i,j),i=1,CP%Transfer%num_redshifts),j=1,CP%InitPower%nn)
!        close(1)


       
      end program viewpower

      subroutine GetPower(minkh,dlnkh,matterpk,cmatterpk,bmatterpk,nmatterpk,Ts)
        use CAMB
        use MyModule
        use transfer
!        use InitialPower
        implicit none
        integer i
    
        type(MatterTransferData) MTrans
        integer itf,in
        integer error
        real minkh,dlnkh
        integer npoints,nmatterpk
        real, dimension(:,:), allocatable :: outpower
        real matterpk(nmatterpk),cmatterpk(nmatterpk),bmatterpk(nmatterpk),Ts(nmatterpk)
        character(LEN=80) fmt
        logical OK
        real(dl) logmink,k,h


        minkh = 5e-5
        dlnkh = 0.02
        Mtrans = InData%MTrans


        itf = 1

!        print *, Mtrans%num_q_trans,MTrans%TransferData(Transfer_kh,Mtrans%num_q_trans,itf)

         npoints = alog(MTrans%TransferData(Transfer_kh,Mtrans%num_q_trans,itf)/minkh)/dlnkh+1
         if(npoints .gt. nmatterpk) then
            print *,'#################################################'
            print *,'#################################################'
            print *,'Small number of power spectrum array ',nmatterpk
            print *,'The npoints in CAMB is ', npoints
            print *,'#################################################'
            print *,'#################################################'
            stop
         else
            nmatterpk = npoints
         endif
         in = 1
         call Transfer_GetMatterPower(Mtrans,matterpk(1),itf,in,minkh,dlnkh,npoints)
         call Transfer_GetSelectedMatterPower(Mtrans,bmatterpk(1),itf,in,minkh,dlnkh,npoints,Transfer_b)
         call Transfer_GetSelectedMatterPower(Mtrans,cmatterpk(1),itf,in,minkh,dlnkh,npoints,Transfer_cdm)
         h = CP%H0/100
         logmink = log(minkh)
         do i = 1, npoints
            k = exp(logmink+dlnkh*(i-1))*h
            Ts(i) = matterpk(i)/ScalarPower(k,in)/(k*pi*twopi*h**3)
         enddo

        end subroutine GetPower

        subroutine Transfer_GetSelectedMatterPower(MTrans,outpower, itf, &
                        in, minkh, dlnkh, npoints,Transfer_matter)
          use CAMB
          use mymodule
          !Allows for non-smooth priordial spectra
          !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
          !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
          !in units of (h Mpc^{-1})^3.   
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PK
        
          integer, intent(in) :: itf, in, npoints
          real, intent(out) :: outpower(npoints)
          real, intent(in) :: minkh, dlnkh
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
          integer ik, llo,il,lhi,lastix
          real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
          real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
          integer Transfer_matter
          

          if (npoints < 2) &
                stop 'Need at least 2 points in Transfer_GetMatterPower'

!         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
!            stop 'Transfer_GetMatterPower: kh out of computed region'
!          end if
          if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
                .and. FeedbackLevel > 0 )  write(*,*) &
             'Warning: extrapolating matter power in Transfer_GetMatterPower'

          
          if (CP%NonLinear/=NonLinear_None) then
           call Transfer_GetMatterPowerData(MTrans, PK, in, itf)
           call NonLinear_GetRatios(PK)
          end if
           
          h = CP%H0/100
          logmink = log(minkh)
          do ik=1,MTrans%num_q_trans
             kh = MTrans%TransferData(Transfer_kh,ik,itf)
             k = kh*h
             kvals(ik) = log(kh)
             atransfer=MTrans%TransferData(Transfer_matter,ik,itf)
             if (CP%NonLinear/=NonLinear_None) &
                 atransfer = atransfer* PK%nonlin_ratio(ik,1) !only one element, this itf
             matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
                 !Put in power spectrum later: transfer functions should be smooth, initial power may not be                
          end do
             
          call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)

            llo=1
            lastix = npoints + 1
            do il=1, npoints
               xi=logmink + dlnkh*(il-1)
               if (xi < kvals(1)) then
                 outpower(il)=-30.
                 cycle
               end if
               do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
                  llo=llo+1
               end do
               if (llo == MTrans%num_q_trans) then
                   lastix = il
                   exit
               end if
               lhi=llo+1
               ho=kvals(lhi)-kvals(llo) 
               a0=(kvals(lhi)-xi)/ho
               b0=(xi-kvals(llo))/ho
              
               outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
                       +(b0**3-b0)*ddmat(lhi))*ho**2/6
              
            end do

            do while (lastix <= npoints)
               !Do linear extrapolation in the log
               !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
               outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
               lastix = lastix+1
            end do

            outpower = exp(max(-30.,outpower))

            do il = 1, npoints
               k = exp(logmink + dlnkh*(il-1))*h
               outpower(il) = outpower(il) * ScalarPower(k,in) 
            end do

          if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_Free(PK)

        end subroutine Transfer_GetSelectedMatterPower


