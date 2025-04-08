 !Interface module for CAMB. Call CAMB_GetResults to do the work.

     module CAMB
         use Precision
         use ModelParams
         use ModelData
         use Transfer
         use GaugeInterface
         use InitialPower
         implicit none

         Type CAMBdata
            Type (ClTransferData) :: ClTransScal,ClTransTens,ClTransVec
            Type (MatterTransferData) :: MTrans
            Type (CAMBparams) :: Params
         end Type CAMBdata

!         public CAMB_GetTransfers, CAMB_GetResults, CAMB_GetCls, CAMB_SetDefParams, & 
!                CAMB_ValidateParams, CAMB_GetAge,CAMB_InitCAMBdata, 
     contains

       subroutine CAMB_GetTransfers(Params, OutData, error)
        use CAMBmain
        use lensing
        type(CAMBparams) :: Params
        type (CAMBdata)  :: OutData
        integer, optional :: error !Zero if OK
  
      !Set internal types from OutData so it always 'owns' the memory, prevent leaks
      
        MT =  OutData%MTrans

        CTransScal = OutData%ClTransScal  
        CTransVec  = OutData%ClTransVec  
        CTransTens = OutData%ClTransTens
  

        call CAMB_GetResults(Params, error)
        OutData%Params = Params
        OutData%MTrans = MT
        OutData%ClTransScal = CTransScal
        OutData%ClTransVec  = CTransVec
        OutData%ClTransTens = CTransTens
  
       end subroutine CAMB_GetTransfers

       subroutine CAMB_InitCAMBdata(Dat)
        type (CAMBdata) :: Dat
     
        call Ranges_Init(Dat%ClTransScal%q)
        call Ranges_Init(Dat%ClTransVec%q)
        call Ranges_Init(Dat%ClTransTens%q)

        nullify(Dat%ClTransScal%Delta_p_l_k)
        nullify(Dat%ClTransVec%Delta_p_l_k)
        nullify(Dat%ClTransTens%Delta_p_l_k)
        nullify(Dat%MTrans%sigma_8,Dat%MTrans%TransferData,Dat%MTrans%q_trans)
             
       end subroutine CAMB_InitCAMBdata


       subroutine CAMB_FreeCAMBdata(Dat)
            type (CAMBdata) :: Dat

            call Free_ClTransfer(Dat%ClTransScal)
            call Free_ClTransfer(Dat%ClTransVec)
            call Free_ClTransfer(Dat%ClTransTens)
            call Transfer_Free(Dat%MTrans)

       end subroutine CAMB_FreeCAMBdata


       subroutine CAMB_TransfersToPowers(CData)
        use CAMBmain
        use lensing
        type (CAMBdata) :: CData

        CP = CData%Params
        call InitializePowers(CP%InitPower,CP%curv)
        if (CData%Params%WantCls) then
          call ClTransferToCl(CData%ClTransScal,CData%ClTransTens, CData%ClTransvec)  
          if (CP%DoLensing) call lens_Cls
          if (CP%OutputNormalization == outCOBE) call COBEnormalize
        end if
        if (CData%Params%WantTransfer) call Transfer_Get_sigma8(Cdata%MTrans,8._dl)
     
       end subroutine CAMB_TransfersToPowers
  


       !Call this routine with a set of parameters to generate the results you want.
       subroutine CAMB_GetResults(Params, error)
        use CAMBmain
        use lensing
        type(CAMBparams) :: Params
        integer, optional :: error !Zero if OK
        type(CAMBparams) P
        logical :: separate = .true. !whether to do P_k in separate call or not

        if ((Params%DoLensing .or. num_redshiftwindows>0) &
           .and. Params%NonLinear==NonLinear_Lens) separate = .false.
        
         if (Params%WantCls .and. Params%WantScalars) then
          P = Params
          if (separate) then
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          end if
          P%WantTensors = .false.
          P%WantVectors = .false.
          if (present(error)) then
               call CAMBParams_Set(P, error) !set other derived variables in ModelParams (modules.f90) 
              if (error /= 0) return
              else
               call CAMBParams_Set(P)  
            end if
          call cmbmain
          call_again = .true.
          !Need to store CP%flat etc, but must keep original P_k settings
          CP%Transfer%high_precision = Params%Transfer%high_precision
          CP%WantTransfer = Params%WantTransfer
          CP%WantTensors = Params%WantTensors
          CP%WantVectors = Params%WantVectors
          CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
          Params = CP            
         end if
 
         if (Params%WantCls .and. Params%WantTensors) then
          P=Params
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          P%WantScalars = .false.
          P%WantVectors = .false.
          if (present(error)) then
               call CAMBParams_Set(P, error) !set other derived variables in ModelParams (modules.f90) 
              if (error /= 0) return
              else
               call CAMBParams_Set(P)  
            end if
           call cmbmain
           call_again = .true.
           CP%Transfer%high_precision = Params%Transfer%high_precision
           CP%WantTransfer = Params%WantTransfer
           CP%WantScalars = Params%WantScalars
           CP%WantVectors = Params%WantVectors
           CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
           Params = CP            
   
         end if
 
         if (Params%WantCls .and. Params%WantVectors) then
          P=Params
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          P%WantScalars = .false.
          P%WantTensors = .false.
          if (present(error)) then
               call CAMBParams_Set(P, error) !set other derived variables in ModelParams (modules.f90) 
              if (error /= 0) return
              else
               call CAMBParams_Set(P)  
            end if
           call cmbmain
           call_again = .true.
          CP%Transfer%high_precision = Params%Transfer%high_precision
          CP%WantTransfer = Params%WantTransfer
          CP%WantTensors = Params%WantTensors
          CP%WantScalars = Params%WantScalars
          CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
          Params = CP            
   
         end if

         if (Params%WantTransfer .and. &
          .not. (Params%WantCls .and. Params%WantScalars .and. .not. separate)) then
          P=Params
          P%WantCls = .false.
          P%WantScalars = .false.
          P%WantTensors = .false.
          P%WantVectors = .false.
          if (present(error)) then
               call CAMBParams_Set(P, error) !set other derived variables in ModelParams (modules.f90) 
              if (error /= 0) return
              else
               call CAMBParams_Set(P)  
            end if
           call cmbmain

          !Need to store num redshifts etc
          CP%WantScalars = Params%WantScalars
          CP%WantCls =  Params%WantCls
          CP%WantTensors = Params%WantTensors
          CP%WantVectors = Params%WantVectors
          Params = CP            

         end if

        call_again = .false.
   

        if (.not. CP%OnlyTransfers) then

         if (CP%WantCls .and. CP%OutputNormalization == outCOBE) call COBEnormalize

         if (CP%DoLensing) then
           call lens_Cls 
         end if

        end if

        end subroutine CAMB_GetResults


        !Return real (NOT double precision) arrays of the computed CMB  Cls
        !Output is l(l+1)C_l/2pi
        !If GC_Conventions = .false. use E-B conventions (as the rest of CAMB does)
        subroutine CAMB_GetCls(Cls, lmax, in, GC_conventions)
          integer, intent(IN) :: lmax, in
          logical, intent(IN) :: GC_conventions
          real, intent(OUT) :: Cls(2:lmax,1:4)
          integer l
      
          Cls = 0
          do l=2, lmax
             if (CP%WantScalars .and. l<= CP%Max_l) then
             Cls(l,1:2) = Cl_scalar(l, in,  C_Temp:C_E)
             Cls(l,4) = Cl_scalar(l, in,  C_Cross)
             end if
             if (CP%WantTensors .and. l <= CP%Max_l_tensor) then
                Cls(l,1:4) = Cls(l,1:4) + Cl_tensor(l, in,  CT_Temp:CT_Cross)
             end if
          end do
          if (GC_conventions) then
             Cls(:,2:3) = Cls(:,2:3)/2
             Cls(:,4)   = Cls(:,4)/sqrt(2.0)             
          end if
 
        end subroutine CAMB_GetCls

        function CAMB_GetAge(P)
           !Return age in gigayears, returns -1 on error
           type(CAMBparams), intent(in) :: P
           real(dl) CAMB_GetAge
           real(dl) atol,a1,a2,rombint, dtda
           real(dl), parameter :: Mpc = 3.085678e22_dl, &
                 c = 2.99792458e8_dl, Gyr=3.1556926e16
           integer error
           external rombint, dtda


           call  CAMBParams_Set(P, error)

           if (error/=0) then
            CAMB_GetAge = -1
           else

           atol = 1d-4
           a1=0
           a2=1
           CAMB_GetAge = rombint(dtda,a1,a2,atol)*Mpc/c/Gyr 
           end if
    
         end function CAMB_GetAge

      
        subroutine CAMB_SetDefParams(P)
            type(CAMBparams), intent(out) :: P

            P%WantTransfer= .false.
            P%WantCls = .true.
            P%Want_CMB = .true.

            P%omegab  = .045
            P%omegac  = 0.255
            P%omegav  = 0.7
            P%omegan  = 0
            P%H0      = 65

            P%TCMB    = 2.726
            P%YHe     = 0.24
            P%Num_Nu_massless =3.04
            P%Num_Nu_massive  =0
            P%Nu_mass_splittings = .false.
            P%Nu_mass_eigenstates = 0
           
            P%Scalar_initial_condition =initial_adiabatic
            P%NonLinear = NonLinear_none
            
            call SetDefPowerParams(P%InitPower)
          
            P%Reionization = .true.
            P%use_optical_depth = .false.
            P%Reion%redshift = 6
            P%Reion%fraction = 1

            P%Transfer%high_precision=.false.
    
            P%OutputNormalization = outNone

            P%WantScalars = .true.
            P%WantVectors = .false.
            P%WantTensors = .false.
            
            P%Max_l=1500
            P%Max_eta_k=3000
            P%Max_l_tensor=400
            P%Max_eta_k_tensor=800
            !Set up transfer just enough to get sigma_8 OK
            P%Transfer%kmax=0.9  
            P%Transfer%k_per_logint=0
            P%Transfer%num_redshifts=1
            P%Transfer%redshifts=0

            P%AccuratePolarization = .true.
            P%AccurateReionization = .false.
            P%AccurateBB = .false.

            P%DoLensing = .false.

            P%MassiveNuMethod = Nu_best
            P%OnlyTransfers = .false.

         end subroutine CAMB_SetDefParams


         !Stop with error is not good
         function CAMB_ValidateParams(P) result(OK)
            type(CAMBparams), intent(in) :: P
            logical OK

             OK = .true.
             if (.not. P%WantTransfer .and. .not. P%WantCls) then
                OK = .false.
                write(*,*) 'There is nothing to do! Do transfer functions or Cls.'
             end if

             if (P%h0 < 20._dl.or.P%h0 > 100._dl) then
               OK = .false.
               write(*,*) '  Warning: H0 has units of km/s/Mpc. You have:', P%h0
            end if
             if (P%tcmb < 2.7d0.or.P%tcmb > 2.8d0) then
                write(*,*) '  Warning: Tcmb has units of K.  Your have:', P%tcmb
             end if

             if (P%yhe < 0.2d0.or.P%yhe > 0.8d0) then
                OK = .false.
                write(*,*) &
                     '  Warning: YHe is the Helium fraction of baryons.', &
                     '  Your have:', P%yhe
             end if
             if (P%Num_Nu_massive < 0.or.P%Num_Nu_massive > 3.1) then
                OK = .false.
                write(*,*) &
                     'Warning: Num_Nu_massive is strange:',P%Num_Nu_massive 
              end if
             if (P%Num_Nu_massless < 0.or.P%Num_Nu_massless > 3.1) then
                OK = .false.
                write(*,*) &
                     'Warning: Num_nu_massless is strange:', P%Num_Nu_massless
              end if
             if (P%Num_Nu_massive < 1 .and. P%omegan > 0.0) then
                OK = .false.
                write(*,*) &
                     'Warning: You have omega_neutrino > 0, but no massive species'
              end if


             if (P%omegab<0.001 .or. P%omegac<0 .or. P%omegab>1 .or. P%omegac>3) then
                OK = .false.
                write(*,*) 'Your matter densities are strange'
             end if

             if (P%WantScalars .and. P%Max_eta_k < P%Max_l .or.  &
                  P%WantTensors .and. P%Max_eta_k_tensor < P%Max_l_tensor) then
                OK = .false.
                write(*,*) 'You need Max_eta_k larger than Max_l to get good results'
             end if
             if (P%Reionization) then
               if (P%use_optical_depth) then
                  if (P%Reion%optical_depth<0 .or. P%Reion%optical_depth > 0.9) then
                    OK = .false.
                    write(*,*) 'Optical depth is strange. You have:',P%Reion%optical_depth 
                  end if
               else
                  if (P%Reion%redshift < 0 .or. P%Reion%Redshift > 30) then
                     OK = .false.
                      write(*,*) 'Reionization redshift strange. You have: ',P%Reion%Redshift
                  end if
                  if (P%Reion%fraction < 0.1 .or. P%Reion%fraction > 1 + P%Yhe) then
                     OK = .false.
                      write(*,*) 'Reionization fraction strange. You have: ',P%Reion%fraction
                  end if
               end if
             end if

             if (P%WantTransfer) then
              if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
                OK = .false.
                write(*,*) 'Maximum ',  max_transfer_redshifts, &
                     'redshifts. You have: ', P%transfer%num_redshifts 
              end if
              if (P%transfer%kmax < 0.01 .or. P%transfer%kmax > 50000 .or. &
                     P%transfer%k_per_logint>0 .and.  P%transfer%k_per_logint <1) then
                 OK = .false.
                 write(*,*) 'Strange transfer function settings.'
              end if
              if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
                OK = .false.
                write(*,*) 'Maximum ',  max_transfer_redshifts, &
                     'redshifts. You have: ', P%transfer%num_redshifts 
              end if


             end if

         end function CAMB_ValidateParams

  end module CAMB


  function dtda(a)
          use Precision
          implicit none
          real(dl) dtda,dtauda,a
          external dtauda
          
          dtda= dtauda(a)*a
  end function

        

