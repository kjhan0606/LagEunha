      subroutine kjhan_Transfer_SaveMatterPower(MTrans, State,FileNames, all21cm)
      use constants
      !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
      Type(MatterTransferData), intent(in) :: MTrans
      Type(CAMBdata) :: State
      character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
      character(LEN=name_tag_len) :: columns(3)
      integer itf, i, unit
      integer points
      real, dimension(:,:), allocatable :: outpower
      real minkh,dlnkh
      Type(MatterPowerData) :: PK_data
      integer ncol
      logical, intent(in), optional :: all21cm
      logical all21
      !JD 08/13 Changes in here to PK arrays and variables
      integer itf_PK

      all21 = DefaultFalse(all21cm)
      if (all21) then
        ncol = 3
      else
        ncol = 1
      end if

      do itf=1, State%CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = State%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
                allocate(outpower(points,ncol))

                !Sources
                if (all21) then
                    call Transfer_Get21cmPowerData(MTrans, State, PK_data, itf_PK)
                else
                    call Transfer_GetMatterPowerData(State, MTrans, PK_data, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) then
                        call State%CP%NonLinearModel%GetNonLinRatios(State, PK_data)
                        PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
                        call MatterPowerdata_getsplines(PK_data)
                    end if
                end if

                outpower(:,1) = exp(PK_data%matpower(:,1))
                !Sources
                if (all21) then
                    outpower(:,3) = exp(PK_data%vvpower(:,1))
                    outpower(:,2) = exp(PK_data%vdpower(:,1))

                    outpower(:,1) = outpower(:,1)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,2) = outpower(:,2)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,3) = outpower(:,3)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                end if

                call MatterPowerdata_Free(PK_Data)
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData(Transfer_kh,i,1),outpower(i,:)
                end do
                close(unit)
            else
                if (all21) stop 'Transfer_SaveMatterPower: if output all assume not interpolated'
                minkh = 1e-4
                dlnkh = 0.02
                points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
                !             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
                allocate(outpower(points,1))
                call Transfer_GetMatterPowerS(State, MTrans, outpower(1,1), itf, minkh,dlnkh, points)

                columns(1) = 'P'
                unit = 1
                open(unit,file=FileNames(itf), form='formatted') 
                pamp = (MTrans%sigma_8(State%Transfer%PK_num_redshifts))**2
                write(unit, 'sigma8@thisRedshift= ', pamp)
                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1)
                end do
                close(unit)
            end if

            deallocate(outpower)
        end if
    end do

    end subroutine kjhan_Transfer_SaveMatterPower

