c    
c    
c    
c    struct of star particle = SN_mark, sw_t_past
c    
c
c
c
        subroutine TEST
        implicit real (a-h,l-z)
        real*8 delta_e
        call readstellarmodel(0)

c       Metal is Z (fractional abundance of metal = x,y,z)
c       Mass unit is Msun, Time unit is Myr
        now_t    =15.0
        del_t    =1.0
        x_ssp    =0.72
        y_ssp    =0.277
        z_ssp    =0.003
        ssp_mass =1e3
        SN_mark = 0

c    calculate SN feedback, r_uniform=random generation
        if(now_t.gt.10.and.now_t.le.60) then 
          call calcsnprob(z_ssp,now_t,del_t,prob)
          r_uniform=0
          if(SN_mark.eq.0.and.r_uniform.le.prob) then
            call calcsnfb(z_ssp,ssp_mass,delta_e,delta_m,
     +                                         delta_x,delta_y,delta_z)
            write(*,*) 'Ejected energy [erg] =', delta_e
            write(*,*) 'Ejected Mass  [Msun] =', delta_m
            write(*,*) 'Ejected H     [Msun] =', delta_x
            write(*,*) 'Ejected He    [Msun] =', delta_y
            write(*,*) 'Ejected Metal [Msun] =', delta_z
            SN_mark  = 1
            sw_t_past=now_t
          endif
        else
c    calculate SW feedback 
          if(SN_mark.eq.1) then 
            sw_t_now = now_t
            call calcswfb(z_ssp,ssp_mass,sw_t_now,sw_t_past,delta_m)
            write(*,*) 'Ejected Mass  [Msun]=',delta_m
            write(*,*) 'Ejected H     [Msun]=',delta_m*x_ssp
            write(*,*) 'Ejected He    [Msun]=',delta_m*y_ssp
            write(*,*) 'Ejected Metal [Msun]=',delta_m*z_ssp
          endif
        endif
        end
c---------------------------------------------------------------------
        subroutine calcswfb(z,m,sw_t_now,sw_t_past,delta_m)
        implicit real (a-h,o-z)
        real m
        parameter(nz=16,nt=64)

        real metal(nz),fit_a(nz),fit_b(nz),fit_t1(nz),fit_t2(nz)
        real del_n(nz),del_m(nz),del_x(nz),del_y(nz),del_z(nz)
        real time(nt),sw_ratio(nz,nt)

        common /sn_r/metal,fit_a,fit_b,fit_t1,fit_t2
        common /sn_fb/del_n,del_m,del_x,del_y,del_z
        common /sw_fb/time,sw_ratio

        diff_log_z=(log10(metal(nz))-log10(metal(1)))/(nz-1)
        if(z .lt. metal(1) ) then
          z = metal(1)
        endif
        if(z .gt. metal(nz) ) then
          z = metal(nz)
        endif
        n_1=((log10(z)-log10(metal(1)))/diff_log_z)+1
        n_2=n_1+1
        if(abs(z-metal(n_1)).lt.abs(z-metal(n_2))) then
          nofz=n_1
        else
          nofz=n_2
        endif

        if(sw_t_past.le.60.and.sw_t_now.gt.60) then
          delta_m   = abs(1.0-del_m(nofz) -sw_ratio(nofz,1))*m
          sw_t_past = sw_t_now
        else if(sw_t_now-sw_t_past.ge.60)then 
          diff_t    = (time(nt)-time(1))/(nt-1)
          noft_past = ((sw_t_past-time(1))/diff_t)+1
          noft_pres = ((sw_t_now -time(1))/diff_t)+1
          past_m    = sw_ratio(nofz,noft_past)
          pres_m    = sw_ratio(nofz,noft_pres)
          delta_m   = abs(past_m-pres_m)*m
          sw_t_past = sw_t_now
        else
          delta_m = 0.0
        endif
        end
c---------------------------------------------------------------------
        subroutine calcSNfb(z,m,delta_e,delta_m,
     +                                          delta_x,delta_y,delta_z)
        implicit real (a-h,o-z)
        real m
        parameter (nz=16)
        real*8 delta_e
        real metal(nz),fit_a(nz),fit_b(nz),fit_t1(nz),fit_t2(nz)
        real del_n(nz),del_m(nz),del_x(nz),del_y(nz),del_z(nz)
        common /sn_r/metal,fit_a,fit_b,fit_t1,fit_t2
        common /sn_fb/del_n,del_m,del_x,del_y,del_z

        diff_log_z=(log10(metal(nz))-log10(metal(1)))/(nz-1)
        if(z .lt. metal(1) ) then
          z = metal(1)
        endif
        if(z .gt. metal(nz)) then
          z = metal(nz)
        endif
        n_1=((log10(z)-log10(metal(1)))/diff_log_z)+1
        n_2=n_1+1
        if(abs(z-metal(n_1)).lt.abs(z-metal(n_2))) then
          n=n_1
        else
          n=n_2
        endif
        delta_e=del_n(n)*m*1.d51
        delta_m=del_m(n)*m
        delta_x=del_x(n)*delta_m
        delta_y=del_y(n)*delta_m
        delta_z=del_z(n)*delta_m
        return
        end
c---------------------------------------------------------------------
        subroutine calcsnprob(z,t,del_t,prob)
        implicit real (a-h,o-z)
        parameter (nz=16)
        real metal(nz),fit_a(nz),fit_b(nz),fit_t1(nz),fit_t2(nz)
        common /sn_r/metal,fit_a,fit_b,fit_t1,fit_t2

        diff_log_z=(log10(metal(nz))-log10(metal(1)))/(nz-1)
        if(z .lt. metal(1)) then
          z = metal(1)
        endif
        if(z. gt. metal(nz)) then
          z = metal(nz)
        endif
        n_1=((log10(z)-log10(metal(1)))/diff_log_z)+1
        n_2=n_1+1
        if(abs(z-metal(n_1)).lt.abs(z-metal(n_2))) then
          n=n_1
        else
          n=n_2
        endif
        if(t.le.fit_t1(n).or.t+del_t.ge.fit_t2(n)) then
          if(t.lt.fit_t1(n)) then
            prob=0.0
          else
            prob=1.0
          endif
          prob_1=0.0
          prob_2=0.0
        else
          prob_1=10.0**(fit_a(n)*(t+del_t)  +fit_b(n))
     +          -10.0**(fit_a(n)*(t      )  +fit_b(n))
          prob_2=10.0**(fit_a(n)*(fit_t2(n))+fit_b(n))
     +          -10.0**(fit_a(n)*(t      )  +fit_b(n))
          prob=prob_1/prob_2
        endif
        return
        end 
c---------------------------------------------------------------------
        subroutine readstellarmodel(myid)
        implicit real (a-h,o-z)
        include 'mpif.h'
        parameter (nz=16,nt=64)
        real metal(nz),fit_a(nz),fit_b(nz),fit_t1(nz),fit_t2(nz)
        real del_n(nz),del_m(nz),del_x(nz),del_y(nz),del_z(nz)
        real time(nt),sw_ratio(nz,nt)
        common /sn_r/metal,fit_a,fit_b,fit_t1,fit_t2
        common /sn_fb/del_n,del_m,del_x,del_y,del_z
        common /sw_fb/time,sw_ratio

        if(myid .eq. 0) then
        open(1,file='SN_rate_kroupa.dat')
        do i=1,nz
          read(1,*) metal(i),fit_a(i),fit_b(i),fit_t1(i),fit_t2(i)
        enddo
        close(1)
        open(1,file='SN_feedback_kroupa.dat')
        do i=1,nz
          read(1,*) del_n(i),del_m(i),del_x(i),del_y(i),del_z(i)
        enddo
        close(1)
        open(1,file='SW_feedback_kroupa.dat')
        read(1,*) (time(i),i=1,nt)
        read(1,*) ((sw_ratio(i,j),j=1,nt),i=1,nz)
        close(1)
        endif
        call MPI_Bcast(metal,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(fit_a,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(fit_b,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(fit_t1,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(fit_t2,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(del_n,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(del_m,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(del_x,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(del_y,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(del_z,nz,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(time,nt,MPI_REAL,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(sw_ratio,nz*nt,MPI_REAL,0,MPI_COMM_WORLD,ierror)

        return
        end
c--------------------------------------------------------------------
