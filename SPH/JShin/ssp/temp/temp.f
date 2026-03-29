c    
c    
c    
c    struct of star particle = SN_mark, sw_t_past
c    
c
c
c
        implicit real (a-h,l-z)
        real*8 delta_e

        call readstellarmodel(0)

        open(1,file='input')
        read(1,*) z, m, sw_t_now,sw_t_past
        close(1)

        call calcswfb(z,m,sw_t_now,sw_t_past,delta_m)

        open(1,file='output')
        write(1,*) sw_t_past,delta_m
        close(1)

        stop
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
        if(z .gt. metal(nz)) then
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
          write(*,*) 'first sw feedback'
        else if(sw_t_now-sw_t_past.ge.60)then 
          diff_t    = (time(nt)-time(1))/(nt-1)
          noft_past = ((sw_t_past-time(1))/diff_t)+1
          noft_pres = ((sw_t_now -time(1))/diff_t)+1
          past_m    = sw_ratio(nofz,noft_past)
          pres_m    = sw_ratio(nofz,noft_pres)
          delta_m   = abs(past_m-pres_m)*m
          sw_t_past = sw_t_now
          write(*,*) 'sw feedback'
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
        if(z .gt. metal(nz) ) then
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
        write(*,*) prob,prob_1,prob_2
        return
        end 
c---------------------------------------------------------------------
        subroutine readstellarmodel(myid)
        implicit real (a-h,o-z)
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
        return
        end
c--------------------------------------------------------------------
