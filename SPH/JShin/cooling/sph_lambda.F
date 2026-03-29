c    Serve heating/cooling arrays for a given t, n, z, r
c    t ; temperature (K)
c       n ; hydrogen number density (cm-3)
c       z ; metallicity scale to solar value (Z_sun)
c    r ; redshift
c
c
c    This is the main function.
c    You need to include some of these common variables.

        subroutine  coolnheat(redshift,gas_temp,gas_dens,gas_metl,
     &             UV_shielding, dTime,Ntot,Temp)
        implicit real*8 (a-h,o-z)
        parameter (nnn=100)
        real*8 kB
        parameter(kB=1.380650424D-16)
        integer UV_shielding
        common /cool_1/cool_UV(nnn,nnn,nnn,nnn),heat_UV(nnn,nnn,nnn,nnn)
        common /cool_2/grid_t(nnn),grid_n(nnn),grid_z(nnn),grid_r(nnn)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO(nnn,nnn),grid_t2(nnn),grid_z2(nnn)
        common /cool_5/c_UV_tnz(nnn,nnn,nnn),h_UV_tnz(nnn,nnn,nnn)
        integer iflag
        save iflag/0/
        real iredshift
        save iredshift /-100/


c       -----
c    You have to call this in the main where redshift is not determined.
c    Common/cool_1~cool_4 contain whole data about cooling/heating.
        if(iflag .eq. 0) then
        call read_cooling_data
        iflag = 1
        endif

c    -----
c    You will extract cooling/heating data of redshift what you want.
c       Common/cool_5 contains cooling/heating data of the redshift.
        if(redshift .ne. iredshift)
        call get_cool_heat_redshift(redshift)
        iredshift = redshift
        endif

c       -----
c    Extract the final cooling/heating rate of T,N,Z.
c    T: temperature(K), N: hydrogen density(cm-3), Z: metallicity (Z_sun)
c    UV_shielding indicates self-shieling of cloud from UV. (1:yes,0:no)
c       gas_temp =1e6
c       gas_dens =10.0
c       gas_metl =2.0
c       UV_shielding=1
        call get_cool_heat_final(gas_temp,gas_dens,gas_metl,
     +                       UV_shielding,cooling_rate,heating_rate)

        dTemp = (heating_rate-cooling_rate)/(1.5*Ntot*kB)
        dTemp = dTemp*dTime
cccccccccccccccccccccccccccccccccccccccccccc
        Temp = Temp + dTemp*dTime
        if(Temp .lt.0.) Temp  = 0.


c       ----
c       write(*,*) 'Cooling [erg s-1 cm-3]=', cooling_rate
c       write(*,*) 'Heating [erg s-1 cm-3]=', heating_rate
c       stop
        return
        end
c---------------------------------------------------------------------
        subroutine get_cool_heat_final(gas_temp,gas_dens,gas_metl,
     +                     UV_shielding,cooling_rate,heating_rate)
        implicit real*8 (a-h,o-z)
        parameter (nnn=100)
        integer UV_shielding
        common /cool_2/grid_t(nnn),grid_n(nnn),grid_z(nnn),grid_r(nnn)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO(nnn,nnn),grid_t2(nnn),grid_z2(nnn)
        common /cool_5/c_UV_tnz(nnn,nnn,nnn),h_UV_tnz(nnn,nnn,nnn)

        if(UV_shielding.eq.1) then
c       When gas cloud is seld-shielding from the external UV source.
          do i=1,num_t2
            if(gas_temp.lt.grid_t2(i))then
              goto 123
            endif
          enddo
123       n_x=i
          do i=1,num_z2
            if(gas_metl.lt.grid_z2(i))then
              goto 124
            endif
          enddo
124       n_y=i
          if(n_x.eq.1.or.nx.eq.num_t2+1) then
            if(n_x.eq.1) then
              n_x_2=1
              n_x_1=2
            else
              n_x_1=num_t2-1
              n_x_2=num_t2
            endif
          else
            n_x_1=n_x-1
            n_x_2=n_x
          endif
          if(n_y.eq.1.or.n_y.eq.num_z2+1) then
            if(n_y.eq.1) then
              n_y_2=1
              n_y_1=2
            else  
              n_y_1=num_z2-1
              n_y_2=num_z2
            endif
          else
            n_y_1=n_y-1
            n_y_2=n_y
          endif
          x1  =log10(grid_t2(n_x_1))
          x2  =log10(grid_t2(n_x_2))
          xx  =log10(gas_temp      )
          y1  =log10(grid_z2(n_y_1))
          y2  =log10(grid_z2(n_y_2))
          yy  =log10(gas_metl      )
          qc11=log10(cool_NO(n_x_1,n_y_1))
          qc12=log10(cool_NO(n_x_1,n_y_2))
          qc21=log10(cool_NO(n_x_2,n_y_1))
          qc22=log10(cool_NO(n_x_2,n_y_2))
          result_c_log=qc11*(x2-xx)*(y2-yy)/((x2-x1)*(y2-y1))+
     +                 qc21*(xx-x1)*(y2-yy)/((x2-x1)*(y2-y1))+
     +                 qc11*(x2-xx)*(yy-y1)/((x2-x1)*(y2-y1))+
     +                 qc22*(xx-x1)*(yy-y1)/((x2-x1)*(y2-y1))
          cooling_rate=10.0**result_c_log*gas_dens**2
          heating_rate=0.0
          return
        else 
c    When gas cloud is exposed to the external UV source & CMB.
          do i=1,num_t
            if(gas_temp.lt.grid_t(i))then
              goto 125
            endif
          enddo
125       n_x=i
          do i=1,num_n
            if(gas_dens.lt.grid_n(i))then
              goto 126
            endif
          enddo
126       n_y=i
          do i=1,num_z
            if(gas_metl.lt.grid_z(i))then
              goto 127
            endif
          enddo
127       n_z=i
          if(n_x.eq.1.or.num_t+1) then
            if(n_x.eq.1) then 
              n_x_2=1
              n_x_1=2
            else 
              n_x_1=num_t-1
              n_x_2=num_t
            endif
          else
            n_x_1=n_x-1
            n_x_2=n_x
          endif
          if(n_y.eq.1.or.num_n+1) then
            if(n_y.eq.1) then
              n_y_2=1
              n_y_1=2
            else
              n_y_1=num_n-1
              n_y_2=num_n
            endif
          else
            n_y_1=n_y-1
            n_y_2=n_y
          endif
          if(n_z.eq.1.or.num_z+1) then
            if(n_z.eq.1) then
              n_z_2=1
              n_z_1=2
            else
              n_z_1=num_z-1
              n_z_2=num_z
            endif
          else
            n_z_1=n_z-1
            n_z_2=n_z
          endif
          xd =log10(gas_temp)-log10(grid_t(n_x_1))
          yd =log10(gas_dens)-log10(grid_n(n_y_1))
          zd =log10(gas_metl)-log10(grid_z(n_z_1))
          ct1=log10(c_UV_tnz(n_x_1,n_y_1,n_z_1))*(1.0-zd)+
     +        log10(c_UV_tnz(n_x_1,n_y_1,n_z_2))*zd
          ct2=log10(c_UV_tnz(n_x_1,n_y_2,n_z_1))*(1.0-zd)+
     +        log10(c_UV_tnz(n_x_1,n_y_2,n_z_2))*zd
          cs1=log10(c_UV_tnz(n_x_2,n_y_1,n_z_1))*(1.0-zd)+
     +        log10(c_UV_tnz(n_x_2,n_y_1,n_z_2))*zd
          cs2=log10(c_UV_tnz(n_x_2,n_y_2,n_z_1))*(1.0-zd)+
     +        log10(c_UV_tnz(n_x_2,n_y_2,n_z_2))*zd
          ht1=log10(h_UV_tnz(n_x_1,n_y_1,n_z_1))*(1.0-zd)+
     +        log10(h_UV_tnz(n_x_1,n_y_1,n_z_2))*zd
          ht2=log10(h_UV_tnz(n_x_1,n_y_2,n_z_1))*(1.0-zd)+
     +        log10(h_UV_tnz(n_x_1,n_y_2,n_z_2))*zd
          hs1=log10(h_UV_tnz(n_x_2,n_y_1,n_z_1))*(1.0-zd)+
     +        log10(h_UV_tnz(n_x_2,n_y_1,n_z_2))*zd
          hs2=log10(h_UV_tnz(n_x_2,n_y_2,n_z_1))*(1.0-zd)+
     +        log10(h_UV_tnz(n_x_2,n_y_2,n_z_2))*zd
          cw1=ct1*(1.0-yd)+ct2*yd
          cw2=cs1*(1.0-yd)+cs2*yd
          hw1=ht1*(1.0-yd)+ht2*yd
          hw2=hs1*(1.0-yd)+hs2*yd
          result_c_log=cw1*(1.0-xd)+cw2*xd
          result_h_log=hw1*(1.0-xd)+hw2*xd
          cooling_rate=10.0**result_c_log
          heating_rate=10.0**result_h_log
          return
        endif
        end
c---------------------------------------------------------------------
ccc    Serve heating/cooling arrays for a given redshift
ccc    Interpolation in redshift space (logscale)
        subroutine get_cool_heat_redshift(redshift)
        implicit real*8 (a-h,o-z)
        parameter (nnn=100)
        common /cool_1/cool_UV(nnn,nnn,nnn,nnn),heat_UV(nnn,nnn,nnn,nnn)
        common /cool_2/grid_t(nnn),grid_n(nnn),grid_z(nnn),grid_r(nnn)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_5/c_UV_tnz(nnn,nnn,nnn),h_UV_tnz(nnn,nnn,nnn)

        if(redshift.eq.grid_r(1).or.redshift.eq.grid_r(n_r)) then
c1       redshift = min(grid_r) or max(grid_r)
          if(redshift.eq.grid_r(1)) then
            num_redshift=1
          else
            num_redshift=num_r
          endif
          do n_z=1,num_z
            do n_n=1,num_n
              do n_t=1,num_t
                c_UV_tnz(n_t,n_n,n_z)=cool_UV(n_t,n_n,n_z,num_redshift)
                h_UV_tnz(n_t,n_n,n_z)=heat_UV(n_t,n_n,n_z,num_redshift)
              enddo
            enddo
          enddo
          return
        else if(redshift.gt.grid_r(num_r)) then
c2    redshift > max(grid_r)
          x1=log10(grid_r(num_r-1)+1.0)         
          x2=log10(grid_r(num_r  )+1.0)         
          xx=log10(redshift+1.0)
          do n_z=1,num_z
            do n_n=1,num_n
              do n_t=1,num_t
                cq1=log10(cool_UV(n_t,n_n,n_z,num_r-1))
                cq2=log10(cool_UV(n_t,n_n,n_z,num_r  ))
                hq1=log10(heat_UV(n_t,n_n,n_z,num_r-1))
                hq2=log10(heat_UV(n_t,n_n,n_z,num_r  ))
                c_UV_tnz(n_t,n_n,n_z)=10.0**(cq1+(cq2-cq1)*(xx-x1)/
     +                                  (x2-x1))
                h_UV_tnz(n_t,n_n,n_z)=10.0**(hq1+(hq2-hq1)*(xx-x1)/
     +                                  (x2-x1))
              enddo
            enddo
          enddo
          return
        else
c3    min(grid_r) < redshift < max(grid_r)
          do n_r=1,num_r
            if(grid_r(n_r).gt.redshift) then
              x1=log10(grid_r(n_r-1)+1.0)
              xx=log10(redshift     +1.0)
              x2=log10(grid_r(n_r)  +1.0)
              do n_z=1,num_z
                do n_n=1,num_n
                  do n_t=1,num_t
                    cq1=log10(cool_UV(n_t,n_n,n_z,n_r-1))
                    cq2=log10(cool_UV(n_t,n_n,n_z,n_r  ))
                    hq1=log10(heat_UV(n_t,n_n,n_z,n_r-1))
                    hq2=log10(heat_UV(n_t,n_n,n_z,n_r  ))
                    c_UV_tnz(n_t,n_n,n_z)=10.0**(cq1+(cq2-
     +                 cq1)*(xx-x1)/(x2-xx))
                    h_UV_tnz(n_t,n_n,n_z)=10.0**(hq1+(hq2-
     +                 hq1)*(xx-x1)/(x2-xx))
                  enddo
                enddo
              enddo
              return
            endif
          enddo
        endif
        end
c---------------------------------------------------------------------
ccc    Read heating/cooing data.
        subroutine read_cooling_data
        implicit real*8 (a-h,o-z)
        parameter(nnn=100)
        common /cool_1/cool_UV(nnn,nnn,nnn,nnn),heat_UV(nnn,nnn,nnn,nnn)
        common /cool_2/grid_t(nnn),grid_n(nnn),grid_z(nnn),grid_r(nnn)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO(nnn,nnn),grid_t2(nnn),grid_z2(nnn)
        open(1,file='cooling_heating_UV.dat')
        read(1,*) num_t,num_n,num_z,num_r
        read(1,*) (grid_t(n_t),n_t=1,num_t)
        read(1,*) (grid_n(n_n),n_n=1,num_n)
        read(1,*) (grid_z(n_z),n_z=1,num_z)
        read(1,*) (grid_r(n_r),n_r=1,num_r)
        read(1,*) ((((cool_UV(n_t,n_n,n_z,n_r),n_t=1,num_t),
     &               n_n=1,num_n),n_z=1,num_z),n_r=1,num_r)
        read(1,*) ((((heat_UV(n_t,n_n,n_z,n_r),n_t=1,num_t),
     &               n_n=1,num_n),n_z=1,num_z),n_r=1,num_r)
        close(1)
        open(1,file='cooling_NO.dat')
        read(1,*) num_t2,num_z2
        read(1,*) (grid_t2(n_t),n_t=1,num_t2)
        read(1,*) (grid_z2(n_z),n_z=1,num_z2)
        read(1,*) ((cool_NO(n_t,n_z),n_t=1,num_t2),n_z=1,num_z2)
        close(1)
        end
c-------------------------------------------------------------------
