c    Serve heating/cooling arrays for a given t, n, z, r
c    t ; temperature (K)
c       n ; hydrogen number density (cm-3)
c       z ; metallicity scale to solar value (Z_sun)
c    r ; redshift
c
c
c    This is the main function.
c    You need to include some of these common variables.
        function coolheat(gas_temp,gas_dens,gas_metl,UV_shield)
        implicit double precision (a-h,o-z)
        double precision coolheat
        parameter (n=100)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)

c       -----
c    You have to call this in the main where redshift is not determined.
c    Common/cool_1~cool_4 contain whole data about cooling/heating.
       
c    -----
c    You will extract cooling/heating data of redshift what you want.
c       Common/cool_5 contains cooling/heating data of the redshift.
c       redshift = 5.0
c       call get_cool_heat_redshift(redshift)

c       -----
c    Extract the final cooling/heating rate of T,N,Z.
c    T: temperature(K), N: hydrogen density(cm-3), Z: metallicity (Z_sun)
c    UV_shielding indicates self-shieling of cloud from UV. (1:yes,0:no)
c       gas_temp =1e4
c       gas_dens =10.0
c       gas_metl =1.0
c       UV_shielding=1
        call getcoolheatfinal(gas_temp,gas_dens,gas_metl,
     +                       UV_shielding,cooling_rate,heating_rate)

        coolheat = cooling_rate + heating_rate
c       ----
c       write(*,*) 'Cooling [erg s-1 cm-3]=', cooling_rate
c       write(*,*) 'Heating [erg s-1 cm-3]=', heating_rate
        return
        end
c---------------------------------------------------------------------
c       subroutine getcoolheatfinal(temp,dens,metal,uv_shielding,
c    &           coolrate,heatrate)
c       implicit real*8 (a-h,o-z)
c       integer uv_shielding
c       call get_cool_heat_final(temp,dens,metal,uv_shielding,coolrate,
c    &        heatrate)
c       return
c       end
c---------------------------------------------------------------------
        subroutine getcoolheatfinal(gas_temp,gas_dens,gas_metl,
     +                     UV_shielding,cooling_rate,heating_rate)
        implicit double precision (a-h,o-z)
        parameter (n=100)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)

        if(gas_temp .le. 0.d0) then
           gas_temp_log = min(g_t_log(1),g_t2_log(1)) - 1
        else
        gas_temp_log=log10(gas_temp)
        endif
        if(gas_metl .le. 0.d0) then
           gas_metl_log = min(g_z_log(1),g_z2_log(1)) - 1
        else
           gas_metl_log=log10(gas_metl)
        endif
        gas_dens_log=log10(gas_dens)

        if(UV_shielding.eq.1) then
c       When gas cloud is seld-shielding from the external UV source.
          if(gas_temp_log.lt.g_t2_log(1)) then 
            n_x_2=1 
            n_x_1=2
          else if(gas_temp_log.gt.g_t2_log(num_t2)) then
            n_x_1=num_t2-1
            n_x_2=num_t2
          else
            diff_t2_log=(g_t2_log(num_t2)-g_t2_log(1))/(num_t2-1)
            n_x=((gas_temp_log-g_t2_log(1))/diff_t2_log)+1
            n_x_1=n_x
            n_x_2=n_x+1
          endif
          if(gas_metl_log.lt.g_z2_log(1)) then
            n_y_2=1
            n_y_1=2
          else if(gas_metl_log.gt.g_z2_log(num_z2)) then
            n_y_1=num_z2-1
            n_y_2=num_z2
          else
            diff_z2_log=(g_z2_log(num_z2)-g_z2_log(1))/(num_z2-1)
            n_y=((gas_metl_log-g_z2_log(1))/diff_z2_log)+1
            n_y_1=n_y
            n_y_2=n_y+1
          endif
          x1  =g_t2_log(n_x_1)
          x2  =g_t2_log(n_x_2)
          xx  =gas_temp_log      
          y1  =g_z2_log(n_y_1)
          y2  =g_z2_log(n_y_2)
          yy  =gas_metl_log      
          qc11=cool_NO_log(n_x_1,n_y_1)
          qc12=cool_NO_log(n_x_1,n_y_2)
          qc21=cool_NO_log(n_x_2,n_y_1)
          qc22=cool_NO_log(n_x_2,n_y_2)
          result_c_log=qc11*(x2-xx)*(y2-yy)/((x2-x1)*(y2-y1))+
     +                 qc21*(xx-x1)*(y2-yy)/((x2-x1)*(y2-y1))+
     +                 qc11*(x2-xx)*(yy-y1)/((x2-x1)*(y2-y1))+
     +                 qc22*(xx-x1)*(yy-y1)/((x2-x1)*(y2-y1))
          cooling_rate=10.0**result_c_log*gas_dens**2
          heating_rate=0.0
          return
        else 
c    When gas cloud is exposed to the external UV source & CMB.
          if(gas_temp_log.lt.g_t_log(1)) then
            n_x_2=1
            n_x_1=2
          else if(gas_temp_log.gt.g_t_log(num_t)) then
            n_x_1=num_t-1
            n_x_2=num_t
          else
            diff_t_log=(g_t_log(num_t)-g_t_log(1))/(num_t-1)
            n_x=((gas_temp_log-g_t_log(1))/diff_t_log)+1
            n_x_1=n_x
            n_x_2=n_x+1
          endif
          if(gas_dens_log.lt.g_n_log(1)) then
            n_y_2=1
            n_y_1=2
          else if(gas_dens_log.gt.g_n_log(num_n)) then
            n_y_1=num_n-1
            n_y_2=num_n
          else
            diff_n_log=(g_n_log(num_n)-g_n_log(1))/(num_n-1)
            n_y=((gas_dens_log-g_n_log(1))/diff_n_log)+1
            n_y_1=n_y
            n_y_2=n_y+1
          endif
          if(gas_metl_log.lt.g_z_log(1)) then
            n_z_2=1
            n_z_1=2
          else if(gas_metl_log.gt.g_z_log(num_z)) then
            n_z_1=num_z-1
            n_z_2=num_z
          else
            diff_z_log=(g_z_log(num_z)-g_z_log(1))/(num_z-1)
            n_z=((gas_metl_log-g_z_log(1))/diff_z_log)+1
            n_z_1=n_z
            n_z_2=n_z+1
          endif
          xd =gas_temp_log-g_t_log(n_x_1)
          yd =gas_dens_log-g_n_log(n_y_1)
          zd =gas_metl_log-g_z_log(n_z_1)
          ct1=(c_UV_tnz_log(n_x_1,n_y_1,n_z_1))*(1.0-zd)+
     +        (c_UV_tnz_log(n_x_1,n_y_1,n_z_2))*zd
          ct2=(c_UV_tnz_log(n_x_1,n_y_2,n_z_1))*(1.0-zd)+
     +        (c_UV_tnz_log(n_x_1,n_y_2,n_z_2))*zd
          cs1=(c_UV_tnz_log(n_x_2,n_y_1,n_z_1))*(1.0-zd)+
     +        (c_UV_tnz_log(n_x_2,n_y_1,n_z_2))*zd
          cs2=(c_UV_tnz_log(n_x_2,n_y_2,n_z_1))*(1.0-zd)+
     +        (c_UV_tnz_log(n_x_2,n_y_2,n_z_2))*zd
          ht1=(h_UV_tnz_log(n_x_1,n_y_1,n_z_1))*(1.0-zd)+
     +        (h_UV_tnz_log(n_x_1,n_y_1,n_z_2))*zd
          ht2=(h_UV_tnz_log(n_x_1,n_y_2,n_z_1))*(1.0-zd)+
     +        (h_UV_tnz_log(n_x_1,n_y_2,n_z_2))*zd
          hs1=(h_UV_tnz_log(n_x_2,n_y_1,n_z_1))*(1.0-zd)+
     +        (h_UV_tnz_log(n_x_2,n_y_1,n_z_2))*zd
          hs2=(h_UV_tnz_log(n_x_2,n_y_2,n_z_1))*(1.0-zd)+
     +        (h_UV_tnz_log(n_x_2,n_y_2,n_z_2))*zd
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
c       subroutine getcoolheat(redshift)
c       implicit real*8 (a-h,o-z)
c       call get_cool_heat_redshift(redshift)
c       return
c       end
c---------------------------------------------------------------------
ccc    Serve heating/cooling arrays for a given redshift
ccc    Interpolation in redshift space (logscale)
        subroutine getcoolheatredshift(redshift)
        implicit double precision (a-h,o-z)
        parameter (n=100)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)

        if(redshift.lt.0) then
          write(*,*) 'Error on get_cool_heat_redshift.'
          write(*,*) 'Redshift can not be less than 0.'
          stop
        endif
        redshift_log=log10(redshift+1.0)

        if(redshift_log.eq.g_r_log(1).or.
     &     redshift_log.eq.g_r_log(num_r)) then
c1       redshift = min(g_r) or max(g_r)
          if(redshift_log.eq.g_r_log(1)) then
            num_redshift=1
          else
            num_redshift=num_r
          endif
          do n_z=1,num_z
            do n_n=1,num_n
              do n_t=1,num_t
                c_UV_tnz_log(n_t,n_n,n_z)=
     &           cool_UV_log(n_t,n_n,n_z,num_redshift)
                h_UV_tnz_log(n_t,n_n,n_z)=
     &           heat_UV_log(n_t,n_n,n_z,num_redshift)
              enddo
            enddo
          enddo
          return
        else if(redshift_log.gt.g_r_log(num_r)) then
c2    redshift > max(g_r)
        x1=g_r_log(num_r-1)
        x2=g_r_log(num_r)
        xx=redshift_log
        do n_z=1,num_z
          do n_n=1,num_n
            do n_t=1,num_t
              cq1=cool_UV_log(n_t,n_n,n_z,num_r-1)
              cq2=cool_UV_log(n_t,n_n,n_z,num_r  )
              hq1=heat_UV_log(n_t,n_n,n_z,num_r-1)
              hq2=heat_UV_log(n_t,n_n,n_z,num_r  )
              c_UV_tnz_log(n_t,n_n,n_z)=(cq1+(cq2-cq1)*(xx-x1)/(x2-x1))
              h_UV_tnz_log(n_t,n_n,n_z)=(hq1+(hq2-hq1)*(xx-x1)/(x2-x1))
            enddo
          enddo
        enddo
        return
        else
c3    min(g_r) < redshift < max(g_r)
        diff_r_log=(g_r_log(num_r)-g_r_log(1))/(num_r-1)
        index_r1  =((redshift_log-g_r_log(1))/diff_r_log)+1
        index_r2  =index_r1+1
        x1        =g_r_log(index_r1)
        x2        =g_r_log(index_r2)
        xx        =redshift_log 
        do n_z=1,num_z
          do n_n=1,num_n
            do n_t=1,num_t
              cq1=cool_UV_log(n_t,n_n,n_z,index_r1)
              cq2=cool_UV_log(n_t,n_n,n_z,index_r2)
              hq1=heat_UV_log(n_t,n_n,n_z,index_r1)
              hq2=heat_UV_log(n_t,n_n,n_z,index_r2)
              c_UV_tnz_log(n_t,n_n,n_z)=(cq1+(cq2-cq1)*(xx-x1)/(x2-xx))
              h_UV_tnz_log(n_t,n_n,n_z)=(hq1+(hq2-hq1)*(xx-x1)/(x2-xx))
            enddo
          enddo
        enddo
        return
        endif
        end
c       subroutine readcoolingtable
c       call read_cooling_data
c       return
c       end
c---------------------------------------------------------------------
ccc    Read heating/cooing data.
        subroutine readcoolingdata(myid)
        implicit double precision (a-h,o-z)
        include "mpif.h"
        parameter(n=100)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)
        if(myid .eq. 0) then
        open(1,file='cooling_heating_UV.dat')
        read(1,*) num_t,num_n,num_z,num_r
        read(1,*) (g_t_log(n_t),n_t=1,num_t)
        read(1,*) (g_n_log(n_n),n_n=1,num_n)
        read(1,*) (g_z_log(n_z),n_z=1,num_z)
        read(1,*) (g_r_log(n_r),n_r=1,num_r)
        read(1,*) ((((cool_UV_log(n_t,n_n,n_z,n_r),n_t=1,num_t),
     &               n_n=1,num_n),n_z=1,num_z),n_r=1,num_r)
        read(1,*) ((((heat_UV_log(n_t,n_n,n_z,n_r),n_t=1,num_t),
     &               n_n=1,num_n),n_z=1,num_z),n_r=1,num_r)
        close(1)
        open(1,file='cooling_NO.dat')
        read(1,*) num_t2,num_z2
        read(1,*) (g_t2_log(n_t),n_t=1,num_t2)
        read(1,*) (g_z2_log(n_z),n_z=1,num_z2)
        read(1,*) ((cool_NO_log(n_t,n_z),n_t=1,num_t2),n_z=1,num_z2)
        close(1)
        endif
        call MPI_Bcast(num_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(num_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(num_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(num_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(num_t2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(num_z2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_t_log,num_t,MPI_DOUBLE_PRECISION,0,
     &                       MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_n_log,num_n,MPI_DOUBLE_PRECISION,0,
     &                       MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_z_log,num_z,MPI_DOUBLE_PRECISION,0,
     &                       MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_r_log,num_r,MPI_DOUBLE_PRECISION,0,
     &                       MPI_COMM_WORLD,ierror)
        call MPI_Bcast(cool_UV_log,n*n*n*n,
     &          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(heat_UV_log,n*n*n*n,
     &          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_t2_log,num_t2,MPI_DOUBLE_PRECISION,0,
     &         MPI_COMM_WORLD,ierror)
        call MPI_Bcast(g_z2_log,num_z2,MPI_DOUBLE_PRECISION,0,
     &         MPI_COMM_WORLD,ierror)
        call MPI_Bcast(cool_NO_log,n*n,MPI_DOUBLE_PRECISION,0,
     &         MPI_COMM_WORLD,ierror)

        return
        end
c-------------------------------------------------------------------
