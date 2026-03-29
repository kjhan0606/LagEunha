c    Serve heating/cooling arrays for a given t, n, z, r
c    t ; temperature (K)
c       n ; hydrogen number density (cm-3)
c       z ; metallicity scale to solar value (Z_sun)
c    r ; redshift
c
c
c    This is the main function.
c    You need to include some of these common variables.
        function heatcool()
        implicit double precision (a-h,o-z)
        include "mpif.h"
        parameter (n=64,nmax=32,nn=16)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)
        common /evol_1_sh0/delta_e0(nn,nn,nn,nmax)
        common /evol_1_sh1/delta_e1(nn,nn,nn,nmax)
        common /evol_2/g_tt_log(nn),g_nn_log(nn),g_zz_log(nn)
        common /evol_3/diff_tt_log_inv,diff_nn_log_inv,diff_zz_log_inv
        common /evol_4/g_tt(nn),g_nn(nn),g_zz(nn)

        integer mpi_status,ierr,myid,nproc
        integer status(mpi_status_size)

        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,myid,ierr)
        call mpi_comm_size(mpi_comm_world,nproc,ierr)

        global_dt=2723229.11051711
        redshift =0.0
        n_shd    =1

        call readcoolingdata(myid)
        call getcoolheatredshift(redshift)
        call calcevoldata(global_dt,n_shd,myid,nproc)
        call mpi_finalize(ierr)
        heatcool=0
        return
        end
c---------------------------------------------------------------------
        subroutine getdeltae(gas_temp,gas_nden,gas_metl,n_step,n_shd,
     +                           final_de)
        implicit double precision (a-h,o-z)
        parameter (n=64,nmax=32,nn=16)
        common /evol_1_sh0/delta_e0(nn,nn,nn,nmax)
        common /evol_1_sh1/delta_e1(nn,nn,nn,nmax)
        common /evol_2/g_tt_log(nn),g_nn_log(nn),g_zz_log(nn)
        common /evol_3/diff_tt_log_inv,diff_nn_log_inv,diff_zz_log_inv
        common /evol_4/g_tt(nn),g_nn(nn),g_zz(nn)

        gas_temp_log=log10(gas_temp)
        gas_nden_log=log10(gas_nden)
        gas_metl_log=log10(gas_metl)

        if(gas_nden_log.lt.g_nn_log(1)) then
          gas_nden_log=g_nn_log(1)
        endif
        if(gas_nden_log.gt.g_nn_log(nn)) then
          gas_nden_log=g_nn_log(nn)
        endif
        
        n_x=(gas_temp_log-g_tt_log(1))*diff_tt_log_inv+1
        n_x_1=n_x
        n_x_2=n_x+1
        if(n_x_1.lt.1) then 
          n_x_1=2
          n_x_2=1
        else if(n_x_2.gt.nn) then
          n_x_1=nn-1
          n_x_2=nn
        endif
        n_y=(gas_nden_log-g_nn_log(1))*diff_nn_log_inv+1
        n_y_1=n_y
        n_y_2=n_y+1
        if(n_y_1.lt.1) then
          n_y_1=2
          n_y_2=1
        else if(n_y_2.gt.nn) then
          n_y_1=nn-1
          n_y_2=nn
        endif
        n_z=(gas_metl_log-g_zz_log(1))*diff_zz_log_inv+1
        n_z_1=n_z
        n_z_2=n_z+1
        if(n_z_1.lt.1) then
          n_z_1=2
          n_z_2=1
        else if(n_z_2.gt.nn) then
          n_z_1=nn-1
          n_z_2=nn
        endif

        if(n_shd.eq.0) then
          de_111=delta_e0(n_x_1,n_y_1,n_z_1,n_step)
          de_112=delta_e0(n_x_1,n_y_1,n_z_2,n_step)
          de_121=delta_e0(n_x_1,n_y_2,n_z_1,n_step)
          de_122=delta_e0(n_x_1,n_y_2,n_z_2,n_step)
          de_211=delta_e0(n_x_2,n_y_1,n_z_1,n_step)
          de_212=delta_e0(n_x_2,n_y_1,n_z_2,n_step)
          de_221=delta_e0(n_x_2,n_y_2,n_z_1,n_step)
          de_222=delta_e0(n_x_2,n_y_2,n_z_2,n_step)
        else 
          de_111=delta_e1(n_x_1,n_y_1,n_z_1,n_step)
          de_112=delta_e1(n_x_1,n_y_1,n_z_2,n_step)
          de_121=delta_e1(n_x_1,n_y_2,n_z_1,n_step)
          de_122=delta_e1(n_x_1,n_y_2,n_z_2,n_step)
          de_211=delta_e1(n_x_2,n_y_1,n_z_1,n_step)
          de_212=delta_e1(n_x_2,n_y_1,n_z_2,n_step)
          de_221=delta_e1(n_x_2,n_y_2,n_z_1,n_step)
          de_222=delta_e1(n_x_2,n_y_2,n_z_2,n_step)
        endif
        xd =(gas_temp-g_tt(n_x_1))/(g_tt(n_x_2)-g_tt(n_x_1))
        yd =(gas_nden-g_nn(n_y_1))/(g_nn(n_y_2)-g_nn(n_y_1))
        zd =(gas_metl-g_zz(n_z_1))/(g_zz(n_z_2)-g_zz(n_z_1))
        t1=(de_111)*(1.0-zd)+(de_112)*zd
        t2=(de_121)*(1.0-zd)+(de_122)*zd
        s1=(de_211)*(1.0-zd)+(de_212)*zd
        s2=(de_221)*(1.0-zd)+(de_222)*zd
        w1=t1*(1.0-yd)+t2*yd
        w2=s1*(1.0-yd)+s2*yd
        final_de=(w1*(1.0-xd)+w2*xd)

        if(isnan(final_de)) then
          print *, 'Error in getdeltae'
          print *, 't, n, z, xd,yd, zd,t1,t2,s1,s2,w1,w2,final_de',
     +    gas_temp,gas_nden,gas_metl,xd,yd,zd,t1,t2,s1,s2,w1,w2,final_de
          stop
        endif
        return
        end
c---------------------------------------------------------------------
        subroutine calcevoldata(global_dt,n_shd,myid,nproc)
        parameter (n=64,nmax=32,nn=16)
        implicit double precision(a-h,o-z)
        include "mpif.h"
        integer mpi_status,myid,nproc,count,req
        integer status(mpi_status_size)
        common /evol_1_sh0/delta_e0(nn,nn,nn,nmax)
        common /evol_1_sh1/delta_e1(nn,nn,nn,nmax)
        common /evol_2/g_tt_log(nn),g_nn_log(nn),g_zz_log(nn)
        common /evol_3/diff_tt_log_inv,diff_nn_log_inv,diff_zz_log_inv
        common /evol_4/g_tt(nn),g_nn(nn),g_zz(nn)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        integer POS(4)
c       -----------------------------
        g_min_nn=1e-6
        g_max_nn=1e3
c       -----------------------------
        diff_tt_log=(g_t_log(num_t)-g_t_log(1))/(nn-1)
        diff_nn_log=(log10(g_max_nn)-log10(g_min_nn))/(nn-1)
        diff_zz_log=(g_z_log(num_z)-g_z_log(1))/(nn-1)
        diff_tt_log_inv=1.0/diff_tt_log
        diff_nn_log_inv=1.0/diff_nn_log
        diff_zz_log_inv=1.0/diff_zz_log
c       -----------------------------
        do i=1,nn
          g_tt_log(i)=g_t_log(1)+diff_tt_log*(i-1.0)
          g_tt(i)    =10**g_tt_log(i)
        enddo
        do i=1,nn
          g_nn_log(i)=log10(g_min_nn)+diff_nn_log*(i-1.0)
          g_nn(i)    =10**g_nn_log(i)
        enddo
        do i=1,nn
          g_zz_log(i)=g_z_log(1)+diff_zz_log*(i-1.0)
          g_zz(i)    =10**g_zz_log(i)
        enddo
c       -----------------------------
        iwork=0
        ifin =1
        idle =2
        irep =3
        iskip=4
        itag3=3

        call real_time(real_time1)
        num_send=0
        num_recv=0
        if(myid.eq.0) then
         do i=1,nmax
         do j=1,nn
         do k=1,nn
         do l=1,nn
           do
             call mpi_probe(mpi_any_source,itag1,mpi_comm_world,
     +                                         status,ierr)
             isrc=status(mpi_source)
             call mpi_recv(iflag,1,mpi_integer,isrc,itag1,
     +                              mpi_comm_world,status,ierr)
             if(iflag.eq.idle) exit
             call mpi_recv(POS,4,mpi_integer,isrc,itag3,
     +                              mpi_comm_world,status,ierr)
             call mpi_recv(rrr,1,mpi_double_precision,isrc,itag3,
     +                            mpi_comm_world,status,ierr)
             if(isnan(rrr)) then
               print *, 'Error in calcevoldata.'
               print *, 'i,j,k,l,r=',POS(1),POS(2),POS(3),POS(4),rrr
               stop
             endif
             if(n_shd.eq.0) then
               delta_e0(POS(2),POS(3),POS(4),POS(1))=rrr
             else    
               delta_e1(POS(2),POS(3),POS(4),POS(1))=rrr
             endif
             num_recv = num_recv + 1
           enddo
           POS(1) = i
           POS(2) = j
           POS(3) = k
           POS(4) = l
           call mpi_send(iwork,1,mpi_integer,isrc,itag1,
     +                              mpi_comm_world,ierr)
           call mpi_send(POS,4,mpi_integer,isrc,itag3,
     +                              mpi_comm_world,ierr)
           num_send=num_send+1
c          print *, 'NOw in send/recv ',num_send,num_recv 
         enddo
         enddo
         enddo
         enddo
         ifinish = 0
         do iii = num_recv+1, num_send, +1
           do 
           call mpi_probe(mpi_any_source,itag1,mpi_comm_world,
     +                                         status,ierr)
           isrc=status(mpi_source)
           call mpi_recv(iflag,1,mpi_integer,isrc,itag1,
     +                              mpi_comm_world,status,ierr)
           if(iflag.eq.irep) then
             call mpi_recv(POS,4,mpi_integer,isrc,itag3,
     +                              mpi_comm_world,status,ierr)
             call mpi_recv(rrr,1,mpi_double_precision,isrc,itag3,
     +                            mpi_comm_world,status,ierr)
             if(isnan(rrr)) then
               print *, 'Error in calcevoldata.'
               print *, 'i,j,k,l,r=',POS(1),POS(2),POS(3),POS(4),rrr
               stop
             endif
             if(n_shd.eq.0) then
               delta_e0(POS(2),POS(3),POS(4),POS(1))=rrr
             else
               delta_e1(POS(2),POS(3),POS(4),POS(1))=rrr
             endif
             num_recv = num_recv + 1
c            print *,'getting report',num_send,num_recv,isrc
             exit
           else
             call mpi_send(ifin,1,mpi_integer,isrc,itag1,
     +                              mpi_comm_world,ierr)
             ifinish = ifinish + 1
c            print *,'sending finish signal',ifinish,isrc
           endif
           enddo
         enddo
         do i=2,nproc-ifinish, +1
           call mpi_probe(mpi_any_source,itag1,mpi_comm_world,
     +                                         status,ierr)
           isrc=status(mpi_source)
            call mpi_recv(iflag,1,mpi_integer,isrc,itag1,
     &                  mpi_comm_world,status,ierr)
            call mpi_send(ifin,1,mpi_integer,isrc,itag1,
     +                              mpi_comm_world,ierr)
         enddo
       else
         do
           call mpi_send(idle,1,mpi_integer,0,itag1,
     +                              mpi_comm_world,ierr)
           call mpi_recv(iflag,1,mpi_integer,0,itag1,
     +                               mpi_comm_world,status,ierr)
           if(iflag.eq.ifin) then
             exit
           else
             call mpi_recv(POS,4,mpi_integer,0,itag3,
     +                               mpi_comm_world,status,ierr)
             call calcft(POS(1),POS(2),POS(3),POS(4),global_dt,n_shd,rr)
             call mpi_send(irep,1,mpi_integer,0,itag1,
     +                               mpi_comm_world,ierr)
             call mpi_send(POS,4,mpi_integer,0,itag3,
     +                               mpi_comm_world,ierr)
             call mpi_send(rr,1,mpi_double_precision,0,itag3,
     +                               mpi_comm_world,ierr)
           endif
         enddo
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        if(n_shd.eq.0) then
           call mpi_bcast(delta_e0,nn*nn*nn*nmax,mpi_double_precision,
     &                   0,mpi_comm_world,ierr)
        else
           call mpi_bcast(delta_e1,nn*nn*nn*nmax,mpi_double_precision,
     &                   0,mpi_comm_world,ierr)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call real_time(real_time2)
        if(myid.eq.0) then
          print *, 'cal.T. for deltae[sec]=', real_time2-real_time1
        endif
c       call MPI_Finalize(ierr)
c       stop
        return
        end
c---------------------------------------------------------------------
        subroutine calcft(i_step,i_t,i_n,i_z,global_dt,n_shd,de)
        implicit double precision(a-h,o-y)
        parameter (n=64,nmax=32,nn=16)
        common /evol_2/g_tt_log(nn),g_nn_log(nn),g_zz_log(nn)    
        common /evol_3/diff_tt_log_inv,diff_nn_log_inv,diff_zz_log_inv
        common /evol_4/g_tt(nn),g_nn(nn),g_zz(n)
c       -----------------------------
        yp         =0.25
        sun_z      =0.02
        wm         =1.0
        delta_t_max=100.0
        criteria   =0.2
        kBoverhm   =8.24971e+07
        hm         =1.67356e-24
        kB         =1.38065e-16
        sectoyr    =3.17098e-08
        yrtosec    =3.15360e+07
c    ----------------------------
        if(n_shd.eq.0)then 
          dt_real_yr=global_dt*2.0**(-i_step+1.0)
          delta_t_min=dt_real_yr*0.01
          vt_now=g_tt(i_t)
          vn_now=g_nn(i_n)
          vz_now=g_zz(i_z)
          rho=vn_now*hm/(1.0-yp-vz_now*sun_z)
          ve_now=1.5*kBoverhm*rho*vt_now
          call getcoolheatfinal(vt_now,vn_now,vz_now,0,c,h)
          temp_dt=ve_now/abs(h-c)*criteria*sectoyr
          delta_t=min(temp_dt,delta_t_min)
          de=0.0
          do time_now=0.0,dt_real_yr,delta_t
            rho=vn_now*hm/(1.0-yp-vz_now*zsum)
            ve_now=1.5*kBoverhm*rho*vt_now
            call getcoolheatfinal(vt_now,vn_now,vz_now,0,c,h)
            temp_dt=ve_now/abs(h-c)*criteria*sectoyr
            delta_t=min(temp_dt,delta_t_min)
            delta_result=(h-c)*delta_t*yrtosec
            de=de+delta_result
            vt_now=vt_now*(1.0+delta_result/ve_now)
          enddo
          if(vt_now.lt.10.0) then
            vt_now=10.0
            de=1.5*kBoverhm*rho*(g_tt(i_t)-10.0)
          endif
        else
          dt_real_yr=global_dt*2.0**(-i_step+1.0)
          delta_t_min=dt_real_yr*0.01
          vt_now=g_tt(i_t)
          vn_now=g_nn(i_n)
          vz_now=g_zz(i_z)
          rho=vn_now*hm/(1.0-yp-vz_now*sun_z)
          ve_now=1.5*kBoverhm*rho*vt_now
          call getcoolheatfinal(vt_now,vn_now,vz_now,1,c,h)
          temp_dt=ve_now/abs(h-c)*criteria*sectoyr
          delta_t=min(temp_dt,delta_t_min)
          de=0.0
          do time_now=0.0,dt_real_yr,delta_t
            rho=vn_now*hm/(1.0-yp-vz_now*zsum)
            ve_now=1.5*kBoverhm*rho*vt_now
            call getcoolheatfinal(vt_now,vn_now,vz_now,1,c,h)
            temp_dt=ve_now/abs(h-c)*criteria*sectoyr
            delta_t=min(temp_dt,delta_t_min)
            delta_result=(h-c)*delta_t*yrtosec
            de=de+delta_result
            vt_now=vt_now*(1.0+delta_result/ve_now)
          enddo
          if(vt_now.lt.10.0) then
            vt_now=10.0
            de=1.5*kBoverhm*rho*(g_tt(i_t)-10.0)
          endif
        endif
        return
        end
c---------------------------------------------------------------------
        subroutine getcoolheatfinal(gas_temp,gas_nden,gas_metl,
     +                     UV_shield,cooling_rate,heating_rate)
        implicit double precision (a-h,o-z)
        parameter (n=64)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)
        integer UV_shield

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
        gas_nden_log=log10(gas_nden)

        if(UV_shield.eq.1) then
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
     +                 qc12*(x2-xx)*(yy-y1)/((x2-x1)*(y2-y1))+
     +                 qc22*(xx-x1)*(yy-y1)/((x2-x1)*(y2-y1))
          cooling_rate=10.0**result_c_log*gas_nden**2
          heating_rate=0.0
          return
        else 
c    When gas cloud is exposed to the external UV source & CMB.
        diff_t_log=(g_t_log(num_t)-g_t_log(1))/(num_t-1)
        diff_n_log=(g_n_log(num_n)-g_n_log(1))/(num_n-1)
        diff_z_log=(g_z_log(num_z)-g_z_log(1))/(num_z-1)

          if(gas_temp_log.lt.g_t_log(1)) then
            n_x_2=1
            n_x_1=2
          else if(gas_temp_log.gt.g_t_log(num_t)) then
            n_x_1=num_t-1
            n_x_2=num_t
          else
            n_x=((gas_temp_log-g_t_log(1))/diff_t_log)+1
            n_x_1=n_x
            n_x_2=n_x+1
          endif
          if(gas_nden_log.lt.g_n_log(1)) then
            n_y_2=1
            n_y_1=2
          else if(gas_nden_log.gt.g_n_log(num_n)) then
            n_y_1=num_n-1
            n_y_2=num_n
          else
            n_y=((gas_nden_log-g_n_log(1))/diff_n_log)+1
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
            n_z=((gas_metl_log-g_z_log(1))/diff_z_log)+1
            n_z_1=n_z
            n_z_2=n_z+1
          endif
          xd =(gas_temp_log-g_t_log(n_x_1))/diff_t_log
          yd =(gas_nden_log-g_n_log(n_y_1))/diff_n_log
          zd =(gas_metl_log-g_z_log(n_z_1))/diff_z_log
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
c---------------------------------------------------------------------
        subroutine getcoolheatredshift(redshift)
        implicit double precision (a-h,o-z)
        parameter (n=64)
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
              c_UV_tnz_log(n_t,n_n,n_z)=(cq1+(cq2-cq1)*(xx-x1)/(x2-x1))
              h_UV_tnz_log(n_t,n_n,n_z)=(hq1+(hq2-hq1)*(xx-x1)/(x2-x1))
            enddo
          enddo
        enddo
        return
        endif
        end
c---------------------------------------------------------------------
ccc    Read heating/cooing data.
        subroutine readcoolingdata(myid)
        implicit double precision (a-h,o-z)
        include "mpif.h"
        parameter(n=64)
        common /cool_1/cool_UV_log(n,n,n,n),heat_UV_log(n,n,n,n)
        common /cool_2/g_t_log(n),g_n_log(n),g_z_log(n),g_r_log(n)
        common /cool_3/num_t,num_n,num_z,num_r,num_t2,num_z2
        common /cool_4/cool_NO_log(n,n),g_t2_log(n),g_z2_log(n)
        common /cool_5/c_UV_tnz_log(n,n,n),h_UV_tnz_log(n,n,n)
        if(myid.eq.0) then
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
c---------------------------------------------------------------------
         subroutine real_time(seconds)
         implicit none
         integer clock_count,clock_max,clock_rate
         real*8 seconds
         call system_clock(clock_count,clock_rate,clock_max)
         seconds = real(clock_count,kind=8)/real(clock_rate,kind=8)
         return
         end
c---------------------------------------------------------------------
