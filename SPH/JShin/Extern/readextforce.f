c234567
      subroutine readextforce(array,nsize)
      real*8 array(*),r


      open(1,file='init_ennha_potential.dat',form='unformatted')

      read(1) ((r,array(i)),i=1,nsize)

      close(1)


      return
      end



      subroutine readstarparticle(mass,x,y,z,vx,vy,vz,nsize)
      integer*8 nsize
      real*8 mass(*),x(*),y(*),z(*)
      real*8 vx(*),vy(*),vz(*)


      open(1,file='init_eunha_disc.dat',form='unformatted')
      read(1) ((mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)),i=1,nsize)
      close(1)
      return
      end

      subroutine readsphparticle(mass,x,y,z,vx,vy,vz,entropy,nsize)
      integer*8 nsize
      real*8 mass(*),x(*),y(*),z(*)
      real*8 vx(*),vy(*),vz(*),entropy(*)


      print *, 'Now reading ',nsize,' gas particles'
      open(1,file='init_eunha_gas.dat',form='unformatted')
      read(1) ((mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),
     &       entropy(i)),i=1,nsize)
      close(1)
      return
      end
