c234567
c     parameter(N=90909)
      parameter(N=909091)
      integer*8 nsize
      real*8 mass(N),x(N),y(N),z(N)
      real*8 vx(N),vy(N),vz(N),entropy(N)


      nsize = N


      print *, 'Now reading ',nsize,' gas particles'
c     open(1,file='init_eunha_gas.dat',form='unformatted')
      open(1,file='init_eunha_disc.dat',form='unformatted')
c     read(1) ((mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),
c    &       entropy(i)),i=1,nsize)
      read(1) ((mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)),
     &        i=1,nsize)
      close(1)
      open(2,file='list.dat')
      do i =1, nsize,+1
      write(2,100) x(i),y(i),z(i),vx(i),vy(i),vz(i)
      enddo
      close(2)
100   format(6(g15.7,1x))
      stop
      end
