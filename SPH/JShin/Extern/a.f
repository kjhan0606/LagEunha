c234567
	  real*8 array(128),r(128)
	  integer nsize

	  nsize = 128


	  open(1,file='init_ennha_potential.dat',form='unformatted')

	  read(1) ((r(i),array(i)) ,i=1,nsize)


	  do i = 1, nsize
	     print *, r(i),array(i)
	  enddo

	  close(1)


	  stop
	  end

