c234567
      real*8 gasdens,gasmetal,netheatcoolrate,gastemp
	  real*8 redshift,coolrate,heatrate



	  redshift = 128




      call readcoolingdata()
	  call getcoolheatredshift(redshift)
	  

	  gasdens = 1.E-20
	  gasmetal = 0.02
	  gastemp = 259


	  call getcoolheatfinal(gastemp,gasdens,gasmetal,0,coolrate,heatrate)



      print *, 'den/temp/metal',gasdens,gasmetal,gastemp
	  print *, coolrate,heatrate

	  stop
	  end



