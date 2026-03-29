c234567
      real red, Temp,inpOmegaB,inpOmegaC,inpOmegaL,inpHOinp,inpTnow
      real inpYp


      red = 47
      inpOmegaB = 0.044
      inpOmegaC = 0.26
      inpOmegaL = 0.74
      inpHOinp = 72
      inpTnow = 2.725
      inpYp = 0.25

      call myrecfast(red,Temp,inpOmegaB,inpOmegaC,inpOmegaL,inpHOinp,
     &    inpTnow, inpYp)
      print *, Temp

      stop
      end
