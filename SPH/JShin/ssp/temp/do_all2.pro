pro test

num_z=10.0

min_z=0.0001
max_z=0.03
diff_z=(alog10(max_z)-alog10(min_z))/(num_z-1.0)
z=10^(findgen(num_z)*diff_z+alog10(min_z))

energyE=dblarr(num_z)
massE  =fltarr(num_z)
xE     =fltarr(num_z)
yE     =fltarr(num_z)
zE     =fltarr(num_z)

temp_arr=dblarr(5)

for k=1, num_z do begin
   openw, 1, 'input'
   printf,1, z(k-1),1e4
   close, 1
   spawn, './temp'
   print, 'ing...',k
   openr, 1, 'output'
   readf, 1, temp_arr
   close, 1
   energyE(k-1) = temp_arr(0)
   massE(k-1)   = float(temp_arr(1))
   xE(k-1)      = float(temp_arr(2))
   yE(k-1)      = float(temp_arr(3))
   zE(k-1)      = float(temp_arr(4))
endfor
save,filename='store2.sav',z,energyE,massE,zE
stop
end
