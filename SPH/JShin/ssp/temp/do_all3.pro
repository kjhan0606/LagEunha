pro test

num_t= 100.0
delta_t = 15

init_t=34.0
swlasttime=init_t
init_m=1e3
mass     =fltarr(num_t)
nowtime  =fltarr(num_t)
swtime   =fltarr(num_t)
for i=1, num_t do begin
   openw, 1, 'input'
   nowtime(i-1) = init_t+delta_t*(i-1)
   printf,1, 0.02,init_m,nowtime(i-1),swlasttime
   close, 1
   spawn, './temp'
   print, 'ing ...', i
   openr, 1, 'output'
   readf, 1, swlasttime,delta_m
   close, 1
   if(i eq 1) then begin
     mass(0) = init_m
   endif else begin
     mass(i-1)   = mass(i-2)-delta_m
   endelse
   swtime(i-1) = swlasttime
endfor

save, filename='store3.sav',nowtime,swtime,mass
stop
end
