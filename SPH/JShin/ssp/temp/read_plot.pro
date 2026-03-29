pro test

readcol, 'SN_rate_kroupa.dat', zzz,fit_a,fit_b,fit_t1,fit_t2,format='f,f,f,f,f'

num_t=20.0&num_z=2.0

min_t=1.0
max_t=100.0
diff_t=(max_t-min_t)/(num_t-1.0)
t=findgen(num_t)*diff_t+min_t

min_z=0.0001
max_z=0.03
diff_z=(alog10(max_z)-alog10(min_z))/(num_z-1.0)
z=10^(findgen(num_z)*diff_z+alog10(min_z))

prob   =fltarr(num_t,num_z)
prob_1 =fltarr(num_t,num_z)
prob_2 =fltarr(num_t,num_z)
del_t=1.0

;----------------------------
for i=1,num_t do begin
  for j=1,num_z do begin
     aaa=interpol(fit_a, zzz,z(j-1))
     bbb=interpol(fit_b, zzz,z(j-1))
     tt1=interpol(fit_t1,zzz,z(j-1))
     tt2=interpol(fit_t2,zzz,z(j-1))
     if(t(i-1) le tt1 or t(i-1)+del_t ge tt2) then begin
       if(t(i-1) le tt1) then begin
         temp_prob=0.0
       endif else begin
         temp_prob=1.0
       endelse
       temp_prob_1=0.0
       temp_prob_2=0.0
     endif else begin
       temp_prob_1=10.0^(aaa*(t(i-1)+del_t)+bbb)-10.0^(aaa*(t(i-1))+bbb)
       temp_prob_2=10.0^(aaa*(tt2         )+bbb)-10.0^(aaa*(t(i-1))+bbb)
       temp_prob=temp_prob_1/temp_prob_2
     endelse
     prob(i-1,j-1)  =temp_prob
     prob_1(i-1,j-1)=temp_prob_1
     prob_2(i-1,j-1)=temp_prob_2
     print, temp_prob,temp_prob_1,temp_prob_2
  endfor
endfor
;----------------------------
stop
end
