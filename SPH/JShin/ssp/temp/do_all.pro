pro test

num_t=20.0&num_z=5.0

min_t=1.0
max_t=100.0
diff_t=(max_t-min_t)/(num_t-1.0)
t=findgen(num_t)*diff_t+min_t

min_z=0.0001
max_z=0.03
diff_z=(alog10(max_z)-alog10(min_z))/(num_z-1.0)
z=10^(findgen(num_z)*diff_z+alog10(min_z))

prob=fltarr(num_t,num_z)

for i=1, num_t do begin
    for k=1, num_z do begin
      openw, 1, 'input'
      printf,1, t(i-1),1.0,1.0-0.25-z(k-1),0.25,z(k-1),1e3
      close, 1
      spawn, './temp'
      openr, 1, 'output'
      readf, 1, temp_prob
      close, 1
      prob(i-1,k-1)=temp_prob
    endfor
endfor
save,filename='store.sav',t,z,prob
stop
end
