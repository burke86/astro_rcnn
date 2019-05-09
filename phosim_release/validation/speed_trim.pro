chipsize=13.3

openw,1,'validation_4A_00_catalog'
printf,1,'Opsim_obshistid 4001'
printf,1,'Opsim_rawseeing 0.65'
printf,1,'SIM_NSNAP 1'
printf,1,'SIM_VISTIME 15.0'
ss='object 0.0 0.0 0.0 35.0 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none'
printf,1,ss
close,1
openw,1,'validation_4A_00_commands'
printf,1,'zenith_v 1000.0'
printf,1,'raydensity 0.0'
close,1

readcol,'mag_dist',mm,point,diffuse
mm=double(mm)
point=point*chipsize/60.*chipsize/60./(!Pi*2.1*2.1)
diffuse=diffuse*chipsize/60.*chipsize/60./(!Pi*2.1*2.1)
openw,1,'validation_4A_01_catalog'
printf,1,'Opsim_obshistid 4000'
printf,1,'Opsim_rawseeing 0.65'
printf,1,'SIM_NSNAP 2'
printf,1,'SIM_VISTIME 33.0'
under=1.0
seed=1000
for ii=0,N_elements(mm)-1 do begin
   nn=long((point(ii))/under)
   if nn eq 0 then goto,jump1a
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*chipsize/60.0-chipsize/2.0/60.0
      dec=randomu(seed,1)*chipsize/60.0-chipsize/2.0/60.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,ss
   endfor
   jump1a:
endfor
for ii=0,N_elements(mm)-1 do begin
   nn=long((diffuse(ii))/under)
   if nn eq 0 then goto,jump1b
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*chipsize/60.0-chipsize/2.0/60.0
      dec=randomu(seed,1)*chipsize/60.0-chipsize/2.0/60.0
      angle=randomu(seed,1)*360.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 sersic2d 1.0 0.8 '+string(angle)+' 1.0 none none'
      printf,1,ss
   endfor
   jump1b:
endfor
close,1


readcol,'mag_dist',mm,point,diffuse
mm=double(mm)
range=3.0
point=point*range/60.*range/60./(!Pi*2.1*2.1)
diffuse=diffuse*range/60.*range/60./(!Pi*2.1*2.1)
openw,1,'../examples/small_catalog'
printf,1,'Opsim_obshistid 99999999'
printf,1,'Opsim_rawseeing 0.65'
printf,1,'SIM_NSNAP 1'
printf,1,'SIM_VISTIME 15.0'
under=1.0
seed=1000
for ii=0,N_elements(mm)-1 do begin
   nn=long((point(ii))/under)
   if nn eq 0 then goto,jump2a
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*range/60.0-range/60.0/2.0
      dec=randomu(seed,1)*range/60.0-range/60.0/2.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,ss
   endfor
   jump2a:
endfor
for ii=0,N_elements(mm)-1 do begin
   nn=long((diffuse(ii))/under)
   if nn eq 0 then goto,jump2b
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*range/60.0-range/60.0/2.0
      dec=randomu(seed,1)*range/60.0-range/60.0/2.0
      angle=randomu(seed,1)*360.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 sersic2d 1.0 0.8 '+string(angle)+' 1.0 none none'
      printf,1,ss
   endfor
   jump2b:
endfor
close,1


readcol,'mag_dist',mm,point,diffuse
mm=double(mm)
range=12.0
point=point*range/60.*range/60./(!Pi*2.1*2.1)
diffuse=diffuse*range/60.*range/60./(!Pi*2.1*2.1)
openw,1,'../examples/large_catalog'
printf,1,'Opsim_obshistid 99999999'
printf,1,'Opsim_rawseeing 0.65'
printf,1,'SIM_NSNAP 1'
printf,1,'SIM_VISTIME 15.0'
under=1.0
seed=1000
for ii=0,N_elements(mm)-1 do begin
   nn=long((point(ii))/under)
   if nn eq 0 then goto,jump3a
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*range/60.0-range/60.0/2.0
      dec=randomu(seed,1)*range/60.0-range/60.0/2.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,ss
   endfor
   jump3a:
endfor
for ii=0,N_elements(mm)-1 do begin
   nn=long((diffuse(ii))/under)
   if nn eq 0 then goto,jump3b
   for jj=0L,nn-1 do begin
      ra=randomu(seed,1)*range/60.0-range/60.0/2.0
      dec=randomu(seed,1)*range/60.0-range/60.0/2.0
      angle=randomu(seed,1)*360.0
      ss='object 0.0 '+string(ra)+' '+string(dec)+' '+string(mm(ii))+' ../sky/sed_flat.txt 0 0 0 0 0 0 sersic2d 1.0 0.8 '+string(angle)+' 1.0 none none'
      printf,1,ss
   endfor
   jump3b:
endfor
close,1


openw,1,'validation_4A_01_commands'
ss='zenith_v '+string(-2.5*alog10(10.^(-0.4*22.08)/under))
printf,1,ss
close,1



end
