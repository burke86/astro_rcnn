PRO CONVERTTIFF

data=mrdfits('/Users/Shared/full_plane/focalplane_0_0.fits')

data = rebin(data,1024,1024)
help,data
scale=1

data=alog(data>1)
data=byte((data/max(data)*255*scale)<255)


loadct,5
tvlct,red,green,blue,/get

newdata=fltarr(3,1024,1024)

for i=0L,1023 do begin
if i mod 100 eq 0 then print,i
for j=0L,1023 do begin
    newdata(0,i,j)=red(data(i,j))
    newdata(1,i,j)=green(data(i,j))
    newdata(2,i,j)=blue(data(i,j))
endfor
endfor


write_tiff,'~/Desktop/focalplane_0_0.tiff',newdata
;print, dat_min
;print, dat_max
END
