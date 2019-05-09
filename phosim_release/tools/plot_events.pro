;window,0,xsize=800,ysize=800
set_plot,'PS'
device,/color
device,xsize=8.5,/inches
device,ysize=8.5,/inches
device,xoffset=0.0,/inches
device,yoffset=0.0,/inches
!x.thick=2
!y.thick=2
!p.charthick=1.4

!p.multi=[0,2,2]
loadct,39
;device,decomposed=0
data=mrdfits('output.fits',1)
zp=median(data(where(data.surface eq 212)).z)
xp=median(data(where(data.surface eq 212)).x)
yp=median(data(where(data.surface eq 212)).y)

phi=0.0 & theta=0.0 & psi=0.0

a11=cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi)
a12=cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi)
a13=sin(psi)*sin(theta)
a21=-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi)
a22=-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi)
a23=cos(psi)*sin(theta)
a31=sin(theta)*sin(phi)
a32=-sin(theta)*cos(phi)
a33=cos(theta)

xp=data.x*a11+data.y*a12+data.z*a13
yp=data.x*a21+data.y*a22+data.z*a23
zp=data.x*a31+data.y*a32+data.z*a33
surf=data.surface

for ii=0,8 do begin

if ii eq 0 then plot,[min(xp),max(xp)],[min(zp),max(zp)],/nodata,/xstyle,/ystyle,title='Atmosphere'
if ii eq 1 then plot,[min(xp),max(xp)],[min(zp(where(surf eq 203))),max(zp(where(surf eq 202)))],/nodata,/xstyle,/ystyle,title='Mirrors'
if ii eq 2 then plot,[min(xp(where(surf eq 204))),max(xp(where(surf eq 204)))],$
[min(zp(where(surf eq 204))),max(zp(where(surf eq 212)))],/nodata,/xstyle,/ystyle,title='Lenses'
if ii eq 3 then plot,[min(xp(where(surf eq 212))),max(xp(where(surf eq 212)))],$
[min(zp(where(surf eq 212))),max(zp(where(surf eq 301)))],/nodata,/xstyle,/ystyle,title='Detector'

if ii eq 4 then plot,[min(xp),max(xp)],[min(yp),max(yp)],/nodata,/xstyle,/ystyle,title='Atmosphere'
if ii eq 5 then plot,[min(xp),max(xp)],[min(yp(where(surf eq 203))),max(yp(where(surf eq 202)))],/nodata,/xstyle,/ystyle,title='Mirrors'
if ii eq 6 then plot,[min(xp(where(surf eq 204))),max(xp(where(surf eq 204)))],$
[min(yp(where(surf eq 204))),max(yp(where(surf eq 204)))],/nodata,/xstyle,/ystyle,title='Lenses'
if ii eq 7 then plot,[min(xp(where(surf eq 212))),max(xp(where(surf eq 212)))],$
[min(yp(where(surf eq 212))),max(yp(where(surf eq 301)))],/nodata,/xstyle,/ystyle,title='Detector'

wavelength=fltarr(N_elements(data))
detected=fltarr(N_elements(data))
detected_postatm=fltarr(N_elements(data))
detected_postaperture=fltarr(N_elements(data))
detected_postfilter=fltarr(N_elements(data))
wc=0L
x1=0 & z1=0
photoncount=0L
for i=0L,N_elements(data)-1 do begin


x0=xp(i)
z0=zp(i)
if ii gt 3 then z0=yp(i)

if data(i).surface eq 0 then begin
    photoncount=photoncount+1
 wavelength(wc)=data(i).z
 wc=wc+1
endif


if data(i).surface eq 301 then detected(wc-1)=1
if data(i).surface eq 105 then detected_postatm(wc-1)=1
if data(i).surface eq 212 then detected_postfilter(wc-1)=1
if data(i).surface eq 201 then detected_postaperture(wc-1)=1


    if photoncount gt 500 and ii lt 8 then goto,skipq
    if photoncount gt 500 then goto,skipp

cc=0
;if data(i).surface eq 100 then cc=200
if data(i).surface eq 101 then cc=160
if data(i).surface eq 102 then cc=200
if data(i).surface eq 103 then cc=160
if data(i).surface eq 104 then cc=200
if data(i).surface eq 105 then cc=160

if data(i).surface eq 201 then cc=200
if data(i).surface eq 202 then cc=160
if data(i).surface eq 203 then cc=200
if data(i).surface eq 204 then cc=40
if data(i).surface eq 205 then cc=120
if data(i).surface eq 206 then cc=40
if data(i).surface eq 207 then cc=120
if data(i).surface eq 208 then cc=40
if data(i).surface eq 209 then cc=120
if data(i).surface eq 210 then cc=40
if data(i).surface eq 211 then cc=120
if data(i).surface eq 212 then cc=40
if data(i).surface eq 300 then cc=120
if data(i).surface eq 301 then cc=250
if data(i).surface eq 213 then cc=200
if data(i).surface eq 214 then cc=200
if data(i).surface eq 302 then cc=0
if (ii eq 2 or ii eq 3 or ii eq 6 or ii eq 7) and cc eq 200 then goto,skipplot
if (ii eq 2 or ii eq 3 or ii eq 6 or ii eq 7) and cc eq 160 then goto,skipplot
if cc eq 0 then goto,skipplot
if ii eq 4 and data(i).surface gt 200 then goto,skipplot

if data(i).surface ne 0 then oplot,[x0,x1],[z0,z1],color=cc

skipplot:
x1=x0
z1=z0

skipp:
endfor
skipq:
endfor



wavelength=wavelength(0:(wc-1))
detected=detected(0:(wc-1))
detected_postatm=detected_postatm(0:(wc-1))
detected_postfilter=detected_postfilter(0:(wc-1))
detected_postaperture=detected_postaperture(0:(wc-1))

!p.multi=0

q=histogram(wavelength,min=0,bin=0.005,max=1.2)
xx=findgen(N_elements(q))*0.005
plot,xx,q,psym=10,xtitle='Wavelength (microns)',ytitle='Photons/bin'

q=histogram(wavelength(where(detected_postatm eq 1)),min=0,bin=0.005,max=1.2)
xx=findgen(N_elements(q))*0.005
oplot,xx,q,psym=10,color=250,linestyle=2

q=histogram(wavelength(where(detected_postaperture eq 1)),min=0,bin=0.005,max=1.2)
xx=findgen(N_elements(q))*0.005
oplot,xx,q,psym=10,color=200

q=histogram(wavelength(where(detected_postfilter eq 1)),min=0,bin=0.005,max=1.2)
xx=findgen(N_elements(q))*0.005
oplot,xx,q,psym=10,color=120,linestyle=2

q=histogram(wavelength(where(detected eq 1)),min=0,bin=0.005,max=1.2)
xx=findgen(N_elements(q))*0.005
oplot,xx,q,psym=10,color=40



data=mrdfits('focalplane.fits',0)
xv=0.0
yv=0.0
for i=0L,4095 do begin
for j=0L,4095 do begin
    xv=xv+i*data(i,j)
    yv=yv+j*data(i,j)
endfor
endfor
xv=xv/total(data)
yv=yv/total(data)
;window,1,xsize=500,ysize=500


!p.multi=0
plot,[0,100],[0,100],/nodata
loadct,27
contour,alog10(data((xv-50):(xv+49),(yv-50):(yv+49))>1),findgen(100),findgen(100),/overplot,/fill,nlevels=11
;tvscl,congrid(alog10(data((xv-50):(xv+50),(yv-50):(yv+50))>1),500,500)
device,/close
set_plot,'X'
END
