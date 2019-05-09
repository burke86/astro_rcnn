pro optics_validation

!EXCEPT=0
sigma = fltarr(2,5)
  for n=1,2 do begin
     for m=1,5 do begin

        temp = strcompress('lsst_e_21'+string(n)+string(m)+'_f2_R22_S11_E000.fits',/remove)

FITS_READ, temp, data, header
length = size(data)


cx = 0
cy = 0
high = 0
for i=0,length[1]-1 do begin
for j=0,length[2]-1 do begin
	if data[i,j] gt high then begin
		cx = i +1
		cy = j +1
		high = data[i,j]
	endif
endfor
endfor
;print, cx, cy, data[cx-1,cy-1], high

data=data-median(data)

x=dblarr(200*200.)
y=dblarr(200*200.)
w=dblarr(200*200.)
k=0L
for i=cx-100,cx+99 do begin
for j=cy-100,cy+99 do begin
        x(k)=i-(length[1]/2.0-0.5)
        y(k)=j-(length[2]/2.0-0.5)
        w(k)=data(i,j)
        k=k+1
endfor
endfor
x=x[0:k-1]
y=y[0:k-1]
w=w[0:k-1]

t1=0.
for i=0L,N_elements(x)-1 do t1=t1+(x(i)*y(i)*w(i))
t2=0.
for i=0L,N_elements(x)-1 do t2=t2+(w(i))
t3=0.
for i=0L,N_elements(x)-1 do t3=t3+(x(i)*w(i))
t4=0.
for i=0L,N_elements(x)-1 do t4=t4+(y(i)*w(i))
t5=0.
for i=0L,N_elements(x)-1 do t5=t5+(x(i)*x(i)*w(i))
t6=0.
for i=0L,N_elements(x)-1 do t6=t6+(y(i)*y(i)*w(i))

resultx=(t5/t2-t3*t3/t2/t2)
resulty=(t6/t2-t4*t4/t2/t2)

r= sqrt(resultx*resulty)
seeing=r
medx=0.
medy=0.

last = 0.
for ty=0,100 do begin
error=0
weight=(exp(-((x-medx)^2+(y-medy)^2)/2.0/seeing/seeing)/((sqrt(2*!Pi*seeing*seeing))^2))

t1=0.
for i=0L,N_elements(x)-1 do t1=t1+(x(i)*y(i)*w(i)*weight(i))
t2=0.
for i=0L,N_elements(x)-1 do t2=t2+(w(i)*weight(i))
t3=0.
for i=0L,N_elements(x)-1 do t3=t3+(x(i)*w(i)*weight(i))
t4=0.
for i=0L,N_elements(x)-1 do t4=t4+(y(i)*w(i)*weight(i))
t5=0.
for i=0L,N_elements(x)-1 do t5=t5+(x(i)*x(i)*w(i)*weight(i))
t6=0.
for i=0L,N_elements(x)-1 do t6=t6+(y(i)*y(i)*w(i)*weight(i))

covxy=(t1/t2-t3*t4/t2/t2)
resultx=(t5/t2-t3*t3/t2/t2)
resulty=(t6/t2-t4*t4/t2/t2)
medx=(t3/t2)
medy=(t4/t2)


if medx lt -length[1]/4 then begin
    medx=-length[1]/4.
    error=error+1
endif
if medx gt length[1]/4 then begin
    medx=length[1]/4.
    error=error+1
endif
if medy lt -length[2]/4 then begin
    medy=-length[2]/4.
    error=error+1
endif
if medy gt length[2]/4 then begin
    medy=length[2]/4.
    error=error+1
endif
if resultx lt 1.0/sqrt(2.0) then begin
    resultx=1.0/sqrt(2.0)
    error=error+1
endif
if resulty lt 1.0/sqrt(2.0) then begin
    resulty=1.0/sqrt(2.0)
    error=error+1
endif

seeing=sqrt(resultx+resulty)
rms=seeing
ellip=((resultx-resulty)^2+(2.0*covxy)^2)/(resultx+resulty)^2
pa=0.5*atan((2.0*covxy),(resultx-resulty))
e1=(resultx-resulty)/(resultx+resulty)
e2=(2.0*covxy)/(resultx+resulty)
ellip=sqrt(ellip)
;if ty mod 10 eq 0 or ty lt 10 then print,rms,ty;,ellip,pa,e1,e2,medx,medy,error,ty
seeing=rms
IF (last eq rms) THEN break
last = rms
endfor

sigma[n-1,m-1]=rms
endfor
endfor


;make contour plots

set_plot,'PS'
device,filename='optics.ps'
device,xsize=10.0,/inches
device,ysize=7.0,/inches
device,/color
!x.thick=2
!y.thick=2
!p.charsize=1.0
!p.charthick=2
!p.multi=[0,5,2,0,0]
!y.margin=[4,3]
!x.margin=[8,3]
loadct,39
readcol,'../raytrace/version',a,b,format='(A,F)',/SILENT
ver='Revision: '+string(fix(b[0]),format='(I6)')

temp = strarr(10)
ii=0
for n=1,2 do begin
     for m=1,5 do begin

        temp(ii) =        strcompress('lsst_e_21'+string(n)+string(m)+'_f2_R22_S11_E000.fits',/remove)
ii++
endfor
endfor
data0=mrdfits(temp(0),0)
data1=mrdfits(temp(1),0)
data2=mrdfits(temp(2),0)
data3=mrdfits(temp(3),0)
data4=mrdfits(temp(4),0)
data5=mrdfits(temp(5),0)
data6=mrdfits(temp(6),0)
data7=mrdfits(temp(7),0)
data8=mrdfits(temp(8),0)
data9=mrdfits(temp(9),0)

image0=data0[1905:2095,1911:2161]
image1=data1[1905:2095,1911:2161]
image2=data2[1905:2095,1911:2161]
image3=data3[1905:2095,1911:2161]
image4=data4[1905:2095,1911:2161]
image5=data5[1905:2095,1911:2161]
image6=data6[1905:2095,1911:2161]
image7=data7[1905:2095,1911:2161]
image8=data8[1905:2095,1911:2161]
image9=data9[1905:2095,1911:2161]


sss0=size(image0)
sss1=size(image1)
sss2=size(image2)
sss3=size(image3)
sss4=size(image4)
sss5=size(image5)
sss6=size(image6)
sss7=size(image7)
sss8=size(image8)
sss9=size(image9)

;print,sss1
xx0=(findgen(sss0[1])-sss0[1]/2)*0.1
yy0=(findgen(sss0[2])-sss0[2]/2)*0.1
xx1=(findgen(sss1[1])-sss1[1]/2)*0.1
yy1=(findgen(sss1[2])-sss1[2]/2)*0.1
xx2=(findgen(sss2[1])-sss2[1]/2)*0.1
yy2=(findgen(sss2[2])-sss2[2]/2)*0.1
xx3=(findgen(sss3[1])-sss3[1]/2)*0.1
yy3=(findgen(sss3[2])-sss3[2]/2)*0.1
xx4=(findgen(sss4[1])-sss4[1]/2)*0.1
yy4=(findgen(sss4[2])-sss4[2]/2)*0.1
xx5=(findgen(sss5[1])-sss5[1]/2)*0.1
yy5=(findgen(sss5[2])-sss5[2]/2)*0.1
xx6=(findgen(sss6[1])-sss6[1]/2)*0.1
yy6=(findgen(sss6[2])-sss6[2]/2)*0.1
xx7=(findgen(sss7[1])-sss7[1]/2)*0.1
yy7=(findgen(sss7[2])-sss7[2]/2)*0.1
xx8=(findgen(sss8[1])-sss8[1]/2)*0.1
yy8=(findgen(sss8[2])-sss8[2]/2)*0.1
xx9=(findgen(sss9[1])-sss9[1]/2)*0.1
yy9=(findgen(sss9[2])-sss9[2]/2)*0.1

contour,alog(image0>1),xx0,yy0,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='optics '+strcompress(sigma[0,0]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image1>1),xx1,yy1,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='optics '+strcompress(sigma[0,1]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image2>1),xx2,yy2,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='optics '+strcompress(sigma[0,2]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image3>1),xx3,yy3,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='optics '+strcompress(sigma[0,3]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image4>1),xx4,yy4,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='optics '+strcompress(sigma[0,4]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image5>1),xx5,yy5,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='tracking '+strcompress(sigma[1,0]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image6>1),xx6,yy6,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='tracking '+strcompress(sigma[1,1]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image7>1),xx7,yy7,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='tracking '+strcompress(sigma[1,2]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image8>1),xx8,yy8,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='tracking '+strcompress(sigma[1,3]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'
contour,alog(image9>1),xx9,yy9,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,subtitle='tracking '+strcompress(sigma[1,4]*0.02*2*sqrt(2*alog(2)),/remove)+' (arcseconds)'


device,/close

end


