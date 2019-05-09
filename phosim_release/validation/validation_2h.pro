

pro validation_2H,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 2H'
readcol,'senM_35_19_50.txt',a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,/silent



sens=dblarr(55,22,35)
sensp=dblarr(55,22,35)
senspn=dblarr(22,35)

for ii=0L,55-1 do begin
for jj=0L,22-1 do begin
for kk=0L,35-1 do begin
   nn=kk*19+(jj-3)
   if jj lt 3 then goto,skip
   if ii eq 5 then sens(ii,jj,kk)=a(nn)
   if ii eq 6 then sens(ii,jj,kk)=b(nn)
   if ii eq 7 then sens(ii,jj,kk)=c(nn)
   if ii eq 8 then sens(ii,jj,kk)=d(nn)
   if ii eq 9 then sens(ii,jj,kk)=e(nn)
   if ii eq 10 then sens(ii,jj,kk)=f(nn)
   if ii eq 11 then sens(ii,jj,kk)=g(nn)
   if ii eq 12 then sens(ii,jj,kk)=h(nn)
   if ii eq 13 then sens(ii,jj,kk)=i(nn)
   if ii eq 14 then sens(ii,jj,kk)=j(nn)
   if ii eq 15 then sens(ii,jj,kk)=k(nn)
   if ii eq 16 then sens(ii,jj,kk)=l(nn)
   if ii eq 17 then sens(ii,jj,kk)=m(nn)
   if ii eq 18 then sens(ii,jj,kk)=n(nn)
   if ii eq 19 then sens(ii,jj,kk)=o(nn)
   if ii eq 20 then sens(ii,jj,kk)=p(nn)
   if ii eq 21 then sens(ii,jj,kk)=q(nn)
   if ii eq 22 then sens(ii,jj,kk)=r(nn)
   if ii eq 23 then sens(ii,jj,kk)=s(nn)
   if ii eq 24 then sens(ii,jj,kk)=t(nn)
   if ii eq 25 then sens(ii,jj,kk)=u(nn)
   if ii eq 26 then sens(ii,jj,kk)=v(nn)
   if ii eq 27 then sens(ii,jj,kk)=w(nn)
   if ii eq 28 then sens(ii,jj,kk)=x(nn)
   if ii eq 29 then sens(ii,jj,kk)=y(nn)
   if ii eq 30 then sens(ii,jj,kk)=z(nn)
   if ii eq 31 then sens(ii,jj,kk)=aa(nn)
   if ii eq 32 then sens(ii,jj,kk)=ab(nn)
   if ii eq 33 then sens(ii,jj,kk)=ac(nn)
   if ii eq 34 then sens(ii,jj,kk)=ad(nn)
   if ii eq 35 then sens(ii,jj,kk)=ae(nn)
   if ii eq 36 then sens(ii,jj,kk)=af(nn)
   if ii eq 37 then sens(ii,jj,kk)=ag(nn)
   if ii eq 38 then sens(ii,jj,kk)=ah(nn)
   if ii eq 39 then sens(ii,jj,kk)=ai(nn)
   if ii eq 40 then sens(ii,jj,kk)=aj(nn)
   if ii eq 41 then sens(ii,jj,kk)=ak(nn)
   if ii eq 42 then sens(ii,jj,kk)=al(nn)
   if ii eq 43 then sens(ii,jj,kk)=am(nn)
   if ii eq 44 then sens(ii,jj,kk)=an(nn)
   if ii eq 45 then sens(ii,jj,kk)=ao(nn)
   if ii eq 46 then sens(ii,jj,kk)=ap(nn)
   if ii eq 47 then sens(ii,jj,kk)=aq(nn)
   if ii eq 48 then sens(ii,jj,kk)=ar(nn)
   if ii eq 49 then sens(ii,jj,kk)=as(nn)
   if ii eq 50 then sens(ii,jj,kk)=at(nn)
   if ii eq 51 then sens(ii,jj,kk)=au(nn)
   if ii eq 52 then sens(ii,jj,kk)=av(nn)
   if ii eq 53 then sens(ii,jj,kk)=aw(nn)
   if ii eq 54 then sens(ii,jj,kk)=ax(nn)
   skip:
endfor
endfor
endfor




xmargo=!x.margin
ymargo=!y.margin

;!x.margin=[0,0]
;!y.margin=[0,0]
r=byte([findgen(128)/127.*255.,fltarr(128)+255.0])
g=byte([findgen(96)/95.*255.,255.0+fltarr(64),255.0-findgen(96)/95.*255.])
b=byte([0.0,fltarr(127)+255.0,255.0-findgen(128)/127.*255.])
tvlct,r,g,b


device,decomposed=0
aa=findgen(33)/32*2*!Pi
usersym,cos(aa),sin(aa)

maxy=max([abs(sens)])*1.1

names=strarr(55)
names(0)='M13 Piston'
names(1)='M13 Decenter x'
names(2)='M13 Decenter y'
names(3)='M13 Tilt x'
names(4)='M13 Tilt y'
names(5)='M2 Piston'
names(6)='M2 Decenter x'
names(7)='M2 Decenter y'
names(8)='M2 Tilt x'
names(9)='M2 Tilty'
names(10)='Camera Piston'
names(11)='Camera Decenter x'
names(12)='Camera Decenter y'
names(13)='Camera Tilt x'
names(14)='Camera Tilt y'
j=1L
for i=15,34 do begin
   ss='M13 Bending Mode'+string(j,format='(I2)')
   names(i)=ss
   j=j+1
endfor
j=1L
for i=35,54 do begin
   ss='M2 Bending Mode'+string(j,format='(I2)')
   names(i)=ss
   j=j+1
endfor

tolerance=dblarr(23)
tolerance(3)=30.7650
tolerance(4)=81.8775
tolerance(5)=81.8775
tolerance(6)=32.1567
tolerance(7)=32.1567
tolerance(8)=60.1262
tolerance(9)=60.1262
tolerance(10)=11.3430
tolerance(11)=25.8841
tolerance(12)=25.8841
tolerance(13)=45.2123
tolerance(14)=45.2123
tolerance(15)=11.8523
tolerance(16)=11.8523
tolerance(17)=22.5635
tolerance(18)=22.5635
tolerance(19)=37.2334
tolerance(20)=37.2334
tolerance(21)=4.1688


dx=dblarr(55)
dx(0)=100.0
dx(1)=500.0
dx(2)=500.0
dx(3)=0.01*3600.0
dx(4)=0.01*3600.0
dx(5)=100.0
dx(6)=500.0
dx(7)=500.0
dx(8)=0.01*3600.0
dx(9)=0.01*3600.0
dx(10)=100.0
dx(11)=500.0
dx(12)=500.0
dx(13)=0.02*3600.0
dx(14)=0.02*3600.0
dx(15:34)=0.5
dx(35:54)=0.25

worsterr=0.0


for ff=0L,35-1 do begin
   obsid='9999'
   if ff lt 10 then filename='opd_'+obsid+'_'+string(ff,format='(I1)')+'.fits.gz'
   if ff ge 10 then filename='opd_'+obsid+'_'+string(ff,format='(I2)')+'.fits.gz'
   ss='python fit_lsst_opd_annuZ.py '+filename+' > temp'
   spawn,ss
   readcol,'temp',aaa,bbb,/silent
   for iii=0L,N_elements(aaa)-1 do senspn(aaa(iii),ff)=bbb(iii)
endfor


toterr=0L
tottest=0L

for dof=0L,55-1 do begin

!p.multi=[0,2,2]

for ff=0L,35-1 do begin
   if dof lt 10 then obsid='90'+string(dof,format='(I1)')+'1'
   if dof ge 10 then obsid='9'+string(dof,format='(I2)')+'1'
   if ff lt 10 then filename='opd_'+obsid+'_'+string(ff,format='(I1)')+'.fits.gz'
   if ff ge 10 then filename='opd_'+obsid+'_'+string(ff,format='(I2)')+'.fits.gz'
   ss='python fit_lsst_opd_annuZ.py '+filename+' > temp'
   spawn,ss
   readcol,'temp',aaa,bbb,/silent
   for iii=0L,N_elements(aaa)-1 do sensp(dof,aaa(iii),ff)=(bbb(iii)-senspn(aaa(iii),ff))/dx(dof)
endfor


maxy=max(abs([sens(dof,*,0)]))
maxw=max(abs([sensp(dof,*,0)]))
if maxw gt maxy then maxy=maxw
;if maxy lt 1.0 then maxy=1.0
maxy=maxy*1.1
maxz=0.0
ss='Sensitivity for '+names(dof)
plot,[0,23],[-maxy,maxy],/nodata,/xstyle,/ystyle,xtitle='Zernike Number',ytitle=ss
nerr=tolerance*sqrt(2.0)*1e-3/dx(dof)
errplot,findgen(22)+1,-nerr+maxy/2.0,nerr+maxy/2.0,color=40

for ff=0L,35-1 do begin
for zz=0L,22-1 do begin
   cc=0L
   if (dof ge 5 and zz ge 3) then begin
      sigma=abs(sensp(dof,zz,ff)-sens(dof,zz,ff))/nerr(zz)
      tottest=tottest+1
      if sigma gt 1 then begin
         toterr=toterr+1
         print,toterr,tottest,dof,zz,ff,sens(dof,zz,ff),sensp(dof,zz,ff),sigma*nerr(zz)*dx(dof)/1e-3
         if sigma*nerr(zz)*dx(dof)/1e-3 gt worsterr then worsterr=sigma*nerr(zz)*dx(dof)/1e-3
      endif
      if sigma gt 1 then cc=240 else cc=0
      if ff eq 0 then oplot,[zz+1],[sens(dof,zz,ff)],psym=8,color=cc,symsize=0.75
   endif
   if ff eq 0 then oplot,[zz+1],[sensp(dof,zz,ff)],psym=7,color=cc,symsize=0.375
   if (abs(sens(dof,zz,ff)-sens(dof,zz,0)) gt maxz) then maxz=abs(sens(dof,zz,ff)-sens(dof,zz,0))
   if (abs(sensp(dof,zz,ff)-sensp(dof,zz,0)) gt maxz) then maxz=abs(sensp(dof,zz,ff)-sensp(dof,zz,0))
endfor

endfor

;if maxz lt 1.0 then maxz=1.0
maxz=maxz*1.1
plot,[0,23],[-maxz,maxz],/nodata,/xstyle,/ystyle,xtitle='Zernike Number',ytitle='Relative Sensitivity to On-axis'
nerr=tolerance*sqrt(2.0)*1e-3/dx(dof)
errplot,findgen(22)+1,-nerr+maxz/2.0,nerr+maxz/2.0,color=40

for ff=0L,35-1 do begin
for zz=0L,22-1 do begin
   cc=0L
   if (dof ge 5 and zz ge 3) then begin
      sigma=abs(sensp(dof,zz,ff)-sens(dof,zz,ff))/nerr(zz)
      if sigma gt 1 then cc=240 else cc=0
      oplot,[zz+1],[sens(dof,zz,ff)-sens(dof,zz,0)],psym=8,color=cc,symsize=0.75
   endif
   oplot,[zz+1],[sensp(dof,zz,ff)-sensp(dof,zz,0)],psym=7,color=cc,symsize=0.375
endfor

endfor

!p.multi=[3,3,2]


;!p.multi=[2,2,2]

   if dof lt 10 then obsid='90'+string(dof,format='(I1)')+'0'
   if dof ge 10 then obsid='9'+string(dof,format='(I2)')+'0'
   filename='opd_'+obsid+'_'+'0.fits.gz'
   data=mrdfits(filename,/silent)
   maxz=max([abs(data)])*1.01
   ll=-maxz+findgen(256)/255.*(maxz*2.0)
   contour,data,/fill,nlevels=256,levels=ll,/xstyle,/ystyle

   if dof lt 10 then obsid='90'+string(dof,format='(I1)')+'1'
   if dof ge 10 then obsid='9'+string(dof,format='(I2)')+'1'
   filename='opd_'+obsid+'_'+'0.fits.gz'
   data=mrdfits(filename,/silent)
   maxz=max([abs(data)])*1.1
   ll=-maxz+findgen(256)/255.*(maxz*2.0)
   contour,data,/fill,nlevels=32,levels=ll,/xstyle,/ystyle

;   if dof lt 10 then obsid='90'+string(dof,format='(I1)')+'0'
;   if dof ge 10 then obsid='9'+string(dof,format='(I2)')+'0'
;   filename='opd_'+obsid+'_'+'34.fits.gz'
;   data=mrdfits(filename,/silent)
;   maxz=max([abs(data)])*1.01
;   ll=-maxz+findgen(256)/255.*(maxz*2.0)
;   contour,data,/fill,nlevels=256,levels=ll,/xstyle,/ystyle

   if dof lt 10 then obsid='90'+string(dof,format='(I1)')+'1'
   if dof ge 10 then obsid='9'+string(dof,format='(I2)')+'1'
   filename='opd_'+obsid+'_'+'34.fits.gz'
   data=mrdfits(filename,/silent)
   maxz=max([abs(data)])*1.1
   ll=-maxz+findgen(256)/255.*(maxz*2.0)
   contour,data,/fill,nlevels=32,levels=ll,/xstyle,/ystyle


endfor

!p.multi=[0,2,1]
loadct,1,/silent


z1=dblarr(220,550)
z2=dblarr(220,550)
z3=dblarr(22,55)
for ff=0,35-1 do begin
   for zz=0,22-1 do begin
      for dof=0,55-1 do begin
         for i=0,9 do begin
            for j=0,9 do begin
               z1(zz*10+i,dof*10+j)=z1(zz*10+i,dof*10+j)+(sens(dof,zz,ff))^2
               z2(zz*10+i,dof*10+j)=z2(zz*10+i,dof*10+j)+(sensp(dof,zz,ff))^2
               z3(zz,dof)=z3(zz,dof)+(sensp(dof,zz,ff))^2
            endfor
         endfor

      endfor
   endfor
endfor

z1=sqrt(z1)
z2=sqrt(z2)
z3=sqrt(z3)
z4=z1
z5=z2
for j=0L,550-1 do begin
max4=1e-30
max5=1e-30
for i=0L,220-1 do begin
   if z1(i,j) gt max4 then max4=z1(i,j)
   if z2(i,j) gt max5 then max5=z2(i,j)
endfor
for i=0L,220-1 do begin
   z4(i,j)=z1(i,j)/max4
   z5(i,j)=z2(i,j)/max5
endfor
endfor


   maxz=max([z1,z2])
   g1=where(z1 ne 0)
   g2=where(z2 ne 0)
   minz=min([z1(g1),z2(g2)])
   minz=maxz/1e4
   z1=256-(alog(z1>minz)-alog(minz))/(alog(maxz)-alog(minz))*256
   z2=256-(alog(z2>minz)-alog(minz))/(alog(maxz)-alog(minz))*256
   x1=findgen(220)/10.
   y1=findgen(550)/10.
   ll=findgen(256)
   contour,z1,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle
   contour,z2,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle


   loadct,3,/silent
   maxz=max([z4,z5])
   minz=min([z4,z5])
   z4=256-(z4-minz)/(maxz-minz)*256
   z5=256-(z5-minz)/(maxz-minz)*256
   x1=findgen(220)/10.
   y1=findgen(550)/10.
   ll=findgen(256)
   contour,z4,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle
   contour,z5,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle



r=byte([findgen(128)/127.*255.,fltarr(128)+255.0])
g=byte([findgen(96)/95.*255.,255.0+fltarr(64),255.0-findgen(96)/95.*255.])
b=byte([0.0,fltarr(127)+255.0,255.0-findgen(128)/127.*255.])
tvlct,r,g,b


for iii=0,1 do begin
z1=dblarr(220,550)
z2=dblarr(220,550)
for ff=0,35-1 do begin
   for zz=0,22-1 do begin
      for dof=0,55-1 do begin
         for i=0,9 do begin
            for j=0,9 do begin
               z1(zz*10+i,dof*10+j)=z1(zz*10+i,dof*10+j)+(sensp(dof,zz,31+iii*2)-sensp(dof,zz,0))
               z2(zz*10+i,dof*10+j)=z2(zz*10+i,dof*10+j)+(sensp(dof,zz,32+iii*2)-sensp(dof,zz,0))
            endfor
         endfor

      endfor
   endfor
endfor
for j=0L,550-1 do begin
max4=0.0
max5=0.0
for i=0L,220-1 do begin
   if abs(z1(i,j)) gt max4 then max4=abs(z1(i,j))
   if abs(z2(i,j)) gt max5 then max5=abs(z2(i,j))
endfor
for i=0L,220-1 do begin
   z4(i,j)=z1(i,j)/max4
   z5(i,j)=z2(i,j)/max5
endfor
endfor
   maxz=max([z4,z5])
   minz=min([z4,z5])
   z4=(z4-minz)/(maxz-minz)*256
   z5=(z5-minz)/(maxz-minz)*256
   x1=findgen(220)/10.
   y1=findgen(550)/10.
   ll=findgen(256)
   contour,z4,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle
   contour,z5,x1,y1,nlevels=256,levels=ll,/cell_fill,/xstyle,/ystyle
endfor



openw,1,'sensitivity_matrix_phosim.txt'
for dof=0,55-1 do begin
   for zz=0,22-1 do begin
      for ff=0,35-1 do begin
         printf,1,dof,zz,ff,sensp(dof,zz,ff)
      endfor
   endfor
endfor

close,1

print,'Worst error = ',worsterr,' nm'
print,'Errors = ',toterr,' out of ',tottest,' matrix elements'

tolerance_low(nnn,0)=0.0
tolerance_high(nnn,0)=1.0
value(nnn,0)=toterr
name(nnn,0)='WF Sensitivity'
unit(nnn,0)=' matrix elements'
comparison(nnn,0)='ZEMAX Code'
task(nnn,0)='2H OPD Sensitivity'
!x.margin=xmargo
!y.margin=ymargo



END
