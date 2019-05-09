;;
;; @package phosim
;; @file validation_2B.pro
;; @brief validation task 2B
;;
;; @brief Created by:
;; @author Nathan Todd (Purdue)
;;
;; @brief Modified by:
;; @author John R. Peterson (Purdue)
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

pro validation_2B,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison


  print,'Task 2B'
  !p.multi=[0,5,3]
  !x.margin=[0,0]
  !y.margin=[0,0]

loadct,5
  !EXCEPT=0
  sigma = fltarr(2,25,3)
  ellip1 = fltarr(2,25,3)
  ellip2 = fltarr(2,25,3)
  centx=fltarr(2,25,3)
  centy=fltarr(2,25,3)
  temperature=fltarr(75)
  tempder=fltarr(75)
  altitude=fltarr(75)

  cc=0L & ct=0L
  for n=0,1 do begin
     if n eq 0 then high=24
     if n eq 1 then high=4

     for m=0,high do begin

        for o=0,2 do begin

           if o eq 0 and cc lt 10 then temp = strcompress('lsst_e_210'+string(cc,format='(I1)')+'_f2_R22_S11_E000.fits.gz',/remove)
           if o eq 0 and cc ge 10 then temp = strcompress('lsst_e_21'+string(cc,format='(I2)')+'_f2_R22_S11_E000.fits.gz',/remove)
           if o eq 1 and cc lt 10 then temp = strcompress('lsst_e_210'+string(cc,format='(I1)')+'_f2_R32_S21_E000.fits.gz',/remove)
           if o eq 1 and cc ge 10 then temp = strcompress('lsst_e_21'+string(cc,format='(I2)')+'_f2_R32_S21_E000.fits.gz',/remove)
           if o eq 2 and cc lt 10 then temp = strcompress('lsst_e_210'+string(cc,format='(I1)')+'_f2_R42_S21_E000.fits.gz',/remove)
           if o eq 2 and cc ge 10 then temp = strcompress('lsst_e_21'+string(cc,format='(I2)')+'_f2_R42_S21_E000.fits.gz',/remove)

           FITS_READ, temp, data, header

           measurepsf,data,rms,e1,e2,medx,medy,flux

           sigma(n,m,o)=rms
           ellip1(n,m,o)=e1
           ellip2(n,m,o)=e2
           centx(n,m,o)=medx
           centy(n,m,o)=medy

           temp=0.0 & tempvar=0.0 & zenith=0.0
           for tt=0,N_elements(header)-1 do begin
              if (strmid(header(tt),0,7) eq 'TEMPERA') then temp=strmid(header(tt),10,21)
              if (strmid(header(tt),0,7) eq 'TEMPVAR') then tempvar=strmid(header(tt),10,21)
              if (strmid(header(tt),0,7) eq 'ZENITH ') then zenith=strmid(header(tt),10,21)
           endfor
           if n eq 0 then temperature(ct)=temp
           if n eq 0 then tempder(ct)=tempvar
           if n eq 0 then altitude(ct)=90.0-zenith

           ct=ct+1L
        endfor
           cc=cc+1L
     endfor
  endfor


;make contour plots

  temp = strarr(90)
  tempb = strarr(90)
  tempc = strarr(90)
  ii=0
  for n=0,5 do begin
  for m=n*5+0,n*5+4 do begin
     if m lt 10 then temp(ii) = strcompress('lsst_e_210'+string(m,format='(I1)')+'_f2_R22_S11_E000.fits.gz',/remove)
     if m ge 10 then temp(ii) = strcompress('lsst_e_21'+string(m,format='(I2)')+'_f2_R22_S11_E000.fits.gz',/remove)
     tempb(ii)=m
     if n eq 5 then tempb(ii)=m-25
     tempc(ii)=0
     ii=ii+1
  endfor
  for m=n*5+0,n*5+4 do begin
     if m lt 10 then temp(ii) = strcompress('lsst_e_210'+string(m,format='(I1)')+'_f2_R32_S21_E000.fits.gz',/remove)
     if m ge 10 then temp(ii) = strcompress('lsst_e_21'+string(m,format='(I2)')+'_f2_R32_S21_E000.fits.gz',/remove)
     tempb(ii)=m
     if n eq 5 then tempb(ii)=m-25
     tempc(ii)=1
     ii=ii+1
  endfor
  for m=n*5+0,n*5+4 do begin
     if m lt 10 then temp(ii) = strcompress('lsst_e_210'+string(m,format='(I1)')+'_f2_R42_S21_E000.fits.gz',/remove)
     if m ge 10 then temp(ii) = strcompress('lsst_e_21'+string(m,format='(I2)')+'_f2_R42_S21_E000.fits.gz',/remove)
     tempb(ii)=m
     if n eq 5 then tempb(ii)=m-25
     tempc(ii)=2
     ii=ii+1
  endfor
endfor


  for k=0L,89 do begin
     data=mrdfits(temp(k),0,/silent)

     if k ge 74 then a=1 else a=0
     b=tempb(k)
     c=tempc(k)
     psf=sigma(a,b,c)*2*0.01*sqrt(2*alog(2))
     if k lt 74 then ss='aber+pert    '
     if k eq 75 or k eq 80 or k eq 85 then ss='aber+track  '
     if k eq 76 or k eq 81 or k eq 86 then ss='aber+wd shk '
     if k eq 77 or k eq 82 or k eq 87 then ss='aber        '
     if k eq 78 or k eq 83 or k eq 88 then ss='aber+dm see '
     if k eq 79 or k eq 84 or k eq 89 then ss='aber+det    '
     mx=0L & my=0L & bx=50 & by=60
     mi=0L & for i=0L,3999 do for j=0L,4071 do if data(i,j) gt mi then begin & mi=data(i,j) & mx=i & my=j & endif
     image=data[(mx-bx):(mx+bx),(my-by):(my+by)]
     sss=size(image)
     xx=(findgen(sss[1])-sss[1]/2+(mx-2000))*0.5
     yy=(findgen(sss[2])-sss[2]/2+(my-2036))*0.5
     mxx=max(alog(image>1))
     ccol=findgen(100)/99.0*255.0
     ccol=reverse(ccol)
     contour,alog(image>1),xx,yy,/fill,/xstyle,/ystyle,nlevels=100,c_colors=ccol
     xyouts,-20.0+(mx-2000)*0.5,20.0+(my-2036)*0.5,ss+string(psf,format='(F4.2)')+'"',/data
endfor


  sigmaopt=sigma*0.5
  centxopt=centx*0.5
  centyopt=centy*0.5
  ellip1opt=ellip1
  ellip2opt=ellip2
  ss='Optics Perturbations & Tracking'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2B; '+vers
  xyouts,0.7,0.98,ss,/normal

  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=0.381

  tolerance_low(nnn,1)=0.0
  tolerance_low(nnn,2)=0.0
  tolerance_low(nnn,3)=0.0
  tolerance_low(nnn,4)=0.0
  tolerance_high(nnn,1)=0.097
  tolerance_high(nnn,2)=0.090
  tolerance_high(nnn,3)=0.259
  tolerance_high(nnn,4)=0.246

  tolerance_high(nnn,5)=0.21
  tolerance_high(nnn,6)=0.21
  tolerance_high(nnn,7)=0.21
  tolerance_high(nnn,8)=0.21
  tolerance_high(nnn,9)=0.21

  tolerance_high(nnn,10)=2.0
  tolerance_high(nnn,11)=2.0
  tolerance_high(nnn,12)=2.0
  tolerance_high(nnn,13)=2.0
  tolerance_high(nnn,14)=2.0

  name(nnn,0)='Total non-atm PSF size'
  name(nnn,1)='  Optics Design'
  name(nnn,2)='  Internal Seeing'
  name(nnn,3)='  Pert/Mis/Track/WdShk'
  name(nnn,4)='  Charge Diffusion'

  name(nnn,5)='PSF Sqrt(ellip)*FWHM'
  name(nnn,6)='  Optics Design'
  name(nnn,7)='  Internal Seeing'
  name(nnn,8)='  Pert/Mis/Track/WdShk'
  name(nnn,9)='  Charge Diffusion'

  name(nnn,10)='PSF Centroid'
  name(nnn,11)='  Optics Design'
  name(nnn,12)='  Internal Seeing'
  name(nnn,13)='  Pert/Mis/Track/WdShk'
  name(nnn,14)='  Charge Diffusion'

  task(nnn,0)='2B Instrumental PSF Properties'

  unit(nnn,0)=' arcsec'
  comparison(nnn,0)='LSST Design'

  unit(nnn,1)=' arcsec'
  comparison(nnn,1)='LSST Design'
  unit(nnn,2)=' arcsec'
  comparison(nnn,2)='LSST Design'
  unit(nnn,3)=' arcsec'
  comparison(nnn,3)='LSST Design'
  unit(nnn,4)=' arcsec'
  comparison(nnn,4)='LSST Design'

  unit(nnn,5)=' arcsec'
  comparison(nnn,5)='LSST Design'
  unit(nnn,6)=' arcsec'
  comparison(nnn,6)='LSST Design'
  unit(nnn,7)=' arcsec'
  comparison(nnn,7)='LSST Design'
  unit(nnn,8)=' arcsec'
  comparison(nnn,8)='LSST Design'
  unit(nnn,9)=' arcsec'
  comparison(nnn,9)='LSST Design'
  unit(nnn,10)=' arcsec'
  comparison(nnn,10)='LSST Design'
  unit(nnn,11)=' arcsec'
  comparison(nnn,11)='LSST Design'
  unit(nnn,12)=' arcsec'
  comparison(nnn,12)='LSST Design'
  unit(nnn,13)=' arcsec'
  comparison(nnn,13)='LSST Design'
  unit(nnn,14)=' arcsec'
  comparison(nnn,14)='LSST Design'

  value(nnn,0)=sqrt((mean(sigmaopt(0,*,*))*2.35)^2+$
  (mean(sigmaopt(1,0,*))*2.35)^2+$
  (mean(sigmaopt(1,1,*))*2.35)^2+$
  (mean(sigmaopt(1,3,*))*2.35)^2+$
  (mean(sigmaopt(1,4,*))*2.35)^2-$
  5.0*(mean(sigmaopt(1,2,*))*2.35)^2)/50.0
  value(nnn,1)=mean(sigmaopt(1,2,*))*2.35/50.0
  value(nnn,2)=sqrt((mean(sigmaopt(1,3,*))*2.35)^2-(mean(sigmaopt(1,2,*))*2.35)^2)/50.0
  value(nnn,3)=sqrt((mean(sigmaopt(0,*,*))*2.35)^2+(mean(sigmaopt(1,0,*))*2.35)^2+(mean(sigmaopt(1,1,*))*2.35)^2-3.0*(mean(sigmaopt(1,2,*))*2.35)^2)/50.0
  value(nnn,4)=sqrt((mean(sigmaopt(1,4,*))*2.35)^2-(mean(sigmaopt(1,2,*))*2.35)^2)/50.0

  value(nnn,6)=mean(((ellip1opt(1,2,*)^2+ellip2opt(1,2,*)^2))^(0.25))*value(nnn,1)
  value(nnn,7)=mean(((ellip1opt(1,3,*)^2+ellip2opt(1,3,*)^2))^(0.25))*value(nnn,2)
  value(nnn,8)=mean(((ellip1opt(0,*,*)^2+ellip2opt(0,*,*)^2))^(0.25))*value(nnn,3)
  value(nnn,9)=mean(((ellip1opt(1,4,*)^2+ellip2opt(1,4,*)^2))^(0.25))*value(nnn,4)
  value(nnn,5)=sqrt(total((value(nnn,6:9))^2))

  medx0=centxopt(1,2,0) & medy0=centyopt(1,2,0)
  medx1=centxopt(1,2,1) & medy1=centyopt(1,2,1)
  medx2=centxopt(1,2,2) & medy2=centyopt(1,2,2)
  value(nnn,11)=0.0
  value(nnn,12)=0.33*sqrt(mean((centxopt(1,3,0)-medx0)^2+(centyopt(1,3,0)-medy0)^2))/50.0+$
                0.33*sqrt(mean((centxopt(1,3,1)-medx1)^2+(centyopt(1,3,1)-medy1)^2))/50.0+$
                0.33*sqrt(mean((centxopt(1,3,2)-medx2)^2+(centyopt(1,3,2)-medy2)^2))/50.0
  value(nnn,13)=0.33*sqrt(mean((centxopt(0,*,0)-medx0)^2+(centyopt(0,*,0)-medy0)^2))/50.0+$
                0.33*sqrt(mean((centxopt(0,*,1)-medx1)^2+(centyopt(0,*,1)-medy1)^2))/50.0+$
                0.33*sqrt(mean((centxopt(0,*,2)-medx2)^2+(centyopt(0,*,2)-medy2)^2))/50.0
  value(nnn,14)=0.33*sqrt(mean((centxopt(1,4,0)-medx0)^2+(centyopt(1,4,0)-medy0)^2))/50.0+$
                0.33*sqrt(mean((centxopt(1,4,1)-medx1)^2+(centyopt(1,4,1)-medy1)^2))/50.0+$
                0.33*sqrt(mean((centxopt(1,4,2)-medx2)^2+(centyopt(1,4,2)-medy2)^2))/50.0
  value(nnn,10)=sqrt(total((value(nnn,11:14))^2))

  !p.multi=0
  !x.margin=[9,2]
  !y.margin=[3,1]

  psf=dblarr(75)
  for i=0,24-1 do for j=0,2 do psf(i*3+j)=sigma(0,i,j)*2*0.01*sqrt(2*alog(2))
  plot,temperature,psf,psym=4,xtitle='Temperature',ytitle='PSF Size'
  result=regress(temperature,psf,con=con)
  xx=findgen(2000)/10.-100.0
  oplot,xx,result(0)*xx+con
  ss='Slope= '+string(result(0))
  xyouts,0.6,0.9,ss,/normal

  plot,tempder,psf,psym=4,xtitle='Temperature Derivative (degrees per hour)',ytitle='PSF Size'

  plot,altitude,psf,psym=4,xtitle='Altitude (degrees)',ytitle='PSF Size'


loadct,39
END
