;;
;; @package phosim
;; @file validation_2D.pro
;; @brief validation task 2D
;;
;; @brief Created by:
;; @author John R. Peterson (Purdue)
;;
;; @brief Modified by:
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

pro validation_2D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 2D'
  !p.multi=[0,2,1]

  data=mrdfits('lsst_e_2300_f2_R22_S11_E000.fits.gz',0,/silent)
  image=data(1950:2050,1966:2106)
  sss=size(image)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*10.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*10.0
  ll=findgen(11)/1.5
  cl=[255,findgen(10)/10.*255]
  contour,alog10(image>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=ll,c_colors=cl
  ss='Spider Diffraction Accuracy'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2D; '+vers
  xyouts,0.7,0.98,ss,/normal



;diffraction calculation
  z=dblarr(sss(1),sss(2))
  readcol,'diffraction_pattern.txt',aa,bb,cc,/silent
  for i=0L,N_elements(aa)-1 do begin
     if (aa(i)-50 ge 0 and aa(i)-50 lt sss(1) and $
         bb(i)-70 ge 0 and bb(i)-70 lt sss(2)) then begin
        z(aa(i)-50,bb(i)-70)=cc(i)
     endif
  endfor
  cimage=z/total(z)*total(image)

;200/1.00000  32.6,21.8,45.6
;200/0.50000  16.1,27.7,56.2
;200/0.25000  71.5,11.4.17.1
;200/0.12500  89.6, 4.1, 6.3
;200/0.06250  92.2, 3.1, 4.7
;200/0.03125  92.2, 3.1, 4.7

;400/1.00000  60.1, 2.6,37.3
;400/0.50000  34.8, 8.4,56.8
;400/0.25000  87.9, 1.7,10.4
;400/0.12500  95.9, 0.6, 3.5
;400/0.06250  (called 5)      97.1, 0.4, 2.5
;400/0.03125

;800/1.00000  60.0, 2.6,37.3
;800/0.50000  34.8, 8.3,56.9
;800/0.25000  87.9, 1.7,10.4
;800/0.12500
;800/0.06250
;800/0.03125

;answer:      95.6, 0.9, 3.5

  range=50.0
  ranges=20.0

  phi=atan(yy,xx)
  rad=sqrt(xx*xx+yy*yy)

  cent=0. & centm=0. & spike=0. & spikem=0.
  for i=0L,sss(1)-1 do begin
     for j=0L,sss(2)-1 do begin
        rad=sqrt(xx(i)*xx(i)+yy(j)*yy(j))
        phi=atan(yy(j),xx(i))

        if rad lt range then begin
           cent=cent+image(i,j)
           centm=centm+cimage(i,j)
        endif

        if ((abs(xx(i)*sin(!Pi/4.)+yy(j)*cos(!Pi/4.)) lt ranges or $
             abs(xx(i)*sin(3.0*!Pi/4.)+yy(j)*cos(3.0*!Pi/4.)) lt ranges) and $
            rad ge range) then begin
           spike=spike+image(i,j)
           spikem=spikem+cimage(i,j)
        endif



     endfor
  endfor




  contour,alog10((round(cimage>1))),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=ll,c_colors=cl
  xyouts,0.11,0.83,'PhoSim Diffraction with flat SED for r band',/normal
  xyouts,0.11,0.8,'Core   '+string(cent/total(image)*100.0,format='(F4.1)')+'%',/normal
  xyouts,0.11,0.77,'Disk   '+string((total(image)-cent-spike)/total(image)*100.0,format='(F4.1)')+'%',/normal
  xyouts,0.11,0.74,'Spikes '+string(spike/total(image)*100.0,format='(F4.1)')+'%',/normal
  xyouts,0.59,0.83,'Fraunhofer Integral at 0.57, 0.62, & 0.67 microns',/normal
  xyouts,0.59,0.8,'Core   '+string(centm/total(cimage)*100.0,format='(F4.1)')+'%',/normal
  xyouts,0.59,0.77,'Disk   '+string((total(cimage)-centm-spikem)/total(cimage)*100.0,format='(F4.1)')+'%',/normal
  xyouts,0.59,0.74,'Spikes '+string(spikem/total(cimage)*100.0,format='(F4.1)')+'%',/normal
  value(nnn,0)=sqrt(1./3.*((cent/total(image)*100.0-centm/total(cimage)*100.0)^2+$
                           (spike/total(image)*100.0-spikem/total(cimage)*100.0)^2+$
                           ((total(image)-spike-cent)/total(image)*100.0-(total(cimage)-spikem-centm)/total(cimage)*100.0)^2))

  nr=floor(sqrt(sss(1)*sss(1)+sss(2)*sss(2))/2.0/3.0)+2
  nphi=round(2*!Pi*100.0/3.0)+2.0

  tolerance_low(nnn,0)=0.
  tolerance_high(nnn,0)=0.5
  name(nnn,0)='Error in Power in wing, core, spikes'

  unit(nnn,0)='%'
  comparison(nnn,0)='Fraunhofer Calculation'
  task(nnn,0)='2D Spider Diffraction'

END
