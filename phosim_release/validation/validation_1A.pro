;;
;; @package phosim
;; @file validation_1A.pro
;; @brief validation task 1A
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

pro validation_1A,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 1A'

  !p.multi=[0,3,2]


  for qqq=0,2 do begin

     e1arr=fltarr(6)
     e2arr=fltarr(6)
     medxarr=fltarr(6)
     medyarr=fltarr(6)
     sigarr=fltarr(6)
     e1arre=fltarr(6)
     e2arre=fltarr(6)
     medxarre=fltarr(6)
     medyarre=fltarr(6)
     sigarre=fltarr(6)

     for nnq=0,11 do begin

        if nnq eq 0 then nn=0
        if nnq eq 1 then nn=1
        if nnq eq 2 then nn=0
        if nnq eq 3 then nn=1
        if nnq eq 4 then nn=2
        if nnq eq 5 then nn=3
        if nnq eq 6 then nn=2
        if nnq eq 7 then nn=3
        if nnq eq 8 then nn=4
        if nnq eq 9 then nn=5
        if nnq eq 10 then nn=4
        if nnq eq 11 then nn=5
        filename='lsst_e_10'+string(nn,format='(I1)')+string(qqq,format='(I1)')+'_f2_R22_S11_E000.fits.gz'
        data=mrdfits(filename,0,/silent)
        data=data(1940:2059,1976:2095)


        if ((nnq mod 4 eq 2) or (nnq mod 4 eq 3)) then begin
           data2=rebin(data,12,12)
           for iiq=0L,119 do begin
              for jjq=0L,119 do begin
                 data(iiq,jjq)=data2(floor(iiq/10),floor(jjq/10))
              endfor
           endfor
        endif
        measurepsf,data,rms,e1,e2,medx,medy,flux
        e1arr(nn)=e1
        e2arr(nn)=e2
        sigarr(nn)=rms*0.02*2.35
        medxarr(nn)=medx*0.02
        medyarr(nn)=medy*0.02
        e1arre(nn)=sqrt(2.0)/sqrt(total(data))
        e2arre(nn)=sqrt(2.0)/sqrt(total(data))
        sigarre(nn)=sqrt((rms*0.02*2.35/sqrt(total(data)))^2+0.02^2)
        medxarre(nn)=sqrt(2.0)*rms*0.02*2.35/sqrt(total(data))
        medyarre(nn)=sqrt(2.0)*rms*0.02*2.35/sqrt(total(data))
        sss=size(data)
        xx=(findgen(sss(1))-sss(1)/2.0)*1.0
        yy=(findgen(sss(2))-sss(2)/2.0)*1.0
        contour,data/mean(data),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,c_colors=findgen(11)*10+30,levels=(findgen(11)+1.0)*1.2

        if nn mod 2 eq 0 then begin
           xyouts,-57.0,52.0,'FWHM: '+string(rms*2.35*0.02,format='(F6.3)')+' +/-'+string(sqrt((sigarr(nn)/sqrt(total(data)))^2+0.02^2),format='(F6.3)')
           xyouts,-57.0,47.0,'e1: '+string(e1,format='(F6.3)')+' +/-'+string((1.0/sqrt(total(data))),format='(F6.3)')
           xyouts,-57.0,42.0,'e2: '+string(e2,format='(F6.3)')+' +/-'+string((1.0/sqrt(total(data))),format='(F6.3)')
           xyouts,-57.0,37.0,'x: '+string(medx*0.02,format='(F6.3)')+' +/-'+string((sigarr(nn)/sqrt(total(data))),format='(F6.3)')
           xyouts,-57.0,32.0,'y: '+string(medy*0.02,format='(F6.3)')+' +/-'+string((sigarr(nn)/sqrt(total(data))),format='(F6.3)')
           xyouts,-57,-50,'Geometric Approach'
           if nn eq 0 then xyouts,-57,-55,'Monochromatic SED; 1.5 ms exp'
           if nn eq 2 then xyouts,-57,-55,'Full SED; 1.5 ms exp'
           if nn eq 4 then xyouts,-57,-55,'Full SED; 15 s exp'
        endif
        if nn mod 2 eq 1 then begin
           xyouts,-57.0,52.0,'FWHM: '+string(rms*2.35*0.02,format='(F6.3)')+' +/-'+string(sqrt((sigarr(nn)/sqrt(total(data)))^2+0.02^2),format='(F6.3)')+' ('+string(sqrt((sigarr(nn)-sigarr(nn-1))^2/(sigarre(nn)^2+sigarre(nn-1)^2)),format='(F4.1)')+' !4r!3)'

           xyouts,-57.0,47.0,'e1: '+string(e1,format='(F6.3)')+' +/-'+string((1.0/sqrt(total(data))),format='(F6.3)')+' ('+string(sqrt((e1arr(nn)-e1arr(nn-1))^2/(e1arre(nn)^2+e1arre(nn-1)^2)),format='(F4.1)')+' !4r!3)'
           xyouts,-57.0,42.0,'e2: '+string(e2,format='(F6.3)')+' +/-'+string((1.0/sqrt(total(data))),format='(F6.3)')+' ('+string(sqrt((e2arr(nn)-e2arr(nn-1))^2/(e2arre(nn)^2+e2arre(nn-1)^2)),format='(F4.1)')+' !4r!3)'
           xyouts,-57.0,37.0,'x: '+string(medx*0.02,format='(F6.3)')+' +/-'+string((sigarr(nn)/sqrt(total(data))),format='(F6.3)')+' ('+string(sqrt((medxarr(nn)-medxarr(nn-1))^2/(medxarre(nn)^2+medxarre(nn-1)^2)),format='(F4.1)')+' !4r!3)'
           xyouts,-57.0,32.0,'y: '+string(medy*0.02,format='(F6.3)')+' +/-'+string((sigarr(nn)/sqrt(total(data))),format='(F6.3)')+' ('+string(sqrt((medyarr(nn)-medyarr(nn-1))^2/(medyarre(nn)^2+medyarre(nn-1)^2)),format='(F4.1)')+' !4r!3)'
           chi2=(sigarr(nn)-sigarr(nn-1))^2/(sigarre(nn)^2+sigarre(nn-1)^2)+$
                (e1arr(nn)-e1arr(nn-1))^2/(e1arre(nn)^2+e1arre(nn-1)^2)+$
                (e2arr(nn)-e2arr(nn-1))^2/(e2arre(nn)^2+e2arre(nn-1)^2)+$
                (medxarr(nn)-medxarr(nn-1))^2/(medxarre(nn)^2+medxarre(nn-1)^2)+$
                (medyarr(nn)-medyarr(nn-1))^2/(medyarre(nn)^2+medyarre(nn-1)^2)
           xyouts,25,-45,'!4V!3!D!4t!3!N!U2!N='+string(chi2/5.0,format='(F5.2)')


           if nnq eq 3 and qqq eq 0 then begin
              value(nnn,5)=abs(sigarr(nn)-sigarr(nn-1))
              tolerance_high(nnn,5)=sqrt(sigarre(nn)^2+sigarre(nn-1)^2)*3.0
              value(nnn,6)=sqrt(sqrt((e1arr(nn)-e1arr(nn-1))^2+(e2arr(nn)-e2arr(nn-1))^2))*value(nnn,5)
              tolerance_high(nnn,6)=sqrt(sqrt(e1arre(nn)^2+e1arre(nn-1)^2+e2arre(nn)^2+e2arre(nn-1)^2))*tolerance_high(nnn,5)
              value(nnn,7)=sqrt((medxarr(nn)-medxarr(nn-1))^2+(medyarr(nn)-medyarr(nn-1))^2)
              tolerance_high(nnn,7)=sqrt(medxarre(nn-1)^2+medxarre(nn)^2+medyarre(nn-1)^2+medyarre(nn)^2)*3.0
           endif
           if nnq eq 7 and qqq eq 0 then begin
              value(nnn,9)=abs(sigarr(nn)-sigarr(nn-1))
              tolerance_high(nnn,9)=sqrt(sigarre(nn)^2+sigarre(nn-1)^2)*3.0
              value(nnn,10)=sqrt(sqrt((e1arr(nn)-e1arr(nn-1))^2+(e2arr(nn)-e2arr(nn-1))^2))*value(nnn,9)
              tolerance_high(nnn,10)=sqrt(sqrt(e1arre(nn)^2+e1arre(nn-1)^2+e2arre(nn)^2+e2arre(nn-1)^2))*tolerance_high(nnn,9)
              value(nnn,11)=sqrt((medxarr(nn)-medxarr(nn-1))^2+(medyarr(nn)-medyarr(nn-1))^2)
              tolerance_high(nnn,11)=sqrt(medxarre(nn-1)^2+medxarre(nn)^2+medyarre(nn-1)^2+medyarre(nn)^2)*3.0
           endif
           if nnq eq 11 and qqq eq 0 then begin
              value(nnn,1)=abs(sigarr(nn)-sigarr(nn-1))
              tolerance_high(nnn,1)=sqrt(sigarre(nn)^2+sigarre(nn-1)^2)*3.0
              value(nnn,2)=sqrt(sqrt((e1arr(nn)-e1arr(nn-1))^2+(e2arr(nn)-e2arr(nn-1))^2))*value(nnn,1)
              tolerance_high(nnn,2)=sqrt(sqrt(e1arre(nn)^2+e1arre(nn-1)^2+e2arre(nn)^2+e2arre(nn-1)^2))*tolerance_high(nnn,1)
              value(nnn,3)=sqrt((medxarr(nn)-medxarr(nn-1))^2+(medyarr(nn)-medyarr(nn-1))^2)
              tolerance_high(nnn,3)=sqrt(medxarre(nn-1)^2+medxarre(nn)^2+medyarre(nn-1)^2+medyarre(nn)^2)*3.0
           endif


           xyouts,-57,-50,'Fourier Approach'
           if nn eq 1 then xyouts,-57,-55,'Monochromatic SED; 1.5 ms exp'
           if nn eq 3 then xyouts,-57,-55,'Full SED; 1.5 ms exp'
           if nn eq 5 then xyouts,-57,-55,'Full SED; 15 s exp'
           if ((nnq mod 4 eq 2) or (nnq mod 4 eq 3)) then begin
              data=data*100
              datapr=datapr*100
           endif
           resid=(datapr-data*total(datapr)/total(data))/sqrt(((datapr>2)+(sqrt(data>2)*total(datapr)/total(data))^2))
           if nnq mod 4 eq 1 then bad=where(datapr le 2 or data le 2)
           if nnq mod 4 eq 3 then bad=where(datapr le 20 or data le 20)
           if N_elements(bad) eq 1 then if bad(0) eq -1 then goto,skipb
           resid(bad)=0
skipb:
           contour,((resid>(-3))<3),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,levels=[-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5],c_colors=[20,40,60,255,200,220,240]
           if nnq mod 4 eq 1 then good=where(datapr gt 2 and data gt 2)
           if nnq mod 4 eq 3 then good=where(datapr gt 20 and data gt 20)
           if N_elements(good) eq 1 then if good(0) eq -1 then goto,skipg
           chi2=total(resid(good)^2)/float(N_elements(good))
skipg:
           xyouts,25,-45,'!4V!3!D!4t!3!N!U2!N='+string(chi2,format='(F5.2)')
           xyouts,-57,52,'Standard Deviations'

        endif

        datapr=data

        if nnq eq 3 or nnq eq 7 or nnq eq 11 then begin
           ss='Diffraction Approximation'
           if qqq eq 1 then ss=ss+' with Low Frequency Screen'
           if qqq eq 2 then ss=ss+' with High Frequency Screen'
           xyouts,0.1,0.98,ss,/normal
           ss='Validation Task 1A; '+vers
           xyouts,0.7,0.98,ss,/normal
        endif

        if nnq eq 3 and qqq eq 0 then value(nnn,4)=chi2
        if nnq eq 7 and qqq eq 0 then value(nnn,8)=chi2
        if nnq eq 11 and qqq eq 0 then value(nnn,0)=chi2
     endfor


  endfor

  tolerance_low(nnn,0)=0.
  tolerance_high(nnn,0)=2.
  tolerance_high(nnn,4)=4.0
  tolerance_high(nnn,8)=4.0


  name(nnn,0)='Integrated PSF Comparison'
  name(nnn,1)=' FWHM error'
  name(nnn,2)=' Sqrt(ellip)*FWHM error'
  name(nnn,3)=' Centroid error'
  name(nnn,4)='Instantaneous Monochr PSF Comparison'
  name(nnn,5)=' FWHM error'
  name(nnn,6)=' Sqrt(ellip)*FWHM error'
  name(nnn,7)=' Centroid error'
  name(nnn,8)='Instantaneous Full SED PSF Comparison'
  name(nnn,9)=' FWHM error'
  name(nnn,10)=' Sqrt(ellip)*FWHM error'
  name(nnn,11)=' Centroid error'

  unit(nnn,0)=' (!4V!3!U2!N/dof)'
  comparison(nnn,0)='PhoSim FFT Calculation'
  unit(nnn,1)=' arcsec'
  comparison(nnn,1)='PhoSim FFT Calculation'
  unit(nnn,2)=' arcsec'
  comparison(nnn,2)='PhoSim FFT Calculation'
  unit(nnn,3)=' arcsec'
  comparison(nnn,3)='PhoSim FFT Calculation'

  unit(nnn,4)=' (!4V!3!U2!N/dof)'
  comparison(nnn,4)='PhoSim FFT Calculation'
  unit(nnn,5)=' arcsec'
  comparison(nnn,5)='PhoSim FFT Calculation'
  unit(nnn,6)=' arcsec'
  comparison(nnn,6)='PhoSim FFT Calculation'
  unit(nnn,7)=' arcsec'
  comparison(nnn,7)='PhoSim FFT Calculation'

  unit(nnn,8)=' (!4V!3!U2!N/dof)'
  comparison(nnn,8)='PhoSim FFT Calculation'
  unit(nnn,9)=' arcsec'
  comparison(nnn,9)='PhoSim FFT Calculation'
  unit(nnn,10)=' arcsec'
  comparison(nnn,10)='PhoSim FFT Calculation'
  unit(nnn,11)=' arcsec'
  comparison(nnn,11)='PhoSim FFT Calculation'

  task(nnn,0)='1A Diffraction Approximation'

END
