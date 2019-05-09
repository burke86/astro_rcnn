;;
;; @package phosim
;; @file validation_4B.pro
;; @brief validation task 4B
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

pro validation_4B,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 4B'
  dytr=fltarr(6,2)

  data=mrdfits('lsst_e_4100_f0_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4101_f0_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(0,0)=approx
  dytr(0,1)=exact

  data=mrdfits('lsst_e_4110_f1_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4111_f1_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(1,0)=approx
  dytr(1,1)=exact

  data=mrdfits('lsst_e_4120_f2_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4121_f2_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(2,0)=approx
  dytr(2,1)=exact

  data=mrdfits('lsst_e_4130_f3_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4131_f3_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(3,0)=approx
  dytr(3,1)=exact

  data=mrdfits('lsst_e_4140_f4_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4141_f4_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(4,0)=approx
  dytr(4,1)=exact

  data=mrdfits('lsst_e_4150_f5_R22_S11_E000.fits.gz',/silent)
  approx=total(data)
  data=mrdfits('lsst_e_4151_f5_R22_S11_E000.fits.gz',/silent)
  exact=total(data)
;print,approx,exact,abs(approx-exact)/exact,'+/-',sqrt(approx+exact)/exact

  dytr(5,0)=approx
  dytr(5,1)=exact

  !p.multi=[0,1,2]
  plot,[-0.5,5.5],[0,1e6],xtitle='Filter',ytitle='Photon Counts',/nodata,/xstyle,/ystyle
  oplot,findgen(6)-0.02,dytr(*,0),psym=4,color=40
  errplot,findgen(6)-0.02,dytr(*,0)-sqrt(dytr(*,0)),dytr(*,0)+sqrt(dytr(*,0)),width=0.0
  oplot,findgen(6)+0.02,dytr(*,1),psym=4,color=250
  errplot,findgen(6)+0.02,dytr(*,0)-sqrt(dytr(*,0)),dytr(*,0)+sqrt(dytr(*,0)),width=0.0

  legend,psym=[4,4],['With dynamic transmission optimization','Exact calculation'],color=[40,250]

  plot,[-0.5,5.5],[-1,1],xtitle='Filter',ytitle='Percent Error',/nodata,/xstyle,/ystyle
  oplot,findgen(6),(dytr(*,0)-dytr(*,1))/dytr(*,1)*100,psym=4,color=40
  errplot,findgen(6),(dytr(*,0)-dytr(*,1)-sqrt(dytr(*,0)+dytr(*,1)))/dytr(*,1)*100.0,(dytr(*,0)-dytr(*,1)+sqrt(dytr(*,0)+dytr(*,1)))/dytr(*,1)*100.0,width=0.0
  oplot,findgen(1000)/1000.*6.0-0.5,fltarr(1000),linestyle=1

  value(nnn,0)=abs(total(dytr(*,0))-total(dytr(*,1)))/total(dytr(*,0))*100.0
  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=sqrt(total(dytr(*,0))+total(dytr(*,1)))/total(dytr(*,0))*100.0*2.0

  ss='Dynamic Transmission Optimization Test'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4B; '+vers
  xyouts,0.7,0.98,ss,/normal

  name(nnn,0)='Color Flux Accuracy'
  unit(nnn,0)='%'
  comparison(nnn,0)='Exact Calculation (2!4r!3)'
  task(nnn,0)='4B Colors w/ Dyn Trans on/off'

END
