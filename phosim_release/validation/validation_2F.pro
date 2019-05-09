;;
;; @package phosim
;; @file validation_2F.pro
;; @brief validation task 2F
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

pro validation_2F,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 2F'
  flux=fltarr(6,2)


  delta=0.0

  data=mrdfits('lsst_e_2500_f0_R22_S11_E000.fits.gz',/silent)
  flux(0,0)=total(data)
  flux(0,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(400.0-320.0+0.1+delta)/6.626e-27

  data=mrdfits('lsst_e_2501_f1_R22_S11_E000.fits.gz',/silent)
  flux(1,0)=total(data)
  flux(1,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(560.0-400.0+delta)/6.626e-27

  data=mrdfits('lsst_e_2502_f2_R22_S11_E000.fits.gz',/silent)
  flux(2,0)=total(data)
  flux(2,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(680.0-560.0+0.1+delta)/6.626e-27

  data=mrdfits('lsst_e_2503_f3_R22_S11_E000.fits.gz',/silent)
  flux(3,0)=total(data)
  flux(3,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(820.0-680.0+0.1+delta)/6.626e-27

  data=mrdfits('lsst_e_2504_f4_R22_S11_E000.fits.gz',/silent)
  flux(4,0)=total(data)
  flux(4,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(920.0-820.0+0.1+delta)/6.626e-27

  data=mrdfits('lsst_e_2505_f5_R22_S11_E000.fits.gz',/silent)
  flux(5,0)=total(data)
  flux(5,1)=10.^((16.+48.6)/(-2.5))*15.0*!Pi*(418.0^2-255.8^2)/500.0*(1020.0-920.0+0.1+delta)/6.626e-27

  !p.multi=[0,1,2]
  plot,[-0.5,5.5],[0,5e6],xtitle='Filter',ytitle='Photon Counts',/nodata,/xstyle,/ystyle
  oplot,findgen(6)-0.02,flux(*,0),psym=4,color=40
  errplot,findgen(6)-0.02,flux(*,0)-sqrt(flux(*,0)),flux(*,0)+sqrt(flux(*,0)),width=0.0
  oplot,findgen(6)+0.02,flux(*,1),psym=4,color=250
  errplot,findgen(6)+0.02,flux(*,0)-sqrt(flux(*,0)),flux(*,0)+sqrt(flux(*,0)),width=0.0

  legend,psym=[4,4],['PhoSim','Exact calculation'],color=[40,250]

  plot,[-0.5,5.5],[-2,2],xtitle='Filter',ytitle='Percent Error',/nodata,/xstyle,/ystyle
  oplot,findgen(6),(flux(*,0)-flux(*,1))/flux(*,1)*100,psym=4,color=40,symsize=0.5
  errplot,findgen(6),(flux(*,0)-flux(*,1)-sqrt(flux(*,0)))/flux(*,1)*100.0,(flux(*,0)-flux(*,1)+sqrt(flux(*,0)))/flux(*,1)*100.0,width=0.0
  oplot,findgen(1000)/1000.*6.0-0.5,fltarr(1000),linestyle=1

  value(nnn,0)=abs(total(flux(*,0))-total(flux(*,1)))/total(flux(*,1))*100.0
  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=sqrt(total(flux(*,0)))/total(flux(*,1))*100.0*3.0

  oplot,findgen(1000)/1000.*6.0-0.5,fltarr(1000)+tolerance_high(nnn,0),linestyle=1
  oplot,findgen(1000)/1000.*6.0-0.5,fltarr(1000)-tolerance_high(nnn,0),linestyle=1


  ss='Throughput Accuracy Test'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2F; '+vers
  xyouts,0.7,0.98,ss,/normal

  name(nnn,0)='Throughput Accuracy'
  unit(nnn,0)='%'
  comparison(nnn,0)='Theoretical Calculation (2!4r!3)'
  task(nnn,0)='2F Throughput Accuracy'

END
