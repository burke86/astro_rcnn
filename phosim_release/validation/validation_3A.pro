;;
;; @package phosim
;; @file validation_3A.pro
;; @brief validation task 3A
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

pro validation_3A,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3A'
  !p.multi=[0,2,2]
  loadct,5,/silent


  ;BNL flat numbers
  tolerance_low(nnn,0)=2.9*0.0
  tolerance_high(nnn,0)=2.9*2.0
  tolerance_low(nnn,1)=13.1*0.0
  tolerance_high(nnn,1)=13.1*2.0
  tolerance_low(nnn,2)=0.8*0.0
  tolerance_high(nnn,2)=0.8*2.0
  tolerance_low(nnn,3)=11.7*0.0
  tolerance_high(nnn,3)=11.7*2.0
  tolerance_low(nnn,4)=1.1*0.0
  tolerance_high(nnn,4)=1.1*2.0
  tolerance_low(nnn,5)=12.0*0.0
  tolerance_high(nnn,5)=12.0*2.0
  tolerance_low(nnn,6)=1.2*0.0
  tolerance_high(nnn,6)=1.2*2.0
  tolerance_low(nnn,7)=9.9*0.0
  tolerance_high(nnn,7)=9.9*2.0


  for i=0L,3 do begin

     if i eq 0 then wave=350.0
     if i eq 1 then wave=650.0
     if i eq 2 then wave=950.0
     if i eq 3 then wave=980.0

     if i eq 0 then file='3000_f0_R22_S11'
     if i eq 1 then file='3001_f2_R22_S11'
     if i eq 2 then file='3002_f5_R22_S11'
     if i eq 3 then file='3003_f5_R22_S11'
     filename='lsst_e_'+file+'_E000.fits.gz'

     data=mrdfits(filename,0,header,/silent)



     md=median(data)
     data2=(data-md)/md
     result=moment(data2)
     g=where(abs(data2) lt 3*sqrt(result(1)))
     clippedresult=moment(data2(g))
     noise=sqrt(md)/md
     noise2=sqrt(md)/(1.7*1000.0+md)

     ll=md+md*(findgen(11)-5.0)/40.0
     loadct,5,/silent
     contour,congrid(smooth(data,5),1024,1024),nlevels=11,/xstyle,/ystyle,/fill,levels=ll
     loadct,39,/silent
     ss=string(wave,format='(I3)')+' nm Monochromatic Flat Before Readout'
     xyouts,0.12,0.9,ss,/normal,color=255

     q=histogram(alog10(data>1),min=0.0,max=5.0,bin=0.1)
     xx=10.^(findgen(N_elements(q))*0.1)

     plot,xx,q,psym=10,/ylog,yr=[1,max(q)],/xlog,xr=[min(xx),max(xx)],xtitle='Electrons',ytitle='Electrons/pixel',/xstyle

     value(nnn,i*2)=sqrt(clippedresult(1)-noise^2)*100.0

     name(nnn,i*2)='PRNU at '+string(wave,format='(I3)')+' nm'
     unit(nnn,i*2)='%'
     comparison(nnn,i*2)='Prototype Devices at BNL'

     ss='PRNU at '+string(wave,format='(I3)')+' nm = '+$
        string(sqrt(result(1)-noise^2)*100.0,format='(F5.1)')+$
        '% / '+string(sqrt(clippedresult(1)-noise^2)*100.0,format='(F5.1)')+'% 3!4r!3 clipped'
     xyouts,2.0,1e7,ss

     ss='BNL measurements:  '+string(0.5*(tolerance_low(nnn,i*2)+tolerance_high(nnn,i*2)),format='(F5.1)')+'% 3!4r!3 clipped'

     xyouts,2.0,1e6,ss


     for j=0L,1 do begin
        for k=0L,7 do begin
           filename='lsst_a_'+file+'_C'+string(j,format='(I1)')+string(k,format='(I1)')+'_E000.fits.gz'
           data=mrdfits(filename,0,/silent)+2147483647L
           if j eq 0 and k eq 0 then begin
              ssize=size(data)
              datatotal=fltarr(ssize(1)*8,ssize(2)*2)
           endif
           for l=0L,ssize(1)-1 do begin
              for m=0L,ssize(2)-1 do begin
                 datatotal(l+k*ssize(1),m+j*ssize(2))=data(l,m)
              endfor
           endfor
        endfor
     endfor

     md=median(datatotal)
     data2=(datatotal-md)/md
     result=moment(data2)

     ll=md+md*(findgen(11)-5.0)/40.0
     loadct,5,/silent
     contour,congrid(smooth(datatotal,5),1024,1024),nlevels=11,/xstyle,/ystyle,/fill,levels=ll
     loadct,39,/silent
     ss=string(wave,format='(I3)')+' nm Monochromatic Flat After Readout'
     xyouts,0.12,0.45,ss,/normal,color=255

     q=histogram(alog10(datatotal>1),min=0.0,max=5.0,bin=0.1)
     xx=10.^(findgen(N_elements(q))*0.1)

     plot,xx,q,psym=10,/ylog,yr=[1,max(q)],/xlog,xr=[min(xx),max(xx)],xtitle='ADC',ytitle='ADC/pixel',/xstyle

     value(nnn,i*2+1)=sqrt(result(1)-noise2^2)*100.0

     name(nnn,i*2+1)='ERNU at '+string(wave,format='(I3)')+' nm'
     unit(nnn,i*2+1)='%'
     comparison(nnn,i*2+1)='Prototype Devices at BNL'

    ss='ERNU at '+string(wave,format='(I3)')+$
       ' nm = '+string(sqrt(result(1)-noise2^2)*100.0,format='(F5.1)')+'%'
     xyouts,2.0,1e7,ss

     ss='BNL measurements:  '+string(0.5*(tolerance_low(nnn,i*2+1)+tolerance_high(nnn,i*2+1)),format='(F5.1)')+'% '
     xyouts,2.0,1e6,ss


     ss='Detector Defects'
     xyouts,0.1,0.98,ss,/normal
     ss='Validation Task 3A; '+vers
     xyouts,0.7,0.98,ss,/normal

   endfor

     task(nnn,0)='3A Detector Defects'
     loadct,39,/silent



END
