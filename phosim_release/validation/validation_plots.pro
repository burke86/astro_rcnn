;;
;; @package phosim
;; @file validation_plots.pro
;; @brief validation plots
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

pro validation_plots,regression

  set_plot,'PS'
  device,filename='validation.ps'
  device,xsize=10.0,/inches
  device,ysize=7.0,/inches
  device,xoffset=0.0,/inches
  device,yoffset=0.0,/inches
  device,/color
  !x.thick=4
  !y.thick=4
  !p.charsize=1.0
  !p.charthick=4
  !p.thick=4
  !x.omargin=[2,2]
  !y.omargin=[2,2]
  xmargo=[9,2]
  ymargo=[3,1]
  !x.margin=xmargo
  !y.margin=ymargo
  loadct,39

  readcol,'../bin/version',a,b,format='(A,A)',/SILENT
  sss=strlen(b(0))
  vers='Commit: ..'+strmid(b(0),sss-7,sss-1)

  tn=30L
  value=dblarr(tn,tn)
  tolerance_low=dblarr(tn,tn)
  tolerance_high=dblarr(tn,tn)
  result_string=strarr(tn,tn)
  comparison=strarr(tn,tn)
  name=strarr(tn,tn)
  task=strarr(tn,tn)
  unit=strarr(tn,tn)

;VALIDATION TASK 1A:  DIFFRACTION APPROXIMATION
  num=0
  validation_1A,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 1B:  SCREEN CONVERGENCE
  num=num+1
  validation_1B,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 1C:  ATMOSPHERE PSF
;VALIDATION TASK 1E:  HIGH AIRMASS ATMOSPHERE PSF
;VALIDATION TASK 1F:  INSTRUMENT+ATMOSPHERE PSF
  num=num+1
  validation_1C,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 1D:  ATMOSPHERE ASTROMETRY
  num=num+1
  validation_1D,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 1G:   DCR
  num=num+1
  validation_1G,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 1H:   STRUCTURE FUNCTION
  num=num+1
  validation_1H,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2A:  SPOT DIAGRAM
  num=num+1
  validation_2A,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2B:  OPTICS & TRACKING
  num=num+1
  validation_2B,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2C:  DESIGN THROUGHPUT
  num=num+1
  validation_2C,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2D:  DIFFRACTION
  num=num+1
  validation_2D,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2F:  BASIC THROUGHPUT
  num=num+1
  validation_2F,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 2H:  SENSITIVITY MATRIX
  num=num+1
  validation_2H,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 3A:  DIFFUSE
  num=num+1
  validation_3A,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 3B:  CHARGE DIFFUSION
  num=num+1
  validation_3B,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 3C:  SIGNAL VS VARIANCE
  num=num+1
  validation_3C,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 3D:  SIZE VS INTENSITY
  num=num+1
  validation_3D,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 4A:  SPEED TEST
  num=num+1
  validation_4A,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 4B:   DYNAMIC TRANS OPT TEST
  num=num+1
  validation_4B,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 4C:   BRIGHT STAR OPT TEST
  num=num+1
  validation_4C,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 4D:   BACKGROUND OPT TEST
  num=num+1
  validation_4D,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;VALIDATION TASK 4E:   BACKGROUND TEST
  num=num+1
  validation_4E,num,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

;;need to fix 1C/1E/1F naming convention

;-------------------------------
;VALIDATION SUMMARY
;-------------------------------
finish:
  print,'Summary'
  !p.multi=0
  !p.charsize=0.6
  page=1
  value_regression=fltarr(tn,tn)
  readcol,'regression.txt',aa,bb,cc,/silent
  for kk=0L,N_elements(aa)-1 do begin
     value_regression(aa(kk),bb(kk))=cc(kk)
  endfor

  readcol,'regression_version.txt',a,b,format='(A,A)',/SILENT
  sss=strlen(b(0))
  vers_regr='Commit: ..'+strmid(b(0),sss-7,sss-1)

  leftm=-0.05
  topm=1.02
  col1=leftm & col2=0.18+leftm & col3=0.38+leftm & col4=0.53+leftm & col5=0.68+leftm & col6=0.87+leftm & col7=0.99+leftm

  intpos=leftm+0.45
  regpos=leftm+0.82

  pos=topm
  linestep=0.0175
  ssq='_____________________________________________________________________________________________________________________________________________________________________________________________________'

  for regr=0,1 do begin

     if regr eq 1 then begin
        value=value_regression
     endif

     for ii=0,19 do begin
        for jj=0,19 do begin
           if jj eq 0 then result_string(ii,jj)=string(value(ii,jj),format='(F7.3)')+unit(ii,jj)
           if jj ne 0 then result_string(ii,jj)=string(value(ii,jj),format='(F7.3)')+unit(ii,jj)
        endfor
     endfor

     if regr eq 0 then begin
        result_string_current=result_string
        value_current=value
     endif else begin
        result_string_regression=result_string
        result_string=result_string_current
        value_regression=value
        value=value_current
     endelse

  endfor


  readcol,'unittest_output',unita1,unitb1,unitc1,unitd1,unite1,unitf1,format='(A,A,F,F,A,A)',delimiter='|',/silent
  readcol,'unittest_output_regression',unita2,unitb2,unitc2,unitd2,unite2,unitf2,format='(A,A,F,F,A,A)',delimiter='|',/silent
  if N_elements(unita1) ne N_elements(unita2) then begin
     spawn,'cp unittest_output unittest_output_regression'
     readcol,'unittest_output_regression',unita2,unitb2,unitc2,unitd2,unite2,unitf2,format='(A,A,F,F,A,A)',delimiter='|',/silent
  endif
  passunit=N_elements(where(unitf1 eq 'Pass'))
  sameunit=0L & diffunit=0L
  totalunit=N_elements(where(unitf1 eq 'Fail' or unitf1 eq 'Pass'))
  good=where(unitf1 eq 'Fail' or unitf1 eq 'Pass')
  xyouts,col1,pos,'UNIT TESTS:'
  pos=pos-linestep
  pos=pos-linestep

  for i=0L,N_elements(good)-1 do begin

        if pos eq topm then begin
           plot,[0,1],[0,1],/nodata,xstyle=4,ystyle=4
           ss='Validation Summary for '+vers+'; Regression Check with '+vers_regr
           xyouts,col1,pos,ss
           ss='Page '+string(page,format='(I1)')+'/4'
           xyouts,col7,pos,ss

           pos=pos-linestep
           xyouts,col1,pos,ssq
           pos=pos-linestep
           ss='Validation Task'
           xyouts,col1,pos,ss
           ss='Validation Metric'
           xyouts,col2,pos,ss
           ss='Accuracy or Value'
           xyouts,col3,pos,ss
           ss='P/F (Tolerance)'
           xyouts,col4,pos,ss
           ss='Comparison Source'
           xyouts,col5,pos,ss
           ss='Previous Value'
           xyouts,col6,pos,ss
           ss='% Difference'
           xyouts,col7,pos,ss

           pos=pos-linestep
           xyouts,col1,pos,ssq
           pos=pos-linestep
           pos=pos-linestep/2.0
        endif

     xyouts,col1,pos,'   '+unita1(good(i))
     if unitf1(good(i)) eq 'Pass' then pass=1 else pass=0
     if pass eq 1 then cc=40 else cc=250
     xyouts,col2,pos,unitb1(good(i))
     if unite1(good(i)) eq 'none' then xyouts,col3,pos,unitc1(good(i)),color=cc
     if unite1(good(i)) ne 'none' then begin
        ss=string(unitc1(good(i)))+' '+unite1(good(i))
        xyouts,col3,pos,unitc1(good(i)),color=cc
     endif
     ss=unitf1(good(i))+' ( = '+string(unitd1(good(i)))+' )'
     xyouts,col4,pos,ss,color=cc
     xyouts,col5,pos,'Analytic Result',color=cc
     if unitf1(good(i)) eq unitf2(good(i)) then begin
        cc=100
        xyouts,col6,pos,unitc2(good(i)),color=cc
        sameunit=sameunit+1L
     endif else begin
        cc=229
        xyouts,col6,pos,unitc2(good(i)),color=cc
        diffunit=diffunit+1L
     endelse
     pcterr=(unitc1(good(i))-unitc2(good(i)))/(abs(unitc2(good(i)))>1e-30)*100.0
     ss=string(pcterr,format='(F8.2)')+'%'
     xyouts,col7,pos,ss,color=cc
     pos=pos-linestep

     if i mod 55 eq 0 then begin
        pos=topm
     endif

  endfor

  pos=pos-linestep
  pos=pos-linestep

  ss='UNIT TEST SUMMARY:   '+string(passunit,format='(I3)')+' Pass and '+string(totalunit-passunit,format='(I3)')+' Fail'
  xyouts,intpos,pos,ss

  ss='REGRESSION CHECK:   '+string(sameunit,format='(I3)')+' Same and '+string(diffunit,format='(I3)')+' Different'
  xyouts,regpos,pos,ss

  pos=pos-linestep
  pos=pos-linestep
  pos=pos-linestep

  ; REMOVE THIS WHEN GET OVER 55 UNIT TESTS
  pos=topm
  page=page+1L

  tpass=0L & tfail=0L
  tsame=0L & tdiff=0L
  for ii=0,19 do begin
     for jj=0,19 do begin

        if pos eq topm then begin
           plot,[0,1],[0,1],/nodata,xstyle=4,ystyle=4
           ss='Validation Summary for '+vers+'; Regression Check with '+vers_regr
           xyouts,col1,pos,ss
           ss='Page '+string(page,format='(I1)')+'/4'
           xyouts,col7,pos,ss

           pos=pos-linestep
           xyouts,col1,pos,ssq
           pos=pos-linestep
           ss='Validation Task'
           xyouts,col1,pos,ss
           ss='Validation Metric'
           xyouts,col2,pos,ss
           ss='Accuracy or Value'
           xyouts,col3,pos,ss
           ss='P/F (Tolerance)'
           xyouts,col4,pos,ss
           ss='Comparison Source'
           xyouts,col5,pos,ss
           ss='Previous Value'
           xyouts,col6,pos,ss
           ss='% Difference'
           xyouts,col7,pos,ss

           pos=pos-linestep
           xyouts,col1,pos,ssq
           pos=pos-linestep
           pos=pos-linestep/2.0
        endif

        if ii eq 0 and jj eq 0 then begin

           xyouts,col1,pos,'INTEGRATION TESTS:'
           pos=pos-linestep
           pos=pos-linestep

        endif

        xyouts,col1,pos,task(ii,jj)
        if (tolerance_low(ii,jj) eq 0 and tolerance_high(ii,jj) eq 0) then goto,skip
        xyouts,col2,pos,name(ii,jj)

        if (value(ii,jj) ge tolerance_low(ii,jj) and value(ii,jj) le tolerance_high(ii,jj)) then pass=1 else pass=0
        if pass eq 1 then ss='Pass ' else ss='Fail '
        if pass eq 1 then cc=40 else cc=250
        if pass eq 1 then tpass=tpass+1
        if pass eq 0 then tfail=tfail+1

        xyouts,col3,pos,result_string(ii,jj),color=cc
        xyouts,col5,pos,comparison(ii,jj),color=cc

        if tolerance_low(ii,jj) eq 0.0 then begin
           ss=ss+' ( < '+string(tolerance_high(ii,jj),format='(F7.3)')+' )'
        endif else if tolerance_high(ii,jj) ge 1e30 then begin
           ss=ss+' ( > '+string(tolerance_low(ii,jj),format='(F7.3)')+' )'
        endif else begin
           ss=ss+' ( '+string(0.5*(tolerance_low(ii,jj)+tolerance_high(ii,jj)),format='(F5.2)')+' +/- '+string(0.5*(tolerance_high(ii,jj)-tolerance_low(ii,jj)),format='(F5.2)')+' )'
        endelse

        xyouts,col4,pos,ss,color=cc

        pcterr=(value(ii,jj)-value_regression(ii,jj))/(abs(value_regression(ii,jj))>1e-30)*100.0
        if (abs(pcterr) gt 0.01) then pass=1 else pass=0
        if pass eq 1 then cc=220 else cc=100
        if pass eq 1 then tdiff=tdiff+1
        if pass eq 0 then tsame=tsame+1

        xyouts,col6,pos,result_string_regression(ii,jj),color=cc
        ss=string(pcterr,format='(F8.2)')+'%'
        xyouts,col7,pos,ss,color=cc

        pos=pos-linestep
skip:
     endfor
     pos=pos-linestep
     if ii eq 5 or ii eq 14 then begin
        pos=topm
        page=page+1L
     endif
  endfor



  pos=pos-linestep

  ss='INTEGRATION TEST SUMMARY:   '+string(tpass,format='(I3)')+' Pass and '+string(tfail,format='(I3)')+' Fail'
  xyouts,intpos,pos,ss
  ss='REGRESSION CHECK:   '+string(tsame,format='(I3)')+' Same and '+string(tdiff,format='(I3)')+' Different'
  xyouts,regpos,pos,ss

  pos=pos-linestep
  pos=pos-linestep
  pos=pos-linestep
  pos=pos-linestep
  xyouts,col1,pos,ssq
  pos=pos-linestep
  pos=pos-linestep
  ss='ALL '+string(tpass+tfail+totalunit,format='(I3)')+' TESTS SUMMARY:   '+$
     string(tpass+passunit,format='(I3)')+' Pass and '+string(tfail+totalunit-passunit,format='(I3)')+' Fail'
  xyouts,intpos,pos,ss
  ss='REGRESSION CHECK:   '+string(tsame+sameunit,format='(I3)')+' Same and '+string(tdiff+diffunit,format='(I3)')+' Different'
  xyouts,regpos,pos,ss


  ; RECENT DEVELOPMENT HISTORY
  spawn,'git log origin/dev --format=%H'' ''%ad'' ''%an'' ''%s > temp.txt'
  readcol,'temp.txt',githash,format='(A)',/silent
  readcol,'../bin/version',a1,b1,format='(A,A)',/SILENT
  readcol,'regression_version.txt',a2,b2,format='(A,A)',/SILENT
  c=strarr(N_elements(githash))
  openr,2,'temp.txt'
  readf,2,c
  close,2
  pos=topm
  plot,[0,1],[0,1],/nodata,xstyle=4,ystyle=4
  xyouts,col1,pos,ssq
  pos=pos-linestep
  xyouts,col1,pos,'Recent Development History'
  pos=pos-linestep
  xyouts,col1,pos,ssq
  pos=pos-linestep
  for kk=0,58 do begin
     cc=0
     if strcmp(githash(kk),b1(0)) gt 0 then cc=40
     if strcmp(githash(kk),b2(0)) gt 0 then cc=100
     xyouts,col1,pos,c(kk),color=cc
     pos=pos-linestep
  endfor

  ; SAVE REGRESSION
  if regression eq 'y' then begin
     openw,1,'regression.txt'
     for kk=0,19 do begin
        for ll=0,19 do begin
           printf,1,kk,ll,value_current(kk,ll)
        endfor
     endfor
     close,1
     spawn,'cp ../bin/version regression_version.txt'
     spawn,'cp unittest_output unittest_output_regression'
  endif

  device,/close

END
