;;
;; @package phosim
;; @file validation_3C.pro
;; @brief validation task 3C
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

pro validation_3C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3C'
  !p.multi=[0,1,2]
  aa=findgen(33)/32.*2*!Pi
  usersym,cos(aa),sin(aa)

  x0=65000.0
  x0=1e5


  tolerance_low(nnn,0)=2.9*0.0
  tolerance_high(nnn,0)=2.9*2.0


  N=15L & Na=3L
  flux=fltarr(N)
  variance=fltarr(N)


  auto=fltarr(Na,Na,N)
  autoc=fltarr(Na,Na,N)

 ;approx DES data values from AR
  auto_dec=fltarr(4,4)
  auto_dec(0,1)=0.018
  auto_dec(1,0)=0.006
  auto_dec(1,1)=0.005
  auto_dec(0,2)=0.0025
  auto_dec(2,0)=0.0021
  auto_dec(2,1)=0.002
  auto_dec(1,2)=0.002
  auto_dec(2,2)=0.0012
  auto_dec(3,0)=0.0011
  auto_dec(0,3)=0.0010
  auto_dec(3,1)=0.0008
  auto_dec(1,3)=0.0007
  auto_dec(3,2)=0.0006
  auto_dec(2,3)=0.0005
  auto_dec(3,3)=0.0002

  prbuffer=0

  for i=0L,N-1 do begin

     if i lt 10 then file='320'+string(i,format='(I1)')+'_f2_R22_S11'
     if i ge 10 then file='32'+string(i,format='(I2)')+'_f2_R22_S11'
     j=i+20
     if j lt 10 then file2='320'+string(j,format='(I1)')+'_f2_R22_S11'
     if j ge 10 then file2='32'+string(j,format='(I2)')+'_f2_R22_S11'

     filename='lsst_e_'+file+'_E000.fits.gz'
     data=mrdfits(filename,0,header,/silent)
     filename='lsst_e_'+file2+'_E000.fits.gz'
     data2=mrdfits(filename,0,header,/silent)



     dataa=data+data2
     datab=data-data2

     column=dblarr(4000)
     for ii=0L,4000-1 do begin
        for jj=0L,4072-1 do begin
           column(ii)=column(ii)+dataa(ii,jj)
        endfor
     endfor
     maxv=max(column)
     minx=4000 & maxx=0
     for ii=0L,4000-1 do begin
        if column(ii) gt maxv*0.9 then begin
           if ii lt minx then minx=ii
           if ii gt maxx then maxx=ii
        endif
     endfor
     buffer=round((maxx-minx)/2)
     if buffer gt 1980 then buffer=1980
     if max(data) eq 100000 then buffer=prbuffer
     prbuffer=buffer




     resulta=moment(dataa((2000-buffer):(2000+buffer),(2036-buffer):(2036+buffer)))
     resultb=moment(datab((2000-buffer):(2000+buffer),(2036-buffer):(2036+buffer)))

     flux(i)=resulta(0)/2.0
     variance(i)=resultb(1)/2.0

     aa=0.0 & bb=0.0 & cc=0.0 & dd=0.0 & ee=0.0
     for mm=0L,Na-1 do begin
        for nn=0L,Na-1 do begin
           shortbuffer=buffer-max([mm,nn])
           for k=2000-shortbuffer,2000+shortbuffer,1 do begin
              for l=2036-shortbuffer,2036+shortbuffer,1 do begin
                 auto(mm,nn,i)=auto(mm,nn,i)+(datab(k,l)*datab(k-mm,l+nn)-resultb(0)*resultb(0)/4.0)
                 auto(mm,nn,i)=auto(mm,nn,i)+(datab(k,l)*datab(k-mm,l-nn)-resultb(0)*resultb(0)/4.0)
                 auto(mm,nn,i)=auto(mm,nn,i)+(datab(k,l)*datab(k+mm,l+nn)-resultb(0)*resultb(0)/4.0)
                 auto(mm,nn,i)=auto(mm,nn,i)+(datab(k,l)*datab(k+mm,l-nn)-resultb(0)*resultb(0)/4.0)
                 autoc(mm,nn,i)=autoc(mm,nn,i)+4.0
              endfor
           endfor
        endfor
     endfor
     for mm=0L,Na-1 do begin
        for nn=0L,Na-1 do begin
           auto(mm,nn,i)=auto(mm,nn,i)/autoc(mm,nn,i)/flux(i)/2.0
        endfor
     endfor
   endfor

  varianceerr=flux/sqrt((2.0*buffer)^2)*sqrt(2.0)

  plot,flux,variance>1,psym=4,xr=[1e1,1e6],yr=[1e1,1e6],xtitle='Signal (electrons)',ytitle='Variance (electrons)',/xstyle,/ystyle,/xlog,/ylog

  errplot,flux,variance-varianceerr,variance+varianceerr,width=0.0

  xxx=10.^(1.0+5.0*findgen(1000)/1000.)
  oplot,xxx,xxx


;  bnlflux=[20000.0,40000.0,60000.0,80000.,100000.0,120000.,140000.,160000.]/150000.*100000.
;  bnlvariance=[20000.0,40000.0,60000.0*32.0/33.0,80000.0*40./45.,100000.*50./55.,100000.*40./55.,100000.*5./55.,100000.*5./55.]/150000.*100000.

  minchi=1e30
  for aa=-3.0,0.0,0.01 do begin
  for bb=4.0,6.0,0.02 do begin
  for cc=3.5,5.0,0.02 do begin
     yyy=(1.0-10.^aa*(flux/1e5)^1)*(1.0-1.0/(1.0+exp(-(flux-10.^bb)/(10.^cc))))
     good=where(variance ne 0.0)
     chi=total((yyy(good)-variance(good)/flux(good))^2)
     if chi lt minchi then begin
        bestaa=aa
        bestbb=bb
        bestcc=cc
        minchi=chi
     endif
  endfor
  endfor
  endfor
;  print,10.^bestaa,10.^bestbb,10.^bestcc
  yyy=(1.0-10.^bestaa*(xxx/1e5)^1)*(1.0-1.0/(1.0+exp(-(xxx-10.^bestbb)/(10.^bestcc))))*xxx
  oplot,xxx,yyy,linestyle=0



;  oplot,bnlflux,bnlvariance,psym=5,color=250
  readcol,'ptc.txt',pta,ptb,ptc,ptd,/silent
  gg=where(pta eq 0)
  oplot,ptc(gg),ptd(gg),psym=8,color=120,symsize=0.5
  gg=where(pta eq 1)
  oplot,ptc(gg),ptd(gg),psym=8,color=80,symsize=0.5

  for ii=0,1 do begin
  minchi=1e30
  for aa=-3.0,0.0,0.01 do begin
  for bb=4.0,6.0,0.02 do begin
  for cc=3.5,5.0,0.02 do begin
     gg=where(pta eq ii)
     yyy=(1.0-10.^aa*(ptc(gg)/1e5)^1)*(1.0-1.0/(1.0+exp(-(ptc(gg)-10.^bb)/(10.^cc))))
     chi=total((yyy-ptd(gg)/ptc(gg))^2)
     if chi lt minchi then begin
        bestdd=aa
        bestee=bb
        bestff=cc
        minchi=chi
     endif
  endfor
  endfor
  endfor
  ;print,10.^bestdd,10.^bestee,10.^bestff
  yyy=(1.0-10.^bestdd*(xxx/1e5)^1)*(1.0-1.0/(1.0+exp(-(xxx-10.^bestee)/(10.^bestff))))*xxx
  if ii eq 0 then oplot,xxx,yyy,color=120,linestyle=0
  if ii eq 1 then oplot,xxx,yyy,color=80,linestyle=0
  if ii eq 0 then begin
     bestdd1=bestdd
     bestee1=bestee
     bestff1=bestff
  endif
  endfor

  oplot,flux,variance,psym=4

  legend,psym=[4,8,8],['PhoSim','ITL','E2V'],color=[0,120,80],symsize=[1,0.5,0.5],/top,/left

  nu=variance/flux
  pcterr=varianceerr/flux

  plot,flux,nu,psym=4,/xlog,xr=[1e1,1e6],xtitle='Signal (electrons)',ytitle='Variance/Signal',yr=[0,1.1],/xstyle

  errplot,flux,nu-pcterr,nu+pcterr,width=0.0

  yyy=(1.0-10.^bestaa*(xxx/1e5)^1)*(1.0-1.0/(1.0+exp(-(xxx-10.^bestbb)/(10.^bestcc))))
  oplot,xxx,yyy,linestyle=0

;   bnlnu=bnlvariance/bnlflux

;  oplot,bnlflux,bnlnu,psym=5,color=250
  gg=where(pta eq 0)
  oplot,ptc(gg),ptd(gg)/ptc(gg),psym=8,color=120,symsize=0.5
  gg=where(pta eq 1)
  oplot,ptc(gg),ptd(gg)/ptc(gg),psym=8,color=80,symsize=0.5

  yyy=(1.0-10.^bestdd1*(xxx/1e5)^1)*(1.0-1.0/(1.0+exp(-(xxx-10.^bestee1)/(10.^bestff1))))
  oplot,xxx,yyy,linestyle=0,color=120
  yyy=(1.0-10.^bestdd*(xxx/1e5)^1)*(1.0-1.0/(1.0+exp(-(xxx-10.^bestee)/(10.^bestff))))
  oplot,xxx,yyy,linestyle=0,color=80

  oplot,xxx,1.0+fltarr(N_elements(xxx)),linestyle=2

  oplot,flux,nu,psym=4

  ss='Variance vs. Signal on Difference Images'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3C; '+vers
  xyouts,0.7,0.98,ss,/normal

  !p.multi=0
  plot,[1e2,2e5],[-0.01,0.10],/nodata,xtitle='Signal',ytitle='Autocovariance or Sub-linear Variance',/xstyle,/ystyle,/xlog
;  plot,[1e2,2e5],[1e-3,1.0],/nodata,xtitle='Signal',ytitle='Autocovariance or Sub-linear Variance',/xstyle,/ystyle,/xlog,/ylog
  ccc=40L
     for mm=0L,Na-1 do begin
        for nn=0L,Na-1 do begin
;           if mm ne 0 or nn ne 0 then begin
              ccc=round(20L+float(sqrt(mm^2+nn^2)/sqrt(2.0)/(float(Na)-1.0))*(254.0-20.0))
;              print,mm,nn,auto(mm,nn,*)
              if mm ne 0 or nn ne 0 then oplot,flux,auto(mm,nn,*),psym=-5,color=ccc,linestyle=1
              if mm eq 0 and nn eq 0 then oplot,flux,(1-auto(mm,nn,*)),psym=-5,color=ccc,linestyle=1
              xx=10.^(findgen(1000)/1000.*5.0+1.0)
              oplot,[x0],[auto_dec(mm,nn)],psym=4,color=ccc
;              print,mm,nn,autoc(mm,nn,0),1.0/sqrt(autoc(mm,nn,0))
              oplot,xx,(xx/x0)^2*auto_dec(mm,nn),color=ccc
;              if mm eq 0 and nn eq 0 then oplot,bnlflux,1.0-bnlnu,psym=4,color=ccc
              if mm eq 0 and nn eq 0 then oplot,ptc,1.0-ptd/ptc,psym=4,color=ccc
              ss=string(mm,format='(I1)')+','+string(nn,format='(I1)')
              NQ=round(3.0*N/4.0)
              if mm ne 0 or nn ne 0 then xyouts,flux(NQ),auto(mm,nn,NQ),ss,/data,color=ccc
              if mm eq 0 and nn eq 0 then xyouts,flux(NQ),1.0-auto(mm,nn,NQ),ss,/data,color=ccc
;           endif
        endfor
     endfor

  ss='Variance vs. Signal on Difference Images'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3C; '+vers
  xyouts,0.7,0.98,ss,/normal

  task(nnn,0)='3C Variance vs. Signal'

END
