;;
;; @package phosim
;; @file validation_2C.pro
;; @brief validation task 2C
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

pro validation_2C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 2C'
  !p.multi=[0,2,2]

  err=1.0
  minph=500

  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,0.5,40.0,'Atmosphere',/data

  ss='System Throughput'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2C; '+vers
  xyouts,0.7,0.98,ss,/normal


  readcol,'baseline/atmos.dat',a,b,/SILENT
  x=a/1000.0 & y=b*100.0
  oplot,x,y

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,c(good)/b(good)*100.0,color=cc,psym=-3,linestyle=lc

     good=where(o gt minph)
     yint=interpol(y,x,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(c(good)/b(good)*100.0-yint)*err,color=20,psym=3

     if ii le 5 then begin
        if ii eq 0 then yintt=(yint-c(good)/b(good)*100.0)^2
        if ii ge 0 then yintt=[yintt,(yint-c(good)/b(good)*100.0)^2]
        value(nnn,7)=sqrt(mean(yintt))*10.0
     endif

  endfor

  xyouts,0.32,122,'baseline',/data
  sstring='airmass=1.2, on-axis, no clouds/aerosols/var'
  xyouts,0.32,115,sstring,color=250,/data
  sstring='airmass=2, on-axis, no clouds/aerosols/var'
  xyouts,0.32,108,sstring,color=50,/data
  sstring='airmass=1.2, 1.6!Uo!N off, no clouds/aerosols/var'
  xyouts,0.32,101,sstring,color=130,/data



  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,1.0,110.0,'Mirror',/data

  sstring='residual'
  xyouts,0.32,11,sstring,color=20,/data

  readcol,'baseline/m1.dat',a,b,/SILENT
  readcol,'baseline/m2.dat',a,c,/SILENT
  readcol,'baseline/m3.dat',a,d,/SILENT
  x=a/1000.0
  y=b*c*d*100.0
  oplot,x,y


  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif 
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,f(good)/c(good)*100.0,color=cc,psym=-3,linestyle=lc
     yint=interpol(y,x,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(f(good)/c(good)*100.0-yint)*err,color=20,psym=3
  endfor



  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,1.0,110.0,'L1+L2',/data

  readcol,'baseline/lens1.dat',a,b,/SILENT
  readcol,'baseline/lens2.dat',a,c,/SILENT
  x=a/1000.0
  y=b*c*100.0
  oplot,x,y

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,j(good)/f(good)*100.0,color=cc,psym=-3,linestyle=lc
     yint=interpol(y,x,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(j(good)/f(good)*100.0-yint)*err,color=20,psym=3
  endfor

  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,1.0,110.0,'Filter',/data

  readcol,'baseline/filter_u.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x0=a/1000.0 & y0=b*100.0
  readcol,'baseline/filter_g.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x1=a/1000.0 & y1=b*100.0
  readcol,'baseline/filter_r.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x2=a/1000.0 & y2=b*100.0
  readcol,'baseline/filter_i.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x3=a/1000.0 & y3=b*100.0
  readcol,'baseline/filter_z.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x4=a/1000.0 & y4=b*100.0
  readcol,'baseline/filter_y.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x5=a/1000.0 & y5=b*100.0

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,l(good)/j(good)*100.0,color=cc,psym=-3,linestyle=lc
     if filter eq 0 then yint=interpol(y0,x0,a(good)/1000.0+0.3)
     if filter eq 1 then yint=interpol(y1,x1,a(good)/1000.0+0.3)
     if filter eq 2 then yint=interpol(y2,x2,a(good)/1000.0+0.3)
     if filter eq 3 then yint=interpol(y3,x3,a(good)/1000.0+0.3)
     if filter eq 4 then yint=interpol(y4,x4,a(good)/1000.0+0.3)
     if filter eq 5 then yint=interpol(y5,x5,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(l(good)/j(good)*100.0-yint)*err,color=20,psym=3
  endfor


  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,0.5,40.0,'L3',/data


  ss='System Throughput'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2C; '+vers
  xyouts,0.7,0.98,ss,/normal


  readcol,'baseline/lens3.dat',a,b,/SILENT
  oplot,a/1000.0,b*100
  x=a/1000.0 
  y=b*100.0

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,n(good)/l(good)*100.0,color=cc,psym=-3,linestyle=lc
     yint=interpol(y,x,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(n(good)/l(good)*100.0-yint)*err,color=20,psym=3
  endfor


  xyouts,0.32,122,'baseline',/data
  sstring='airmass=1.2, on-axis, no clouds/aerosols/var'
  xyouts,0.32,115,sstring,color=250,/data
  sstring='airmass=2, on-axis, no clouds/aerosols/var'
  xyouts,0.32,108,sstring,color=50,/data
  sstring='airmass=1.2, 1.6!Uo!N off, no clouds/aerosols/var'
  xyouts,0.32,101,sstring,color=130,/data


  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,1.0,110.0,'Detector',/data

  sstring='residual'
  xyouts,0.32,11,sstring,color=20,/data

  readcol,'baseline/detector.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  x=a/1000.0
  y=b*100.0

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,o(good)/n(good)*100.0,color=cc,psym=-3,linestyle=lc
     yint=interpol(y,x,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(o(good)/n(good)*100.0-yint)*err,color=20,psym=3
  endfor


  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,0.35,110.0,'Total Throughput',/data

  readcol,'baseline/atmos.dat',aa,bb,/SILENT

  readcol,'baseline/m1.dat',aa,bb,/SILENT
  bua=bb & bga=bb & bra=bb & bia=bb & bza=bb & bya=bb
  readcol,'baseline/m2.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb
  readcol,'baseline/m3.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb
  readcol,'baseline/lens1.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb
  readcol,'baseline/lens2.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb
  readcol,'baseline/filter_u.dat',aa,bb,/SILENT
  bua=bua*bb
  readcol,'baseline/filter_g.dat',aa,bb,/SILENT
  bga=bga*bb
  readcol,'baseline/filter_r.dat',aa,bb,/SILENT
  bra=bra*bb
  readcol,'baseline/filter_i.dat',aa,bb,/SILENT
  bia=bia*bb
  readcol,'baseline/filter_z.dat',aa,bb,/SILENT
  bza=bza*bb
  readcol,'baseline/filter_y.dat',aa,bb,/SILENT
  bya=bya*bb
  readcol,'baseline/lens3.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb
  readcol,'baseline/detector.dat',aa,bb,/SILENT
  bua=bua*bb & bga=bga*bb & bra=bra*bb & bia=bia*bb & bza=bza*bb & bya=bya*bb

  baa=aa

  readcol,'baseline/total_u.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  bu=b
  readcol,'baseline/total_g.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  bg=b
  readcol,'baseline/total_r.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  br=b
  readcol,'baseline/total_i.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  bi=b
  readcol,'baseline/total_z.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  bz=b
  readcol,'baseline/total_y.dat',a,b,/SILENT
  oplot,a/1000.0,b*100.0
  by=b

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,o(good)/b(good)*100.0,color=cc,psym=-3,linestyle=lc
     if ii lt 10 then openw,1,'totalthroughput_'+string(ii,format='(I1)')+'.txt'
     if ii ge 10 then openw,1,'totalthroughput_'+string(ii,format='(I2)')+'.txt'
     for jj=0L,N_elements(good)-1 do begin
        printf,1,a(good(jj))/1000.0+0.3,o(good(jj))/b(good(jj))*100.0
     endfor
     close,1
     if ii eq 0 then du=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bu)/10.))
     if ii eq 1 then dg=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bg)/10.))
     if ii eq 2 then dr=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(br)/10.))
     if ii eq 3 then di=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bi)/10.))
     if ii eq 4 then dz=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bz)/10.))
     if ii eq 5 then dy=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(by)/10.))
     if ii eq 6 then eu=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bu)/10.))
     if ii eq 7 then eg=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bg)/10.))
     if ii eq 8 then er=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(br)/10.))
     if ii eq 9 then ei=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bi)/10.))
     if ii eq 10 then ez=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(bz)/10.))
     if ii eq 11 then ey=1000.*(2.5*alog10(total(o(good)/b(good)))-2.5*alog10(total(by)/10.))
     mlambda=0.0 & mtotal=0.0
     for jj=0L,N_elements(good)-1 do begin
        mlambda=mlambda+(a(good(jj))/1000.0+0.3)*(o(good(jj))/b(good(jj)))
        mtotal=mtotal+(o(good(jj))/b(good(jj)))
     endfor
     mlambda=mlambda/mtotal
     if ii eq 0 then mu=mlambda
     if ii eq 1 then mg=mlambda
     if ii eq 2 then mr=mlambda
     if ii eq 3 then mi=mlambda
     if ii eq 4 then mz=mlambda
     if ii eq 5 then my=mlambda
  endfor

  ss='u: '+string(eu,format='(I4)')+' to '+string(du,format='(I4)')+' millimags'
;xyouts,0.1,0.42,ss,/normal
  ss='g: '+string(eg,format='(I4)')+' to '+string(dg,format='(I4)')+' millimags'
;xyouts,0.1,0.40,ss,/normal
  ss='r: '+string(er,format='(I4)')+' to '+string(dr,format='(I4)')+' millimags'
;xyouts,0.1,0.38,ss,/normal
  ss='i: '+string(ei,format='(I4)')+' to '+string(di,format='(I4)')+' millimags'
;xyouts,0.1,0.36,ss,/normal
  ss='z: '+string(ez,format='(I4)')+' to '+string(dz,format='(I4)')+' millimags'
;xyouts,0.1,0.34,ss,/normal
  ss='y: '+string(ey,format='(I4)')+' to '+string(dy,format='(I4)')+' millimags'
;xyouts,0.1,0.32,ss,/normal

  ss='<!4k!3!Du!N> = '+string(mu,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.46,ss,/normal
  ss='<!4k!3!Dg!N> = '+string(mg,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.44,ss,/normal
  ss='<!4k!3!Dr!N> = '+string(mr,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.42,ss,/normal
  ss='<!4k!3!Di!N> = '+string(mi,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.40,ss,/normal
  ss='<!4k!3!Dz!N> = '+string(mz,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.38,ss,/normal
  ss='<!4k!3!Dy!N> = '+string(my,format='(F4.2)')+' !4l!3m'
  xyouts,0.32,0.36,ss,/normal



  plot,[0.3,1.2],[-10,130],/nodata,/xstyle,xtitle='Wavelength (!4l!3m)',ytitle='Monochromatic Efficiency (%)',/ystyle
  xyouts,0.75,110.0,'Instrumental Throughput',/data

  oplot,aa/1000.,bua*100.
  x0=aa(0:(N_elements(bua)-1))/1000. & y0=bua*100.
  oplot,aa/1000.,bga*100.
  x1=aa(0:(N_elements(bga)-1))/1000. & y1=bga*100.
  oplot,aa/1000.,bra*100.
  x2=aa(0:(N_elements(bra)-1))/1000. & y2=bra*100.
  oplot,aa/1000.,bia*100.
  x3=aa(0:(N_elements(bia)-1))/1000. & y3=bia*100.
  oplot,aa/1000.,bza*100.
  x4=aa(0:(N_elements(bza)-1))/1000. & y4=bza*100.
  oplot,aa/1000.,bya*100.
  x5=aa(0:(N_elements(bya)-1))/1000. & y5=bya*100.

  for ii=0L,17 do begin
     filter=(ii mod 6)
     if ii le 5 then begin
        tname='throughput_lsst_e_220'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=250 & lc=1
     endif
     if ii ge 6 and ii le 11 then begin
        tname='throughput_lsst_e_221'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R22_S11_E00'
        cc=50 & lc=1
     endif
     if ii ge 12 and ii le 17 then begin
        tname='throughput_lsst_e_222'+string(filter,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R24_S12_E00'
        cc=130 & lc=1
     endif
     tname1=tname+'0.txt'
     tname2=tname+'1.txt'
     readcol,tname1,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,/SILENT
     readcol,tname2,a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,/SILENT
     a=a1 & b=b1+b2 & c=c1+c2 & d=d1+d2 & e=e1+e2 & f=f1+f2 & g=g1+g2
     h=h1+h2 & i=i1+i2 & j=j1+j2 & k=k1+k2 & l=l1+l2 & m=m1+m2 & n=n1+n2 & o=o1+o2
     good=where(o gt minph)
     oplot,a(good)/1000.0+0.3,o(good)/c(good)*100.0,color=cc,psym=-3,linestyle=lc
     if filter eq 0 then yint=interpol(y0,x0,a(good)/1000.0+0.3)
     if filter eq 1 then yint=interpol(y1,x1,a(good)/1000.0+0.3)
     if filter eq 2 then yint=interpol(y2,x2,a(good)/1000.0+0.3)
     if filter eq 3 then yint=interpol(y3,x3,a(good)/1000.0+0.3)
     if filter eq 4 then yint=interpol(y4,x4,a(good)/1000.0+0.3)
     if filter eq 5 then yint=interpol(y5,x5,a(good)/1000.0+0.3)
     if ii le 5 then oplot,a(good)/1000.0+0.3,(o(good)/c(good)*100.0-yint)*err,color=20,psym=3
     if ii eq 0 then du=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bua)))
     if ii eq 1 then dg=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bga)))
     if ii eq 2 then dr=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bra)))
     if ii eq 3 then di=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bia)))
     if ii eq 4 then dz=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bza)))
     if ii eq 5 then dy=1000.*(2.5*alog10(total(o(good)/c(good)))-2.5*alog10(total(bya)))
     if ii eq 0 then begin
        zu=total(o(good)/c(good))
        zut=total(bua)
     endif
     if ii eq 1 then begin
        zg=total(o(good)/c(good))
        zgt=total(bga)
     endif
     if ii eq 2 then begin
        zr=total(o(good)/c(good))
        zrt=total(bra)
     endif
     if ii eq 3 then begin
        zi=total(o(good)/c(good))
        zit=total(bia)
     endif
     if ii eq 4 then begin
        zz=total(o(good)/c(good))
        zzt=total(bza)
     endif
     if ii eq 5 then begin
        zy=total(o(good)/c(good))
        zyt=total(bya)
     endif
  endfor

  ss='u   '+string(du,format='(I4)')+' millimags'
  xyouts,0.60,0.46,ss,/normal
  ss='g   '+string(dg,format='(I4)')+' millimags'
  xyouts,0.60,0.44,ss,/normal
  ss='r   '+string(dr,format='(I4)')+' millimags'
  xyouts,0.60,0.42,ss,/normal
  ss='i   '+string(di,format='(I4)')+' millimags'
  xyouts,0.60,0.40,ss,/normal
  ss='z   '+string(dz,format='(I4)')+' millimags'
  xyouts,0.60,0.38,ss,/normal
  ss='y   '+string(dy,format='(I4)')+' millimags'
  xyouts,0.60,0.36,ss,/normal
  ss='rms '+string(sqrt(1./6.*(du*du+dg*dg+dr*dr+di*di+dz*dz+dy*dy)),format='(I4)')+' millimags'
  xyouts,0.60,0.34,ss,/normal

  tolerance_low(nnn,0)=0.
  tolerance_high(nnn,0)=100.
  tolerance_low(nnn,7)=0.0
  tolerance_high(nnn,7)=50.0

  value(nnn,0)=sqrt(1./6.*(du*du+dg*dg+dr*dr+di*di+dz*dz+dy*dy))
  value(nnn,1)=-48.6-2.5*alog10(6.626e-27*500.0/zu/(!Pi*(418.0^2-255.8^2)))
  value(nnn,2)=-48.6-2.5*alog10(6.626e-27*500.0/zg/(!Pi*(418.0^2-255.8^2)))
  value(nnn,3)=-48.6-2.5*alog10(6.626e-27*500.0/zr/(!Pi*(418.0^2-255.8^2)))
  value(nnn,4)=-48.6-2.5*alog10(6.626e-27*500.0/zi/(!Pi*(418.0^2-255.8^2)))
  value(nnn,5)=-48.6-2.5*alog10(6.626e-27*500.0/zz/(!Pi*(418.0^2-255.8^2)))
  value(nnn,6)=-48.6-2.5*alog10(6.626e-27*500.0/zy/(!Pi*(418.0^2-255.8^2)))
  tolerance_low(nnn,1)=-48.6-2.5*alog10(6.626e-27*500.0/zut/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_low(nnn,2)=-48.6-2.5*alog10(6.626e-27*500.0/zgt/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_low(nnn,3)=-48.6-2.5*alog10(6.626e-27*500.0/zrt/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_low(nnn,4)=-48.6-2.5*alog10(6.626e-27*500.0/zit/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_low(nnn,5)=-48.6-2.5*alog10(6.626e-27*500.0/zzt/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_low(nnn,6)=-48.6-2.5*alog10(6.626e-27*500.0/zyt/(!Pi*(418.0^2-255.8^2)))-0.2
  tolerance_high(nnn,1)=-48.6-2.5*alog10(6.626e-27*500.0/zut/(!Pi*(418.0^2-255.8^2)))+0.2
  tolerance_high(nnn,2)=-48.6-2.5*alog10(6.626e-27*500.0/zgt/(!Pi*(418.0^2-255.8^2)))+0.2
  tolerance_high(nnn,3)=-48.6-2.5*alog10(6.626e-27*500.0/zrt/(!Pi*(418.0^2-255.8^2)))+0.2
  tolerance_high(nnn,4)=-48.6-2.5*alog10(6.626e-27*500.0/zit/(!Pi*(418.0^2-255.8^2)))+0.2
  tolerance_high(nnn,5)=-48.6-2.5*alog10(6.626e-27*500.0/zzt/(!Pi*(418.0^2-255.8^2)))+0.2
  tolerance_high(nnn,6)=-48.6-2.5*alog10(6.626e-27*500.0/zyt/(!Pi*(418.0^2-255.8^2)))+0.2


  name(nnn,0)='Avg band flat AB!D500nm!N zeropoint'
  name(nnn,1)='  u flat AB!D500nm!N zeropoint'
  name(nnn,2)='  g flat AB!D500nm!N zeropoint'
  name(nnn,3)='  r flat AB!D500nm!N zeropoint'
  name(nnn,4)='  i flat AB!D500nm!N zeropoint'
  name(nnn,5)='  z flat AB!D500nm!N zeropoint'
  name(nnn,6)='  y flat AB!D500nm!N zeropoint'
  name(nnn,7)='RMS Atm Transmission Error'

  unit(nnn,0)=' millimags'
  comparison(nnn,0)='LSST Design'
  unit(nnn,1)=' mags'
  comparison(nnn,1)='LSST Design'
  unit(nnn,2)=' mags'
  comparison(nnn,2)='LSST Design'
  unit(nnn,3)=' mags'
  comparison(nnn,3)='LSST Design'
  unit(nnn,4)=' mags'
  comparison(nnn,4)='LSST Design'
  unit(nnn,5)=' mags'
  comparison(nnn,5)='LSST Design'
  unit(nnn,6)=' mags'
  comparison(nnn,6)='LSST Design'
  unit(nnn,7)=' millimags'
  comparison(nnn,7)='MODTRAN Code'

  task(nnn,0)='2C Throughput'

END
