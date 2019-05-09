;;
;; @package phosim
;; @file validation_1H.pro
;; @brief validation task 1H
;;
;; @brief Created by:
;; @author En-Hsin Peng (Purdue)
;;
;; @brief Modified by:
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;


FUNCTION calculate_structure_function, img, pixelsize
  seed=100L
  dim=size(img)
  N=dim(1)
  Dstr=DBLARR(N)
  npx=DBLARR(N)
  maxr=FLOOR(ALOG(N)/ALOG(sqrt(2.0)))-1
  NLoops=1000
  FOR rd=0L,maxr DO BEGIN
     d=floor(sqrt(2)^(rd))
     FOR l=0L, Nloops-1 DO BEGIN
        ix=FLOOR(RANDOMU(seed)*N)
        iy=FLOOR(RANDOMU(seed)*N)
        phi=randomu(seed)*2*!Pi
        ixn=floor(ix+d*cos(phi))
        iyn=floor(iy+d*sin(phi))
        if (ixn ge 0 and ixn lt N and iyn ge 0 and iyn lt N) then begin
           Dstr[d]=Dstr[d]+(img[ix,iy]-img[ixn,iyn])^2
           npx[d]=npx[d]+1.0
        endif
     endfor
  endfor
  posid=WHERE(npx GT 0)
  r=FINDGEN(N)*pixelsize
  RETURN, [[r[posid]], [Dstr[posid]/npx[posid]]]
end



PRO validation_1H,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison
  print,'Task 1H'

  obsids=['1991','1992','1993','1994','1995','1996','1997','1998','1999']
  mm=5
  mm_all=N_ELEMENTS(obsids)
  mm=mm_all
  savec=fltarr(mm)
  ;Adams & Skrutskie 2MASS
  adams_fig4=[[0.09, 0.10, 0.10, 0.11, 0.12, 0.13, 0.13, 0.14, 0.15, 0.17, 0.18, 0.19, 0.21, 0.23, 0.25, 0.25, 0.27, 0.29, 0.29, 0.32, 0.32, 0.35, 0.38, 0.40, 0.42, 0.45, 0.47, 0.52, 0.55, 0.56, 0.63, 0.72, 0.82, 0.92, 0.97, 1.22, 1.33, 1.80, 2.08, 2.70, 3.70, 5.39], [1.64, 1.77, 1.70, 1.75, 1.83, 1.90, 1.79, 1.97, 2.15, 2.15, 2.25, 2.39, 2.48, 2.51, 2.39, 2.73, 2.45, 2.91, 3.16, 2.42, 2.67, 2.57, 2.70, 2.98, 2.80, 3.70, 2.87, 3.05, 2.84, 3.89, 2.84, 3.20, 2.87, 3.09, 3.09, 3.89, 3.66, 5.46, 4.89, 7.94, 11.71, 44.40]] 

  ;Ivezic et al. 2007, AJ 134, 3, 973.
  ivezic_fig22_cir=[[0.42,0.84,1.26,1.67,2.09],[0.019,0.025,0.030,0.035,0.043]]
  ivezic_fig22_tri=[[0.33,0.66,1.00,1.34,1.66,1.99,2.33],[0.007,0.014,0.022,0.030,0.039,0.047,0.053]]
  ivezic_fig22_sqr=[[0.41,0.84,1.26,1.68,2.09],[0.012,0.020,0.026,0.032,0.039]]

  !p.multi=[0,3,2]
  ;airglow screen
  ;print, ' airglow'
  strfun=DBLARR(mm,2)
  FOR m=0, mm-1 DO BEGIN
     obsid=obsids[m]
     if m gt 3 then imagepr4=imagepr3
     if m gt 2 then imagepr3=imagepr2
     if m gt 1 then imagepr2=imagepr
     if m gt 0 then imagepr=image
     image=mrdfits('airglowscreen_'+obsid+'.fits.gz',0,/silent)

    pixelsize=15.0/3600  ; in degrees
    D_airglow=calculate_structure_function(image,pixelsize)
    D_airglow[*,1]=SQRT(D_airglow[*,1])*100
    IF m EQ 0 THEN begin
       PLOT, adams_fig4[*,0], adams_fig4[*,1], LINESTYLE=0, ytitle='sqrt[D(r)] (%)', xtitle='r (deg)', xr=[0.01,10], yr=[0.2,100], /xstyle,/ystyle, /XLOG, /YLOG
       xyouts,0.1,0.6,'Airglow Screens',/norm
    endif

    OPLOT, D_airglow[*,0], D_airglow[*,1], LINESTYLE=1,color=250
    strfun[m,0]=INTERPOL(D_airglow[*,1],D_airglow[*,0],0.1)
    strfun[m,1]=INTERPOL(D_airglow[*,1],D_airglow[*,0],1.0)
  ENDFOR
  legend,linestyle=[1,0],['PhoSim','Adams & Skrutskie 1996'],box=0,charsize=0.8,/top,/left,color=[250,0]
  ss='Structure Functions'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1H; '+vers
  XYOUTS,0.7,0.98,ss,/normal

  value(nnn,0)=MEAN(strfun[*,0])
  model=INTERPOL(adams_fig4[*,1],adams_fig4[*,0],0.1)
  tolerance_high(nnn,0)=model+0.1
  tolerance_low(nnn,0)=model-0.1
  unit(nnn,0)=' %'
  name(nnn,0)='Airglow Fluctuation at 0.1 deg'
  comparison(nnn,0)='Adams & Skrutskie'

  value(nnn,1)=MEAN(strfun[*,1])
  model=INTERPOL(adams_fig4[*,1],adams_fig4[*,0],1.0)
  tolerance_high(nnn,1)=model+1.0
  tolerance_low(nnn,1)=model-1.0
  unit(nnn,1)=' %'
  name(nnn,1)='Airglow Fluctuation at 1.0 deg'
  comparison(nnn,1)='Adams & Skrutskie'

  ximage=findgen(1024)*15./3600.
  yimage=findgen(1024)*15./3600.
  loadct,3,/silent
  contour,image,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  contour,imagepr,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  contour,imagepr2,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  contour,imagepr3,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  contour,imagepr4,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  loadct,39,/silent

  ;cloud screen
  N=1024
  v=20.0  ; wind velocity (m/s) ~ 3-15 deg/min
  ap=2.5  ; SDSS aperture size (m)
  expt=55 ; exposure time (s)
  pixelscale=256.0 ; (cm)
  wx=ROUND(SQRT(ap*(ap+v*expt))*100/pixelscale)
  wy=wx
  totalcloudmean=0.0
  FOR m=0, mm_all-1 DO BEGIN
    obsid=obsids[m]
    if m gt 3 then cloudmeanpr4=cloudmeanpr3
    if m gt 2 then cloudmeanpr3=cloudmeanpr2
    if m gt 1 then cloudmeanpr2=cloudmeanpr
    if m gt 0 then cloudmeanpr=cloudmean
    SPAWN, "grep cloudmean atmosphere_"+obsid+".pars| awk '{print $3}'", cloudmean
    totalcloudmean=totalcloudmean+cloudmean[0]+cloudmean[1]
  ENDFOR
  totalcloudmean=totalcloudmean/mm_all

  strfun=DBLARR(mm_all,2)
  FOR m=0, mm_all-1 DO BEGIN
     obsid=obsids[m]
     SPAWN, "grep height atmosphere_"+obsid+".pars| awk '{print $3}'", height
     SPAWN, "grep cloudvary atmosphere_"+obsid+".pars| awk '{print $3}'", cloudvary

    pixelsize=pixelscale/((height[1]+height[2])/2*1e5)/!PI*180  ; in degrees
    if m gt 3 then begin
       image1pr4=image1pr3
       image2pr4=image2pr3
       cloudvarypr4=cloudvarypr3
    endif
    if m gt 2 then begin
       image1pr3=image1pr2
       image2pr3=image2pr2
       cloudvarypr3=cloudvarypr2
    endif
    if m gt 1 then begin
       image1pr2=image1
       image2pr2=image2
       cloudvarypr2=cloudvary
    endif
    if m gt 0 then begin
       image1pr=image1
       image2pr=image2
       cloudvarypr=cloudvary
    endif

    image=mrdfits('cloudscreen_'+obsid+'_1.fits.gz',0,/silent)
    image1=[[image,image,image],[image,image,image],[image,image,image]]
    image=mrdfits('cloudscreen_'+obsid+'_2.fits.gz',0,/silent)
    image2=[[image,image,image],[image,image,image],[image,image,image]]
    simage=SMOOTH(cloudvary[0]*image1+cloudvary[1]*image2,[wx,wy])
    D_cloud=calculate_structure_function(simage[N:2*N-1,N:2*N-1],pixelsize)
    D_cloud[*,1]=SQRT(D_cloud[*,1])

    IF m EQ 0 THEN begin
       PLOT, ivezic_fig22_cir[*,0], (totalcloudmean/1.3*ivezic_fig22_cir[*,1])*0.5+ivezic_fig22_cir[*,1]*0.5, PSYM=4, ytitle='sqrt[D(r)] (mag)', xtitle='r (deg)', xr=[0.03,10], yr=[1e-4,0.2], /xstyle,/ystyle, /XLOG, /YLOG
       errplot, ivezic_fig22_cir[*,0], ivezic_fig22_cir[*,1],totalcloudmean/1.3*ivezic_fig22_cir[*,1],width=0.0
       xyouts,0.1,0.6,'Cloud Screens',/norm
    endif

    OPLOT, D_cloud[*,0], D_cloud[*,1], LINESTYLE=1,color=250

    strfun[m,0]=INTERPOL(D_cloud[*,1],D_cloud[*,0],ivezic_fig22_cir[0,0])
    strfun[m,1]=INTERPOL(D_cloud[*,1],D_cloud[*,0],ivezic_fig22_cir[4,0])
  ENDFOR
  ss='Structure Functions'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1H; '+vers
  XYOUTS,0.7,0.98,ss,/normal
  legend,linestyle=[1],['PhoSim'],box=0,charsize=0.8,color=[250],position=[0.04,0.11]
  legend,psym=[4],['Ivezic et al. 2007'],box=0,charsize=0.8,position=[0.04,0.08]

  value(nnn,2)=MEAN(strfun[*,1])
  model=ivezic_fig22_cir[4,1]*totalcloudmean/1.3
  model2=ivezic_fig22_cir[4,1]
  tolerance_high(nnn,2)=model2
  tolerance_low(nnn,2)=model
  unit(nnn,2)=' mag'
  name(nnn,2)='Cloud SF at 2 deg'
  comparison(nnn,2)='Ivezic et al. 2007'

  cloudimage=(((10.^(-0.4*(cloudvary[0]*image1+cloudvary[1]*image2+cloudmean[0]+cloudmean[1])))<1.0)>0.0)
  cloudimagepr=(((10.^(-0.4*(cloudvarypr[0]*image1pr+cloudvarypr[1]*image2pr+cloudmeanpr[0]+cloudmeanpr[1])))<1.0)>0.0)
  cloudimagepr2=(((10.^(-0.4*(cloudvarypr2[0]*image1pr2+cloudvarypr2[1]*image2pr2+cloudmeanpr2[0]+cloudmeanpr2[1])))<1.0)>0.0)
  cloudimagepr3=(((10.^(-0.4*(cloudvarypr3[0]*image1pr3+cloudvarypr3[1]*image2pr3+cloudmeanpr3[0]+cloudmeanpr3[1])))<1.0)>0.0)
  cloudimagepr4=(((10.^(-0.4*(cloudvarypr4[0]*image1pr4+cloudvarypr4[1]*image2pr4+cloudmeanpr4[0]+cloudmeanpr4[1])))<1.0)>0.0)
  ximage=findgen(1024)*2.56
  yimage=findgen(1024)*2.56
  r=byte(255*(1-findgen(256)/255.)+135*(findgen(256)/255.))
  g=byte(255*(1-findgen(256)/255.)+206*(findgen(256)/255.))
  b=byte(255*(1-findgen(256)/255.)+235*(findgen(256)/255.))
  r(0)=0 & g(0)=0 & b(0)=0
  tvlct,r,g,b
  lll=(findgen(11)/10.0)
  contour,cloudimage(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)',levels=lll
  contour,cloudimagepr(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)',levels=lll
  contour,cloudimagepr2(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)',levels=lll
  contour,cloudimagepr3(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)',levels=lll
  contour,cloudimagepr4(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)',levels=lll
  loadct,39,/silent

  ;turbulence screen
  ARCSEC=0.000004841369
  natmospherefile=7
  totalseeing=0.67
  wavelength=0.62
  alt=89.0
  zenith=(90.0-alt)*!pi/180
  zenithfactor=(1/cos(zenith))^0.6
  wavelengthFactor=(wavelength/0.5)^(-0.2)
  seeing=(totalseeing+1e-6)*zenithfactor*wavelengthFactor
  strfun=DBLARR(mm,2)
  FOR m=0, mm-1 DO BEGIN
    obsid=obsids[m]
    SPAWN, "grep seeing atmosphere_"+obsid+".pars | grep -v total | awk '{print $3}'", seefactor
    SPAWN, "grep outer atmosphere_"+obsid+".pars | awk '{print $3}'", outerscales
    seefactor=seefactor*zenithfactor
    density_norm=DBLARR(natmospherefile)
    see_norm=DBLARR(natmospherefile)
    F=5 ; 5x larger than the fine screen
    N=1024*F
    image=DBLARR(N,N)
;    mx=64 & my=64 & nm=64
;    cx=256 & cy=256 & nc=4

    mx=0 & my=0 & nm=128
    cx=0 & cy=0 & nc=16
    lx=0 & ly=0 & nl=2

    outerScaleCorrection=0.0
    FOR i=0L, natmospherefile-1 DO BEGIN
      imagep=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_largep.fits.gz',/remove),0,/silent) ; 512cm
      imagec=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_coarsep.fits.gz',/remove),0,headerc,/silent) ; 64cm
      imagem=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_mediump.fits.gz',/remove),0,/silent) ; 8cm
      imagef=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_finep.fits.gz',/remove),0,/silent)   ; 1cm
      imagetmp=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_coarsex.fits.gz',/remove),0,headerx,/silent)
      density_norm[i]=FXPAR(headerc,'NORM')
      see_norm[i]=FXPAR(headerx,'NORM')
;      print,'a',see_norm[i],density_norm[i]
      outerScaleCorrection=outerScaleCorrection+(seefactor[i]/(zenithfactor*(totalseeing/2.35+1e-6)))^(5./3.)*density_norm[i]/1.0e-6
      image1=REBIN(imagem[mx:mx+nm*F-1,my:my+nm*F-1],N,N)+REBIN(imagec[cx:cx+nc*F-1,cy:cy+nc*F-1],N,N)+REBIN(imagep[lx:lx+nl*F-1,ly:ly+nl*F-1],N,N)
      FOR j=0L, F-1 DO BEGIN
        FOR q=0L, F-1 DO BEGIN
          image1[j*1024:(j+1)*1024-1,q*1024:(q+1)*1024-1]=image1[j*1024:(j+1)*1024-1,q*1024:(q+1)*1024-1]+imagef
        ENDFOR
      ENDFOR
      image=image+image1*density_norm[i]*seefactor[i]/(totalseeing*zenithfactor/2.35)
   ENDFOR
;    outerScaleCorrection=3e-3
;    outerscaleCorrection=4.0
;    print,'b',outerScaleCorrection
;    outerScaleCorrection=1.0/(outerScaleCorrection);
    image=image*sqrt(0.0229*(0.98*(1e-4*wavelength)/(seeing*ARCSEC))^(-5./3.))*4.0e8
    pixelsize=1/100.0  ; in meters
    D_turb=calculate_structure_function(image,pixelsize)
    if m eq 0 then turbimage1=image
    if m eq 1 then turbimage2=image
    if m eq 2 then turbimage3=image
    if m eq 3 then turbimage4=image
    if m eq 4 then turbimage5=image


    IF m EQ 0 THEN BEGIN
      rr=0.01*EXP(ALOG(30.0/0.01)/20.0)^(FINDGEN(20)+1)
      r0=202140*wavelength/1e6/totalseeing
;      print,'r0 ',r0
      y1=6.88*(rr/r0)^(5./3.)
      plot,rr,y1, LINESTYLE=0, ytitle='D(r)', xtitle='r (m)', xr=[0.02,40], yr=[1e-2,1e5], /xstyle,/ystyle, /XLOG, /YLOG
      xyouts,0.1,0.6,'Turbulence Screens',/norm

   endif
    outerscale=total(seefactor^(5./3.)*outerscales)/total(seefactor^(5./3.))
    outerscale=10.0^(1.0+float(m)/(float(mm-1))*1.0)
    k0=(2*!Pi)/outerscale
    u=findgen(100000)/100000.*1e3
    y2=fltarr(N_elements(rr))
    for i=0L,N_elements(rr)-1 do begin
       integrand=u*(1-beselj(u,0))/((u*u+k0*k0*rr(i)*rr(i))^(11./6.))
       f=int_tabulated(u,integrand)
       y2(i)=0.173*(outerscale/r0)^(5./3.)*5./3.*(k0*rr(i))^(5./3.)*f
    endfor
    ;Lucke & Young 2007, AO 256.
    oplot,rr,y2,linestyle=2
    OPLOT, D_turb[*,0], D_turb[*,1],color=250,linestyle=1
    savec(m)=D_turb[0,1]
    strfun[m,0]=INTERPOL(D_turb[*,1],D_turb[*,0],0.02)
    strfun[m,1]=INTERPOL(D_turb[*,1],D_turb[*,0],0.1)
  ENDFOR
  legend,linestyle=[1,2,0],['PhoSim','von Karman (Lucke 2007)','Kolmogorov (Fried 1965)'],box=0,charsize=0.8,/top,/left,color=[250,0,0]
  k=3
  value(nnn,k)=MEAN(strfun[*,1])
  model=6.88*(0.1/r0)^(5.0/3)
  tolerance_high(nnn,k)=model*1.2
  tolerance_low(nnn,k)=model*0.8
  unit(nnn,k)=' rad!U2!N'
  name(nnn,k)='Phase Structure Function at 0.1 m'
  comparison(nnn,k)='Kolmogorov Spectrum'

  k=4
  value(nnn,k)=ALOG(MEAN(strfun[*,1])/MEAN(strfun[*,0]))/ALOG(0.1/0.02)
  tolerance_high(nnn,k)=5.0/3*1.1
  tolerance_low(nnn,k)=5.0/3*0.9
  unit(nnn,k)=' '
  name(nnn,k)='Phase Structure Function Slope Index'
  comparison(nnn,k)='Kolmogorov Spectrum'

  sss=size(turbimage1)
  xxx=findgen(sss(1))*pixelsize
  yyy=findgen(sss(2))*pixelsize
  loadct,5,/silent
  contour,turbimage1,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage2,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage3,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage4,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage5,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  loadct,39,/silent

  ss='Structure Functions'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1H; '+vers
  XYOUTS,0.7,0.98,ss,/normal

  task(nnn,0)='1H Structure Functions'

END
