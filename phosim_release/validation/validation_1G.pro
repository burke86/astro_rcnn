;;
;; @package phosim
;; @file validation_1G.pro
;; @brief validation task 1G
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

FUNCTION star_centroid, image
  dim=SIZE(image)
  maximg=MAX(image,I)
  IX=I MOD dim(1)
  IY=I/dim(1)
  xlo=IX-100
  xhi=IX+100
  if xlo LT 0 then xlo=0
  if xhi GE dim(1) then xhi=dim(1)-1
  ylo=IY-100
  yhi=IY+100
  if ylo LT 0 then ylo=0
  if yhi GE dim(2) then yhi=dim(2)-1
  subdata=image[xlo:xhi,ylo:yhi]
  measurepsf,subdata,rms,e1,e2,medx,medy,flux
  medx=medx+(xlo+xhi)/2.0
  medy=medy+(ylo+yhi)/2.0
  RETURN, [medx,medy]
end

FUNCTION airRefraction, wavelength, temperature, pressure, water_pressure
  n=64.328+29498.1/(146-1/wavelength/wavelength)+255.4/(41-1/wavelength/wavelength)
  n=n*pressure*(1+(1.049-0.0157*temperature)*1e-6*pressure)/720.883/(1+0.003661*temperature)
  n=n-((0.0624-0.000680/wavelength/wavelength)/(1+0.003661*temperature)*water_pressure)
  n=1e-6*n+1
  RETURN, n
end

FUNCTION airRefraction_owen, wavelength, temperature, pressure, water_pressure
  sigma=1.0/wavelength
  temp=temperature+273.15
  ps=pressure/760.00*1013.25
  pw=water_pressure/760.00*1013.25
  dw=(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/temp-710.792/temp/temp+7.75141e4/temp/temp/temp))*pw/temp
  ds=(1+ps*(57.90e-8-9.325e-4/temp+0.25844/temp/temp))*ps/temp
  n=(2371.34+683939.7/(130.0-sigma^2)+4547.3/(38.9-sigma^2))*ds
  n=n+(6478.31-58.058*sigma^2-0.71150*sigma^4+0.08851*sigma^6)*dw
  n=1e-8*n+1
  RETURN, n
end

PRO validation_1G,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison
  print,'Task 1G'
  !p.multi=0
  chipids=['R02_S11','R20_S11','R24_S11','R42_S11','R22_S11']
  alt=89-findgen(9)*10.3
  zen=90-alt
  dx=DBLARR(5,9)
  dy=DBLARR(5,9)
  for j=0,8 do begin
     for i=0, 4 do begin
        print,i,j
        if alt(j) lt 10 then ss='0'+string(floor(alt(j)))
        if alt(j) ge 10 then ss=string(floor(alt(j)))
      image=mrdfits(strcompress('lsst_e_18'+ss+'_f1_'+chipids[i]+'_E000.fits.gz',/remove),0,/silent)
      xy_red=star_centroid(image)
      image=mrdfits(strcompress('lsst_e_19'+ss+'_f1_'+chipids[i]+'_E000.fits.gz',/remove),0,/silent)
      xy_blue=star_centroid(image)
      dx[i,j]=xy_blue[0]-xy_red[0]
      dy[i,j]=xy_blue[1]-xy_red[1]
      print,xy_blue
      print,xy_red
    endfor
  endfor
  dr=SQRT(dx^2+dy^2)*0.2  ; in arcsec
  PLOT, asin(sqrt((sin(zen*!pi/180.))^2+(sin(1.41*!Pi/180.))^2))*180./!Pi, dr[0,*], psym=1, ytitle='Separation (arcsec)', xtitle='Zenith angle (deg)', yr=[-1,10.0], xr=[-1,90],/xstyle,/ystyle
  OPLOT, abs(zen+1.41), dr[1,*], psym=2
  OPLOT, abs(zen-1.41), dr[2,*], psym=4
  OPLOT, asin(sqrt((sin(zen*!pi/180.))^2+(sin(1.41*!Pi/180.))^2))*180./!Pi, dr[3,*], psym=5
  OPLOT, zen, dr[4,*], psym=6
  print,dx
  print,dy


  water_pressure=8
  pressure=520
  temperature=20
  z=FINDGEN(90)
  model=TAN(z*!pi/180)*206265*(airRefraction(0.42,temperature,pressure,water_pressure)-airRefraction(0.52,temperature,pressure,water_pressure))
  OPLOT, z, model, linestyle=0

  model=TAN(z*!pi/180)*206265*(airRefraction_owen(0.42,temperature,pressure,water_pressure)-airRefraction_owen(0.52,temperature,pressure,water_pressure))
  OPLOT, z, model, linestyle=1,color=250

  height=2600.0
  latitude=-30.66*!Pi/180.0
  beta=0.001254*(273.15+temperature)/273.15
  kappa=1.0+0.005302*(sin(latitude))^2-0.00000583*(sin(2.0*latitude))^2-0.000000315*height
  gamma1=airRefraction_owen(0.42,temperature,pressure,water_pressure)-1.0
  gamma2=airRefraction_owen(0.52,temperature,pressure,water_pressure)-1.0
  model=(kappa*gamma1*(1-beta)*TAN(z*!pi/180)-kappa*gamma1*(beta-gamma1/2.0)*(tan(z*!pi/180.))^3)*206265.0-$
        (kappa*gamma2*(1-beta)*TAN(z*!pi/180)-kappa*gamma2*(beta-gamma2/2.0)*(tan(z*!pi/180.))^3)*206265.0
  OPLOT, z, model, linestyle=2,color=250

  XYOUTS, 2, 4.0, 'Differential Positions for Monochromatic Sources at 420nm and 520nm'
  legend,linestyle=[0,1,2], ['Filippenko 1982: Plane-parallel approx. + n(!4k!3) of Barrell 1951 ',$
                             'Stone 1996: Plane-parallel approx. + n(!4k!3) of Owens 1967 ',$
                             'Stone 1996: Curved Earth 2-term atm approx. model + n(!4k!3) of Owens 1967 '],color=[0,250,250]
  legend,psym=[6,1,5,2,4], ['PhoSim, at field center','PhoSim, +1.41 perp to zen','PhoSim, -1.41 perp to zen','PhoSim, +1.41 par to zen','PhoSim, -1.41 par to zen'],/bottom,/right


  ss='Differential Chromatic Refraction'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1G; '+vers
  XYOUTS,0.7,0.98,ss,/normal

  k=0
  for j=1, 3, 2 do begin
    for i=0, 4 do begin
      value(nnn,k)=dr[i,j]
      if i eq 0 then z=asin(sqrt((sin(zen(j)*!pi/180.))^2+(sin(1.41*!Pi/180.))^2))*180./!Pi
      if i eq 1 then z=abs(zen[j]+1.41)
      if i eq 2 then z=abs(zen[j]-1.41)
      if i eq 3 then z=asin(sqrt((sin(zen(j)*!pi/180.))^2+(sin(1.41*!Pi/180.))^2))*180./!Pi
      if i eq 4 then z=zen[j]
      model=(kappa*gamma1*(1-beta)*TAN(z*!pi/180)-kappa*gamma1*(beta-gamma1/2.0)*(tan(z*!pi/180.))^3)*206265.0-$
        (kappa*gamma2*(1-beta)*TAN(z*!pi/180)-kappa*gamma2*(beta-gamma2/2.0)*(tan(z*!pi/180.))^3)*206265.0
      tolerance_high(nnn,k)=model*1.1
      tolerance_low(nnn,k)=model*0.9
      unit(nnn,k)=' arcsec'
      name(nnn,k)='Differential Refraction'
      comparison(nnn,k)='Stone 1996'
      k=k+1
    endfor
  endfor

  !p.multi=0
  chipids=['R02_S11','R20_S11','R24_S11','R42_S11']
     plot,[-2,2],[-2,2],/nodata
  for j=0, 8 do begin
        ss=0.5/max(sqrt(dx(*)*dx(*)+dy(*)*dy(*)))
        arrow,-1.41,0,-1.41+ss*dx(0,j),0+ss*dy(0,j),/data,hsize=-0.2
        arrow,0,-1.41,0.0+ss*dx(1,j),-1.41+ss*dy(1,j),/data,hsize=-0.2
        arrow,0,1.41,0.0+ss*dx(2,j),1.41+ss*dy(2,j),/data,hsize=-0.2
        arrow,1.41,0,1.41+ss*dx(3,j),0+ss*dy(3,j),/data,hsize=-0.2
        arrow,0.0,0,0.0+ss*dx(4,j),0+ss*dy(4,j),/data,hsize=-0.2
  endfor
  !p.multi=0

  task(nnn,0)='1G Differential Chromatic Refraction'
end
