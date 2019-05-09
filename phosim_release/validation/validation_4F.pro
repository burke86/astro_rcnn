;;@package phosim
;; @file validation_temp.pro
;; @brief validation 4E, F
;;
;; @brief Created by:
;; @author John R. Peterson (Purdue)
;;
;; @brief Modified by:
;; @author Caleb Remocaldo (Purdue)
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.

pro validation_4f, nnn, vers, value, tolerance_low, tolerance_high, task, name, unit, comparison
print, 'Task 4F'
  
!p.multi = [0, 2, 2]
aa = findgen(33)/32.*2*!Pi
usersym, cos(aa), sin(aa), /fill

spawn, '/bin/ls wiyn_odi_e_* > temp.txt'
;spawn,  '/bin/ls lsst_e_* > temp.txt'
readcol, 'temp.txt', fileList, format = '(A)', /silent

;Array initialization
alt = dblarr(N_elements(fileList))
salt = dblarr(N_elements(fileList))
dist2moon = dblarr(N_elements(fileList))
mphase = dblarr(N_elements(fileList))
intensity = dblarr(N_elements(fileList))
filter = dblarr(N_elements(fileList))
mjd = dblarr(N_elements(fileList))
malt = dblarr(N_elements(fileList))
zen = dblarr(N_elements(fileList))

prff = -1
for i = 0L, N_elements(fileList)-1 do begin
   data = mrdfits(fileList(i), /silent)
   fitsHeader = headfits(fileList(i), /silent)
   catNum = strmid(fileList(i),  11,  3)
   ;catNum = strmid(fileList(i),  7,  3) ;Catalog number
   
   zen(i) = sxpar(fitsHeader, 'ZENITH')
   alt(i) = 90.00 - zen(i)
   salt(i) = sxpar(fitsHeader, 'SUNALT') * (180/!Pi)
   dist2moon(i) = sxpar(fitsHeader, 'MOONDIS') * (180/!Pi)
   mphase(i) = sxpar(fitsHeader, 'MOONPHS')
   filter(i) = sxpar(fitsHeader, 'FILTNM')
   malt(i) = sxpar(fitsHeader, 'MOONALT') * (180/!Pi)
   mjd(i) = sxpar(fitsHeader, 'MJD-OBS')

   ;should be counts per second per pixel
   intensity(i) = ((mean(data(70:3950, 30:3610))*100.0/15.0)>1.0)*15.0
endfor

xLabel = ['Altitude', 'Sun Altitude', 'Distance to Moon', 'Lunar Phase']
for i = 0,  3 do begin
if i eq 0 then x = alt
if i eq 1 then x = salt
if i eq 2 then x = dist2moon
if i eq 3 then x = mphase

plot, x, intensity, /ylog, psym = 8, xtitle = xLabel(i), ytitle = 'Counts/exp/pixel', /nodata, xr = [min(x)-(max(x)-min(x))*0.1, max(x)+(max(x)-min(x))*0.1], /xstyle
for f = 0,  5 do begin
   g = where(filter eq f)
   if N_elements(g) gt 1 then begin
      oplot, x(g), intensity(g), psym = 8, color = 35+f*42, symsize = 0.5
   endif else begin
      if g(0) ne -1 then begin
         oplot, [x(g)], [intensity(g)], psym = 8, color = 35+f*42, symsize = 0.5
      endif
   endelse
endfor
endfor

readcol, 'yoachim.txt', yalt, yaz, ymjd, yu, yg, yr, yi, yz, yyy, format = '(D, D, D, D, D, D, D, D, D)'
lookingUp = where(yalt eq 90.0)
ymjd = ymjd(lookingUp)
yu = yu(lookingUp)
yg = yg(lookingUp)
yr = yr(lookingUp)
yi = yi(lookingUp)
yz = yz(lookingUp)
yyy = yyy(lookingUp)
yu = 1.0/25.0*0.0574*0.5*5455.*10.^(-0.4*(yu-25.0))
yg = 1.0/25.0*0.1735*0.5*5455.*10.^(-0.4*(yg-25.0))
yr = 1.0/25.0*0.1502*0.5*5455.*10.^(-0.4*(yr-25.0))
yi = 1.0/25.0*0.1272*0.5*5455.*10.^(-0.4*(yi-25.0))
yz = 1.0/25.0*0.0872*0.5*5455.*10.^(-0.4*(yz-25.0))
yyy = 1.0/25.0*0.0199*0.5*5455.*10.^(-0.4*(yyy-25.0))

!p.multi = 0

readcol, 'canon.txt', ca, cb, cc, cd, ce, format = '(D, F, F, F, F)', /silent
readcol, 'photodiode.txt', pa, pb, pc, pd, format = '(D, D, D, D)', /silent
minpa = floor(min(pa))
pa = pa-minpa
ca = ca-minpa

goto, skipcalib

skipcalib:

miny = min([intensity])/10.0
maxy = max([intensity])*10.0

;determine nights
sortedPa = sort(pa)
pas = pa(sortedPa)
nightArr = pas(0)
night = 1L
for i = 0L, N_elements(pas)-2 do begin
   if pas(i+1)-pas(i) gt 0.25 then begin
      nightArr = [nightArr, pas(i+1)]
      if night eq 1 then dayArr = pas(i) else dayArr = [dayArr, pas(i)]
      night = night+1L
   endif
endfor
dayArr = [dayArr, max(pas)]


;fit to data
for jj = 0, 6 do begin
   if jj eq 0 then begin & ff = 0 & ww = cb & uu = ca & endif
   if jj eq 1 then begin & ff = 1 & ww = cc & uu = ca & endif
   if jj eq 2 then begin & ff = 2 & ww = cd & uu = ca & endif
   if jj eq 3 then begin & ff = 1 & ww = ce & uu = ca & endif
   if jj eq 4 then begin & ff = 2 & ww = pb & uu = pa & endif
   if jj eq 5 then begin & ff = 5 & ww = pc & uu = pa & endif
   if jj eq 6 then begin & ff = 4 & ww = pd & uu = pa & endif

      zz = dblarr(N_elements(uu))
      yy = dblarr(N_elements(uu))
      xx = dblarr(N_elements(uu))

for ii = 0, N_elements(xx)-1 do begin
   xx(ii) = uu(ii)
   yy(ii) = ww(ii)
   g1 = where(filter eq ff)
   g = where(abs(uu(ii)-(mjd(g1)-minpa)) eq min(abs([uu(ii)-(mjd(g1)-minpa)])))
   zz(ii) = intensity(g1(g(0)))*0.5
   g = where(abs(uu(ii)-(ymjd-minpa)) eq min(abs([uu(ii)-(ymjd-minpa)])))
   if ff eq 0 then zz(ii) = zz(ii)+0.5*yu(g(0))
   if ff eq 1 then zz(ii) = zz(ii)+0.5*yg(g(0))
   if ff eq 2 then zz(ii) = zz(ii)+0.5*yr(g(0))
   if ff eq 3 then zz(ii) = zz(ii)+0.5*yi(g(0))
   if ff eq 4 then zz(ii) = zz(ii)+0.5*yz(g(0))
   if ff eq 5 then zz(ii) = zz(ii)+0.5*yyy(g(0))

endfor
plot, xx, yy, psym = 5, /ylog, yr = [miny, maxy]
oplot, xx, zz, psym = 4

for ii = 0, N_elements(nightArr)-1 do begin
   g = where(uu ge nightArr(ii) and uu le dayArr(ii))
   g2 = where(uu ge round(nightArr(ii)) and uu le (round(dayArr(ii)-0.4)+0.4))
   frac = 1.0
   vv = regress(zz(g2), yy(g2), const = const, chisq = chi)
   aaa = vv(0) & bbb = const
   goto, skiptrim

   redo:
   frac = frac-0.01
   resid = abs(yy(g2)-(vv(0)*zz(g2)+const))
   ss = sort(resid)
   ss = ss(0:round((N_elements(ss)-1)*frac))
   vv = regress(zz(g2(ss)), yy(g2(ss)), const = const, chisq = chi)
   if chi/N_elements(g) gt 2.0 then goto, redo

   skiptrim:
     ; print, frac, ii, const, 1.0/vv(0), chi/N_elements(g)

    plot, zz(g), ww(g), psym = 4
    xxx = findgen(1000)/1000.*(max(zz(g))-min(zz(g)))+min(zz(g))
    oplot, xxx, const+vv(0)*xxx
    oplot, xxx, bbb+aaa*xxx, linestyle = 1
    ww(g) = (ww(g)-const)/vv(0)

endfor

   if jj eq 0 then begin & cb = ww & endif
   if jj eq 1 then begin & cc = ww & endif
   if jj eq 2 then begin & cd = ww & endif
   if jj eq 3 then begin & ce = ww & endif
   if jj eq 4 then begin & pb = ww & endif
   if jj eq 5 then begin & pc = ww & endif
   if jj eq 6 then begin & pd = ww & endif

endfor

!p.multi = [0, 1, 2]

for ii = 0, N_elements(dayArr)-1 do begin
    g = where(mjd-minpa ge nightArr(ii) and mjd-minpa le dayArr(ii))
    miny = min([intensity(g)])/3.0
    maxy = max([intensity(g)])*3.0
;print, minpa
plot, [nightArr(ii), dayArr(ii)], [10, 5e4], /nodata, /ylog, xtitle = '!4D!3MJD', ytitle = 'Counts/exp/pixel', /xstyle, /ystyle

g = where(pa ge nightArr(ii) and pa le dayArr(ii))
if n_elements(g) gt 1 then begin
   oplot, pa(g), pb(g), color = 130, psym = 3
   oplot, pa(g), pc(g), color = 250, psym = 3
   oplot, pa(g), pd(g), color = 210, psym = 3
endif

g = where(ca ge nightArr(ii) and ca le dayArr(ii))
if n_elements(g) gt 1 then begin
   oplot, ca(g), cb(g), color = 50, psym = 3
   oplot, ca(g), cc(g), color = 110, psym = 3
   oplot, ca(g), cd(g), color = 130, psym = 3
   oplot, ca(g), ce(g), color = 90, psym = 3
endif

terr = 0.0
terrn = 0.0
terrl = fltarr(3)
terrln = fltarr(3)
terrw = fltarr(6)
terrwn = fltarr(6)

for ff = 0, 5 do begin
   g = where(filter eq ff and mjd-minpa ge nightArr(ii) and mjd-minpa le dayArr(ii))
   ss = sort(mjd(g))
   g = g(ss)
   if n_elements(g) gt 1 then oplot, mjd(g)-minpa, intensity(g), color = 50+ff*40, psym = 8, symsize = 0.75
   if ff eq 0 then model = yu
   if ff eq 1 then model = yg
   if ff eq 2 then model = yr
   if ff eq 3 then model = yi
   if ff eq 4 then model = yz
   if ff eq 5 then model = yyy
   oplot, ymjd-minpa, model, color = 50+ff*40, psym = 4, symsize = 0.75
   mint = round(min(mjd(g)-minpa)*10.0)/10.0
   maxt = round(max(mjd(g)-minpa)*10.0)/10.0
   for t0 = mint+0.1, maxt, 0.1 do begin
      g1 = where(mjd(g)-minpa gt t0-0.05 and mjd(g)-minpa le t0+0.05)
      g2 = where(ymjd-minpa gt t0-0.05 and ymjd-minpa le t0+0.05)
      error = round((-2.5*alog10(mean(intensity(g(g1))))+2.5*alog10(mean(model(g2))))*100.0)/100.0
      terrl(terrn mod 3) = terrl(terrn mod 3)+error
      terrln(terrn mod 3) = terrln(terrn mod 3)+1.0
      terr = terr+error
      terrn = terrn+1.0
      terrw(ff) = terrw(ff)+error
      terrwn(ff) = terrwn(ff)+1.0
      print, t0, ff, round(mean(intensity(g(g1)))), round(mean(model(g2))), error
   endfor

if n_elements(g) eq 1 and g(0) ne -1 then oplot, [mjd(g)-minpa], [intensity(g)], color = 50+ff*40, psym = 8, symsize = 0.5
endfor
print, 'Total          =  ', terr/terrn
print, 'Total by time  =  ', terrl/terrln
print, 'Total by band  =  ', terrw/terrwn

endfor
END
