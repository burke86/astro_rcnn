;;
;; @package phosim
;; @file validation_3D.pro
;; @brief validation task 3D
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

pro validation_3D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3D'
  !p.multi=[0,1,2]


  tolerance_low(nnn,0)=2.9*0.0
  tolerance_high(nnn,0)=2.9*2.0




  nn=20
  w=fltarr(20)
  x=fltarr(20)
  y=fltarr(20)
  z=fltarr(20)
  v=fltarr(20,nn)
  sigmax=fltarr(20)
  sigmay=fltarr(20)
  er=fltarr(20)

  data=mrdfits('lsst_e_3300_f2_R22_S11_E000.fits.gz')

  mag=14.5
  for ii=0,19 do begin
     mag=mag+0.5
     for jj=0,nn-1 do begin

        x0=296+180*jj
        y0=3742-180*ii

        datas=data((x0-50):(x0+50),(y0-50):(y0+50))

        measurepsf,datas,a,b,c,d,e,f

        x(ii)=mag
        y(ii)=y(ii)+a*2.35*0.2
;        print,a
        z(ii)=z(ii)+max(datas)
        w(ii)=w(ii)+f
        v(ii,jj)=a*2.35*0.2
        sigmax(ii)=sigmax(ii)+a*2.35*0.2*(sqrt((1.0+b)/2.0))*sqrt(2.0)
        sigmay(ii)=sigmay(ii)+a*2.35*0.2*(sqrt((1.0-b)/2.0))*sqrt(2.0)
     endfor
     result=moment(v(ii,*))
     er(ii)=sqrt(result(1))
  endfor


        g=where(z/float(nn) lt 5000.0)
        y1=median(y(g))/float(nn)
        xx=findgen(1000)/1000.*1e5
        scale=3.3
        range=0.05

        plot,x,y/float(nn),psym=4,/xstyle,/ystyle,yr=[y1-range,y1+range],xtitle='Magnitude',ytitle='FWHM (arcseconds)'
        errplot,x,y/float(nn)-er,y/float(nn)+er,width=0.0

        oplot,x,sigmax/float(nn),psym=5,color=80,symsize=0.5
        oplot,x,sigmay/float(nn),psym=5,color=250,symsize=0.5

        mm=findgen(1000)/1000.*10.0+15.

        vv=z(where(x eq 20))/nn
        pp=vv(0)*10.^(0.4*(20.0-mm))


        oplot,mm,y1*(1.0+0.035/scale*(pp/1e5)),linestyle=2
        oplot,mm,y1*(1.0+0.03/scale*(pp/1e5)),linestyle=2,color=80
        oplot,mm,y1*(1.0+0.04/scale*(pp/1e5)),linestyle=2,color=250
        oplot,mm,y1*(1.0+0.0*(pp/1e5)),linestyle=2


        plot,z/float(nn),y/float(nn),psym=4,/ystyle,yr=[y1-range,y1+range],xtitle='Maximum electrons per pixel',ytitle='FWHM (arcseconds)'
        errplot,z/float(nn),y/float(nn)-er,y/float(nn)+er,width=0.0

        oplot,z/float(nn),sigmax/float(nn),psym=5,color=80,symsize=0.5
        oplot,z/float(nn),sigmay/float(nn),psym=5,color=250,symsize=0.5

        oplot,xx,y1*(1.0+0.035/scale*(xx/1e5)),linestyle=2
        oplot,xx,y1*(1.0+0.03/scale*(xx/1e5)),linestyle=2,color=80
        oplot,xx,y1*(1.0+0.04/scale*(xx/1e5)),linestyle=2,color=250
        oplot,xx,y1*(1.0+0.0*(xx/1e5)),linestyle=2



  ss='Size vs. Intensity'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3D; '+vers
  xyouts,0.7,0.98,ss,/normal


  task(nnn,0)='3D PSF size vs. intensity'
;  print,flux
;  print,variance


END
