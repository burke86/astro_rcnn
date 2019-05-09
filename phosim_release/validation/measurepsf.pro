pro measurepsf,data,rms,e1,e2,medx,medy,flux


;chuck's definition
; the FWHM is this function * sqrt(2.0*alog(2.0))/sqrt(!Pi)
; but then the sigma divides by 2.0*sqrt(2.0*alog(2.0))
;rms2=sqrt((total(data))^2/total(data*data))*0.5/sqrt(!Pi)

data2=double(data)  
length = size(data2)
good=where(data2 ne 0)

x=dblarr(N_elements(good))
y=dblarr(N_elements(good))
w=dblarr(N_elements(good))
for i=0L,N_elements(good)-1 do begin
   xx=good(i) mod length(1)
   yy=good(i) / length(1)
   x(i)=xx-(length[1]/2.0-0.5)
   y(i)=yy-(length[2]/2.0-0.5)
   w(i)=data2(xx,yy)
endfor


t1=double(0.0) & for i=0L,N_elements(x)-1 do t1=t1+(x(i)*y(i)*w(i))
t2=double(0.0) & for i=0L,N_elements(x)-1 do t2=t2+(w(i))
t3=double(0.0) & for i=0L,N_elements(x)-1 do t3=t3+(x(i)*w(i))
t4=double(0.0) & for i=0L,N_elements(x)-1 do t4=t4+(y(i)*w(i))
t5=double(0.0) & for i=0L,N_elements(x)-1 do t5=t5+(x(i)*x(i)*w(i))
t6=double(0.0) & for i=0L,N_elements(x)-1 do t6=t6+(y(i)*y(i)*w(i))

covxy=(t1/t2-t3*t4/t2/t2)
resultx=(t5/t2-t3*t3/t2/t2)
resulty=(t6/t2-t4*t4/t2/t2)

medx=median(x)
medy=median(y)
alphax=sqrt(resultx*2.0)
alphay=sqrt(resulty*2.0)
alphaxy=covxy*sqrt(2.0)

for ty=0,100 do begin
weight=(exp(-((x-medx)^2/alphax/alphax-$
             2.0*alphaxy/alphax/alphax/alphay/alphay*(x-medx)*(y-medy)+$
             (y-medy)^2/alphay/alphay)/2.0/(1-(alphaxy/alphax/alphay)^2)))/$
       (2*!Pi*alphax*alphay*sqrt(1-(alphaxy/alphax/alphay)^2))

t1=0. & for i=0L,N_elements(x)-1 do t1=t1+(x(i)*y(i)*w(i)*weight(i))
t2=0. & for i=0L,N_elements(x)-1 do t2=t2+(w(i)*weight(i))
t3=0. & for i=0L,N_elements(x)-1 do t3=t3+(x(i)*w(i)*weight(i))
t4=0. & for i=0L,N_elements(x)-1 do t4=t4+(y(i)*w(i)*weight(i))
t5=0. & for i=0L,N_elements(x)-1 do t5=t5+(x(i)*x(i)*w(i)*weight(i))
t6=0. & for i=0L,N_elements(x)-1 do t6=t6+(y(i)*y(i)*w(i)*weight(i))
t7=0. & for i=0L,N_elements(x)-1 do t7=t7+(w(i))
t8=0. & for i=0L,N_elements(x)-1 do t8=t8+(weight(i)*weight(i))

covxy=(t1/t2-t3*t4/t2/t2)
resultx=(t5/t2-t3*t3/t2/t2)
resulty=(t6/t2-t4*t4/t2/t2)
medx=(t3/t2)
medy=(t4/t2)
flux=t7

rms=sqrt(resultx+resulty)

e1=(resultx-resulty)/(resultx+resulty)
e2=(2.0*covxy)/(resultx+resulty)

ellip=sqrt(e1^2+e2^2)
pa=0.5*atan(e2,e1)

;print,ellip,pa*180/!Pi,e1,e2,medx,medy,rms,alphax,alphay
;print,ty,abs(alphax-sqrt(2.0*resultx)),abs(alphay-sqrt(2.0*resulty)),abs(alphaxy-2.0*covxy)

IF (abs(alphax-sqrt(2.0*resultx)) lt 1e-6 and abs(alphay-sqrt(2.0*resulty)) lt 1e-6 and abs(alphaxy-2.0*covxy) lt 1e-6) THEN break
     alphax=sqrt((resultx*2.0)>1.0)
     alphay=sqrt((resulty*2.0)>1.0)
     alphaxy=covxy*2.0
endfor

END
