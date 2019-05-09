PRO MAKE_BIG_TIFF

;+++++++++++++++++++++++++
;
; NAME: MAKE_BIG_TIFF
;
; PURPOSE: Combines all the focalplane files
;          to make a big tiff of the 
;          focalplane
;
; PROCEDURE CALLS: num2str
;
; INPUTS: The files that you are joining
;
; OUTPUTS: One big tiff called focalplane.tiff 
;
; WRITTEN BY: Alan Meert
;
; FOR: John Peterson
;      Purdue University
;
; Date: July 30 2008
;
;--------------------------

size = 1024 ; The size to reduce the fits images to. Can't be much bigger or will crash program
scale = 1

second_picture = fltarr((size+1)*15, (size+1)*15)
final_picture = bytarr(3,(size+1)*15, (size+1)*15)

y = 0
x = 0

FOR y = -7, 7 DO BEGIN
    FOR x = -7, 7 DO BEGIN
        IF ~( ( (y EQ -7) || (y EQ 7)) && ((x LT -4) || (x GT 4))) THEN BEGIN
            IF ~( ( (y EQ -6) || (y EQ 6)) && ((x LT -5) || (x GT 5))) THEN BEGIN
                IF ~(((y EQ -5) || (y EQ 5)) && ((x EQ -7) || (x EQ 7))) THEN BEGIN
                    
                    in_file = 'focalplane_' + num2str(x) + '_' + num2str(y) + '.fits'

                    IF FILE_TEST(in_file) THEN BEGIN
                        
                        print, 'Opening ', in_file
                        data = mrdfits(in_file)
                        data = rebin(data, size,size)
                        
                        
                        x_offset = (x+7)*(size+1)
                        y_offset = abs(y-8)*(size+1);makes the picture be correct rather than flipped
                        
                        FOR countx = 0, (size - 1) DO BEGIN
                            FOR county = 0, (size - 1) DO BEGIN
 second_picture[countx + x_offset, y_offset-(county+1)]=data[countx, county] 
                            ENDFOR
                        ENDFOR
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
    ENDFOR
ENDFOR

second_picture = alog(second_picture>1) ;puts picture into log scale

second_picture=byte((second_picture/max(second_picture)*255*scale)<255) ;normalizes picture to brightest pixel

loadct,5 ;normal black blue red colortable
tvlct,red,green,blue,/get

;saving new picture as tiff

for i=0L,15*(size+1)-1 do begin
if i mod 100 eq 0 then print,i
for j=0L,15*(size+1)-1 do begin
    final_picture(0,i,j)=red(second_picture(i,j))
    final_picture(1,i,j)=green(second_picture(i,j))
    final_picture(2,i,j)=blue(second_picture(i,j)) 
endfor
endfor

;Makes corners white rather than black
FOR y=-7, 7 DO BEGIN
    FOR x = -7, 7 DO BEGIN
        IF (((y EQ -7) || (y EQ 7)) && ((x LT -4) || (x GT 4))) $
          || (((y EQ -6) || (y EQ 6)) && ((x LT -5) || (x GT 5))) $
          || (((y EQ -5) || (y EQ 5)) && ((x EQ -7) || (x EQ 7))) THEN BEGIN
            
            x_offset = (x+7)*(size+1)
            y_offset = (y+7)*(size+1)
            
            FOR countx = 0, (size - 1) DO BEGIN
                FOR county = 0, (size - 1) DO BEGIN
                    final_picture(0,countx+x_offset,county+y_offset)=red(255)
                    final_picture(1,countx+x_offset,county+y_offset)=green(255)
                    final_picture(2,countx+x_offset,county+ y_offset)=blue(255)
                ENDFOR
            ENDFOR
        ENDIF
    ENDFOR
ENDFOR




write_tiff,'focalplane.tiff',final_picture

END


