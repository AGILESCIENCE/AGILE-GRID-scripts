pro getsource_name , l , b 
resolve_all
; to convert from Galactic to RA-DEC
euler , l , b , ra , dec , 2

; to convert from degrees to hms
radec , ra, dec, rah, ram, ras , dech , decm, decs

if ras gt 30 then ram += 1
if decs gt 30 then decm += 1
if dech lt 0. then sim='-' else sim='+'
;help,rah,ram,dech , decm
str = string(rah , ram , sim , abs(dech) , decm , format='(I02,I02,a,I02,I02)')
print
print,str
forprint , str , /nocomment , TEXTOUT = 'source_name' , /silent
;stop

end

