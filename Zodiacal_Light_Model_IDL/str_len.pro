function str_len,str

tag  = tag_names(str)
ntag = n_elements(tag)

nbytes = 0
for i=0,ntag-1 do begin
   status = execute('s = size(str(0).'+tag(i)+')')
   n = n_elements(s)-1
   type = s(n-1)
   nelm = s(n)
   
   case type of
     0 : nbytes = 0
     1 : nbytes = nbytes + nelm*1
     2 : nbytes = nbytes + nelm*2
     3 : nbytes = nbytes + nelm*4
     4 : nbytes = nbytes + nelm*4
     5 : nbytes = nbytes + nelm*8
     6 : nbytes = nbytes + nelm*8
   endcase

endfor

if (nbytes eq 0) then message,'Zero length, nested, or otherwise ambiguous structure.',/info 

return,nbytes
end
