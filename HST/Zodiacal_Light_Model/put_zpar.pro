function put_zpar, aend_in, PF1_C0, PF1_C1, PF1_Ka, A1, det1=det1, hg3=hg3

; creates new aend=zpars arrat with new phase function and albedo for channel 1

; INPUT
; aend_in = original aend array (e.g. from zpars.xdr)
; PF1_C0 = phase function parameter 1
; PF1_C1 = phase function parameter 2
; PF1_Ka = phase function parameter 3
; A1     = albedo (used for all geometric components)
; det1 = detector number to replace (default: det = 0, otherwise det=1)
; hg3 = [g1,g2,g3,w1,w2,w3] for weighted sum of HG phase functions

if (keyword_set(det1)) then offset = 1 else offset=0
aend = aend_in
aend[1+3*offset:3+3*offset] = [PF1_C0,PF1_C1,PF1_Ka]
aend[[33,59,85,111,137,170]+offset] = A1

if n_elements(hg3) ge 4 then aend[183:183+n_elements(hg3)-1] = hg3

return,aend
end
