function zl_coeffs_str, Nstruct
 
  common zl_coeffs_str, defined
 
  if N_elements( defined ) NE 1 then begin
 
       zl_coeffs = { zl_coeffs                             , $
                    PIXEL_NO: 0L                           , $
                    ZL_COEFFS: fltarr( 10, 11 )                }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {zl_coeffs}, Nstruct )
end
