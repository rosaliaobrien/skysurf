; ======================================================================
;
; Function ConfInt
;
; Calculate uncertainty on model parameters which were fit by
; minimization of chisquare.
;
; Written By:  BA Franz, ARC, 2/95
;
; Calling Sequence:
;     sigma = ConfInt(covar,nfree,chisqr,confidence=confidence) 
;
; Inputs:
;     covar  - m x m covariance matrix (or hessian) of fitted parameters 
;     nfree  - number of degrees of freedom (e.g., n - m)
;     chisqr - reduced chisquare of the fit (i.e., per degree of freedom)
;
; Optional Keyword Inputs:
;     confidence - desired confidence level (def=68.0)
;     hessian    - set this to 1 if hessian was supplied in place 
;                  the covariance matrix
;
; Output:
;     sigma - m-element vector of parameter uncertainties
;
; Description of Algorithm (by T. J. Sodroski)
; The parameter uncertainties are obtained from an m-dimensional
; ellipsoidal confidence region in m-parameter space.  The length
; of each principal axis i of the confidence region is given by
; sqrt(2E/Li), where E is the largest allowable deviation of
; chisquare from its optimal value for a given level of significance
; (as derived from the F-distribution), and Li is the associated 
; eigenvalue of the hessian matrix for parameter i.  The components
; of each principal axis, Pi, in parameter space are then obtained
; by: [Pi] = [U] * [Wi], where [U] is the unitary transformation
; matrix whose columns are the normalized eigenvectors of the
; Hessian, and [Wi] is the principal axis vector in the coordinate
; system defined by the axes of the ellipsoid.  The uncertainty in
; each parameter is then given by the maximum component of the 
; principal axis vectors along the parameter axis.
;
; References: 
; Bard, Y.  1974, Nonlinear Parameter Estimation, (Orlando: Academic)
; Sodroski, T. J.  1988, Ph.D. thesis, University of Maryland
;
; ======================================================================
function ConfInt,covar,nfree,chisqr,confidence=conf,hessian=hessian 

if (n_elements(conf) eq 0) then conf=68.0

;
; Get number of parameters from size of COVAR
;
s = size(covar)
if ( (s(0) ne 2) or (s(1) ne s(2)) ) then begin
    message,'Covariance matrix must be a square array',/info
   return,-1
endif else nterms = s(1)
diag = indgen(nterms)*(nterms+1)         ; Subscripts of diagonal elements

;
; Get Hessian from COVAR unless Hessian was passed instead
;
if (not keyword_set(Hessian)) then begin
    hess = invert(covar,status)
    case (status) of
      0:
      1: message,'Inversion of covariance matrix failed',/info
      2: message,'Inversion of covariance matrix inaccurate',/info
    endcase
endif else hess = covar

;
; Use f-distribution to determine ep, the largest allowable deviation
; of chisquare from its optimal value within the confidence region.
;
fd = f_test(1.-conf/100.,nterms,nfree)
ep = 2.0 * nterms * chisqr * fd

;
; Compute confidence intervals for each paramter versus other parameters
;
temp1  = hess                             ; The unnormalized hessian
nr_tred2,temp1,temp2,temp3,/double        ; Reduce to tridiagonal form
nr_tqli,temp2,temp3,temp1,/double         ; Compute:
eig    = temp2                            ;   Eigenevalues
eigvec = temp1                            ;   Eigenvectors
paxis_len = dblarr(nterms,nterms)         ; Length of principal axes of
paxis_len(diag) = sqrt((ep/eig) > 1e-30)  ;   confidence region
theta  = abs(paxis_len#transpose(eigvec)) ; Projection of principal axes of 
                                          ;   conf regions onto param axes
;
; Compute parameter uncertainties as max components of principal axes
; in parameter space.
;
sigma = fltarr(nterms)                    
for i=0,nterms-1 do sigma(i) = max(theta(*,i))            

return,sigma
end
