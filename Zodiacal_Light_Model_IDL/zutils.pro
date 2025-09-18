function strnumber,st,val
;+
; NAME:
;   STRNUMBER
; PURPOSE:
;   Function to determine if a string is a valid numeric value.
; CALLING SEQUENCE:
;	result = strnumber( st, [val] )
; INPUTS:
;	st - string
; OUTPUTS:
;	1 is returned as the function value if the string st has a
;	valid numeric value, otherwise, 0 is returned.
; OPTIONAL OUTPUT:
;	val - (optional) value of the string.  real*8
; WARNING:
;    (1)   Note that (as of Version 2.2.2) a blank string (e.g. " ") is not
;          a valid numeric value, although an empty string ("") is.
;    (2)   As of V2.2.2 there is a bug in the IDL ON_IOERROR procedure that
;          will cause the following statement to hang up IDL
;
;            IDL> print,'' + string( strnumber('xxx') )
; HISTORY:
;	version 1  By D. Lindler Aug. 1987
;-
if N_params() EQ 0 then begin
     print,'Syntax - result = strnumber( st, [val] )
     return, 0
endif

On_IOerror,L1			;Go to L1 if conversion error occurs
val = double(st)
return, 1			;No conversion error
L1: return, 0			;Conversion error occured
end
;
FUNCTION REPCHR, OLD, C1, C2
;+NAME/ONE LINE DESCRIPTION: 
;    REPCHR replaces one character with another in a text string.
;
; PURPOSE: 
;   Replace all occurences of one character with another
;   in a text string.   (Use the procedure REPSTR to replace
;   more than one character.)
;
; CALLING SEQUENCE: 
;   NEW = REPCHR(OLD, C1, [C2])
;
; INPUTS:
;   OLD = text string to edit, scalar or vector
;   C1 = character to replace.
; OPTIONAL INPUTS:
;   C2 = character to insert (def = ' ' = space).
;
; OUTPUTS:
;   NEW = edited string.
;
; EXAMPLE:
;   If OLD = 'THIS_IS_THE_TEXT', C1 = '_'
;   then NEW = REPCHR(OLD,C1) ==> NEW = 'THIS IS THE TEXT'
; MODIFICATION HISTORY: 
;   R. Sterner.  28 Oct, 1986.
; Vidya Sagar (ARC) Aug 1991
; Improved Help documentation 
;
;.TITLE
;Routine REPCHR
;-
B = BYTE(OLD)			   ; convert string to a byte array.
CB1 = BYTE(C1)			   ; convert char 1 to byte.
W = WHERE(B EQ CB1(0),NFOUND)	   ; find occurrences of char 1.
IF NFOUND EQ 0 THEN RETURN, OLD	   ; if none, return old string.
IF N_PARAMS(0) LT 3 THEN C2 = ' '  ; default char 2 is space.
CB2 = BYTE(C2)			   ; convert char 2 to byte.
B(W) = CB2(0)			   ; replace char 1 by char 2.
RETURN, STRING(B)		   ; return new string.
END
;
pro remchar,st,char	;Remove character
;+
; NAME:
;    REMCHAR removes all appearances of a character from a string
;
; PURPOSE:
;    Remove all appearances of character (char) from string (st)
;
; CALLING SEQUENCE:
;    REMCHAR, ST, CHAR
;
; INPUTS:
;    ST  - String from which character will be removed.  
;    CHAR- Character to be removed from string. 
;
; EXAMPLE:
;    If a = 'a,b,c,d,e,f,g' then IDL> remchar,a, ','
;      will give a = 'abcdefg'
;
; REVISIONS HISTORY
;    Written D. Lindler October 1986
;    Test if empty string needs to be returned   W. Landsman  Feb 1991
;
;spr 10448
;.title
;Routine REMCHAR
;-
 bst = byte(st)                                 ;Convert string to byte

 bchar = byte(char) & bchar = bchar(0)          ;Convert character to byte

 good = where( bst NE bchar, Ngood)
 if Ngood GT 0 then st = string(bst(good)) else st = ''

 return
 end
;
FUNCTION GETTOK,ST,CHAR
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     GETTOK retrieves a segment of a character string.
;
;DESCRIPTION:
;     IDL function to retrieve the first part of a character string
;     until the character CHAR is encountered.
;
;CALLING SEQUENCE:
;     Y = GETTOK(ST,CHAR)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     Y             O   str        output character string which
;                                  terminates at CHAR.
;
;     ST            I   str        input character string
;
;     CHAR          I   str        character at which output character
;                                  string is to be terminated.
;
;WARNINGS:
;     NONE
;
;EXAMPLE:
;     ST='ABC=999'
;     Y=GETTOK(ST,'=')
;     PRINT,Y
;       ABC
;
;#
;COMMON BLOCKS:
;     NONE
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     GETTOK uses the IDL function STRPOS to determine the position of
;     CHAR in the character string ST, and the IDL funtion STRMID to
;     extract the character string segment which occurs before CHAR.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     NONE
;
;MODIFICATION HISTORY
;     Written by D. Lindler, April 1986
;     Revised by Vidya Sagar, Applied Research Corp.   Aug 1991
;     Improved documentation
;
; SPR 9616
;.TITLE
;ROUTINE GETTOK
;
;-
; if char is a blank treat tabs as blanks
;
	tab='	'
	while strpos(st,tab) ge 0 do begin
		pos=strpos(st,tab)
		strput,st,' ',pos
	end

	;
	; find character in string
	;
	pos=strpos(st,char)
	if pos eq -1 then begin	;char not found?
		token=st
		st=''
		return,token
	endif

	;
	; extract token
	;
	token=strmid(st,0,pos)
	len=strlen(st)
	if pos eq (len-1) then st='' else st=strmid(st,pos+1,len-pos-1)

	;
	;  Return the result.
	;
	return,token
	end
;
PRO ZPARCHECK,PROGNAME,PARAMETER,PARNUM,TYPES,DIMENS,MESSAGE
;+NAME:
;     ZPARCHECK verifies data type of parameters in a procedure call.
;
; PURPOSE:
;	Routine to check user parameters to a procedure and prints
;       a message/returns if the parameters don't match
;
; CALLING SEQUENCE:
;	zparcheck,progname,parameter,parnum,types,dimens,message
;
; ARGUMENTS:
;	progname  I str name of     -  calling procedure
;	parameter I any-variable    -  parameter passed to the routine
;	parnum    I int             -  parameter number for error
;                                       message output
;	types     I int vec or sclr -  valid data types for parameter
;		 	1 - byte        2 - integer  3 - int*4
;			4 - real*4      5 - real*8   6 - complex
;			7 - string      8 - structure
;	ndimens   I int sclr or vec -  number of dimensions
;					of allowed dimensions.
;	message -[I] str            -  string message to be printed
;				      	if an error is found).
;
; WARNINGS:
;	Do not misuse government software.
; EXAMPLES:
;  ;;;In this example, the routine simply returns since
;  ;;;the number '10' is an integer scaler
;  UIDL> zparcheck,'DOGANDPONY',10,7,2,0,'Error you dummy'
;  UIDL>
;  ;;;Here we get an error since 10 is not a long scaler
;  UIDL> zparcheck,'DOGANDPONY',10,7,3,0,'Error you dummy'
;  Parameter 7 (Error you dummy)  of routine DOGANDPONY is an invalid
;                                                             data type
;  Valid dimensions are: scaler
;  Valid types are:  longword
;  UIDL>
;#
;
;       
; COMMON BLOCKS:
;        None
;
; PROCEDURES (AND OTHER PROGRAMMING NOTES):
;	If there is an error in the parameters then a
;	a RETALL issued
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;	None
;
; MODIFICATION HISTORY:
;
;	version 1  D. Lindler  Dec. 86
;	documentation updated.  M. Greason, May 1990.
;.TITLE
;Routine ZPARCHECK
;-
;
; convert types and ndimens to vectors if scalars supplied
;
vtypes=intarr(n_elements(types)) + types
vndims=intarr(n_elements(dimens)) + dimens
;
; get type and size of parameter
;
s=size(parameter)
ndim=s(0)
type=s(ndim+1)
;
; check if parameter defined.
;
if type eq 0 then begin
	err=' is undefined.'
	goto,abort
endif
;
; check for valid dimensions
;
valid=where(ndim eq vndims,nvalid)
if nvalid lt 1 then begin
	err='has wrong number of dimensions'
	goto,abort
endif
;
; check for valid type
;
valid=where(type eq vtypes,ngood)
if ngood lt 1 then begin
	err='is an invalid data type'
	goto,abort
endif
;
return
;
; bad parameter
;
abort:
mess=' '
if n_params(0) lt 6 then message=''
if message ne '' then mess=' ('+message+') '
print,string(7b) + 'Parameter '+strtrim(parnum,2)+mess,$
	' of routine ',strupcase(progname)+' ',err
sdim=' '
for i=0,n_elements(vndims)-1 do begin
	if vndims(i) eq 0 then sdim=sdim+'scalar' $
			  else sdim=sdim+string(vndims(i),'(i3)')
end
print,'Valid dimensions are:'+sdim
;
stype=' '
for i=0,n_elements(vtypes)-1 do begin
	case vtypes(i) of
		1: stype=stype+' byte'
		2: stype=stype+' integer'
		3: stype=stype+' longword'
		4: stype=stype+' real*4'
		5: stype=stype+' real*8'
		6: stype=stype+' complex'
		7: stype=stype+' string'
                8: stype=stype+' structure'
	endcase
endfor
print,'Valid types are:'+stype
;if !debug then stop
retall  ; zparcheck
end
;
pro fdecomp,filename,disk,dir,name,qual,version
;+ NAME/ONE LINE DESCRIPTION OF ROUTINE:
;      FDECOMP decomposes a file name into its components.
;
;  DESCRIPTION:
;      IDL procedure to decompose a file name into its
;      constituent components (disk, directory, file name,
;      file name extension, and version number).
;
;  CALLING SEQENCE:
;      FDECOMP,FILENAME,DISK,DIR,NAME,QUAL,VERSION
;
;  ARGUMENTS (I = input, O = output, [] = optional):
;      FILENAME       I      str      File name to be decomposed
;      DISK           O      str      Disk name
;      DIR            O      str      Directory name
;      NAME           O      str      File name
;      QUAL           O      str      File name extension
;      VERSION        O      str      Version number
;
;  WARNINGS:
;      None
;
;  EXAMPLE:
;      To decompose the file name
;      FIRCOADD:[FIRAS.SKYMAP]FCF_SKY_LLSS.ED_8934301_9026410;2:
;
;      filename =
;            'FIRCOADD:[FIRAS.SKYMAP]FCF_SKY_LLSS.ED_8934301_9026410;2'
;      fdecomp,filename,disk,dir,name,qual,ver
;
;      FDECOMP returns these strings for the file name components:
;            disk = 'fircoadd'
;            dir  = '[firas.skymap]'
;            name = 'fcf_sky_llss'
;            qual = 'ed_8934301_9026410'
;            ver  = '2'
;#
;  COMMON BLOCKS:
;      None
;
;  PROCEDURE (AND OTHER PROGRAMMING NOTES):
;      The input and output variables to this procedure are all string
;      arrays.
;
;  PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;      GETTOK
;
;  MODIFICATION HISTORY:
;      version 1  D. Lindler  Oct 1986
;
; SPR 9616
;.TITLE
;Routine FDECOMP
;-
;
st=filename
;
; get disk
;
if strpos(st,':') ge 0 then disk=gettok(st,':')+':' else disk=''
;
; get dir
;
if strpos(st,']') ge 0 then dir=gettok(st,']')+']' else dir=''
;
; get name
;
name=gettok(st,'.')
;
; get qualifier
;
qual=gettok(st,';')
;
; get version
;
version=st
return
end
;
;
;---------------------------------------------------------------------------
;   I deleted procedure findpro from here.  It was nothing but trouble.
;---------------------------------------------------------------------------
;
;----------------------------------------------------------------------------
; IDL Binary Search Routine
;
; Written By: BA Franz, Applied Research Corp., 12/92
;
;----------------------------------------------------------------------------
function binsrch,array,x

n = n_elements(array)
if (n eq 1) then return,0

found = 0
lo = long(0)
hi = long(n-1)

while (lo le hi) and (not found) do begin
    mid = long((lo+hi)/2)
    if (array(mid) eq x) then $
        found = 1             $
    else begin
        if (x lt array(mid)) then  $
             hi = mid-1            $
        else                       $
             lo = mid+1
    endelse 
endwhile

if (found) then   $
    return,mid    $
else              $
    return,-1

end



;------------------------------------------------------------------
; Function ANYWHERE
;
; IDL function which returns the indices of the elements of vector a
; which are also elements of vector b.
;
; i.e., if     indx = anywhere(a,b)
;              c = a(indx)  
;       then   c is the intersection of sets a and b
; 
; Written By: BA Franz, ARC, 6/93
;
;------------------------------------------------------------------
function anywhere,a,b

aa = a(*)
n = n_elements(aa)
index = lonarr(n)*0-1

bb = b(sort(b))
unique_list,bb,bb

j = 0L
for i=0L,n-1 do begin
    ipos = binsrch(bb,aa(i))
    if (ipos ne -1) then begin
        index(j) = i
        j = j+1
    endif
endfor

if (j gt 0) then return,index(0:j-1) else return,-1
end



;------------------------------------------------------------------
; Function NOWHERE
;
; IDL function which returns the indices of the elements of vector a
; which are not elements of vector b.
;
; i.e., if     indx = nowhere(a,b)
;              c = a(indx)  
;       then   c is the difference of sets a and b
; 
; Written By: BA Franz, ARC, 6/93
;
;------------------------------------------------------------------
function nowhere,a,b

n = n_elements(a)
index = lonarr(n)*0-1

bb = b(sort(b))

j = 0L
for i=0L,n-1 do begin
    ipos = binsrch(bb,a(i))
    if (ipos eq -1) then begin
        index(j) = i
        j = j+1
    endif
endfor

if (j gt 0) then return,index(0:j-1) else return,-1
end
;
;
;
pro readcol,name,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            FORMAT = fmt, DEBUG=debug, SILENT=silent, SKIPLINE = skipline, $
            NUMLINE = numline
;+
; NAME:
;    READCOL reads a free-format ASCII data file with columns of data.
;
; PURPOSE:
;    Read a free-format ASCII data file with columns of data into IDL variables.
;    Lines of data not meeting the specified format (e.g. comments) are
;    ignored.  Columns may be separated by commas or spaces.   Use READFMT
;    to read a fixed-format ASCII file.
;
; CALLING SEQUENCE:
;    readcol, name, v1, [ v2, v3, v4, v5,v6, v7, v8,v9,v10,v11,v12,v13,v14,v15
;             FORMAT = , DEBUG = , SILENT = , SKIPLINE = , NUMLINE = ]
; INPUTS:
;    NAME - Name of ASCII data file, scalar string.  In VMS, an extension of 
;          .DAT is assumed, if not supplied.
;
; OPTIONAL INPUT KEYWORDS:
;    FORMAT - scalar string containing a letter specifying an IDL type
;        for each column of data to be read.  Allowed letters are 
;        A - string data, B - byte, D - double precision, F- floating 
;        point, I - integer, L - longword, and X - skip a column.
;
;        Columns without a specified format are assumed to be floating 
;        point.  Examples of valid values of FMT are
;
;        'A,B,I'        ;First column to read as 6 character string, then 
;                        1 column of byte data, 1 column integer data
;        'L,L,L,L'       ;Four columns will be read as longword arrays.
;        ' '             ;All columns are floating point
;
;        If a FORMAT keyword string is not supplied, then all columns are 
;        assumed to be floating point.
;
;   SILENT - Normally, READCOL will display each line that it skips over.
;          If SILENT is set and non-zero then these messages will be suppressed.
;   DEBUG - If this keyword is non-zero, then additional information is printed
;          as READCOL attempts to read and interpret the file.
;   SKIPLINE - Scalar specifying number of lines to skip at the top of file
;          before reading.   Default is to start at the first line.
;   NUMLINE - Scalar specifying number of lines in the file to read.   Default
;             is to read the entire file
;
; OUTPUTS:
;    V1,V2,V3,...V15 - IDL vectors to contain columns of data.
;         Up to 15 columns may be read.  The type of the output vectors
;         are as specified by FORMAT.
;
; EXAMPLES:
;    Each row in a file POSITION.DAT contains a star name and 6 columns
;    of data giving an RA and Dec in sexigesimal format.   Read into IDL 
;    variables.     (NOTE: The star names must not contain internal spaces.)
;
;        IDL> FMT = 'A,I,I,F,I,I,F'
;        IDL> READCOL,'POSITION',F=FMT,name,hr,min,sec,deg,dmin,dsec  
;
;    The HR,MIN,DEG, and DMIN variables will be integer vectors.
;
;    Alternatively, all except the first column could be specified as
;    floating point.
;        IDL> READCOL,'POSITION',F='A',name,hr,min,sec,deg,dmin,dsec 
;
;    To read just the variables HR,MIN,SEC
;        IDL> READCOL,'POSITION',F='X,I,I,F',HR,MIN,SEC
;
; RESTRICTIONS:
;    This procedure is designed for generality and not for speed.
;    If a large ASCII file is to be read repeatedly, it may be worth
;    writing a specialized reader.
;
;    Columns to be read as strings must not contain spaces or commas, 
;    since these are interpreted as column delimiters.    Use READFMT
;    to read such files.
;
;    Numeric values are converted to specified format.  For example,
;    the value 0.13 read with an 'I' format will be converted to 0.
;
; PROCEDURES CALLED
;     GETTOK, SPEC_DIR, REPCHR, STRNUMBER
; REVISION HISTORY:
;    Written         W. Landsman                 November, 1988
;    Modified	     J. Bloch 			June, 1991
;	(Fixed problem with over allocation of logical units.)
;    Added SKIPLINE and NUMLINE keywords  W. Landsman    March 92
;
;spr 10448
;.title
;Routine READCOL
;-
  On_error,2                           ;Return to caller

  if N_params() lt 2 then begin
     print,'Syntax - readcol, name, v1, [ v2, v3,...v15, '
     print,'        FORMAT= ,SILENT = ,SKIPLINE =, NUMLINE = , /DEBUG]'
     return
  endif

  ncol = N_params() - 1           ;Number of columns of data expected
  vv = 'v' + strtrim( indgen(ncol)+1, 2)
  nskip = 0

  if N_elements(fmt) GT 0 then begin    ;FORMAT string supplied?

    zparcheck, 'READCOL', fmt, 2, 7, 0, 'FORMAT string'
;   Remove blanks from format string
    frmt = strupcase(strcompress(fmt,/REMOVE))   
    remchar, frmt, '('                  ;Remove parenthesis from format
    remchar, frmt, ')'           

;   Determine number of columns to skip ('X' format)
    pos = strpos(frmt, 'X', 0)

    while pos NE -1 do begin
        pos = strpos( frmt, 'X', pos+1)
        nskip = nskip + 1
    endwhile

  endif else begin                     ;Read everything as floating point

    frmt = 'F'
    if ncol GT 1 then for i = 1,ncol-1 do frmt = frmt + ',F'
    if not keyword_set( SILENT ) then message, $
      'Format keyword not supplied - All columns assumed floating point',/INF

  endelse

  nfmt = ncol + nskip
  idltype = intarr(nfmt)
  openr, lun, name, ERROR=err, /GET_LUN
  if err LT 0 then $ 
     message,'Unable to open file ' + spec_dir( name, 'DAT')

; Get number of lines in file

   nlines = 0L
   temp = ' '
   while not EOF(lun) do begin
      readf, lun, temp
      nlines = nlines+1
   endwhile

   if keyword_set(DEBUG) then $
      message,strupcase(name)+' contains ' + strtrim(nlines,2) + ' lines',/INF

   if not keyword_set( SKIPLINE ) then skipline = 0
   nlines = nlines - skipline
   if keyword_set( NUMLINE) then nlines = numline < nlines

; Create output arrays according to specified formats

   k = 0L                                     ;Loop over output columns
   for i = 0L, nfmt-1 do begin

       fmt1 = gettok( frmt, ',' )
       if fmt1 EQ '' then fmt1 = 'F'         ;Default is F format
       case strmid(fmt1,0,1) of 
          'A':  idltype(i) = 7          
          'D':  idltype(i) = 5
          'F':  idltype(i) = 4
          'I':  idltype(i) = 2
          'B':  idltype(i) = 1
          'L':  idltype(i) = 3
          'X':  idltype(i) = 0               ;IDL type of 0 ==> to skip column
         ELSE:  message,'Illegal format ' + fmt1 + ' in field ' + strtrim(i,2)
      endcase

; Define output arrays

      if idltype(i) NE 0 then begin
          st = vv(k) + '= make_array(nlines,TYPE = idltype(i) )'  
          tst = execute(st)
          k = k+1
      endif

   endfor

   free_lun,lun
   openr, lun, name, /get_lun
   ngood = 0L

   if skipline GT 0 then $
       for i = 0, skipline-1 do readf, lun, temp        ;Skip any lines

   for j = 0L, nlines-1 do begin

      readf, lun, temp
      if strlen(temp) LT ncol then begin    ;Need at least 1 chr per output line
          ngood = ngood-1
          if not keyword_set(SILENT) then $
                       message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF
          goto, BADLINE 
       endif
    temp = repchr(temp,',','  ')        ;Replace comma delimiters by spaces
    k = 0

    for i = 0L,nfmt-1 do begin

       temp = strtrim(temp,1)                  ;Remove leading spaces
       var = gettok(temp,' ')                  ;Read next field
       if ( idltype(i) NE 0 ) then begin       ;Expecting data?

          if ( idltype(i) NE 7 ) then begin    ;Check for valid numeric data
             tst = strnumber(var,val)          ;Valid number?
             if tst EQ 0 then begin            ;If not, skip this line
                 if not keyword_set(SILENT) then $ 
                      message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF 
                 ngood = ngood-1
                 goto, BADLINE 
             endif
          st = vv(k) + '(ngood) = val'     

         endif else $
           st = vv(k) + '(ngood) = strtrim(var,2)'

      tst = execute(st)
      k = k+1

    endif  

  endfor

BADLINE:  ngood = ngood+1

   endfor

  free_lun,lun
  if ngood EQ 0 then message,'ERROR - No valid lines found for specified format'

  if not keyword_set(SILENT) then $
        message,strtrim(ngood,2) + ' valid lines read', /INFORM  

; Compress arrays to match actual number of valid lines

  if ngood EQ 0 then return

  for i = 0,ncol-1 do begin 
      tst = execute(vv(i) + '='+ vv(i)+ '(0:ngood-1)')
  endfor

  return
  end


