FUNCTION idl_legal, str
; Convert the variable name into a legal IDL variable name
;--------------------------------------------------

name = str                      ;
status = CALL_EXTERNAL('pdb2idl.so', 'idl_legalvar', name, /i_value);
RETURN, name
END

; primitive PDB function wrappers

FUNCTION pdb_open, filename, write=write
   filename = STRING(filename) ; just make sure it's a string
   IF KEYWORD_SET(write) THEN BEGIN
       status = CALL_EXTERNAL('pdb2idl.so', 'pdb_open', filename, 1, /i_value); /* If has 2 or more arguments-> write */
   ENDIF ELSE BEGIN
       status = CALL_EXTERNAL('pdb2idl.so', 'pdb_open', filename, /i_value)
   ENDELSE
   RETURN, status
END

PRO pdb_close
   status = CALL_EXTERNAL('pdb2idl.so', 'pdb_close', /i_value)
END


; returns an array of the variable names in the already open file
FUNCTION pdb_get_list 
   nvars = CALL_EXTERNAL('pdb2idl.so', 'pdb_list_vars', /i_value)
   IF nvars LE 0 THEN RETURN, 0
   vars = STRARR(nvars) ; the array of variable names
   FOR i=0, nvars-1 DO BEGIN
       str = CALL_EXTERNAL('pdb2idl.so', 'pdb_list_varname', FIX(i), /s_value)
       vars[i] = str
   ENDFOR
   ; status = CALL_EXTERNAL('pdb2idl.so', 'pdb_list_free', /i_value)
   RETURN, vars
END

; returns a long array similar to size(var) containing
; type as number: 0 - unknown 1 - integer 2 - float 3 -  double
; number of dimensions
; size of dimensions (one entry per dimension)
; total amount of data
FUNCTION pdb_query_var, varname
   varname = STRING(varname)
   ndims = CALL_EXTERNAL('pdb2idl.so', 'pdb_query_ndim', varname, /i_value);
   query = LONARR(2*ndims + 3)
   
   typestr = CALL_EXTERNAL('pdb2idl.so', 'pdb_query_type', varname, /s_value);
   type = 0
   IF typestr EQ "integer" THEN type = 1
   IF typestr EQ "float" THEN type = 2
   IF typestr EQ "double" THEN type = 3
   IF typestr EQ "char" THEN type=4

   query[0] = LONG(type)
   query[1] = LONG(ndims)
   
   dims = LONARR(2*ndims + 1)
   
   status = CALL_EXTERNAL('pdb2idl.so', 'pdb_query_dims', varname, dims, /i_value)
   FOR i=0, 2*ndims DO BEGIN
       query[i+2] = dims[i]
   ENDFOR
   
   RETURN, query
END

FUNCTION pdb_read_var, name, query
   name = STRING(name)

   IF((query[0] EQ 0) OR (query[N_ELEMENTS(query)-1] EQ 0)) THEN BEGIN
       PRINT, "Could not read variable "+name
       RETURN, 0
   ENDIF  

   type = query[0]
   ndims = query[1]

   IF type EQ 1 THEN BEGIN
       ; read integers
       var = 0
       CASE ndims OF
           0: var = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_0d_int', name)
           1: BEGIN
               xsize = query[3] - query[2] + 1
               var = LONARR(xsize)
               inds = LONARR(3)
               inds[0] = LONG(query[2])
               inds[1] = LONG(query[3])
               inds[2] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_1d_int', name, inds, var, /i_value)
              END
           ELSE: PRINT, "Cannot read in "+STRTRIM(STRING(ndims),2)+" dimensional integer"
       ENDCASE
   ENDIF ELSE IF (type EQ 2) OR (type EQ 3) THEN BEGIN
       ; read floats or doubles (returned as floats)
       var = 0.0
       CASE ndims OF
           0: var = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_0d_float', name, /f_value)
           1: BEGIN
               xsize = query[3] - query[2] + 1
               var = FLTARR(xsize)
               inds = LONARR(3)
               inds[0] = LONG(query[2])
               inds[1] = LONG(query[3])
               inds[2] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_1d_float', name, inds, var, /i_value)
              END
           2: BEGIN
               xsize = query[3] - query[2] + 1
               ysize = query[5] - query[4] + 1
               var = FLTARR(xsize, ysize)
               inds = LONARR(6)
               inds[0] = LONG(query[2])
               inds[1] = LONG(query[3])
               inds[2] = 1L
               inds[3] = LONG(query[4])
               inds[4] = LONG(query[5])
               inds[5] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_2d_float', name, inds, var, /i_value)
              END
           3: BEGIN
               xsize = query[3] - query[2] + 1
               ysize = query[5] - query[4] + 1
               zsize = query[7] - query[6] + 1
               var = FLTARR(xsize, ysize, zsize)
               inds = LONARR(9)
               inds[0] = LONG(query[2])
               inds[1] = LONG(query[3])
               inds[2] = 1L
               inds[3] = LONG(query[4])
               inds[4] = LONG(query[5])
               inds[5] = 1L
               inds[6] = LONG(query[6])
               inds[7] = LONG(query[7])
               inds[8] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_3d_float', name, inds, var, /i_value)
              END
           4: BEGIN
               xsize = query[3] - query[2] + 1
               ysize = query[5] - query[4] + 1
               zsize = query[7] - query[6] + 1
               tsize = query[9] - query[8] + 1
               var = FLTARR(xsize, ysize, zsize, tsize)
               inds = LONARR(12)
               inds[0] = LONG(query[2])
               inds[1] = LONG(query[3])
               inds[2] = 1L
               inds[3] = LONG(query[4])
               inds[4] = LONG(query[5])
               inds[5] = 1L
               inds[6] = LONG(query[6])
               inds[7] = LONG(query[7])
               inds[8] = 1L
               inds[9] =  LONG(query[8])
               inds[10] = LONG(query[9])
               inds[11] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_4d_float', name, inds, var, /i_value)
              END
           ELSE: PRINT, "Cannot read in "+STRTRIM(STRING(ndims),2)+" dimensional float/double"
       ENDCASE
   ENDIF ELSE BEGIN
       ; read in a string
       CASE ndims OF
           1: BEGIN
               xsize=query[3] - query[2] + 1
               var = CALL_EXTERNAL('pdb2idl.so', 'pdb_read_string', name, LONG(xsize), /s_value)
           END
           ELSE: PRINT, "Cannot read in "+STRTRIM(STRING(ndims),2)+" dimensional string"
       ENDCASE
   ENDELSE
   RETURN, var
END

FUNCTION pdb_write_var, name, var, lowercase=lowercase
   name = STRING(name)
   if keyword_set(LOWERCASE) then name=strlowcase(name)

   ; get info about variable
   query = SIZE(var, /structure)


   type = query.type
   ndims = query.n_dimensions

   IF (type EQ 3) OR (type EQ 2) THEN BEGIN
       ; write integers or longs (written as longs)
       var = LONG(var)
       CASE ndims OF
           0: status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_0d_int', name, var, /i_value)
           1: BEGIN
               xsize = query.n_elements
               inds = LONARR(3)
               inds[0] = LONG(0)
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_1d_int', name, inds, var, /i_value)
              END
           ELSE: PRINT, "Cannot write "+STRTRIM(STRING(ndims),2)+" dimensional ", query.type_name
       ENDCASE
   ENDIF ELSE IF type EQ 4 THEN BEGIN
       ; write floats (written as floats)
       CASE ndims OF 
           0:status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_0d_float', name, var, /i_value)
           1: BEGIN
               xsize = query.n_elements
               inds = LONARR(3)
               inds[0] = LONG(0)
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_1d_float', name, inds, var, /i_value)
           END
           2: BEGIN
               xsize = query.dimensions[0]
               ysize = query.dimensions[1]
               inds = LONARR(6)
               inds[0] = 0L
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               inds[3] = 0L
               inds[4] = LONG(ysize-1)
               inds[5] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_2d_float', name, inds, var, /i_value)
           END
           3: BEGIN
               xsize=query.dimensions[0]
               ysize=query.dimensions[1]
               zsize=query.dimensions[2]

               inds = LONARR(9)
               inds[0] = 0L
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               inds[3] = 0L
               inds[4] = LONG(ysize-1)
               inds[5] = 1L
               inds[6] = 0L
               inds[7] = LONG(zsize-1)
               inds[8] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_3d_float', name, inds, var, /i_value)
           END
           4: BEGIN
               xsize=query.dimensions[0]
               ysize=query.dimensions[1]
               zsize=query.dimensions[2]
               tsize=query.dimensions[3]

               inds = LONARR(12)
               inds[0] = 0L
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               inds[3] = 0L
               inds[4] = LONG(ysize-1)
               inds[5] = 1L
               inds[6] = 0L
               inds[7] = LONG(zsize-1)
               inds[8] = 1L
               inds[9] = 0L
               inds[10] = LONG(tsize-1)
               inds[11] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_4d_float', name, inds, var, /i_value)
           END
           ELSE: PRINT, "Cannot write "+STRTRIM(STRING(ndims),2)+" dimensional ", query.type_name
       ENDCASE
   ENDIF ELSE IF type EQ 5 THEN BEGIN
       ; write doubles (written as doubles)
       CASE ndims OF 
           0:status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_0d_double', name, var, /i_value)
           1: BEGIN
               xsize = query.n_elements
               inds = LONARR(3)
               inds[0] = LONG(0)
               inds[1] = LONG(xsize-1)
               inds[2] = 1L
               status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_1d_double', name, inds, var, /i_value)
              END
           ELSE: BEGIN
               PRINT, "Cannot write "+STRTRIM(STRING(ndims),2)+" dimensional double "+name+". Writing as float"
               status = pdb_write_var(name, FLOAT(var))
           END
       ENDCASE
   ENDIF ELSE IF type EQ 7 THEN BEGIN
       ; write a string
       status = CALL_EXTERNAL('pdb2idl.so', 'pdb_write_string', name, var, /i_value)
   ENDIF ELSE BEGIN
       PRINT, "Sorry, don't know how to write variables of type ", query.type_name
       status = 0
   ENDELSE
   RETURN, status
END

; FUNCTIONS TO BE CALLED BY USER (CRASH-PROOFED I HOPE)

PRO PD_ls, filename
   s = size(filename, /structure)
   IF s.type NE 7 THEN BEGIN
       PRINT, "Useage: PD_ls, 'filename'"
       RETURN
   ENDIF
   status = pdb_open(filename)
   IF status THEN BEGIN
       ls = pdb_get_list()
       pdb_close
       PRINT, "Contents of file "+filename+" : "
       PRINT, " "
       PRINT, ls
       PRINT, " "
   ENDIF ELSE PRINT, "could not open file "+filename
END

PRO PD_desc, filename, name
   s = size(filename, /structure)
   IF s.type NE 7 THEN BEGIN
       PRINT, "Useage: PD_desc, 'filename' [,'variable name']"
       RETURN
   ENDIF
   
   status = pdb_open(filename)
   IF status NE 1 THEN BEGIN
       PRINT, "Could not open file "+filename
       RETURN
   ENDIF
   s = size(name)
   IF s[1] NE 7 THEN BEGIN ; not supplied name - describe all
       varlist = pdb_get_list()
       n = N_ELEMENTS(varlist)
   ENDIF ELSE BEGIN
       n = 1
       varlist = STRARR(1)
       varlist[0] = name
   ENDELSE

   FOR j=0, n-1 DO BEGIN
       name = varlist[j]
       query = pdb_query_var(name)
       str = name + ": "
       CASE query[0] OF
           0: str = str + "unknown type"
           1: str = str + "integer"
           2: str = str + "float"
           3: str = str + "double"
           4: str = str + "char"
       ENDCASE

       FOR i=0, query[1] - 1 DO BEGIN
           str = str + "["+STRTRIM(STRING(query[2*i + 2]),2)+":"+$
             STRTRIM(STRING(query[2*i + 3]),2)+"]"
       ENDFOR
       IF query[2*query[1] + 2] GT 1 THEN BEGIN
           str = str + ". Elements: "+ STRTRIM(STRING(query[2*query[1] + 2]),2)
       ENDIF
       PRINT, str
   ENDFOR
   pdb_close
END

; reads in a single variable from the file
FUNCTION PD_read, filename, name
   s = size(filename, /structure)
   t = size(name, /structure)
   IF (s.type NE 7) OR (t.type NE 7) THEN BEGIN
       PRINT, "Useage: var = PD_read('filename', 'variable name')"
       RETURN, 0
   ENDIF

   status = pdb_open(filename)
   var = 0
   IF status THEN BEGIN
       query = pdb_query_var(name) ; get information about this variable
       var = pdb_read_var(name, query) ; read the variable
       pdb_close
   ENDIF ELSE PRINT, "could not open file "+filename
   RETURN, var
END

; reads in the entire file and returns a structure
FUNCTION PD_import, filename
   s = size(filename, /structure)
   IF s.type NE 7 THEN BEGIN
       PRINT, "Useage: data = PD_import('filename')"
       retstr = {pd_error: 1, errstr: 'Incorrect useage'}
       RETURN, retstr
   ENDIF
   status = pdb_open(filename)
   IF NOT status THEN BEGIN
       ; failed to open
       retstr = {pd_error: 2, errstr: 'Could not open file'} ; return an error code and string
       RETURN, retstr
   ENDIF

   ; get a list of variables in the file
   list = pdb_get_list()
   s = size(list)
   IF s[N_ELEMENTS(s)-2] NE 7 THEN BEGIN ; check if returned value is an array of strings
       ; an error getting the list
       pdb_close
       retstr = {pd_error: 3, errstr: 'Could not list file contents'}
       RETURN, retstr
   ENDIF
   
   nvars = s[1] ; number of variables

   retstr = {pd_error: 0}
   FOR i=0, nvars-1 DO BEGIN ; read in the variables
       query = pdb_query_var(list[i]) ; get information about this variable
       var = pdb_read_var(list[i], query) ; read the variable
       retstr = CREATE_STRUCT(retstr, idl_legal(list[i]), var) ; add to the structure
   ENDFOR

   RETURN, retstr
END

; WRITES A SINGLE VARIABLE TO FILE <FILENAME> WITH NAME <NAME>
; status = 1 if successful, 0 otherwise
PRO PD_write, filename, name, var, status=status
   s = size(filename, /structure)
   t = size(name, /structure)
   IF (s.type NE 7) OR (t.type NE 7) THEN BEGIN
       PRINT, "Useage: var = PD_write('filename', 'variable name', variable)"
       status = 0
       RETURN
   ENDIF
   status = pdb_open(filename, /write)
   IF status EQ 1 THEN BEGIN
       status = pdb_write_var(name, var)
       pdb_close
   ENDIF
END

; WRITES A STRUCTURE TO A PDB FILE (INVERSE OF PD_IMPORT FUNCTION)
; status = 1 if successful, 0 otherwise
PRO PD_export, filename, data, lowercase=lowercase, status=status
   s = size(filename, /structure)
   t = size(data, /structure)
   IF (s.type NE 7) OR (t.type NE 8) THEN BEGIN
       PRINT, "Useage: PD_export, 'filename', structure"
       status = 0
       RETURN
   ENDIF
   
   status = pdb_open(filename, /write)
   IF status NE 1 THEN RETURN

   names = TAG_NAMES(data)
   ; loop through each variable and write to file
   status = 1
   FOR i=0, N_TAGS(data)-1 DO BEGIN
       var = data.(i)
       vn = STRUPCASE(names[i])
       IF vn NE "PD_ERROR" THEN BEGIN
           s = pdb_write_var(names[i], var, lowercase=lowercase)
       ENDIF ELSE s = 1
       IF s NE 1 THEN status = 0
   ENDFOR
   pdb_close
END

PRO pdb2idl
  ; dummy procedure - does nothing
  PRINT, "PDB2IDL Library"
END
