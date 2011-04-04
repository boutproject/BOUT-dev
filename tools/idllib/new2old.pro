; reformats the output from BOUT++ (txyz) to the output from collect_* (xyzt)
FUNCTION new2old, d, ng=ng
   names = TAG_NAMES(d)
   nn = N_ELEMENTS(names)
   
   IF NOT KEYWORD_SET(ng) THEN ng = 2

   result = {t_array:d.t_array, $
             wci:d.wci, rho_s:d.rho_s, $
             zmax:d.zmax, zmin:d.zmin, $
             ngz:d.mz}

   w = WHERE(STRLOWCASE(names) EQ "ni_x", count)
   IF count EQ 1 THEN result = CREATE_STRUCT(result, "ni_x", d.ni_x)
   w = WHERE(STRLOWCASE(names) EQ "te_x", count)
   IF count EQ 1 THEN result = CREATE_STRUCT(result, "te_x", d.te_x)
   w = WHERE(STRLOWCASE(names) EQ "ti_x", count)
   IF count EQ 1 THEN result = CREATE_STRUCT(result, "ti_x", d.ti_x)

   FOR i=0, nn-1 DO BEGIN
       s = SIZE(d.(i))
       IF s[0] EQ 4 THEN BEGIN
           PRINT, "Reformatting "+names[i]
           nt = s[1]
           nx = s[2]
           ny = s[3]
           nz = s[4]
           
           varname = names[i] + "_xyzt"
           var = FLTARR(nx, ny-(2*ng), nz, nt)
           
           orig = d.(i)

           FOR x=0,nx-1 DO BEGIN
               FOR y=0,ny-(1+2*ng) DO BEGIN
                   FOR z=0,nz-1 DO BEGIN
                       FOR t=0,nt-1 DO BEGIN
                           var[x,y,z,t] = orig[t,x,y+ng,z]
                       ENDFOR
                   ENDFOR
               ENDFOR
           ENDFOR
           result = CREATE_STRUCT(varname, var, result)
       ENDIF
   ENDFOR

   result = CREATE_STRUCT("lxmin", 0, "lxmax", nx-1, $
                           "jpmin", 0, "jpmax", ny-(1+2*ng), $
                           "trange", nt, $
                           result)
   
   RETURN, result
END
