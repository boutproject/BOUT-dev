; Collects data together into a single file

zinds = [0,3]

f = file_open("result.nc", /create)

list = file_list("data/BOUT.dmp.0.nc")

; Select only the variables (NYPE, MYSUB etc differ)
w = WHERE(STRCMP(list,'var', 3, /fold_case))
list = list[w]

FOR i=0, N_ELEMENTS(list)-1 DO BEGIN & PRINT, "Reading "+list[i] & var = collect(path="data", var=list[i], zind=zinds) & status = file_write(f, list[i], var) & ENDFOR

file_close, f

exit
