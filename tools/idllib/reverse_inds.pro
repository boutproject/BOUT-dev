; Reverse array indices
; arr[t,z,y,x] -> arr[x,y,z,t]
; Preserves the type of the input variable
; (Float -> Float, Double -> Double etc.)

FUNCTION reverse_inds, var
  nd = SIZE(var, /N_DIM)
  dims = SIZE(var, /DIMENSIONS)
  
  CASE nd OF
    0: data = var
    1: data = var
    2: BEGIN
      data = REFORM(var, dims[1], dims[0])
      FOR i=0, dims[0]-1 DO BEGIN
        FOR j=0, dims[1]-1 DO BEGIN
          data[j,i] = var[i,j]
        ENDFOR
      ENDFOR
    END
    3: BEGIN
      data = REFORM(var, dims[2], dims[1], dims[0])
      FOR i=0, dims[0]-1 DO BEGIN
        FOR j=0, dims[1]-1 DO BEGIN
          FOR k=0, dims[2]-1 DO BEGIN
            data[k,j,i] = var[i,j,k]
          ENDFOR
        ENDFOR
      ENDFOR
    END
    4: BEGIN
      data = REFORM(var, dims[3], dims[2], dims[1], dims[0])
      FOR i=0, dims[0]-1 DO BEGIN
        FOR j=0, dims[1]-1 DO BEGIN
          FOR k=0, dims[2]-1 DO BEGIN
            FOR l=0, dims[3]-1 DO BEGIN
              data[l,k,j,i] = var[i,j,k,l]
            ENDFOR
          ENDFOR
        ENDFOR
      ENDFOR
    END
    5: BEGIN
      data = REFORM(var, dims[4], dims[3], dims[2], dims[1], dims[0])
      FOR i=0, dims[0]-1 DO BEGIN
        FOR j=0, dims[1]-1 DO BEGIN
          FOR k=0, dims[2]-1 DO BEGIN
            FOR l=0, dims[3]-1 DO BEGIN
              FOR m=0, dims[4]-1 DO BEGIN
                data[m,l,k,j,i] = var[i,j,k,l,m]
              ENDFOR
            ENDFOR
          ENDFOR
        ENDFOR
      ENDFOR
    END
    ELSE: BEGIN
      PRINT, "Sorry: Reverse_inds can't handle this many dimensions"
      RETURN, 0
    END
  ENDCASE

  RETURN, data
END
