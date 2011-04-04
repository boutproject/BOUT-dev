
FUNCTION diff, arr
  r = arr[1:*] - arr[0:(N_ELEMENTS(arr)-2)]
  RETURN, [r[0], r]
END
