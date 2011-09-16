pro hlmk_grids
  restore,'helimak.sav'
  plot,(shot_data.set3.vfloat)[*,0]
  oplot,(shot_data.set3.vfloat)[*,1]
  
  oplot,(shot_data.set3.vfloat)[*,2]
  oplot,(shot_data.set3.vfloat)[*,3]
  oplot,(shot_data.set3.vfloat)[*,4]
  oplot,(shot_data.set3.vfloat)[*,5]
  
  oplot,(shot_data.set3.vfloat)[*,6]
end
