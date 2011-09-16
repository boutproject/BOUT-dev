pro read_csv_helimak
  ;restore,'helimak_template.sav'  ;retreive the base template, by itself its useless
  
  ;data = read_ascii('probecurve.csv',template = helimak_template)
  data = rd_tfile('probecurve.csv',58,delim = ',') ; skips empty lines

  R = double(data[2:57,0])
  help,R
  note = data[*,1]

  probe_name = data[1:57,2]

  temp = fltarr(56)
  stemp = strarr(56)
  
  substrt = {name:stemp,density:temp,Te1:temp,Te2:temp,vfloat:temp}
  substrt2 = {density:temp,Te:temp,vfloat:temp}
  
  shot_data = {shot:0.0,set1:substrt,set2:substrt,set3:substrt2,omega:'0',L_c:0.0}

  
  shot_data = replicate(shot_data,9)


  ;count the number of shots
  N_shots = N_ELEMENTS(wc_where(data[1,*],'Te1'))/2


  shot_i = wc_where(data[1,*],'9*')
  Te1_i = wc_where(data[1,*],'Te1*')
  Te2_i = wc_where(data[1,*],'Te2*')
  Te_i = wc_where(data[1,*],'Te')
  den_i = wc_where(data[1,*],'Den*')
  v_i = wc_where(data[1,*],'Den*')
  
  print,N_shots
  ; loop over the fields
  for i = 0, N_shots-1 do begin 
     print,i
     shot_data[i].shot = data[1,shot_i[2*i]]
   
     shot_data[i].set1.density = data[2:57,den_i[3*i]]
     shot_data[i].set2.density = data[2:57,den_i[3*i+1]]
     shot_data[i].set3.density = data[2:57,den_i[3*i+2]]

     shot_data[i].set1.Te1 = data[2:57,Te1_i[2*i]]
     shot_data[i].set2.Te1 = data[2:57,Te1_i[2*i+1]]
     
     shot_data[i].set1.Te2 = data[2:57,Te2_i[2*i]]
     shot_data[i].set2.Te2 = data[2:57,Te2_i[2*i+1]]

     shot_data[i].set3.Te = data[2:57,Te_i[i]]
     
     shot_data[i].set1.vfloat = data[2:57,v_i[3*i]]
     shot_data[i].set2.vfloat = data[2:57,v_i[3*i+1]]
     shot_data[i].set3.vfloat = data[2:57,v_i[3*i+2]]

     shot_data[i].L_c = data[0,Te_i[i]]
     shot_data[i].omega = data[0,Te_i[i]+1]



  endfor

  plot,R,shot_data[0].set1.density
  print,R
  
  save,shot_data,filename= 'helimak.sav'

;;   help,shot_data.set1,/str
;;   help,shot_data[0].set1,/str
  
end
