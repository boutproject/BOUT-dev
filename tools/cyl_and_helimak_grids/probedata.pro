pro read_ken
  restore,'helimak_template.sav'  ;retreive the base template, by itself its useless
  
  data = read_ascii('probecurve.csv',template = helimak_template)
  help,data,/str
  
end
