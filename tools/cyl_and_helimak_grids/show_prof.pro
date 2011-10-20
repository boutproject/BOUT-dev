pro show_prof

!p.multi = [0,1,1]
!x.thick=3
!y.thick=3
!p.thick=3
!p.charthick=2
!p.charsize =2

save_name = 'profile.ps'
pdf_save_name = 'profile.pdf'

set_plot,'ps'
;device,filename=save_name,/color,YSIZE=25,YOFFSET =1,landscape =0
device,filename=save_name,/color,landscape =1



restore,'helimak.sav'
grid = "./Ln/Helimak_32x16_0.10_small_y.nc"
d = file_import(grid)

d.Ni0 = d.Ni0*1e14
 

plot,d.Rxy[*,d.ny/2]*1e2,d.Ni0[*,d.ny/2],$
        TITLE = "Ni0",xtitle = "radius [cm]",ytitle = textoidl('cm^{-3}'),$
        ystyle = 2

oplot,1e2*(shot_data.set3.R)[*,5],(shot_data.set3.density)[*,5]*1e-6,psym = 2

polyfill,[100,115,115,100],$
         [0,0,5e10,5e10],/data,/line_fill,linestyle =1,color= 75;,spacing = 1

device,/close
spawn,strjoin(["ps2pdf ",save_name, " ",pdf_save_name])
set_plot,'X'


end
