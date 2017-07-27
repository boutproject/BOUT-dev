pro show_one, v,r,z, alog=alog, fill=fill, title=title
;
;

nlev=30
contour, v, r, z, /iso, nlev=nlev,fill=fill,title=title, c_col=findgen(nlev)*255/nlev

oplot, r, z, psym=3
xyouts, chars=1.5, 0.7*(min(r)+max(r))/2, (min(z)+max(z))/2, min(v)
xyouts, chars=1.5, 0.7*(min(r)+max(r))/2, 0.9*(min(z)+max(z))/2, max(v)

;;stop
;
;
end


pro show_two, v1,r1,z1, v2,r2,z2, alog=alog, fill=fill, title1=title1, title2=title2
;
;
;

!p.multi=[0,2,1,0,0]
Show_One, v1,r1,z1, alog=alog, fill=fill, title=title1
Show_One, v2,r2,z2, alog=alog, fill=fill, title=title2
!p.multi=0

;
;
;
end




pro show_nphi, arms=arms, du=du, alog=alog, fill=fill, it=it
;
;
;

if not keyword_set(IT) then begin
    it=0
    read, it, prompt='it='
endif
!p.subtitle="it="+STRTRIM(STRING(it),2)


rb06=du.rxy
zb06=du.zxy
rbpp=du.rxy
zbpp=du.zxy


!p.multi=[0,2,1,0,0]

        vb06=REFORM(arms.ni_xyzt[*,*,it])
        title="RMS<Ni>"
        if keyword_set(ALOG) then begin
            vb06=alog(vb06+1e-20)
            title="ALOG("+title+")"
        endif
        SHOW_ONE, vb06,rb06,zb06, alog=alog, fill=fill, title=title

        
        vb06=REFORM(arms.phi_xyzt[*,*,it])
        title="RMS<Phi>"
        if keyword_set(ALOG) then begin
            vb06=alog(vb06+1e-20)
            title="ALOG("+title+")"
        endif
        SHOW_ONE, vb06,rb06,zb06, alog=alog, fill=fill, title=title

!p.multi=0
!p.subtitle=''
;
;
;
end
