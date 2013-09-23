safe_colors

!p.multi=[0,3,2]

plot,gr,xtitle='time (!7s!3!iA!n)',ytitle='!7c!3/!7x!3!iA!n',chars=2,thick=3

plot,psn,rmsni[*,32,-1]/max(rmsni[*,32,-1]),xtitle='!7W!3!iN!n',ytitle='N!ii!n',chars=2,thick=3

plot,psn,rmsti[*,32,-1]/max(rmsti[*,32,-1]),xtitle='!7W!3!iN!n',ytitle='T!ii!n',chars=2,thick=3

plot,psn,rmste[*,32,-1]/max(rmste[*,32,-1]),xtitle='!7W!3!iN!n',ytitle='T!ie!n',chars=2,thick=3

;plot,psn,rmsp[*,32,-1],xtitle='!7W!3!iN!n',ytitle='p',chars=2,thick=3

plot,psn,rmsps[*,32,-1]/max(rmsps[*,32,-1]),xtitle='!7W!3!iN!n',ytitle='!7w!3',chars=2,thick=3

plot,psn,rmsvp[*,32,-1]/max(rmsvp[*,32,-1]),xtitle='!7W!3!iN!n',ytitle='V!ii||!n',chars=2,thick=3

tek_color
