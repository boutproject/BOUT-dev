PRO Poincare_Main, zArr,xArr,cArr,tindex=tindex,time=time,rational=rational,$
             nstart=nstart, nmap=nmap, nsubcycle=nsubcycle, $
             xmin=xmin, xmax=xmax, savefile=savefile, debug=debug
;

COMMON griddata, g, deltaZtor, Ntor

;None of these should be needed since they are all set in Poincare.pro
;IF NOT KEYWORD_SET(nstart) THEN nstart=10
;IF NOT KEYWORD_SET(nmap) THEN nmap=50
;IF NOT KEYWORD_SET(nsubcycle) THEN nsubcycle=1
;IF NOT KEYWORD_SET(savefile) THEN savefile='puncture_plot.idl.dat'
;IF NOT KEYWORD_SET(xmin) THEN xmin=MIN(g.psixy)
;IF NOT KEYWORD_SET(xmax) THEN xmax=MAX(g.psixy)

zmin=0.0
zmax=deltaZtor ;;2*!PI/period

iColor=0
zArr=[0.0]
xArr=[0.0]
cArr=[0]

for istart=1, nStart-2 do begin
    ;;-select flux surface
    if keyword_set(RATIONAL) then begin
        qval=(18+iStart)/15.
        xStart=INTERPOL(xxarr,qqarr,qval)
    endif else begin
        xStart = xmin+(xmax-xmin)*DOUBLE(istart)/DOUBLE(nStart-1)
    endelse

    for isub=1,Nsubcycle do begin

        ;;-select random z-coordinate
        zStart=zMax*RANDOMU(seed)

        iColor=iColor+1
        color=(iColor mod 13) + 1

        x=xStart+0.01*(isub-2)
        z=zStart

        for imap=0,Nmap-1 do begin
            FieldTrace, xIn=x, zIn=z, xout=xNew, zout=zNew, debug=debug
            zNew = zNew MOD zmax ; Periodic
            x=xNew
            z=zNew

            ;;-store data
            zArr=[zArr,zNew]
            xArr=[xArr,xNew]
            cArr=[cArr,color]

        endfor ;-imap
    endfor ;-isub
endfor ;-istart

allPts={x:xArr,z:zArr,c:cArr}
xminplot=xmin
xmaxplot=xmax
save, allPts,zmin,zmax,xminplot,xmaxplot,time,tindex, f=savefile
;
END
