PRO Poincare_Main, zArr,xArr,cArr,tindex=tindex,time=time,rational=rational,$
             nstart=nstart, nmap=nmap, nsubcycle=nsubcycle, $
             xmin=xmin, xmax=xmax, savefile=savefile, debug=debug
;

COMMON griddata, g, deltaZtor, Ntor
COMMON flags, flag, mc_flag

zmin=0.0
zmax=deltaZtor ;;2*!PI/period

iColor=0
zArr=[0.0]
xArr=[0.0]
cArr=[0]

for istart=0, nStart-1 do begin
    ;;-select flux surface
    if keyword_set(RATIONAL) then begin
        qval=(18+iStart)/15.
        xStart=INTERPOL(xxarr,qqarr,qval)
    endif else begin
        xStart = xmin+(xmax-xmin)*DOUBLE(istart)/DOUBLE(nStart-1)
    endelse

;(JPS) changed istart=0,nstart-1 from istart=1,nstart-2
;but need to make sure that xstart doesn't exit computational domain  
    IF istart EQ 0 THEN xStart=xmin+(xmax-xmin)/(2*DOUBLE(nstart-1))
    IF istart EQ nstart-1 THEN xStart=xmax-(xmax-xmin)/(2*DOUBLE(nstart-1))

    for isub=1,Nsubcycle do begin

        ;;-select random z-coordinate
;        zStart=zMax*RANDOMU(seed)
        zStart=0.25*zMax

        iColor=iColor+1
        color=(iColor mod 13) + 1

        x=xStart+0.01*(isub-1)   ; (JPS) Can give x values below xmin for large nstart (CORRECTED)
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
save, allPts,zmin,zmax,xminplot,xmaxplot,time,tindex,flag,mc_flag, f=savefile
;
END
