; IDL routine to read an 'a' file aeqdsk
;
; Example:
;   a = read_aeqdsk("aeqdsk")
;
; Format of A-EQDSK file is specified here:
;   https://fusion.gat.com/THEORY/efit/a_eqdsk.html
; 
; Ben Dudson, University of York, Feb 2010
;

FUNCTION read_aeqdsk, file
  OPENR, fid, file, /GET_LUN, error=errid
  IF errid THEN BEGIN
    PRINT, "Error whilst reading '"+file+"' :" + !ERROR_STATE.MSG
    RETURN, 0
  ENDIF
  
  line=' '
  READF, fid, line ; first line with dates
  PRINT, line
  
  READF, fid, line ; shot number
  eshot = 0L
  READS, line, eshot

  READF, fid, line ; Time
  etime = 0.0d
  READS, line, etime
  PRINT, "Shot number: ", eshot, " @", etime
  
  READF, fid, line
  s = STRSPLIT(line, '* ',/extract)
  IF N_ELEMENTS(s) NE 9 THEN BEGIN
    PRINT, "ERROR: Expecting 9 fields on line 4"
    RETURN, 0
  ENDIF
  READS, s[0], etime
  jflag = 0 & READS, s[1], jflag
  lflag = 0 & READS, s[2], lflag
  limloc = s[3]
  SWITCH limloc OF
    "SNT": PRINT, "Equilibrium is Single Null Top"
    "SNB": PRINT, "Equilibrium is Single Null Bottom"
    "DN": PRINT, "Equilibrium is Double Null"
    ELSE: PRINT, "Equilibrium has unknown code: "+limloc
  ENDSWITCH
  mco2v  = 0 & READS, s[4], mco2v
  mco2r  = 0 & READS, s[5], mco2r
  ; s[6] = 'CLC' ?
  ;qmflag = 0 & READS, s[8], qmflag
  
  PRINT, mco2v, mco2r
  
  ; Set up generator
  status = next_double(fid=fid)
  
  tsaisq = next_double()
  rcencm = next_double()
  bcentr = next_double()
  pasmat = next_double()
  cpasma = next_double()
  rout   = next_double()
  zout   = next_double()
  aout   = next_double()
  eout   = next_double()
  doutu  = next_double()
  doutl  = next_double()
  vout   = next_double()
  rcurrt = next_double()
  zcurrt = next_double()
  qsta   = next_double()
  betat  = next_double()
  betap  = next_double()
  ali    = next_double()
  oleft  = next_double()
  oright = next_double()
  otop   = next_double()
  obott  = next_double()
  qpsi95 = next_double()
  vertn  = next_double()
  
  rco2v  = read_1d(mco2v)
  dco2v  = read_1d(mco2v)
  rco2r  = read_1d(mco2r)
  dco2r  = read_1d(mco2r)
  
  shearb = next_double()
  bpolav = next_double()
  s1     = next_double()
  s2     = next_double()
  s3     = next_double()
  qout   = next_double()
  olefs  = next_double()
  orighs = next_double()
  otops  = next_double()
  sibdry = next_double()
  areao  = next_double()
  wplasm = next_double()
  terror = next_double()
  elongm = next_double()
  qqmagx = next_double()
  cdflux = next_double()
  alpha  = next_double()
  rttt   = next_double()
  psiref = next_double()
  xndnt  = next_double()
  rseps1 = next_double()
  zseps1 = next_double()
  rseps2 = next_double()
  zseps2 = next_double()
  
  PRINT, "rseps1 = ", rseps1, "zseps1 = ", zseps1, "rseps2 = ", rseps2, "zseps1 = ", zseps2
  
  sepexp = next_double()
  obots  = next_double()
  btaxp  = next_double()
  btaxv  = next_double()
  aaq1   = next_double()
  aaq2   = next_double()
  aaq3   = next_double()
  seplim = next_double()
  rmagx  = next_double()
  zmagx  = next_double()
  simagx = next_double()
  taumhd = next_double()
  betapd = next_double()
  betatd = next_double()
  wplasmd= next_double()
  fluxx  = next_double()
  vloopt = next_double()
  taudia = next_double()
  qmerci = next_double()
  tavem  = next_double()
  
  nsilop = LONG(next_double())
  magpri = LONG(next_double())
  nfcoil = LONG(next_double())
  nesum  = LONG(next_double())
  
  PRINT, "nsilop =", nsilop, " magpri =",magpri, " nfcoil=",nfcoil, " nesum=", nesum

  csilop = read_1d(nsilop)
  cmpr2  = read_1d(magpri)
  ccbrsp = read_1d(nfcoil)
  eccurt = read_1d(nesum)
  
  pbinj  = next_double()
  rvsin  = next_double()
  zvsin  = next_double()
  rvsout = next_double()
  zvsout = next_double()
  vsurfa = next_double()
  wpdot  = next_double()
  wbdot  = next_double()
  slantu = next_double()
  slantl = next_double()
  zuperts = next_double()
  chipre = next_double()
  cjor95 = next_double()
  pp95   = next_double()
  ssep   = next_double()
  yyy2   = next_double()
  xnnc   = next_double()
  cprof  = next_double()
  oring  = next_double()
  cjor0  = next_double()
  
  fexpan = next_double()
  qqmin  = next_double()
  chigamt= next_double()
  ssi01  = next_double()
  fexpvs = next_double()
  sepnose= next_double()
  ssi95  = next_double()
  rqqmin = next_double()
  cjor99 = next_double()
  cj1ave = next_double()
  rmidin = next_double()
  rmidout= next_double()
  psurfa = next_double()
  
  FREE_LUN, fid
  
  ; Put data into structure
  
  result = {shot:eshot, time:etime, $
            mco2v:mco2v, mco2r:mco2r, $
            tsaisq:tsaisq, rcencm:rcencm, bcentr:bcentr, pasmat:pasmat, $
            cpasma:cpasma, rout:rout, zout:zout, aout:aout, eout:eout, $
            doutu:doutu, doutl:doutl, vout:vout, rcurrt:rcurrt, zcurrt:zcurrt, $
            qsta:qsta, betat:betat, betap:betap, ali:ali, oleft:oleft, $
            oright:oright, otop:otop, obott:obott, qpsi95:qpsi95, vertn:vertn, $
            rco2v:rco2v, dco2v:dco2v, rco2r:rco2r, dco2r:dco2r, $
            shearb:shearb, bpolav:bpolav, s1:s1, s2:s2, s3:s3, qout:qout, $
            olefs:olefs, orighs:orighs, otops:otops, sibdry:sibdry, $
            areao:areao, wplasm:wplasm, terror:terror, elongm:elongm, $
            qqmagx:qqmagx, cdflux:cdflux, alpha:alpha, rttt:rttt, $
            psiref:psiref, xndnt:xndnt, rseps1:rseps1, rseps2:rseps2, $
            zseps2:zseps2, sepexp:sepexp, obots:obots, btaxp:btaxp, $
            btaxv:btaxv, aaq1:aaq1, aaq2:aaq2, aaq3:aaq3, seplim:seplim, $
            rmagx:rmagx, zmagx:zmagx, simagx:simagx, taumhd:taumhd, $
            betapd:betapd, betatd:betatd, wplasmd:wplasmd, fluxx:fluxx, $
            vloopt:vloopt, taudia:taudia, qmerci:qmerci, tavem:tavem, $
            nsilop:nsilop, magpri:magpri, nfcoil:nfcoil, nesum:nesum, $
            csilop:csilop, cmpr2:cmpr2, ccbrsp:ccbrsp, eccurt:eccurt, $
            pbinj:pbinj, rvsin:rvsin, zvsin:zvsin, rvsout:rvsout, $
            zvsout:zvsout, vsurfa:vsurfa, wpdot:wpdot, wbdot:wbdot, $
            slantu:slantu, slantl:slantl, zuperts:zuperts, chipre:chipre, $
            cjor95:cjor95, pp95:pp95, ssep:ssep, yyy2:yyy2, xnnc:xnnc, $
            cprof:cprof, oring:oring, cjor0:cjor0, fexpan:fexpan, qqmin:qqmin, $
            chigamt:chigamt, ssi01:ssi01, fexpvs:fexpvs, sepnose:sepnose, $
            ssi95:ssi95, rqqmin:rqqmin, cjor99:cjor99, cj1ave:cj1ave, $
            rmidin:rmidin, rmidout:rmidout, psurfa:psurfa}
  RETURN, result
END
