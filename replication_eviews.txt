```eviews
'===============================================================
' EViews 13 — Final Replication: PRQI/NRI, ARDL–Bounds, ECM
' Study: Do Resource–Energy Channels Moderate Unemployment under Economic Shocks? Evidence from Iran
'===============================================================

'--- 0) Workfile setup
close @all
wfcreate(wf=iran_ue, a) 1991 2024
setmaxerrs 1

' Keep Excel next to this .prg file
string progdir = @runpath
string xlsfile = progdir + "Iran's_Macro_Data.xlsx"

'--- 1) Import data (no cleaning; expects headers on first row)
' Expected columns: YEAR, U, CPI, OUTPUTGAP, REER, NETMIGk, NRD, AGRILAND, WATERSTRESSINDEX, CF, LOSSES, CAPPERPERSON, SHOCK
import(type=xlsx, mode=update, namepos=firstrow) {xlsfile}

' Align sample window
smpl 1991 2024

' Alias mapping (non-destructive)
series AGRI_LAND = AGRILAND
series WSI       = WATERSTRESSINDEX
series CAP_PC    = CAPPERPERSON

'==========================
' 2) Composite indices via PCA (correlation matrix)
'==========================

' 2.1) Z-scores
series z_cf       = (CF        - @mean(CF))        / @stdev(CF)         ' desirable
series z_cap_pc   = (CAP_PC    - @mean(CAP_PC))    / @stdev(CAP_PC)     ' desirable
series z_losses   = (LOSSES    - @mean(LOSSES))    / @stdev(LOSSES)     ' undesirable

series z_agri     = (AGRI_LAND - @mean(AGRI_LAND)) / @stdev(AGRI_LAND)  ' desirable
series z_wsi      = (WSI       - @mean(WSI))       / @stdev(WSI)        ' undesirable
series z_nrd      = (NRD       - @mean(NRD))       / @stdev(NRD)        ' undesirable

' Flip undesirable directions
series z_lossesInv = -z_losses
series z_wsiInv    = -z_wsi
series z_nrdInv    = -z_nrd

' 2.2) PCA — PRQI = PC1 of [z_cf, z_cap_pc, z_lossesInv]
group g_prqi z_cf z_cap_pc z_lossesInv
g_prqi.pcomp(ncomp=1, cor)
delete(noerr) PRQI_PC1
if @isobject("PC1") then
  series PRQI_PC1 = PC1
  delete(noerr) PC1
endif

' 2.3) PCA — NRI = PC1 of [z_agri, z_wsiInv, z_nrdInv]
group g_nri z_agri z_wsiInv z_nrdInv
g_nri.pcomp(ncomp=1, cor)
delete(noerr) NRI_PC1
if @isobject("PC1") then
  series NRI_PC1 = PC1
  delete(noerr) PC1
endif

' 2.4) Orient signs (anchors: PRQI↔z_cf, NRI↔z_agri)
scalar sgn_prqi = @sign(@cor(PRQI_PC1, z_cf))
series PRQI     = sgn_prqi*PRQI_PC1
scalar sgn_nri  = @sign(@cor(NRI_PC1,  z_agri))
series NRI      = sgn_nri*NRI_PC1

' Standardized indices for interactions
series zPRQI = (PRQI - @mean(PRQI)) / @stdev(PRQI)
series zNRI  = (NRI  - @mean(NRI))  / @stdev(NRI)

' Optional: factor "loadings" (post-orientation; correlations with index)
scalar L_PRQI_CF       = @cor(z_cf,        PRQI)
scalar L_PRQI_CAP      = @cor(z_cap_pc,    PRQI)
scalar L_PRQI_LOSSES   = @cor(z_lossesInv, PRQI)
scalar L_NRI_AGL       = @cor(z_agri,      NRI)
scalar L_NRI_WSIinv    = @cor(z_wsiInv,    NRI)
scalar L_NRI_NRDinv    = @cor(z_nrdInv,    NRI)

'==========================
' 3) Mean-centered sanctions & interactions
'==========================
series ShockC   = SHOCK - @mean(SHOCK)
series INT_PRQI = ShockC * zPRQI
series INT_NRI  = ShockC * zNRI

'==========================
' 4) Unit-root checks (views)
'==========================
show U.adf
show U.kpss
show CPI.adf
show CPI.kpss
show OUTPUTGAP.adf
show OUTPUTGAP.kpss
show REER.adf
show REER.kpss
show NETMIGk.adf
show NETMIGk.kpss
show NRD.adf
show NRD.kpss
show AGRI_LAND.adf
show AGRI_LAND.kpss
show WSI.adf
show WSI.kpss

'==========================
' 5) ARDL (levels) & Pesaran Bounds
'==========================
equation eq_ardl.ardl U PRQI NRI CPI OUTPUTGAP REER NETMIGk ShockC INT_PRQI INT_NRI, aic maxlag(1)
eq_ardl.bounds
freeze(tab_bounds) eq_ardl.bounds
show tab_bounds

'==========================
' 6) Long-run LS (HAC) and ECM
'==========================
equation eq_lr.ls(cov=hac, kernel=bartlett, lag=2) U c PRQI NRI CPI OUTPUTGAP REER NETMIGk ShockC INT_PRQI INT_NRI
freeze(tab_longrun) eq_lr.output
show tab_longrun

series EC = eq_lr.@resid

equation eq_ecm.ls(cov=hac, kernel=bartlett, lag=2) d(U) c EC(-1) d(U(-1)) _
  d(CPI) d(CPI(-1)) d(OUTPUTGAP) d(OUTPUTGAP(-1)) d(REER) d(REER(-1)) _
  d(NETMIGk) d(NETMIGk(-1)) d(PRQI) d(PRQI(-1)) d(NRI) d(NRI(-1)) _
  d(ShockC) d(ShockC(-1)) d(INT_PRQI) d(INT_PRQI(-1)) d(INT_NRI)
freeze(tab_ecm) eq_ecm.output
show tab_ecm

'==========================
' 7) Long-run sanction effect scenarios (from eq_lr)
'==========================
scalar b7 = eq_lr.@coefs("ShockC")
scalar b8 = eq_lr.@coefs("INT_PRQI")
scalar b9 = eq_lr.@coefs("INT_NRI")

scalar s1 = b7 + b8*(-1) + b9*(-1)   ' PRQI −1 sd, NRI −1 sd
scalar s2 = b7 + b8*(-1) + b9*(+1)   ' PRQI −1 sd, NRI +1 sd
scalar s3 = b7 + b8*(+1) + b9*(-1)   ' PRQI +1 sd, NRI −1 sd
scalar s4 = b7 + b8*(+1) + b9*(+1)   ' PRQI +1 sd, NRI +1 sd
scalar s5 = b7                        ' Indices at mean (0,0)

table tab_scen 6 2
tab_scen(1,1) = "Scenario"                       : tab_scen(1,2) = "Long-run effect of sanctions on U"
tab_scen(2,1) = "PRQI -1 sd, NRI -1 sd"         : tab_scen(2,2) = @str(s1, "f(10.3)")
tab_scen(3,1) = "PRQI -1 sd, NRI +1 sd"         : tab_scen(3,2) = @str(s2, "f(10.3)")
tab_scen(4,1) = "PRQI +1 sd, NRI -1 sd"         : tab_scen(4,2) = @str(s3, "f(10.3)")
tab_scen(5,1) = "PRQI +1 sd, NRI +1 sd"         : tab_scen(5,2) = @str(s4, "f(10.3)")
tab_scen(6,1) = "Indices at mean (0,0)"         : tab_scen(6,2) = @str(s5, "f(10.3)")
show tab_scen

'==========================
' 8) Diagnostics & stability (Tables 7–8; Figure 2)
'==========================
show eq_ecm.residtestserial(2)     ' Breusch–Godfrey LM
show eq_ecm.white                  ' White heteroskedasticity
show eq_ecm.autoregcondhet(2)      ' ARCH–LM
show eq_ecm.normtest               ' Jarque–Bera
show eq_ecm.reset                  ' Ramsey RESET
show eq_ecm.qstat(2)               ' Ljung–Box Q
show eq_ecm.stability(cusum)       ' CUSUM
show eq_ecm.stability(cusumsq)     ' CUSUMSQ

'==========================
' 9) Save workfile with frozen tables
'==========================
save(t=all) iran_ue.wf1

'==========================
' End
'==========================
```
