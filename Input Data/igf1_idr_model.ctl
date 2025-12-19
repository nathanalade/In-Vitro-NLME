$PROBLEM 2-Compartment Oral PK with Type I Indirect Response Model
; Model: Inhibition of IGF-1 Production
; Response: IGF-1 concentration (ng/mL)
; PK: 2-compartment with first-order absorption
; MU-referencing implemented for efficient SAEM/IMP estimation
;
; Literature basis for IGF-1 turnover:
;   - Circulating IGF-1 t1/2 = 12-16h (ternary complex bound)
;   - Baseline = 150-200 ng/mL in healthy adults
;   - Kout derived from t1/2 = ln(2)/t1/2 ≈ 0.05-0.12 h^-1
;   - Zhou et al. 2000, Thorsted et al. 2016

$INPUT ID TIME AMT RATE CMT EVID MDV DV WT DOSE

$DATA igf1_pkpd_data.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9

$MODEL
  NCOMP=4
  COMP=(DEPOT)      ; CMT 1: Absorption compartment
  COMP=(CENTRAL)    ; CMT 2: Central PK compartment
  COMP=(PERIPH)     ; CMT 3: Peripheral PK compartment  
  COMP=(IGF1)       ; CMT 4: IGF-1 response compartment

$PK
;-----------------------------------------------------------------
; MU-Referenced Parameters for Efficient Estimation
;-----------------------------------------------------------------

; --- PK Parameters ---
  MU_1 = THETA(1)                    ; Log Ka
  KA = EXP(MU_1 + ETA(1))
  
  MU_2 = THETA(2)                    ; Log CL
  TVCL = EXP(MU_2)
  CL = TVCL * (WT/70)**0.75 * EXP(ETA(2))
  
  MU_3 = THETA(3)                    ; Log V1 (central volume)
  TVV1 = EXP(MU_3)
  V1 = TVV1 * (WT/70) * EXP(ETA(3))
  
  MU_4 = THETA(4)                    ; Log Q
  TVQ = EXP(MU_4)
  Q = TVQ * (WT/70)**0.75 * EXP(ETA(4))
  
  MU_5 = THETA(5)                    ; Log V2 (peripheral volume)
  TVV2 = EXP(MU_5)
  V2 = TVV2 * (WT/70) * EXP(ETA(5))
  
; --- PD Parameters (Indirect Response) ---
  MU_6 = THETA(6)                    ; Log Kout
  KOUT = EXP(MU_6 + ETA(6))
  
  MU_7 = THETA(7)                    ; Log Baseline IGF-1
  BASE = EXP(MU_7 + ETA(7))
  
  ; Kin derived from steady-state: Kin = Kout * Baseline
  KIN = KOUT * BASE
  
  MU_8 = THETA(8)                    ; Logit Imax (constrained 0-1)
  IMAX = EXP(MU_8 + ETA(8))/(1 + EXP(MU_8 + ETA(8)))
  
  MU_9 = THETA(9)                    ; Log IC50
  IC50 = EXP(MU_9 + ETA(9))

; --- Micro-rate constants ---
  S2 = V1                            ; Scaling for central compartment
  K20 = CL/V1                        ; Elimination rate constant
  K23 = Q/V1                         ; Central to peripheral
  K32 = Q/V2                         ; Peripheral to central

; --- Initialize IGF-1 compartment at baseline ---
  A_0(4) = BASE

$DES
;-----------------------------------------------------------------
; Differential Equations
;-----------------------------------------------------------------

; Depot compartment (first-order absorption)
  DADT(1) = -KA*A(1)

; Central compartment
  DADT(2) = KA*A(1) - K20*A(2) - K23*A(2) + K32*A(3)

; Peripheral compartment  
  DADT(3) = K23*A(2) - K32*A(3)

; Drug concentration in central compartment (ng/mL)
; Note: A(2) in mg, V1 in L -> mg/L = ug/mL; *1000 -> ng/mL
  CP = A(2)/V1 * 1000
  IF(CP.LT.0) CP = 0

; Type I Indirect Response: Inhibition of production
; dR/dt = Kin*(1 - Imax*Cp/(IC50+Cp)) - Kout*R
  IF(CP.GT.0) THEN
    INH = IMAX*CP/(IC50 + CP)
  ELSE
    INH = 0
  ENDIF
  DADT(4) = KIN*(1 - INH) - KOUT*A(4)

$ERROR (OBSERVATION ONLY)
;-----------------------------------------------------------------
; Error Model
;-----------------------------------------------------------------

; PK predictions (ng/mL)
  CPRED = A(2)/V1 * 1000
  IF(CPRED.LT.0) CPRED = 0
  
; PD predictions (IGF-1 ng/mL)
  RPRED = A(4)
  IF(RPRED.LT.0) RPRED = 0

; Combined proportional + additive error for PK
  IF(CMT.EQ.2) THEN
    IPRED = CPRED
    W = SQRT(THETA(10)**2 * IPRED**2 + THETA(11)**2)
    IF(W.LT.0.001) W = 0.001
    IRES = DV - IPRED
    IWRES = IRES/W
    Y = IPRED + W*EPS(1)
  ENDIF

; Combined proportional + additive error for PD
  IF(CMT.EQ.4) THEN
    IPRED = RPRED
    W = SQRT(THETA(12)**2 * IPRED**2 + THETA(13)**2)
    IF(W.LT.0.001) W = 0.001
    IRES = DV - IPRED
    IWRES = IRES/W
    Y = IPRED + W*EPS(2)
  ENDIF

$THETA
;-----------------------------------------------------------------
; Initial Estimates - Based on Simulation Truth
;-----------------------------------------------------------------
; PK Parameters (log-transformed)
  (-0.22)    ; 1 - Log Ka (h^-1), TV = 0.8
  (1.61)     ; 2 - Log CL (L/h), TV = 5.0
  (3.40)     ; 3 - Log V1 (L), TV = 30
  (1.10)     ; 4 - Log Q (L/h), TV = 3.0
  (4.09)     ; 5 - Log V2 (L), TV = 60

; PD Parameters
  (-1.90)    ; 6 - Log Kout (h^-1), TV = 0.15 (t1/2 ~ 4.6h)
  (5.16)     ; 7 - Log Baseline (ng/mL), TV = 175
  (2.94)     ; 8 - Logit Imax, TV = 0.95
  (6.40)     ; 9 - Log IC50 (ng/mL), TV = 600

; Residual Error Parameters
  (0.10)     ; 10 - Proportional error PK (CV ~10%)
  (0.5)      ; 11 - Additive error PK (ng/mL)
  (0.08)     ; 12 - Proportional error PD (CV ~8%)
  (3)        ; 13 - Additive error PD (ng/mL)

$OMEGA BLOCK(5)
;-----------------------------------------------------------------
; Between-Subject Variability - PK (Correlated)
;-----------------------------------------------------------------
  0.0625                             ; ETA(1) Ka, CV ~25%
  0.01 0.0625                        ; ETA(2) CL, CV ~25%
  0.01 0.02 0.0625                   ; ETA(3) V1, CV ~25%
  0.005 0.005 0.005 0.04             ; ETA(4) Q, CV ~20%
  0.005 0.005 0.005 0.01 0.04        ; ETA(5) V2, CV ~20%

$OMEGA 
;-----------------------------------------------------------------
; Between-Subject Variability - PD (Diagonal)
;-----------------------------------------------------------------
  0.0225                             ; ETA(6) Kout, CV ~15%
  0.04                               ; ETA(7) Baseline, CV ~20%
  0.0025                             ; ETA(8) Imax, CV ~5%
  0.09                               ; ETA(9) IC50, CV ~30%

$SIGMA
;-----------------------------------------------------------------
; Residual Variability
;-----------------------------------------------------------------
  1 FIX                              ; EPS(1) for PK (scaled in $ERROR)
  1 FIX                              ; EPS(2) for PD (scaled in $ERROR)

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT 
            SIGL=9 NSIG=3 MSFO=igf1_idr.msf

$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE ID TIME AMT CMT EVID MDV DV WT DOSE
       KA CL V1 Q V2 KOUT BASE KIN IMAX IC50 CP
       IPRED CWRES IWRES 
       ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9
       NOPRINT ONEHEADER FILE=igf1_sdtab.csv

$TABLE ID KA CL V1 Q V2 KOUT BASE KIN IMAX IC50 WT DOSE
       ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9
       NOPRINT ONEHEADER FIRSTONLY FILE=igf1_patab.csv
