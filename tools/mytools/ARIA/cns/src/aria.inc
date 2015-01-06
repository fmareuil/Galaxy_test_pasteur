C =======
C aria.inc
C =======
C
C  data structure (COMMON block) for ARIA 
C  note: most parameters are still in noe.inc
C
C Some physical constants   
C =======================
C Note: constants are in SI units, exponents separately
C
C     GYR     gyromagnetic ration of the proton
C     HBAR    hbar
C
      DOUBLE PRECISION MGYR
      INTEGER EGYR
      PARAMETER (MGYR = 2.6751965D0, EGYR=8)  
C
      DOUBLE PRECISION MHBAR
      INTEGER EHBAR
      PARAMETER (MHBAR = 1.0545919D0, EHBAR=-34)
C
C for conjuate gradient minimzer (not used at present)
      INTEGER CGSTEP, CICYCL, MCICYC
      DOUBLE PRECISION TOLER, MCALTR, STEP
C
      INTEGER NMRATM
      LOGICAL ARINIT
C
C ambiguous assignment statistics
      INTEGER HPNHGL, PNATOM, HPNMAT, HPNMA2, MAXSEL, MATDIM,
     &        AMBMAX, AMBMIN, AMBSTR, RNGMIN, RNGMAX
      DOUBLE PRECISION AMBLEV, AMBCUT
      CHARACTER*4 AMBMOD,AVEMOD,AMBFRM
C
C array dimensions
      INTEGER DIMMAT,DIMNHG,DIMATM
C
C calibration matrix
      INTEGER CALCNT, MAXCAL, HPNCAL
      PARAMETER (MAXCAL=20)
      DOUBLE PRECISION VREFER(MAXCAL,MAXCAL), DREFER(MAXCAL,MAXCAL)
      DOUBLE PRECISION VCOUNT(MAXCAL,MAXCAL), DCOUNT(MAXCAL,MAXCAL)
      DOUBLE PRECISION CALFAC(MAXCAL,MAXCAL)
      DOUBLE PRECISION CALEXP, CALDIS
C
      INTEGER HPDMAT,HPIMAT,HPPOPU,HPIGR,HPRATE,HPHHPR,HPNPR
C
C pseudo atom corrections (not used)
      INTEGER HPNPSE 
C hsqc intensities to correct backcalculated noe
      INTEGER HPNDIA
C violation statistics
      DOUBLE PRECISION VIOTHR, VIOEXC
      INTEGER VIOSTR
C array dimensions (number of spins etc)
      INTEGER NHG, NHY


C
C COMMON BLOCK
C
C     double precision block
      COMMON /ARIA/  TOLER, MCALTR, STEP,     
     &               AMBLEV, AMBCUT,
     &              VREFER, DREFER, CALEXP, CALDIS,
     &              VCOUNT, DCOUNT, CALFAC, VIOTHR, VIOEXC
C     integer block
      COMMON /IARIA/  
     &    CGSTEP, CICYCL, MCICYC, 
     &    NMRATM, MAXSEL, MATDIM, NHG, NHY,
     &    AMBMIN, AMBMAX, AMBSTR,CALCNT, VIOSTR,
     &    RNGMIN, RNGMAX, 
     &    HPNHGL, HPNCAL, HPNDIA, HPNPSE,  HPPOPU, HPIGR, HPNPR,
     &    HPDMAT, HPIMAT, HPRATE, HPHHPR,
     &    HPNMAT, HPNMA2, DIMMAT, DIMNHG, DIMATM

C
C     character string block
      COMMON /CARIA/ AMBMOD, AVEMOD, AMBFRM
C
      logical block
      COMMON /LARIA/ ARINIT
C
C
      SAVE /ARIA/
      SAVE /IARIA/
      SAVE /CARIA/
      SAVE /LARIA/
