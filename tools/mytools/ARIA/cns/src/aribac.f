C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RXPAIR3
     $     (DIM, NLAYERS, HGLNK, DUMMY, NPAIRS, HHPAIR, CUTOFF,
     $     DMAT, RXAVEX, NJMAX)
C     
C     calculate list of J spins that interact with I spins.
C     we do this by calculating several layers around each atom:
C     first, get all JG that interact with IG.
C     then, get all KG that interact with all JG etc
C     in contrast to RXPAIRS, this routine uses an average DMAT distance 
C     matrix accumulated from an ensemble of structures by AMBACC
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
C     
      INTEGER DIM, NLAYERS, DUMMY(*), HGLNK(*), 
     $     NPAIRS(*), HHPAIR(DIM,*), NJMAX
      DOUBLE PRECISION DMAT(DIM,*), RXAVEX, CUTOFF
C     
      INTEGER IG, JG, KG, ILAYER, ILIST, JLIST, NIJPR, IAT, JAT,
     $     NIJPR2, NIJPR1
      DOUBLE PRECISION DISIJ, MAXDIS
      LOGICAL WARN, FOUND
C     
C     begin
C     
C     
C     the routine first sets up a square mask of 1 and 0
C     and then condenses the rows
C     
C     initialize the pair array
C     
      DO IG=1,DIM
         DO JG=1,DIM
            HHPAIR(IG,JG)=0
         END DO
         HHPAIR(IG,IG)=1
      END DO
C     
C     first, calculate the first layer of neighbours
C     
      NJMAX=0
      WARN=.FALSE.
      DO IG=1,DIM-1
         DO JG=IG+1,DIM
            DISIJ=DMAT(IG,JG)
            IF (DISIJ.LT.CUTOFF) THEN
               HHPAIR(IG,JG)=1
               HHPAIR(JG,IG)=1
            END IF
         END DO
      END DO
C     
C     loop over layers
C     
      IF (NLAYERS.GT.1) THEN
         DO ILAYER=2,NLAYERS
            DO IG=1,DIM
               DO KG=1,DIM
                  DUMMY(KG)=0
               END DO
               DO JG=IG,DIM
                  IF (HHPAIR(IG,JG).GT.0) THEN
                     DO KG=1,DIM
                        DUMMY(KG)=DUMMY(KG)+HHPAIR(JG,KG)
                     END DO
                  END IF
               END DO
               DO KG=1,DIM
                  HHPAIR(IG,KG)=HHPAIR(IG,KG)+DUMMY(KG)
               END DO
            END DO
         END DO
      END IF
C     
C     now condense HHPAIR
C     
      MAXDIS = ZERO
      DO IG=1,DIM
         NIJPR=0
         DO KG=1,DIM
            DUMMY(KG)=HHPAIR(IG,KG)
            HHPAIR(IG,KG)=0
         END DO
         DO JG=1,DIM
            IF (DUMMY(JG).GT.0) THEN
               NIJPR=NIJPR+1
               HHPAIR(IG,NIJPR)=JG
               IF (DMAT(IG,JG).GT.MAXDIS) THEN
                  MAXDIS = DMAT(IG,JG)
                  WRITE(6,'(2I5,E15.5)') IG,JG,MAXDIS
               END IF                  
            END IF
         END DO
         NPAIRS(IG)=NIJPR
         IF (NIJPR.LE.1) WARN=.TRUE.
         NJMAX=MAX(NIJPR,NJMAX)
      END DO
C     
      WRITE(6,'(A,I5,E15.5)') ' RXPAIR2: maximum Npair, distance ',
     $     NJMAX,MAXDIS
      IF (WARN) THEN
         WRITE(6,'(A)') ' RXPAIR-WRN: no interactions for some spins'
      END IF
C     
      RETURN
      END
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RXPAIRS
     $     (DIM, NATOM, NLAYERS, HGLNK, IGRAT, NPAIRS, HHPAIR, CUTOFF,
     $     X, Y, Z, DMAT, RXAVEX, NJMAX)
C     
C     calculate list of J spins that interact with I spins.
C     we do this by calculating several layers around each atom:
C     first, get all JG that interact with IG.
C     then, get all KG that interact with all JG etc
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
C     
      INTEGER DIM, NATOM, NLAYERS, IGRAT(*), HGLNK(*), 
     $     NPAIRS(*), HHPAIR(DIM,*), NJMAX
      DOUBLE PRECISION DMAT(DIM,*), RXAVEX, X(*), Y(*), Z(*), CUTOFF
C     
      INTEGER IG, JG, ILAYER, ILIST, JLIST, NIJPR, IAT, JAT, NIJPR2,
     $     NIJPR1
      DOUBLE PRECISION DISIJ, MAXDIS
      LOGICAL WARN, FOUND
C     
C     begin
C     
C     
      NJMAX=0
      WARN=.FALSE.
      DO IG=1,DIM
         NIJPR2=1
         NIJPR1=1
         NIJPR =1
C     enter atom itself on list
         HHPAIR(IG,NIJPR)=IG
         MAXDIS = ZERO
C     loop over layers
         DO ILAYER=1,NLAYERS
            NIJPR1=NIJPR2
            NIJPR2=NIJPR
            DO ILIST=NIJPR1,NIJPR2
               IAT=IGRAT(HHPAIR(IG,ILIST))
               DO JG=1,DIM
                  JAT=IGRAT(JG)
                  IAT=IGRAT(IG)
C     first calculate the distance between IG and JG
                  CALL RM6DIS(IAT,JAT,HGLNK,X,Y,Z,RXAVEX,DISIJ)
                  DMAT(IG,JG)=DISIJ
C     then the distances to JG from the neighbours of IG
                  IAT=IGRAT(HHPAIR(IG,ILIST))
                  CALL RM6DIS(IAT,JAT,HGLNK,X,Y,Z,RXAVEX,DISIJ)
                  IF (DISIJ.LT.CUTOFF) THEN
C     check if the atom is already on the list
                     FOUND=.FALSE.
                     DO JLIST=1,NIJPR
                        IF (HHPAIR(IG,JLIST).EQ.JG) FOUND = .TRUE.
                     END DO
C     if not yet present, enter on the list
                     IF (.NOT. FOUND) THEN
                        NIJPR=NIJPR+1
                        HHPAIR(IG,NIJPR)=JG
                        MAXDIS = MAX(MAXDIS, DMAT(IG,JG))       
                     END IF
                  END IF
               END DO
            END DO
         END DO
         NPAIRS(IG)=NIJPR
         IF (NIJPR.LE.1) WARN=.TRUE.
         NJMAX=MAX(NIJPR,NJMAX)
C         WRITE (6, '(A, 2I5, E15.5)') ' RXPAIR: IG, NIJPR, DISIJ ',
C     $        IG, NIJPR, MAXDIS
      END DO
      WRITE(6,'(A,I5,E15.5)') ' RXPAIR: maximum Npair, distance ',
     $     NJMAX,MAXDIS
      IF (WARN) THEN
         WRITE(6,'(A)') ' RXPAIR-WRN: no interactions for some spins'
      END IF
C     
      RETURN
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RXPAIR2
     $     (DIM, NLAYERS, HGLNK, IGRAT, NPAIRS, HHPAIR, CUTOFF,
     $     DMAT, RXAVEX, NJMAX)
C     
C     calculate list of J spins that interact with I spins.
C     we do this by calculating several layers around each atom:
C     first, get all JG that interact with IG.
C     then, get all KG that interact with all JG etc
C     in contrast to RXPAIRS, this routine uses an average DMAT distance 
C     matrix accumulated from an ensemble of structures by AMBACC
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
C     
      INTEGER DIM, NLAYERS, IGRAT(*), HGLNK(*), 
     $     NPAIRS(*), HHPAIR(DIM,*), NJMAX
      DOUBLE PRECISION DMAT(DIM,*), RXAVEX, CUTOFF
C     
      INTEGER IG, JG, ILAYER, ILIST, JLIST, NIJPR, IAT, JAT, NIJPR2,
     $     NIJPR1
      DOUBLE PRECISION DISIJ, MAXDIS
      LOGICAL WARN, FOUND
C     
C     begin
C     
C     
      NJMAX=0
      WARN=.FALSE.
      DO IG=1,DIM
         NIJPR2=1
         NIJPR1=1
         NIJPR =1
C     enter atom itself on list
         HHPAIR(IG,NIJPR)=IG
         MAXDIS = ZERO
C     loop over layers
         DO ILAYER=1,NLAYERS
            NIJPR1=NIJPR2
            NIJPR2=NIJPR
            DO ILIST=NIJPR1,NIJPR2
               IAT=IGRAT(HHPAIR(IG,ILIST))
               DO JG=1,DIM
                  DISIJ=DMAT(HHPAIR(IG,ILIST),JG)
                  IF (DISIJ.LT.CUTOFF) THEN
C     found one - check if the atom is already on the list
                     FOUND=.FALSE.
                     DO JLIST=1,NIJPR
                        IF (HHPAIR(IG,JLIST).EQ.JG) FOUND = .TRUE.
                     END DO
C     ok, enter on the list
                     IF (.NOT. FOUND) THEN
                        NIJPR=NIJPR+1
                        HHPAIR(IG,NIJPR)=JG
                        MAXDIS = MAX(MAXDIS, DMAT(IG,JG))                  
                     END IF
                  END IF
               END DO
            END DO
         END DO
         NPAIRS(IG)=NIJPR
         IF (NIJPR.LE.1) WARN=.TRUE.
         NJMAX=MAX(NIJPR,NJMAX)
      END DO
      WRITE(6,'(A,I5,E15.5)') ' RXPAIR2: maximum Npair, distance ',
     $     NJMAX,MAXDIS
      IF (WARN) THEN
         WRITE(6,'(A)') ' RXPAIR-WRN: no interactions for some spins'
      END IF
C     
      RETURN
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72

      SUBROUTINE CNTGRP(NATOM,HGLNK,NHG2,NHY2)
C
C     subroutine counts the number of groups and hydrogens in the system
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C
C     global
C      INCLUDE 'aria.inc'
C     
C     input/output
      INTEGER NATOM,HGLNK(*),NHG2,NHY2
C     
C     local
      INTEGER I,I1,POPI
      LOGICAL QNEXT
C     
C     begin
C     
C     count number of groups present and maximum population
      NHG2=0
      NHY2=0
      DO I=1,NATOM
         IF (HGLNK(I).GE.I) THEN
            NHG2=NHG2+1
            I1=I
            POPI=0
            QNEXT=.TRUE.
            DO WHILE (QNEXT)
               POPI=POPI+1
               CALL NEXTHY(I1, HGLNK, QNEXT)
            END DO
C$$$  EQVMAX=MAX(EQVMAX,POPI)
         END IF
         IF (HGLNK(I).GT.0) NHY2=NHY2+1
      END DO
      WRITE(6,'(A,2I5)') ' CNTGRP: NHG2, NHY2 ', NHG2, NHY2
C     
      RETURN
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE GRFPOP(NATOM,HGLNK,GRPOP)
C     
C     returns sqrt(number of hydros) in each group
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C     
C     global
      INCLUDE 'numbers.inc'
C     
C     input/ output
      INTEGER NATOM, HGLNK(*)
      DOUBLE PRECISION GRPOP(*)
C     
C     local
      INTEGER I,I1,NGRP
      DOUBLE PRECISION POPI
C     
C     begin
C     
      NGRP=0
      DO I=1,NATOM
         IF (HGLNK(I).GE.I) THEN
            NGRP=NGRP+1
            I1=I
C first atom in group
            POPI=ONE
C other atoms
            DO WHILE (HGLNK(I1).GT.I1)
               POPI=POPI+ONE
               I1=HGLNK(I1)     
            END DO
            GRPOP(NGRP)=POPI
         END IF
      END DO
C     
      RETURN
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE GRAT(NATOM,HGLNK,IGRAT)
C     
C     sets up pointer array  igrf --> ihy
C     (<=> HGLNK(IAT) >= I)
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
C     input/ output
      INTEGER NATOM, HGLNK(*), IGRAT(*)
C     
C     local
      INTEGER I,NGRP
C     
C     begin
C     
      NGRP=0
      DO I=1,NATOM
         IF (HGLNK(I).GE.I) THEN
            NGRP=NGRP+1
            IGRAT(NGRP)=I
         END IF
      END DO
C     
      RETURN
      END
C
C
      SUBROUTINE RM6DIS(I,J,HGLNK,X,Y,Z,EXPO,DSIX)
C     
C     calculates the weighted distance between two groups of atoms
C     if the exponent is 0, the distance between the geometric
C     averages of the groups is calculated
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C
C     global
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
C     
C     input/ output
      INTEGER I,J,HGLNK(*)
      DOUBLE PRECISION X(*),Y(*),Z(*),EXPO,DSIX
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C     
C     local
      INTEGER II,JJ,L,LZ,LI,LJ
      DOUBLE PRECISION XIJ,YIJ,ZIJ,RIJ,SIJ
      DOUBLE PRECISION XIAVE,YIAVE,ZIAVE,XJAVE,YJAVE,ZJAVE
      DOUBLE PRECISION MINDIS
      LOGICAL QNEXTI, QNEXTJ
C     
C     begin
C     
C set the minimum distance to 1.0 (to the inverse sixth power)
      MINDIS=1.0

      IF ((I.EQ.0).OR.(J.EQ.0)
     $     .OR.(.NOT.FRSTEL(I,HGLNK))
     $     .OR.(.NOT.FRSTEL(J,HGLNK))) THEN
         WRITE(6, '(A,2I8)') ' %RM6DIS: error in atom number',
     $        I, J
         DSIX=9999.0D0
C     
      ELSE
         IF (EXPO.EQ.ZERO) THEN
C     generate average coordinates for the pseudo atoms
            II=I
            JJ=J
            XIAVE=ZERO
            YIAVE=ZERO
            ZIAVE=ZERO
            LI=0
            QNEXTI=.TRUE.
            DO WHILE (QNEXTI)
               XIAVE=XIAVE+X(II)
               YIAVE=YIAVE+Y(II)
               ZIAVE=ZIAVE+Z(II)
               LI=LI+1
            END DO
            XIAVE=XIAVE/FLOAT(LI)
            YIAVE=YIAVE/FLOAT(LI)
            ZIAVE=ZIAVE/FLOAT(LI)
C     
            XJAVE=ZERO
            YJAVE=ZERO
            ZJAVE=ZERO
            LJ=0
            QNEXTJ=.TRUE.
            DO WHILE (QNEXTJ)
               XJAVE=XJAVE+X(JJ)
               YJAVE=YJAVE+Y(JJ)
               ZJAVE=ZJAVE+Z(JJ)
               LJ=LJ+1
            END DO
            XJAVE=XJAVE/FLOAT(LJ)
            YJAVE=YJAVE/FLOAT(LJ)
            ZJAVE=ZJAVE/FLOAT(LJ)
C     
            DSIX=SQRT((XIAVE-XJAVE)**2
     $           +(YIAVE-YJAVE)**2+(ZIAVE-ZJAVE)**2)
C     
         ELSE
C     
            DSIX=ZERO
            II=I
            JJ=J
            L=0
            LZ=0
C     loop over all pairs of atoms
            QNEXTI=.TRUE.
            DO WHILE (QNEXTI)
               QNEXTJ=.TRUE.
               DO WHILE (QNEXTJ)
                  XIJ = X(II)-X(JJ)
                  YIJ = Y(II)-Y(JJ)
                  ZIJ = Z(II)-Z(JJ)
                  SIJ = XIJ**2+YIJ**2+ZIJ**2
                  IF (SIJ.GT.RSMALL) THEN
                     RIJ = SIJ**(-EXPO/TWO)
                     DSIX=DSIX+RIJ
                     L=L+1
                  ELSE
                     RIJ=-ONE
                     LZ=LZ+1
                     DSIX = MINDIS
                  END IF     
                  CALL NEXTHY(JJ, HGLNK, QNEXTJ)
               END DO
               CALL NEXTHY(II, HGLNK, QNEXTI)
            END DO
C     
            IF (L.GT.0) THEN
               DSIX=DSIX/FLOAT(L) ! inverse distances
               DSIX=DSIX**(ONE/EXPO)
            END IF
C     
            IF (DSIX.GT.RSMALL) THEN
               DSIX=ONE/DSIX
            ELSEIF (L.EQ.0) THEN
               DSIX=ZERO
            ELSE
               DSIX=9999.0D0
               WRITE(6,'(A)')
     $              ' %RM6DIS-ERR: infinite distance, set to 9999'
            END IF
         END IF
      END IF
C     
      RETURN
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RMATRIX
     $     (DIM,NATOM, HGLNK,GRPOP,
     $     SUMTRP, RATE,DMAT,HHPAIR,
     $     NPAIR, OMEGA,ZLEAK,MTAU)
C     
C     Calculates the complete rate matrix
C     Ref.: Keepers, J.W., and James, T.L. (1984), J. Magn. Reson. 57, 404.
C     
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C
      IMPLICIT NONE
C
C     global variables
      INCLUDE 'cns.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'aria.inc'
C     
C     input/ output
      INTEGER NATOM, DIM
      INTEGER HGLNK(*)
      DOUBLE PRECISION GRPOP(*)
      DOUBLE PRECISION RATE(DIM,*), DMAT(DIM,*)
      INTEGER HHPAIR(DIM,*),NPAIR(*)
      DOUBLE PRECISION OMEGA,ZLEAK,MTAU
C     
C     workspace
      DOUBLE PRECISION SUMTRP(*)
C     
C     local variables
      INTEGER I,IG,IH,JH,JG,IJ,L
      INTEGER II,JJ,LI,LJ,LR
      DOUBLE PRECISION TAU, TSWS, TS2WS
      DOUBLE PRECISION POPISQ, POPJSQ, POPIJ
      DOUBLE PRECISION XIJ,YIJ,ZIJ,SIJ,RIJ,DIJ,Q,DSIX,DDSIX
      DOUBLE PRECISION OMSQ, SORD
      DOUBLE PRECISION S0,S1,S2
      DOUBLE PRECISION CON0,CON1,CON2,TRPI,TRPJ,DRT,DTRI,DTRJ
      DOUBLE PRECISION CRAT
      INTEGER ETAU,EOM,EDIST,EMU
      LOGICAL DISJUNCT
      CHARACTER*6 PROCNAME
      INTEGER NI,JJG
C     
C     begin
C     
      SORD=ONE
C

      CON0=MGYR**4*MHBAR**2
      CON1=3.0D0*CON0
      CON2=6.0D0*CON0
      OMSQ=2.0D0*PI*OMEGA*1.0D-9
      OMSQ=OMSQ*OMSQ
      WRITE(6,'(A,2E15.5)') ' RMATRIX: tau, omega: ', MTAU,OMEGA
C     
C     correction factors
      EOM=9                     !GHz-->Hz
      ETAU=-9                   !ns-->s
      EDIST=-10                 !A -->m
      EMU=-14                   !Mu0**2/16pi**2
      CRAT=1.0D0                !10**(4*EGYR+2*EHBAR+ETAU-6*EDIST+EMU-1)
C     
C     initialise relaxation matrix etc.
      DO IG=1,DIM
         DO JG=1,DIM
            RATE(JG,IG) = ZERO
         END DO
      END DO
C     
      DO IG=1,DIM
         SUMTRP(IG) = ZERO
      END DO
C     
      IJ=0
      DO IG=1,DIM
         DO JG=1,DIM
C
C get spectral densities
c$$$            CALL TAUGET(IG,JG,HGLNK,MTAU,SORD)
            TAU=MTAU*1.0D9
            TSWS=TAU*TAU*OMSQ
            TS2WS=FOUR*TSWS
            S0 = SORD*CON0*TAU
            S1 = SORD*CON1*TAU/(ONE+TSWS)
            S2 = SORD*CON2*TAU/(ONE+TS2WS)
C     
C     relaxation rate and transition probabilities
C     diagonal elements
            IF (IG.EQ.JG) THEN
               IF (GRPOP(IG).GT.ONE) THEN 
                  Q=DMAT(IG,IG)**(-SIX)
                  DRT=TWO*(GRPOP(IG)-ONE)*(HALF*S1+S2)*CRAT
                  RATE(IG,IG)=Q*DRT
                  TRPI=ZERO
                  TRPJ=ZERO
               ELSE
                  RATE(IG,IG)=ZERO
                  TRPI=ZERO
                  TRPJ=ZERO
               END IF 
C     
C     off-diagonal elements and accumulate transition probabilities
            ELSE
               Q=DMAT(IG,JG)**(-SIX)
               DRT=(S2-S0)*GRPOP(IG)*CRAT
               RATE(IG,JG)=Q*DRT
               DTRI=(S1+S2+S0)*GRPOP(JG)*CRAT
               TRPI=Q*DTRI
               SUMTRP(IG)=SUMTRP(IG)+TRPI
            ENDIF
         END DO
      END DO
C     
C     add z leakage and transition probabilities to diagonal terms
      DO I=1,DIM
         RATE(I,I)=RATE(I,I)+SUMTRP(I)+ZLEAK
      END DO
C
      WRITE(6, '(A)') ' RMATRIX: constructed relaxation matrix '
C     
      RETURN                         
      END
C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RXIBLQ
     $     (NHG,NJMAX,NTAUP1,ABSZIS,GRPOP,DMAT,RTMAT,CALVOL,
     $     INTMAT,HHPAIR,NPAIR,TAUMIX)
C     
C     calculate a sparse spectrum using numerical integration
C     of a sparse relaxation matrix.
C     the time intervals are squared
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
C     global variables
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
C     
C     input/ output  
C     
      INTEGER NHG,NJMAX,NTAUP1
      DOUBLE PRECISION RTMAT(NHG,*), DMAT(NHG,*), INTMAT(NHG,*), 
     $     CALVOL(NHG,NHG,NTAUP1)
      INTEGER  HHPAIR(NHG,*), NPAIR(*)
      DOUBLE PRECISION ABSZIS(*),GRPOP(*)
      DOUBLE PRECISION TAUMIX
C     
C     local variables               
      INTEGER           I,J,NI,NJ,JJ,LL,L,T,NTAUP12
      DOUBLE PRECISION  DELTAU, INTBEG, SUMDIS, SUMVOL,
     $     MAXVOL, MINDIS
      CHARACTER*6       PROCNAME
C     
C     begin
C     
      DO T=1,NTAUP1
         DO J=1,NHG
            DO I=1,NHG
               CALVOL(I,J,T)=ZERO
            END DO
         END DO
      END DO
C     
      ABSZIS(1)=ZERO
      ABSZIS(2)=ONE
      DO T=3,NTAUP1
         ABSZIS(T)=ABSZIS(T-1)*2
      END DO
      DO T=1,NTAUP1
         ABSZIS(T)=ABSZIS(T)/ABSZIS(NTAUP1)*TAUMIX
      END DO
C     
      DELTAU=ABSZIS(2)-ABSZIS(1)
      WRITE(6,'(A,E15.5)') ' RXIBLQ: integration interval, '
     $     , DELTAU
C     
C     rewrite RTMAT to be (I - dtR)
      DO J=1,NHG
         DO I=1,NHG
            RTMAT(I,J)=-DELTAU*RTMAT(I,J)
         END DO
      END DO
      DO I=1,NHG
         RTMAT(I,I)=RTMAT(I,I)+ONE
      END DO
C     
C     first integration step
      T=2
      SUMVOL=ZERO
      DO I=1,NHG
         DO JJ=1,NPAIR(I)
            J=HHPAIR(I,JJ)
            CALVOL(I,J,T)=RTMAT(I,J)
            IF (I.NE.J) THEN
               SUMVOL=SUMVOL+CALVOL(I,J,T)
            END IF
         END DO
      END DO
      WRITE(6, '(A,I5,2E15.5)')
     $     ' RXIBLQ: step, mixing time, sum of volumes: ',
     $     T, ABSZIS(T), SUMVOL
C     
      DO T=3,NTAUP1
         SUMVOL=ZERO
         SUMDIS=ZERO
         DO I=1,NHG
            DO JJ=1,NPAIR(I)
               J=HHPAIR(I,JJ)
               DO LL=1,NPAIR(I)
                  L=HHPAIR(I,LL)
                  CALVOL(I,J,T)= CALVOL(I,J,T)
     $                 +CALVOL(I,L,T-1)*CALVOL(L,J,T-1)
               END DO
               IF (I.NE.J) THEN
                  SUMDIS=SUMDIS+DMAT(I,J)**(-SIX)
                  SUMVOL=SUMVOL+CALVOL(I,J,T)
               END IF
            END DO
         END DO
         WRITE(6, '(A,I5,2E15.5)')
     $        ' RXIBLQ: step, mixing time, sum of volumes: ',
     $        T, ABSZIS(T), SUMVOL
      END DO
C     
C     multiply with the populations (i.e. the magnetization at zero)
      DO T=1,NTAUP1
         SUMVOL=ZERO
         SUMDIS=ZERO
         MINDIS=ZERO
         MAXVOL=-R4BIG
         DO I=1,NHG
            DO JJ=1,NPAIR(I)
               J=HHPAIR(I,JJ)
               CALVOL(I,J,T)=
     $              CALVOL(I,J,T)*GRPOP(J)
               IF (I.NE.J) THEN
                  SUMDIS=SUMDIS+GRPOP(I)*GRPOP(J)*DMAT(I,J)**(-SIX)
                  SUMVOL=SUMVOL+CALVOL(I,J,T)
                  MINDIS=MAX(MINDIS,
     $                 GRPOP(I)*GRPOP(J)*DMAT(I,J)**(-SIX))
                  MAXVOL=MAX(MAXVOL,CALVOL(I,J,T))
               END IF
            END DO
         END DO
      END DO
C     
      DO I=1,NHG
         DO J=1,NHG
            INTMAT(J,I)=ZERO
         END DO
      END DO
C     scale intensities such that they come out ~ d^-6
      DO I=1, NHG
         DO JJ=1,NPAIR(I)
            J=HHPAIR(I,JJ)
            IF (J.GT.NHG) THEN
               WRITE(6,'(A)') ' RXIBLQ: programming error, nhg'
            ELSE
               INTMAT(I,J)=MINDIS/MAXVOL*CALVOL(I,J,NTAUP1)
            END IF
         END DO
      END DO
C     
      WRITE(6, '(A,2E15.5)')
     $     ' RXIBLQ: sum, max of calculated volumes: ',
     $     SUMVOL, MAXVOL
      WRITE(6, '(A,2E15.5)')
     $     ' RXIBLQ: sum, max of calculated r^-6: ',
     $     SUMDIS, MINDIS
C     
      RETURN
      END



C=======================================================================
C==== 6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE RXIBL2
     $     (NHG,NJMAX,NTAUP1,ABSZIS,GRPOP,DMAT,RTMAT,CALVOL,
     $     INTMAT,HHPAIR,NPAIR,TAUMIX)
C     
C     calculate a sparse spectrum using numerical integration
C     of a sparse relaxation matrix.
C     the time intervals are squared
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
C     global variables
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
C     
C     input/ output  
C     
      INTEGER NHG,NJMAX,NTAUP1
      DOUBLE PRECISION RTMAT(NHG,*), DMAT(NHG,*), INTMAT(NHG,*), 
     $     CALVOL(NHG,NJMAX,NTAUP1)
      INTEGER  HHPAIR(NHG,*), NPAIR(*)
      DOUBLE PRECISION ABSZIS(*),GRPOP(*)
      DOUBLE PRECISION TAUMIX
C     
C     local variables               
      INTEGER           I,J,NI,NJ,JJ,LL,L,T,NTAUP12
      DOUBLE PRECISION  DELTAU, INTBEG, SUMDIS, SUMVOL,
     $     MAXVOL, MINDIS
      CHARACTER*6       PROCNAME
C     
C     begin
C     
      DO T=1,NTAUP1
         DO J=1,NJMAX
            DO I=1,NHG
               CALVOL(I,J,T)=ZERO
            END DO
         END DO
      END DO
C     
      ABSZIS(1)=ZERO
      ABSZIS(2)=ONE
      DO T=3,NTAUP1
         ABSZIS(T)=ABSZIS(T-1)*2
      END DO
      DO T=1,NTAUP1
         ABSZIS(T)=ABSZIS(T)/ABSZIS(NTAUP1)*TAUMIX
      END DO
C     
      DELTAU=ABSZIS(2)-ABSZIS(1)
      WRITE(6,'(A,E15.5)') ' RXIBLQ: integration interval, '
     $     , DELTAU
C     
C     rewrite RTMAT to be (I - dtR)
      DO J=1,NHG
         DO I=1,NHG
            RTMAT(I,J)=-DELTAU*RTMAT(I,J)
         END DO
      END DO
      DO I=1,NHG
         RTMAT(I,I)=RTMAT(I,I)+ONE
      END DO
C     
C     first integration step
      T=2
      SUMVOL=ZERO
      DO I=1,NHG
         DO J=1,NHG
            CALVOL(I,J,T)=RTMAT(I,J)
            IF (I.NE.J) THEN
               SUMVOL=SUMVOL+CALVOL(I,J,T)
            END IF
         END DO
      END DO
      WRITE(6, '(A,I5,2E15.5)')
     $     ' RXIBLQ: step, mixing time, sum of volumes: ',
     $     T, ABSZIS(T), SUMVOL
C     
      DO T=3,NTAUP1
         SUMVOL=ZERO
         SUMDIS=ZERO
         DO I=1,NHG
            DO J=1,NHG
               DO L=1,NHG
                  CALVOL(I,J,T)= CALVOL(I,J,T)
     $                 +CALVOL(I,L,T-1)*CALVOL(L,J,T-1)
               END DO
               IF (I.NE.J) THEN
                  SUMDIS=SUMDIS+DMAT(I,J)**(-SIX)
                  SUMVOL=SUMVOL+CALVOL(I,J,T)
               END IF
            END DO
         END DO
         WRITE(6, '(A,I5,2E15.5)')
     $        ' RXIBLQ: step, mixing time, sum of volumes: ',
     $        T, ABSZIS(T), SUMVOL
      END DO
C     
C     multiply with the populations (i.e. the magnetization at zero)
      DO T=1,NTAUP1
         SUMVOL=ZERO
         SUMDIS=ZERO
         MINDIS=R4BIG
         MAXVOL=-R4BIG
         DO I=1,NHG
            DO J=1,NHG
               CALVOL(I,J,T)=
     $              CALVOL(I,J,T)*GRPOP(J)
               IF (I.NE.J) THEN
                  SUMDIS=SUMDIS+DMAT(I,J)**(-SIX)
                  SUMVOL=SUMVOL+CALVOL(I,J,T)
                  MINDIS=MIN(MINDIS,DMAT(I,J))
                  MAXVOL=MAX(MAXVOL,CALVOL(I,J,T))
               END IF
            END DO
         END DO
      END DO
      MINDIS=MINDIS**(-SIX)
C     
      DO I=1,NHG
         DO J=1,NHG
            INTMAT(J,I)=ZERO
         END DO
      END DO
C     scale intensities such that they come out ~ d^-6
      DO I=1, NHG
         DO J=1,NHG
            INTMAT(I,J)=MINDIS/MAXVOL*CALVOL(I,J,NTAUP1)
         END DO
      END DO
C     
      WRITE(6, '(A,2E15.5)')
     $     ' RXIBLQ: sum, max of calculated volumes: ',
     $     SUMVOL, MAXVOL
      WRITE(6, '(A,2E15.5)')
     $     ' RXIBLQ: sum, max of calculated r^-6: ',
     $     SUMDIS, MINDIS
C     
      RETURN
      END
