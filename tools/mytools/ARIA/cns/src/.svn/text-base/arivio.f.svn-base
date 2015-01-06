C     
C=====================================================================
C     
      SUBROUTINE NOEEXC
     1     (ICL,NOENUM,VIOEXC,NSTRUC,MODE,VIOLOL,VIOUPL,
     1     NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     1     NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     1     NOEPP1,NOEPP2,NOEHP1,NOEHP2,NOEWGH,NOEVIO,NUMNOE,
     1     NOESTP,NOECV,NOEXCL,NOEPID,NOECDI)
C     
C     plan: get frequency differences into the scheme
C     1) calculate all exp(frequency differences) and store
C     (e.g. NOEWGT = NOEWGT * ...)
C     this might be done in ariass somewhere.
C     need to use standard distance-weighted frequency of assignment
C     2) normalize such that the smallest difference gets weight 1
C     3) exclude if probability is > 0.5
C     or cumulative figure of merit < 0.5
C     P(Nviol,DelFreq) = Nviol/Nstruct  / (F(DelFreq) * InitWeight)
C     
C     excludes datapoints for which NOEVIO exceeds threshold
C     VIOEXC*NSTRUC by setting NOEXCL to 0
C     
C     Author: Michael Nilges
C     
      IMPLICIT NONE
C     I/O
C     
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
C     
      INTEGER NOENUM,ICL,NSTRUC
      DOUBLE PRECISION VIOEXC,VIOLOL,VIOUPL
      CHARACTER*4 MODE
C     
C     global NOE arrays on HEAP:
C     restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C     classes
      INTEGER NOECND(*)
C     averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C     target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C     chemical shifts
      DOUBLE PRECISION NOEPP1(*),NOEPP2(*),NOEHP1(*),NOEHP2(*)
C     volume, weights
      DOUBLE PRECISION NOEVOL(*),NOEWGH(*)
C     number of violations
      INTEGER NOEVIO(*),NUMNOE(*)
C     time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEXCL(*), NOEPID(*)
      DOUBLE PRECISION NOECDI(*)
C     
C     local
      INTEGER N, NEXC
      DOUBLE PRECISION AWEIGH,AFIGME,FIGMER
      INTEGER NAVE, NFIX
      CHARACTER*4 CEXCL
C     begin
C     
      IF (NOENUM.EQ.0) THEN
         WRITE(6,'(A)') ' ARIVIO-ERR: no restraints'
         RETURN
      END IF
C     
      IF (MODE.EQ.'INIT') THEN
         DO N=1,NOENUM
            NOEXCL(N)=1
            NOEVIO(N)=0
            NUMNOE(N)=0
            NOERRV(N)=ZERO
            NOERAV(N)=ZERO
         END DO
         WRITE (6,'(A)') ' ARIVIO: database initialized'
         RETURN
      END IF
C     
C     else calculate a figure of merit for the peak from the initial weight,
C     the frequency figure of merit, and the violation statistic.
      AWEIGH=ZERO
      AFIGME=ZERO
      NAVE=0
      NFIX=0
      DO N=1,NOENUM
         IF (NOECND(N).EQ.ICL) THEN
            IF (NOEWGH(N).LE.ONE) THEN
               NAVE=NAVE+1
               AWEIGH=AWEIGH+NOEWGH(N)
            ELSE
               NFIX=NFIX+1
            END IF
         END IF
      END DO
C
      IF (NAVE.EQ.0) THEN     
         AWEIGH = ONE
      ELSE
         AWEIGH = AWEIGH/NAVE
      END IF
C     
      WRITE(6,'(A)') 
     &     ' COUNtviolations: average restraint properties '
      WRITE(6,'(A,I5)') 
     &     '     number of fixed restraints: ', NFIX
      WRITE(6,'(A,I5)') 
     &     '     number of soft restraints: ', NAVE
         WRITE(6,'(A,F15.3)')
     &     '     average initial weight: ', AWEIGH
      NAVE=0
      DO N=1,NOENUM
         IF (NOECND(N).EQ.ICL) THEN
            IF (NOERAV(N).GT.ZERO) THEN
               NAVE=NAVE+1
               AFIGME=AFIGME+NOERAV(N)
            END IF
         END IF
      END DO
      IF (NAVE.EQ.0) THEN
         AFIGME = ONE
      ELSE
         AFIGME=AFIGME/NAVE
      END IF
      WRITE(6,'(A,F15.3)')
     &     '     average figure of merit: ', AFIGME
C     
      NEXC=0
      DO N=1,NOENUM
         IF (NOECND(N).EQ.ICL) THEN
            FIGMER=ONE-FLOAT(NOEVIO(N))/MAX(1,NUMNOE(N))
            IF (NOERAV(N).GT.ZERO) THEN
               FIGMER=FIGMER*NOERAV(N) /AFIGME
            END IF
            NOEXCL(N) = 1
            IF (NOEWGH(N).LE.ONE) THEN
               FIGMER=FIGMER*NOEWGH(N)/AWEIGH
            END IF
            IF (FIGMER.LE.(ONE-VIOEXC)) THEN
               IF (MODE.EQ.'EXCL') THEN
                  IF (NOEWGH(N).LE.ONE) THEN
                     NEXC=NEXC+1
                     NOEXCL(N) = 0
                     CEXCL = 'EXCL'
                  ELSE
                     CEXCL = 'KEPT'
                  END IF
                  WRITE(6,'(A,I5,3E15.5,X,I5,X,A)')
     $                 ' ARIEXC: peak, lower, upper, distance',
     $                 NOEPID(N),NOEDIS(N)-NOELOW(N),NOEDIS(N)
     $                 +NOEHIG(N),
     $                 NOECDI(N)**(-SIXTH),NOEVIO(N),CEXCL
               ELSEIF (MODE.EQ.'SET ') THEN
                  NEXC=NEXC+1
                  NOEXCL(N)=1
C     
C     if TOKEN was specified on input, VIOxxL are -1
                  IF (VIOUPL.GE.ZERO) THEN
                     NOEHIG(N) = MAX(ZERO,(VIOUPL-NOEDIS(N)))
                  END IF
                  IF (VIOLOL.GE.ZERO) THEN
                     NOELOW(N) = MAX(ZERO,(NOEDIS(N)-VIOLOL))
                  END IF
                  WRITE(6,'(A,I6, F15.3,F15.3)')
     &                 ' COUNtviolations: peak, new error estimates:',
     $                 NOEPID(N), NOELOW(N), NOEHIG(N)
               ELSEIF (MODE.EQ.'MULT') THEN
                  NEXC=NEXC+1
                  NOEXCL(N)=1
C     
C     if TOKEN was specified on input, VIOxxL are -1
                  IF (VIOUPL.GE.ZERO) THEN
                     NOEHIG(N) = NOEHIG(N)*VIOUPL
                  END IF
                  IF (VIOLOL.GE.ZERO) THEN
                     NOELOW(N) = NOELOW(N)*VIOLOL
                  END IF
                  WRITE(6,'(A,I6, F15.3,F15.3)')
     &                 ' COUNtviolations: peak, new error estimates:',
     $                 NOEPID(N), NOELOW(N), NOEHIG(N)
C     
               ELSEIF (MODE.EQ.'WEIG') THEN
                  NOEXCL(N)=1
                  NOEWGH(N)=EXP(-FIGMER**2)
               ELSE
                  WRITE(6,'(A)')
     &                 ' %COUNTV-ERR: unknown mode'
               END IF
            END IF
         END IF
      END DO
C     
      IF (MODE.EQ.'EXCL') THEN
         WRITE (6,'(A,I5,A)')
     &        ' COUNtviolations: ', NEXC, ' restraints excluded '
      ELSE
         WRITE (6,'(A,I5,A)')
     &        ' COUNtviolations: ', NEXC, ' limits reset '
      END IF
C     
      RETURN
      END
C     
C=====================================================================
C     
      SUBROUTINE ARIVIO
     1     (NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     1     NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     1     NOEPP1,NOEPP2,NOEHP1,NOEHP2,NOEWGH,NOEVIO,NUMNOE,
     1     NOESTP,NOECV,NOEXCL,NOEPID)
C     
C     
C     subroutine counts noe violations and calculates average violation
C     
C     Author: Michael Nilges
C     
      IMPLICIT NONE
C     I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'aria.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
C     global NOE arrays on HEAP
C     restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C     classes
      INTEGER NOECND(*)
C     averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C     target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C     volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
      DOUBLE PRECISION NOEHP1(*),NOEHP2(*)
C     weights
      DOUBLE PRECISION NOEWGH(*)
C     number of violations
      INTEGER NOEVIO(*),NUMNOE(*)
C     time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEXCL(*), NOEPID(*)
C     
C     local
      DOUBLE PRECISION EN, SAVEXP, RN, RN2, LW, DDF, DEVI, WENS
      DOUBLE PRECISION XIJ, YIJ, ZIJ, RIJ, SIJ, XI, YI, ZI, XJ, YJ, ZJ
      INTEGER N, NVIO, NCV, L1, L2, INORR, NORR
      INTEGER I, J, II, JJ, LI, LJ
C     
C     begin
C     
      NVIO=0
      IF (NOENUM.EQ.0) THEN
         RETURN
      END IF
C     
      DO N=1,NOENUM
         IF (NOECV(N).NE.NOEICV) THEN
            NUMNOE(N)=NUMNOE(N)+1
            L1=0
            RN=ZERO
            INORR=0
            DO NORR=NOEORR(N)+1, NOEORR(N+1)
               INORR=INORR+1
C     
               IF (NOEAVE(NOECND(N)).EQ.NOEM6
     &              .OR.(NOEAVE(NOECND(N)).EQ.NOEM3)
     &              .OR.(NOEAVE(NOECND(N)).EQ.NOESUM)) THEN
                  IF (NOEAVE(NOECND(N)).EQ.NOEM3) THEN
                     SAVEXP=THREE
                  ELSEIF (NOEAVE(NOECND(N)).EQ.NOEM6) THEN
                     SAVEXP=SIX
                  ELSE
                     SAVEXP=NAVEXP(NOECND(N))
                  END IF
C     
                  L2=0
                  LW=ZERO
                  RN2=ZERO
C     
C     loop over all pairs of atoms belonging to restraint N
                  DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
                     DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                        I=NOEILS(II)
                        J=NOEJLS(JJ)
                        XIJ=X(I)-X(J)
                        YIJ=Y(I)-Y(J)
                        ZIJ=Z(I)-Z(J)
                        SIJ=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
                        RIJ=ONE/(SIJ**(SAVEXP/TWO))
                        RN2=RN2+RIJ
                        LW=LW+ONE
                        L1=L1+1
                        L2=L2+1
                        RIJ=RIJ/SIJ
                     END DO
                  END DO
C     
                  IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
                     IF (SAVEXP.EQ.OREXP(NOECND(N))) THEN
                        RN2=MAX(RN2/NMONO(NOECND(N)),RSMALL)
                     ELSE
                        RN2=MAX((NMONO(NOECND(N))/RN2)**(ONE/SAVEXP),
     $                       RSMALL)
                     END IF
                  ELSE
                     RN2=MAX((LW/RN2)**(ONE/SAVEXP),RSMALL)
                  END IF
C     
               ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
                  LI=0
                  XI=ZERO
                  YI=ZERO
                  ZI=ZERO
                  DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
                     I=NOEILS(II)
                     XI=XI+X(I)
                     YI=YI+Y(I)
                     ZI=ZI+Z(I)
                     LI=LI+1
                  END DO
                  XI=XI/LI
                  YI=YI/LI
                  ZI=ZI/LI
                  LJ=0
                  XJ=ZERO
                  YJ=ZERO
                  ZJ=ZERO
                  DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                     J=NOEJLS(JJ)
                     XJ=XJ+X(J)
                     YJ=YJ+Y(J)
                     ZJ=ZJ+Z(J)
                     LJ=LJ+1
                  END DO
                  XJ=XJ/LJ
                  YJ=YJ/LJ
                  ZJ=ZJ/LJ
                  XIJ=XI-XJ
                  YIJ=YI-YJ
                  ZIJ=ZI-ZJ
                  RN2=MAX(SQRT(XIJ**2+YIJ**2+ZIJ**2),RSMALL)
               END IF
C     
               IF (INORR.GT.1) THEN
                  IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
                     IF (SAVEXP.EQ.OREXP(NOECND(N))) THEN
                        RN=RN+RN2
                     ELSE
                        RN=RN+RN2**(-OREXP(NOECND(N)))
                     END IF
                  ELSE
                     RN=RN+RN2**(-OREXP(NOECND(N)))
                  END IF
               ELSE
                  RN=RN2
               END IF
            END DO
C     
            RN2=MAX((ONE/RN)**(ONE/OREXP(NOECND(N))),RSMALL)
            RN=RN2
C     
            DEVI=MAX(ZERO,RN-NOEDIS(N)-
     &           MAX(ZERO,NOEHIG(N)-NOEOFF(NOECND(N))))+
     &           MIN(ZERO,RN-NOEDIS(N)+
     &           MAX(ZERO,NOELOW(N)-NOEMOF(NOECND(N))))
            NOERRV(N)=NOERRV(N)+DEVI**2
C     
            IF (ABS(DEVI).GT.VIOTHR) THEN
               WRITE(6,'(A,I5,3E15.5)')
     $              ' ARIEXC: peak, lower, upper, distance',
     $              NOEPID(N), NOEDIS(N)-NOELOW(N),
     $              NOEDIS(N)+NOEHIG(N),
     $              RN
               NOEVIO(N)=NOEVIO(N)+1
               NVIO=NVIO+1
            END IF
         END IF
      END DO
      WRITE(6,'(A,I5,A,f15.5)')
     &     ' COUNTV:',NVIO,' violations above ', VIOTHR
C     
      RETURN
      END
C     
C=====================================================================
C     
      SUBROUTINE NVIOLS
     1     (NOENUM,VIOEXC,NSTRUC,
     1     NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND, 
     1     NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL, 
     1     NOEPP1,NOEPP2,NOEHP1,NOEHP2,NOEWGH,NOEVIO,NUMNOE,
     1     NOESTP,NOECV,NOEXCL,NOEPID)
C     
C     lists datapoints for which NOEVIO exceeds threshold
C     VIOEXC*NSTRUC
C     Author: Michael Nilges
C     
      IMPLICIT NONE
C     I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      INTEGER NOENUM,NSTRUC
      DOUBLE PRECISION VIOEXC
C     
C     global NOE arrays on HEAP
C     restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C     classes
      INTEGER NOECND(*)
C     averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C     target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C     volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
      DOUBLE PRECISION NOEHP1(*),NOEHP2(*)
C     weights
      DOUBLE PRECISION NOEWGH(*)
C     number of violations
      INTEGER NOEVIO(*),NUMNOE(*)
C     time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEXCL(*), NOEPID(*)
C     
C     local
      INTEGER I,J,K,N, NORR
C     
C     begin
C     
      IF (NOENUM.GT.0) THEN
         DO N=1,NOENUM
            IF (NOEXCL(N).EQ.0) THEN
               WRITE(PUNIT,'(A,I5,A)')
     &              ' ========== restraint ',NOEPID(N),' =========='
               DO NORR=NOEORR(N)+1, NOEORR(N+1)
                  WRITE(PUNIT,'(A)') ' set-i-atoms'
                  DO I=NOEIPR(NORR)+1,NOEIPR(NORR+1)
                     K=NOEILS(I)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),
     $                    RES(K),TYPE(K)
                  END DO
                  WRITE(PUNIT,'(A)') ' set-j-atoms'
                  DO J=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                     K=NOEJLS(J)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),
     $                    RES(K),TYPE(K)
                  END DO
               END DO
C     
               WRITE(PUNIT,'(I5,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)')
     &             NOEVIO(N),
     $              ' Violations Above Threshold, Rms Violation=',
     &              Sqrt(NOERRV(N)/MAX(NOEVIO(N),1)),
     &              ' Target distance=',NOEDIS(N),
     &              ' (-',NOELOW(N),'/+',NOEHIG(N),')'
            END IF
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
