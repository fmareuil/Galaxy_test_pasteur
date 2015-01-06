      SUBROUTINE ARISET
C     
C     Sets up parameters for ARIA assignment, calibration and
C     violation analysis.
C     Nilges, M., Macias, M., O'Donoghue, S., & Oschkinat, M. (1997).
C     Automated NOESY Interpretation with Ambiguous Distance Restraints
C     The Refined Solution Structure of the Pleckstrin Homology Domain
C     from beta-Spectrin.
C     J. Mol. Biol. 269, 408-422.
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
      INCLUDE 'cns.inc'
      INCLUDE 'aria.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'numbers.inc'
C     local
      INTEGER ISLCT, JSLCT, NISLCT
C     begin
C     
      IF (ARINIT) THEN
         CALL ARIHP(0,0,NATOM)
         NMRATM=NATOM
         CALL FILL4(HEAP(HPNHGL),NMRATM,1)
         CALL FILL4(HEAP(HPNPSE),NMRATM,1)
         CALL MAKIND(HEAP(HPNPSE),NMRATM,NISLCT)
         WRITE (6,'(A)') ' ARIA arrays initialized'
         ARINIT = .FALSE.
C         CALL MAKIND(NOEHGL,NMRATM,ISLCT)
      END IF
C     
      ISLCT=ALLHP(INTEG4(NATOM))
      JSLCT=ALLHP(INTEG4(NATOM))
      CALL ARISE2(HEAP(ISLCT),HEAP(JSLCT),HEAP(HPNHGL),HEAP(HPNPSE))
      CALL FREHP(JSLCT,INTEG4(NATOM))
      CALL FREHP(ISLCT,INTEG4(NATOM))
C     
      RETURN
      END
C====================================================================
      SUBROUTINE ARISE2(ISLCT,JSLCT,NOEHGL,NOEPSE)
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
C     I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'aria.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER ISLCT(*), JSLCT(*), NOEHGL(*), NOEPSE(*)
C     local
      INTEGER I, J, NISLCT, NJSLCT, NMISS, IPEAK, MATDIM2, NHG2
      INTEGER IAT, JAT, II, JJ, NIAT, NJAT
      CHARACTER*4 CLASS, ACTION
      DOUBLE PRECISION VIOUPL, VIOLOL, TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION OLDNOE, NEWNOE, FSWAP, RNUM
      LOGICAL MATCH, QNEXT
C     parameter
      DOUBLE PRECISION ZERO, ONE, SIX
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, SIX=6.0D0)
C     function
      INTEGER BONDED
      LOGICAL FRSTEL
      EXTERNAL FRSTEL, BONDED
C     
C     begin
C     
C     all defaults are set in NOERES!
      CALL PUSEND('ARIA>')
      DO WHILE (.NOT.DONE)
         CALL NEXTWD('ARIA>')
         CALL MISCOM('ARIA>',USED)
         IF (.NOT.USED) THEN
C     
            IF (WD(1:4).EQ.'HELP') THEN
               CALL CNSHELP('cns-aria')
C     
C====================================================================
            ELSE IF (WD(1:2).EQ.'DO') THEN
               CALL ARIMAN(HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &              HEAP(HPNWGH),HEAP(HPNVOL),HEAP(HPNCV))
C====================================================================
            ELSE IF (WD(1:4).EQ.'FLOA') THEN
               CALL PUSEND('FLOAting-assignment>')
               DO WHILE (.NOT.DONE)
                  CALL NEXTWD('FLOAting-assignment>')
C     
                  IF (WD(1:4).EQ.'HELP') THEN
                     CALL CNSHELP('cns-aria-floatingassignment')
C     
                  ELSE IF (WD(1:4).EQ.'DEFI') THEN
                     CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
                     IF (NISLCT.GT.0) THEN
                        NMISS=0
                        DO I=1,NATOM
                           IF ((ISLCT(I).GT.0).AND.
     $                          (.NOT.FRSTEL(I,NOEHGL)))
     $                          THEN
                              ISLCT(I)=0
                              NISLCT=NISLCT-1
                           END IF
                        END DO
                     END IF
                     IF (NISLCT.EQ.2) THEN
                        CALL MAKIND(ISLCT,NATOM,NISLCT)    
                        IAT=ISLCT(1)
                        JAT=ISLCT(2)
                        NOEPSE(IAT)=JAT
                     ELSE
                        WRITE(6,'(A,A)') ' ARIA-FLOAT-ERR: exactly',
     $                       ' two atom groups have to be selected '
                     END IF
C     
                  ELSE IF (WD(1:4).EQ.'SWAP') THEN
                     CALL NEXTF('SWAP_limit',FSWAP)
                     DO I=1,NATOM
                        CALL GGUBFS(RNUM)
                        IF ((NOEPSE(I).NE.I).AND.(RNUM.LT.FSWAP)) THEN
                           JJ=ABS(NOEPSE(I))
                           CALL ENOE(OLDNOE,'NOAN',0)
                           CALL ARIFLI(IAT,JAT,ISLCT,JSLCT,
     $                          HEAP(HPNILS), HEAP(HPNJLS), NOEHGL)
                           CALL ENOE(NEWNOE,'NOAN',0)
                           WRITE(6,'(A,2F15.5)')
     $                          ' ARIFLI energies: ', OLDNOE, NEWNOE
C     we store the information on the flipping in NOEPSE
                           IF (NEWNOE.GT.OLDNOE) THEN
                              CALL ARIFLI(IAT,JAT,ISLCT,JSLCT,
     $                             HEAP(HPNILS), HEAP(HPNJLS),NOEHGL)
                           ELSE
                              IF (NOEPSE(IAT).EQ.IAT) THEN
                                 NOEPSE(IAT)=-JAT
                              ELSE
                                 NOEPSE(IAT)=-NOEPSE(IAT)
                              END IF
                           END IF
                        END IF
                     END DO
C     
                  ELSE IF (WD(1:4).EQ.'STOR') THEN
C     
C     this writes the state of the assignments to a print file.
                     DO I=1,NATOM
                        IF (NOEPSE(I).LT.0) THEN      
                           II=I
                           NOEPSE(I)=-NOEPSE(I)
                           JJ=NOEPSE(I)
                           NOEPSE(JJ)=I
C     
                           WRITE(PUNIT,'(14A)')
     $                          ' REVE ((segid "',SEGID(II),
     &                          '" and resid ',RESID(II),
     &                          ' and name ',TYPE(II),')',
     $                          ' OR   ( segid "',SEGID(JJ),
     &                          '" and resid ',RESID(JJ),
     &                          ' and name ',TYPE(JJ),'))'
                        END IF
                     END DO
                  ELSE
                     CALL CHKEND('FLOAting-assignment>',DONE)
                  END IF
               END DO
               DONE=.FALSE.
C====================================================================
            ELSE IF (WD(1:4).EQ.'CFLI') THEN
               CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
               IF (NISLCT.GT.0) THEN
                  NMISS=0
                  DO I=1,NATOM
                     IF ((ISLCT(I).GT.0).AND.
     $                    (.NOT.FRSTEL(I,NOEHGL))
     $                    .OR..NOT.INITIA(I,X,Y,Z)) THEN
                        ISLCT(I)=0
                        NMISS=NMISS+1
                        NISLCT=NISLCT-1
                     END IF
                  END DO
               END IF
               IF (NISLCT.EQ.2) THEN
                  CALL MAKIND(ISLCT,NATOM,NISLCT)
                  II=ISLCT(1)
                  JJ=ISLCT(2)      
                  NIAT=0
                  QNEXT=.TRUE.
                  DO WHILE (QNEXT)
                     NIAT=NIAT+1
                     ISLCT(NIAT)=II
                     CALL NEXTHY(II, NOEHGL,QNEXT)
                  END DO
                  NJAT=0
                  QNEXT=.TRUE.
                  DO WHILE (QNEXT)
                     NJAT=NJAT+1
                     JSLCT(NJAT)=JJ
                     CALL NEXTHY(JJ, NOEHGL,QNEXT)
                  END DO
C     
                  IF (NIAT.EQ.NJAT) THEN
                     DO II=1,NIAT
                        IAT=ISLCT(II)
                        JAT=JSLCT(II)
                        TEMPX=X(IAT)
                        TEMPY=Y(IAT)
                        TEMPZ=Z(IAT)
                        X(IAT)=X(JAT)
                        Y(IAT)=Y(JAT)
                        Z(IAT)=Z(JAT)
                        X(JAT)=TEMPX
                        Y(JAT)=TEMPY
                        Z(JAT)=TEMPZ
                     END DO
C     this is maybe a little risky
C     the heavy atom attached to the first 
C     of the hydrogens is swapped, too
                     IF (NIAT.GT.1) THEN
                        IAT=BONDED(ISLCT(1))
                        JAT=BONDED(JSLCT(1))
                        TEMPX=X(IAT)
                        TEMPY=Y(IAT)
                        TEMPZ=Z(IAT)
                        X(IAT)=X(JAT)
                        Y(IAT)=Y(JAT)
                        Z(IAT)=Z(JAT)
                        X(JAT)=TEMPX
                        Y(JAT)=TEMPY
                        Z(JAT)=TEMPZ
                     END IF
                  END IF
               END IF
C     
            ELSE IF (WD(1:4).EQ.'VFLI') THEN
C     
               IF (NISLCT.GT.0) THEN
                  NMISS=0
                  DO I=1,NATOM
                     IF (ISLCT(I).NE.1) THEN
                     ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
                        ISLCT(I)=0
                        NMISS=NMISS+1
                        NISLCT=NISLCT-1
                     END IF
                  END DO
               END IF
               IF (NISLCT.EQ.2) THEN
                  CALL MAKIND(ISLCT,NATOM,NISLCT)
                  TEMPX=XV(ISLCT(1))
                  TEMPY=YV(ISLCT(1))
                  TEMPZ=ZV(ISLCT(1))
                  XV(ISLCT(1))=XV(ISLCT(2))
                  YV(ISLCT(1))=YV(ISLCT(2))
                  ZV(ISLCT(1))=ZV(ISLCT(2))
                  XV(ISLCT(2))=TEMPX
                  YV(ISLCT(2))=TEMPY
                  ZV(ISLCT(2))=TEMPZ
               END IF
C     
C====================================================================
            ELSE IF (WD(1:4).EQ.'FLIP') THEN
C     
               CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
               IF (NISLCT.GT.0) THEN
                  NMISS=0
                  DO I=1,NATOM
                     IF ((ISLCT(I).GT.0).AND.
     $                    (.NOT.FRSTEL(I,NOEHGL)))
     $                    THEN
                        ISLCT(I)=0
                        NISLCT=NISLCT-1
                     END IF
                  END DO
               END IF
C     
               WRITE (6, '(A,I6,A)') ' ARIFLI: ', NISLCT,
     $              ' groups selected'
               IF (NISLCT.EQ.2) THEN
                  CALL MAKIND(ISLCT,NATOM,NISLCT)    
                  IAT=ISLCT(1)
                  JAT=ISLCT(2)
                  CALL ENOE(OLDNOE,'NOAN',0)
                  CALL ARIFLI(IAT,JAT,ISLCT,JSLCT,HEAP(HPNILS),
     $                 HEAP(HPNJLS), NOEHGL)
                  CALL ENOE(NEWNOE,'NOAN',0)
                  WRITE(6,'(A,2F15.5)')
     $                 ' ARIFLI energies: ', OLDNOE, NEWNOE
C     we store the information on the flipping in NOEPSE
                  IF (NEWNOE.GT.OLDNOE) THEN
                     CALL ARIFLI(IAT,JAT,ISLCT,JSLCT,
     $                    HEAP(HPNILS), HEAP(HPNJLS),NOEHGL)
                  ELSE
                     IF (NOEPSE(IAT).EQ.IAT) THEN
                        NOEPSE(IAT)=-JAT
                     ELSE
                        NOEPSE(IAT)=-NOEPSE(IAT)
                     END IF
                  END IF
               END IF
C     
C     
C====================================================================
            ELSE IF (WD(1:4).EQ.'KEEP') THEN
C     do nothing - leave group as it is - just parse
               CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
C     
C====================================================================
            ELSE IF (WD(1:4).EQ.'REVE') THEN
C     
               CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
               IF (NISLCT.GT.0) THEN
                  NMISS=0
                  DO I=1,NATOM
                     IF ((ISLCT(I).GT.0).
     $                    AND.(.NOT.FRSTEL(I,NOEHGL))) THEN
                        ISLCT(I)=0
                        NISLCT=NISLCT-1
                     END IF
                  END DO
               END IF
C     
ccc   WRITE (6, '(A,I6,A)') ' ARIFLI: ', NISLCT, ' groups selected'
               IF (NISLCT.EQ.2) THEN
                  CALL MAKIND(ISLCT,NATOM,NISLCT)    
                  IAT=ISLCT(1)
                  JAT=ISLCT(2)
                  CALL ARIFLI(IAT,JAT,ISLCT,JSLCT,HEAP(HPNILS),
     $                 HEAP(HPNJLS),NOEHGL)
                  IF (NOEPSE(IAT).EQ.IAT) THEN
                     NOEPSE(IAT)=-JAT
                  ELSE
                     NOEPSE(IAT)=-NOEPSE(IAT)
                  END IF
               END IF
C     
C====================================================================
            ELSE IF (WD(1:4).EQ.'MODI') THEN
C     get peakid, then use truncated form of noeass
C     e.g. modify 333 test=0 end
               CALL ARIMOD(HEAP(HPNORR),
     &              HEAP(HPNIPR),
     &              HEAP(HPNILS),HEAP(HPNJPR),HEAP(HPNJLS),
     $              HEAP(HPNCND),
     &              HEAP(HPNRAV),HEAP(HPNRRV),HEAP(HPNNSP),
     $              HEAP(HPNDIS),
     &              HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNWGH),
     $              HEAP(HPNCV),
     &              HEAP(HPNVIO),
     &              HEAP(HPNPID),
     &              HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNHP1),
     $              HEAP(HPNHP2),
     &              HEAP(HPNVOL),NATOM,ISLCT,NISLCT,
     &              JSLCT,NJSLCT,X,Y,Z)
C===================================================================
            ELSE IF (WD(1:4).EQ.'?   ') THEN
               WRITE(6,'(A,I6,A,I6,A,/,A,F8.3,A,I6)')
     &              ' NOE: total number of restraints:',NOENUM,
     &              ' partitioned into ',NOECCN,' classes',
     &              ' NOE: ceiling=',NOECEI
     $              ,' current allocation=',NOEMAX
               IF (NOEICV.GT.0) THEN
                  WRITE(6,'(2A,/,A,I6)')
     &                 ' NOE: data are partitioned into working set',
     $                 ' and test set.',
     &                 ' NOE: test set number=',NOEICV
               END IF
C===================================================================
            ELSE IF (WD(1:4).EQ.'COUN') THEN
C     
C     get the class name:
               CALL PUSEND('COUNtviolations>')
               DO WHILE (.NOT.DONE)
                  CALL NEXTWD('COUNtviolations>')
C     
                  IF (WD(1:4).EQ.'HELP') THEN
                     CALL CNSHELP('cns-aria-countviolations')
C     
                  ELSE IF (WD(1:4).EQ.'RESE') THEN
                     VIOSTR = 0
                     VIOEXC = ZERO
                     VIOTHR = ZERO
                     CALL NOEEXC(0,NOENUM,VIOEXC,VIOSTR,'INIT',
     $                    ZERO,ZERO,
     1                    HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS
     $                    ),HEAP(HPNJPR),HEAP(HPNJLS)
     $                    ,HEAP(HPNCND),HEAP(HPNRAV)
     $                    ,HEAP(HPNRRV),HEAP(HPNDIS)
     $                    ,HEAP(HPNLOW),HEAP(HPNHIG)
     $                    ,HEAP(HPNVOL),HEAP(HPNPP1)
     $                    ,HEAP(HPNPP2),HEAP(HPNHP1)
     $                    ,HEAP(HPNHP2),HEAP(HPNWGH)
     $                    ,HEAP(HPNVIO),HEAP(HPNNUM)
     $                    ,HEAP(HPNNSP),HEAP(HPNCV),HEAP(HPNEXC
     $                    ),HEAP(HPNPID),HEAP(HPCDIS))
C     
                  ELSE IF (WD(1:4).EQ.'THRE') THEN
                     CALL NEXTF('THREshold=',VIOTHR)
                     VIOSTR=VIOSTR+1
                     CALL ARIVIO(HEAP(HPNORR),HEAP(HPNIPR),
     $                    HEAP(HPNILS),HEAP(HPNJPR),
     $                    HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV
     $                    ), HEAP(HPNRRV),HEAP(HPNDIS)
     $                    ,HEAP(HPNLOW),HEAP(HPNHIG)
     $                    ,HEAP(HPNVOL), HEAP(HPNPP1)
     $                    ,HEAP(HPNPP2),HEAP(HPNHP1),
     $                    HEAP(HPNHP2), HEAP(HPNWGH)
     $                    ,HEAP(HPNVIO), HEAP(HPNNUM),
     $                    HEAP(HPNNSP),HEAP(HPNCV),HEAP(HPNEXC)
     $                    , HEAP(HPNPID))
C     
                  ELSE IF (WD(1:4).EQ.'EXCL') THEN
                     CALL NEXTA4('class-name=',CLASS)
                     CALL NEXTF('EXCLude=',VIOEXC)
                     DO I=1,NOECCN
                        CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1
     $                       ,MATCH)
                        IF (MATCH) THEN
                           CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR
     $                          ,'EXCL',ZERO,ZERO,HEAP(HPNORR)
     $                          ,HEAP(HPNIPR),HEAP(HPNILS)
     $                          ,HEAP(HPNJPR),HEAP(HPNJLS)
     $                          ,HEAP(HPNCND),HEAP(HPNRAV)
     $                          ,HEAP(HPNRRV),HEAP(HPNDIS)
     $                          ,HEAP(HPNLOW),HEAP(HPNHIG)
     $                          ,HEAP(HPNVOL),HEAP(HPNPP1)
     $                          ,HEAP(HPNPP2),HEAP(HPNHP1)
     $                          ,HEAP(HPNHP2),HEAP(HPNWGH)
     $                          ,HEAP(HPNVIO),HEAP(HPNNUM)
     $                          ,HEAP(HPNNSP),HEAP(HPNCV)
     $                          ,HEAP(HPNEXC),HEAP(HPNPID)
     $                          ,HEAP(HPCDIS))
                        END IF
                     END DO
C     
                  ELSE IF (WD(1:3).EQ.'SET') THEN
                     CALL NEXTA4('class-name=',CLASS)
                     CALL NEXTF('ratio=',VIOEXC)
C     get lower limit. check for TOKEN.
                     CALL NEXTA4('lower_limit=',ACTION)
                     IF (ACTION.EQ.'TOKE') THEN
                        VIOLOL = -ONE
                     ELSE
                        CALL SAVEWD
                        CALL NEXTF('lower_limit=',VIOLOL)
                     END IF
C     get upper limit. check for TOKEN.
                     CALL NEXTA4('upper_limit=',ACTION)
                     IF (ACTION.EQ.'TOKE') THEN
                        VIOUPL = -ONE
                     ELSE
                        CALL SAVEWD
                        CALL NEXTF('upper_limit=',VIOUPL)
                     END IF
C     
                     DO I=1,NOECCN
                        CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1
     $                       ,MATCH)
                        IF (MATCH) THEN
                           CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR
     $                          ,'SET ',VIOLOL,VIOUPL
     $                          ,HEAP(HPNORR),HEAP(HPNIPR)
     $                          ,HEAP(HPNILS),HEAP(HPNJPR)
     $                          ,HEAP(HPNJLS),HEAP(HPNCND)
     $                          ,HEAP(HPNRAV),HEAP(HPNRRV)
     $                          ,HEAP(HPNDIS),HEAP(HPNLOW)
     $                          ,HEAP(HPNHIG),HEAP(HPNVOL)
     $                          ,HEAP(HPNPP1),HEAP(HPNPP2)
     $                          ,HEAP(HPNHP1),HEAP(HPNHP2)
     $                          ,HEAP(HPNWGH),HEAP(HPNVIO)
     $                          ,HEAP(HPNNUM),HEAP(HPNNSP)
     $                          ,HEAP(HPNCV),HEAP(HPNEXC)
     $                          ,HEAP(HPNPID),HEAP(HPCDIS))
                        END IF
                     END DO
C     
                  ELSE IF (WD(1:4).EQ.'MULT') THEN
                     CALL NEXTA4('class-name=',CLASS)
                     CALL NEXTF('ratio=',VIOEXC)
C     get lower limit. check for TOKEN.
                     CALL NEXTA4('lower_limit=',ACTION)
                     IF (ACTION.EQ.'TOKE') THEN
                        VIOLOL = -ONE
                     ELSE
                        CALL SAVEWD
                        CALL NEXTF('lower_limit=',VIOLOL)
                     END IF
C     get upper limit. check for TOKEN.
                     CALL NEXTA4('upper_limit=',ACTION)
                     IF (ACTION.EQ.'TOKE') THEN
                        VIOUPL = -ONE
                     ELSE
                        CALL SAVEWD
                        CALL NEXTF('upper_limit=',VIOUPL)
                     END IF
C     
                     DO I=1,NOECCN
                        CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1
     $                       ,MATCH)
                        IF (MATCH) THEN
                           CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR
     $                          ,'MULT',VIOLOL,VIOUPL
     $                          ,HEAP(HPNORR),HEAP(HPNIPR)
     $                          ,HEAP(HPNILS),HEAP(HPNJPR)
     $                          ,HEAP(HPNJLS),HEAP(HPNCND)
     $                          ,HEAP(HPNRAV),HEAP(HPNRRV)
     $                          ,HEAP(HPNDIS),HEAP(HPNLOW)
     $                          ,HEAP(HPNHIG),HEAP(HPNVOL)
     $                          ,HEAP(HPNPP1),HEAP(HPNPP2)
     $                          ,HEAP(HPNHP1),HEAP(HPNHP2)
     $                          ,HEAP(HPNWGH),HEAP(HPNVIO)
     $                          ,HEAP(HPNNUM),HEAP(HPNNSP)
     $                          ,HEAP(HPNCV),HEAP(HPNEXC)
     $                          ,HEAP(HPNPID),HEAP(HPCDIS))
                        END IF
                     END DO
C     
                  ELSE IF (WD(1:4).EQ.'WEIG') THEN
                     CALL NEXTA4('class-name=',CLASS)
                     DO I=1,NOECCN
                        CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1
     $                       ,MATCH)
                        IF (MATCH) THEN
                           CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR
     $                          ,'WEIG',VIOLOL,VIOUPL
     $                          ,HEAP(HPNORR),HEAP(HPNIPR)
     $                          ,HEAP(HPNILS),HEAP(HPNJPR)
     $                          ,HEAP(HPNJLS),HEAP(HPNCND)
     $                          ,HEAP(HPNRAV),HEAP(HPNRRV)
     $                          ,HEAP(HPNDIS),HEAP(HPNLOW)
     $                          ,HEAP(HPNHIG),HEAP(HPNVOL)
     $                          ,HEAP(HPNPP1),HEAP(HPNPP2)
     $                          ,HEAP(HPNHP1),HEAP(HPNHP2)
     $                          ,HEAP(HPNWGH),HEAP(HPNVIO)
     $                          ,HEAP(HPNNUM),HEAP(HPNNSP)
     $                          ,HEAP(HPNCV),HEAP(HPNEXC)
     $                          ,HEAP(HPNPID))
                        END IF
                     END DO
C     
                  ELSE IF (WD(1:4).EQ.'LIST') THEN
                     CALL NEXTA4('class-name=',CLASS)
                     CALL NEXTF('EXCLude=',VIOEXC)
                     DO I=1,NOECCN
                        CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1
     $                       ,MATCH)
                        IF (MATCH) THEN
                           CALL NVIOLS(NOENUM,VIOEXC,VIOSTR,
     1                          HEAP(HPNORR),HEAP(HPNIPR),
     $                          HEAP(HPNILS),HEAP(HPNJPR),
     1                          HEAP(HPNJLS),HEAP(HPNCND),
     $                          HEAP(HPNRAV),HEAP(HPNRRV),
     1                          HEAP(HPNDIS),HEAP(HPNLOW),
     $                          HEAP(HPNHIG),HEAP(HPNVOL),
     1                          HEAP(HPNPP1),HEAP(HPNPP2),
     $                          HEAP(HPNHP1),HEAP(HPNHP2),
     1                          HEAP(HPNWGH),HEAP(HPNVIO),
     $                          HEAP(HPNNUM),
     1                          HEAP(HPNNSP),HEAP(HPNCV)
     $                          ,HEAP(HPNEXC),HEAP(HPNPID))
                        END IF
                     END DO
C     
                  ELSE
                     CALL CHKEND('COUNtviolations>',DONE)
                  END IF
               END DO
               DONE=.FALSE.
C====================================================================
            ELSE IF (WD(1:4).EQ.'CALI') THEN
               HPNDIA=ALLHP(IREAL8(NATOM))
               CALL FILLR8(HEAP(HPNDIA),NATOM,1.0D0)
               HPNCAL=ALLHP(INTEG4(NATOM))
               CALL FILL4(HEAP(HPNCAL),NATOM,0)
               CALL CALIBR(ISLCT,HEAP(HPNCAL),HEAP(HPNDIA))
               CALL FREHP(HPNCAL,INTEG4(NATOM))
               CALL FREHP(HPNDIA,IREAL8(NATOM))
C====================================================================
            ELSE IF (WD(1:4).EQ.'ANAL') THEN
               CALL ANALRS(NOEHGL,ISLCT,JSLCT)
C====================================================================
c            ELSE IF (WD(1:4).EQ.'RESE') THEN
c               ARINIT=.TRUE.
c               CALL ARIINI
C====================================================================
            ELSE
               CALL CHKEND('ARIA>',DONE)
            END IF
         END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C     
C====================================================================
C     
      SUBROUTINE ARIFLI(I,J,ISLCT,JSLCT,NOEILS,NOEJLS,NOEHGL)
C     
C     flips the assignments of two atoms in the restraint list
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
C     -I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INTEGER NOEILS(*),NOEJLS(*)
      INTEGER NOEHGL(*),ISLCT(*),JSLCT(*),I,J
C     -local
      INTEGER II,JJ,NIAT,NJAT,IINOE
      LOGICAL QNEXT
C     
C     -begin
C     
      NIAT=0
      II=I
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
         NIAT=NIAT+1
         ISLCT(NIAT)=II
         CALL NEXTHY(II, NOEHGL,QNEXT)
      END DO
      NJAT=0
      JJ=J
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
         NJAT=NJAT+1
         JSLCT(NJAT)=JJ
         CALL NEXTHY(JJ, NOEHGL,QNEXT)
      END DO
C      WRITE(6,'(A,2I5)') ' ARIFLI: atoms in group: ', NIAT, NJAT
      IF (NIAT.EQ.NJAT) THEN
         DO IINOE=1,NOEMAX
            DO II=1,NIAT
               IF (NOEILS(IINOE).EQ.ISLCT(II)) THEN
                  NOEILS(IINOE)=JSLCT(II)
               ELSEIF (NOEILS(IINOE).EQ.JSLCT(II)) THEN
                  NOEILS(IINOE)=ISLCT(II)
               END IF
               IF (NOEJLS(IINOE).EQ.ISLCT(II)) THEN
                  NOEJLS(IINOE)=JSLCT(II)
               ELSEIF (NOEJLS(IINOE).EQ.JSLCT(II)) THEN
                  NOEJLS(IINOE)=ISLCT(II)
               END IF
            END DO
         END DO
      END IF      
C     
      RETURN
      END
C     
C====================================================================
C     
      SUBROUTINE ARIMOD(NOEORR,
     &     NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,NOERAV,NOERRV,
     &     NOESTP,NOEDIS,NOELOW,NOEHIG,NOEWGH,NOECV,NOEVIO,
     &     NOEPID,NOEPP1,NOEPP2,NOEHP1,NOEHP2,NOEVOL,
     &     NATOM,ISLCT,NISLCT,
     &     JSLCT,NJSLCT,X,Y,Z)
C     
C     Subroutine modifies any of the values of the specified restraint. 
C     e.g., MODIFY H200 1111 TEST=0 END
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
C     I/O
      INTEGER NOEORR(*)
      INTEGER NOEIPR(*), NOEILS(*), NOEJPR(*), NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOERAV(*), NOERRV(*), NOEDIS(*)
      DOUBLE PRECISION NOELOW(*), NOEHIG(*), NOEWGH(*)
      INTEGER NOECV(*), NOEVIO(*), NOESTP(*), NOEPID(*)
      DOUBLE PRECISION NOEPP1(*),NOEPP2(*),NOEHP1(*),NOEHP2(*)
     $     ,NOEVOL(*)
      INTEGER NATOM, ISLCT(*), NISLCT, JSLCT(*), NJSLCT
      DOUBLE PRECISION X(*), Y(*), Z(*)
C     local
      INTEGER I, II, JJ, J, TISLCT, TJSLCT, IINOE, PEAKID
      LOGICAL SUCCES,OK,FOUND
      CHARACTER*4 CLASS
C     
C     begin
C     
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('PEAKid=',PEAKID)
      IINOE=0
      FOUND=.FALSE.
      DO WHILE (.NOT.FOUND.AND.IINOE.LT.NOENUM)
         IINOE=IINOE+1
         FOUND = (NOEPID(IINOE).EQ.PEAKID)
     $        .AND.(CLASS.EQ.NOECNM(NOECND(IINOE)))
      END DO
      IF (.NOT.FOUND) THEN
         WRITE(6,'(A)')
     &        ' %ARIMOD-ERR: peak not found'
         IINOE=NOEMAX
      END IF

C     
         CALL PUSEND('ARIMOD>')
         DO WHILE (.NOT. DONE)
            CALL NEXTWD('ARIMOD>')
            CALL MISCOM('ARIMOD>',USED)
            IF (.NOT.USED) THEN
C     
               IF (WD(1:4).EQ.'HELP') THEN
                  CALL CNSHELP('cns-noe')
               ELSEIF (WD(1:4).EQ.'DIST') THEN
                  CALL NEXTF('Distance=',NOEDIS(IINOE))
               ELSEIF (WD(1:4).EQ.'LOW ') THEN
                  CALL NEXTF('lower-error=',NOELOW(IINOE))
               ELSEIF (WD(1:4).EQ.'HIGH') THEN
                  CALL NEXTF('higher-error=',NOEHIG(IINOE))
               ELSE IF (WD(1:4).EQ.'PPM1') THEN
                  CALL NEXTF('PPM(F1)=',NOEPP1(IINOE))
               ELSE IF (WD(1:4).EQ.'PPM2') THEN
                  CALL NEXTF('PPM(F2)=',NOEPP2(IINOE))
               ELSE IF (WD(1:4).EQ.'HPM1') THEN
                  CALL NEXTF('PPM(F1)=',NOEHP1(IINOE))
               ELSE IF (WD(1:4).EQ.'HPM2') THEN
                  CALL NEXTF('PPM(F2)=',NOEHP2(IINOE))
               ELSE IF (WD(1:4).EQ.'VOLU') THEN
                  CALL NEXTF('VOLUME=',NOEVOL(IINOE))
               ELSE IF (WD(1:4).EQ.'WEIG') THEN
                  CALL NEXTF('WEIGHT=',NOEWGH(IINOE))
               ELSE IF (WD(1:4).EQ.'CV  ') THEN
                  CALL NEXTI('CV=',NOECV(IINOE))
               ELSE
                  CALL CHKEND('ARIMOD>',DONE)
               END IF
            END IF
         END DO
         DONE=.FALSE.
      RETURN
      END
C     
C=======================================================================
      SUBROUTINE ATGR(NATOM,HGLNK,IATGR)
C     
C     sets up pointer array  iat --> igra
C     (<=> HGLNK(IAT) >= I; then enter all atoms in group)
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
      INTEGER NATOM, HGLNK(*), IATGR(*)
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C     
C     local
      INTEGER I,I1,NGRP
      LOGICAL QNEXT
C     
C     begin
C     
      DO I=1,NATOM
         IATGR(I)=0
      END DO
C     
      NGRP=0
      DO I=1,NATOM
         IF (FRSTEL(I,HGLNK)) THEN
            NGRP=NGRP+1
            I1=I
            QNEXT=.TRUE.
            DO WHILE (QNEXT)
               CALL NEXTHY(I1,HGLNK,QNEXT)
               IATGR(I1)=NGRP
            END DO
         END IF
      END DO
C     
      RETURN
      END
C====================================================================
C     
      SUBROUTINE NEXTHY(IAT, HGLNK, QNEXT)
C     
C     -returns next element in a group in IAT (cyclically).
C     -QNEXT true if there are elements left in group  
C     
      IMPLICIT NONE
C     
C     -global
C     
C     -input/ output
      INTEGER IAT, HGLNK(*)
      LOGICAL QNEXT
C     
C     -local
      INTEGER JAT
C     
C     -begin
C     
      JAT=HGLNK(IAT)
      QNEXT=(HGLNK(JAT).LT.JAT)
      IAT=JAT
C     
      END
C     
C======================================================================
C     
      SUBROUTINE EQVGRP(NATOM,HGLNK)
C     
C     group hydrogen atoms into groups of unresolved hydrogens.
C     
C     Author and copyright: Michael Nilges, EMBL
C     No warranty implied or expressed
C     All rights reserved
C     
      IMPLICIT NONE
C     global
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'coord.inc'
C     input/output
      INTEGER NATOM,HGLNK(*)
C     local
      INTEGER IISLCT,INHYBO,IPATBO,NSLCT,NHYBMX,IHGR
C     
C     begin
C     
      NHYBMX=4
C     
C     initialise group link array
C     
C     allocate space for the bond lists
      IISLCT=ALLHP(INTEG4(NATOM))
      INHYBO=ALLHP(INTEG4(NATOM))
      IPATBO=ALLHP(INTEG4(NHYBMX*NATOM))
      IHGR=ALLHP(INTEG4(NATOM))
C     
C     -parsing
      CALL PUSEND('EQUIvalent>')
      DO WHILE (.NOT. DONE)
         CALL NEXTWD('EQUIvalent>')
         CALL MISCOM('EQUIvalent>',USED)
         IF (.NOT.USED) THEN
C     
            IF (WD(1:4).EQ.'HELP') THEN
               CALL CNSHELP('cns-aria-equivalent')
C=====================================================================
            ELSE IF (WD(1:4).EQ.'METH') THEN
               CALL ALLGRP(3,.FALSE.,HGLNK,HEAP(INHYBO),4,
     @              HEAP(IPATBO),HEAP(IISLCT),HEAP(IHGR))
C=====================================================================
            ELSE IF (WD(1:4).EQ.'1-2 ') THEN
               CALL ALLGRP(2,.FALSE.,HGLNK,HEAP(INHYBO),4,
     @              HEAP(IPATBO),HEAP(IISLCT),HEAP(IHGR))
C=====================================================================
            ELSE IF ((WD(1:4).EQ.'NONE').OR.(WD(1:4).EQ.'INIT')) THEN
               CALL HGRINI(NATOM,HGLNK)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'SELE') THEN
               CALL SELCTA(HEAP(IISLCT),NSLCT,X,Y,Z,.TRUE.)
               CALL NEWGRP(HEAP(IISLCT),NSLCT,.TRUE.,HGLNK,HEAP(IHGR))
C=====================================================================
            ELSE
               CALL CHKEND('EQUIvalent>',DONE)
C     
            END IF
         END IF
      END DO
      DONE=.FALSE.
C     
      CALL FREHP(IHGR,INTEG4(NATOM))
      CALL FREHP(IPATBO,INTEG4(NHYBMX*NATOM))
      CALL FREHP(INHYBO,INTEG4(NATOM))
      CALL FREHP(IISLCT,INTEG4(NATOM))
C     
      RETURN
      END
C     
C=======================================================================
C     
      SUBROUTINE HGRINI(NATOM,HGLNK)
C     
C     -initialises the hydrogen group array.
C     -Checks if atoms are known and interact.
C     -All subsequent hydrogen tests are based on the hydrogen group array.
C     
C     -Author: Michael Nilges, HHMI and Yale University
C     
      IMPLICIT NONE
C     -global
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
C     -input/output
      INTEGER NATOM,HGLNK(*)
C     -local
      INTEGER I,NH
C     
C     -begin
C     
      NH=0
      DO I=1,NATOM
         HGLNK(I)=0
      END DO
      DO I=1,NATOM
         IF (HYDROG(I)) THEN
            HGLNK(I)=I
            NH=NH+1
         END IF
      END DO
C     
      RETURN
      END
C     
C=======================================================================
C     
      SUBROUTINE ALLGRP(HMIN,OVRWRT,HGLNK,NHYBON,NHYBMX,
     &     PATBON,ISLCT,LHGR)
C     
C     defines hydrogen groups by connectivity and number.
C     if OVRWRT any existing overlapping groups are overwritten.
C     
C     Author: Michael Nilges, HHMI and Yale University
C     
      IMPLICIT NONE
C     global
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C     input/output
      INTEGER HMIN,HGLNK(*),NHYBON(*)
      INTEGER NHYBMX,PATBON(NHYBMX,*),ISLCT(*),LHGR(*)
      LOGICAL OVRWRT
C     local
      INTEGER NHYDRO,I,J,IBT,JBT,I1,J1
C     function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C     
C     begin
C     
C     set up list of bonds to hydrogens
      DO I=1,NATOM
         NHYBON(I)=0
      END DO
C     
      DO I=1,NBOND
         IBT=IB(I)
         JBT=JB(I)
         IF (IBT.GT.0 .AND. JBT.GT.0) THEN
            IF (HGLNK(JBT).GT.0) THEN
               NHYBON(IBT)=NHYBON(IBT)+1
               PATBON(NHYBON(IBT),IBT)=JBT
            END IF
            IF (HGLNK(IBT).GT.0) THEN
               NHYBON(JBT)=NHYBON(JBT)+1
               PATBON(NHYBON(JBT),JBT)=IBT
            END IF
         END IF
      END DO
C     
C     insert groups
      DO I=1,NATOM
C     
         IF (HGLNK(I).EQ.0) THEN
            NHYDRO=NHYBON(I)
            DO I1=1,NATOM
               ISLCT(I1)=0
            END DO
C     
            IF (NHYDRO.GE.HMIN) THEN
C     found a group
               DO J=1,NHYDRO
                  J1=PATBON(J,I)
                  ISLCT(J1)=1
               END DO
               CALL NEWGRP(ISLCT,NHYDRO,OVRWRT,HGLNK,LHGR)
C     
            ELSEIF (NHYDRO.GT.0) THEN
C     one group for each h-atom
               DO J=1,NHYDRO
                  J1=PATBON(J,I)
                  ISLCT(J1)=1
                  CALL NEWGRP(ISLCT,1,OVRWRT,HGLNK,LHGR)
                  ISLCT(J1)=0
               END DO
C     
            END IF
         END IF
      END DO
C     
      RETURN
      END
C     
C=======================================================================
C     
      SUBROUTINE NEWGRP(ISLCT,NSLCT,OVRWRT,HGLNK,LGRP)
C     
C     groups are defined by array ISELCT.
C     non-hydrogen atoms are removed.
C     the groups are connected by cyclic pointers:
C     HGLNK(I1ST) = ILAST; => HGLNK(I1ST) > I1ST
C     HGLNK(I2ND) = I1ST;  => HGLNK(I2ND) < I2ND
C     HGLNK(I3RD) = I2ND;  ...
C     ...
C     HGLNK(I) = I  <=>  I is the only atom in group
C     HGLNK(I) = 0  <=>  I is not a hydrogen
C     if OVRWRT any existing overlapping groups are overwritten.
C     
C     Author and copyright: Michael Nilges, EMBL
C     No warranty implied or expressed
C     All rights reserved
C     
C     
      IMPLICIT NONE
C     

C     global variables
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C     input/output
      INTEGER ISLCT(*),NSLCT,HGLNK(*),LGRP(*)
      LOGICAL OVRWRT
C     local variables
      INTEGER NHYDRO,I,J,J1,JTEMP
      LOGICAL CLEAR
C     function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C     
C     begin
C     
C     remove non-hydrogen and unknown or non-interacting atoms in selection
      DO I=1,NATOM
         IF (HGLNK(I).GT.0) THEN
            ISLCT(I)=ISLCT(I)
         ELSE
            ISLCT(I)=0
         END IF
      END DO
C     gather selected group indices
      NHYDRO=0
      DO I=1,NATOM
         IF ((ISLCT(I)).GT.0) THEN
            NHYDRO=NHYDRO+1
            LGRP(NHYDRO)=I
         END IF
      END DO
C     
      CLEAR=.TRUE.
C     
C     check if any groups exist and delete if OVRWRT
      DO J=1,NHYDRO
         J1=HGLNK(LGRP(J))
         IF (HGLNK(J1).NE.J1) THEN
            IF (OVRWRT) THEN
C     remove group
               DO WHILE ((HGLNK(J1).NE.J1))
                  JTEMP=J1
                  HGLNK(J1)=J1
                  J1=HGLNK(JTEMP)
               END DO
            ELSE
               CLEAR=.FALSE.
            END IF
         END IF
      END DO
C     
      IF (CLEAR) THEN
C     
C     add a new group
         IF (NHYDRO.GT.0) THEN
            DO J=1,NHYDRO
               IF (J.GT.1) THEN
                  J1=J-1
               ELSE
                  J1=NHYDRO
               END IF
               HGLNK(LGRP(J))=LGRP(J1)
            END DO
         END IF
      END IF
C     
      RETURN
      END
C     
C====================================================================
C     
      LOGICAL FUNCTION FRSTEL(IAT, HGLNK)
C     
C     true if IAT is first element on input. 
C     
      IMPLICIT NONE
C     
C     input/ output
      INTEGER IAT, HGLNK(*)
C     
C     begin
C     
      FRSTEL=((HGLNK(IAT).GE.IAT).AND.(HGLNK(IAT).GT.0))
C     
      END
C     
C====================================================================
      SUBROUTINE ARIINI
C     
C     subroutine initializes the ARIA list
C     
C     Author and copyright: Michael Nilges, EMBL
C     No warranty implied or expressed
C     All rights reserved
C     
      IMPLICIT NONE
C     
C     I/O
      INCLUDE 'cns.inc'
      INCLUDE 'aria.inc'
C     begin
C     
C     set heap pointers of dynamic data structure to zero
      HPNHGL=0
      HPNMAT=0
      HPNMA2=0
      HPDMAT=0
      HPIMAT=0
      HPNCAL=0
      HPNPSE=0
      HPNDIA=0
      HPHHPR=0
      HPNPR =0
      HPRATE=0
C     
C     initialize array sizes
      MATDIM=0
      DIMMAT=0
      NHG=0
      DIMNHG=0
      NMRATM=0
      DIMATM=0
C     
C     reset other variables
      CALL ARIRES
      RETURN
      END

C====================================================================
      SUBROUTINE ARIRES
C     
C     subroutine resets all ARIA variables
C     
C     Author and copyright: Michael Nilges, EMBL
C     No warranty implied or expressed
C     All rights reserved
C     
C     
      IMPLICIT NONE
C     
C     I/O
      INCLUDE 'cns.inc'
      INCLUDE 'aria.inc'
      INCLUDE 'consta.inc'
C begin
      CALL ARIHP(0,0,0)
C
      ARINIT = .TRUE.
      VIOSTR = 0
      VIOEXC = 0.0D0
      VIOTHR = 0.0D0
      AMBLEV = 1.0D0
      AMBMAX = 1
      AMBMIN = 1
      AMBSTR = 0
      AMBCUT = R4BIG
      AVEMOD = 'AVER'
      AMBMOD = 'ALL '
      AMBFRM = 'OR-R'
      CALCNT=0
      CALDIS=6.0D0
C
      RETURN
      END
C     
C====================================================================
      SUBROUTINE ARIHP(NMATDIM,NNHG,NNATOM)
C     
C     **********************************************
C     * Author and copyright: Michael Nilges, EMBL *
C     * No warranty implied or expressed           *
C     * All rights reserved                        *
C     **********************************************
C     
      IMPLICIT NONE
C     
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'aria.inc'
      INTEGER NMATDIM,NNHG,NNATOM
C     
C     free old space if things have changed
C     assign new space or initialize variables
      IF (NMATDIM.NE.DIMMAT) THEN
         IF (HPNMA2.NE.0)   CALL FREHP(HPNMA2,INTEG4(DIMMAT))
         IF (HPNMAT.NE.0)   CALL FREHP(HPNMAT,IREAL8(DIMMAT))
         IF (NMATDIM .GT. 0) THEN
            HPNMAT=ALLHP(IREAL8(NMATDIM))
            HPNMA2=ALLHP(INTEG4(NMATDIM))
         ELSE
            HPNMA2=0
            HPNMAT=0
         END IF
      END IF
      IF (NNHG.NE.DIMNHG) THEN
         IF (HPIMAT.NE.0)   CALL FREHP(HPIMAT,IREAL8(DIMNHG*DIMNHG))
         IF (HPDMAT.NE.0)   CALL FREHP(HPDMAT,IREAL8(DIMNHG*DIMNHG))
         IF (HPHHPR.NE.0)   CALL FREHP(HPHHPR,INTEG4(DIMNHG*DIMNHG))
         IF (HPRATE.NE.0)   CALL FREHP(HPRATE,IREAL8(DIMNHG*DIMNHG))
         IF (NNHG .GT. 0) THEN
            HPRATE = ALLHP(IREAL8(NNHG*NNHG))
            HPHHPR = ALLHP(INTEG4(NNHG*NNHG))
            HPDMAT = ALLHP(IREAL8(NNHG*NNHG))
            HPIMAT = ALLHP(IREAL8(NNHG*NNHG))
         ELSE
            HPDMAT=0
            HPHHPR=0
            HPIMAT=0
            HPRATE=0
         END IF
      END IF
      IF (NNATOM.NE.DIMATM) THEN
         IF (HPPOPU.NE.0)   CALL FREHP(HPPOPU,IREAL8(DIMATM))
         IF (HPNPR.NE.0)    CALL FREHP(HPNPR,INTEG4(DIMATM))
         IF (HPIGR.NE.0)    CALL FREHP(HPIGR,INTEG4(DIMATM))
         IF (HPNCAL.NE.0)   CALL FREHP(HPNCAL,INTEG4(DIMATM))
         IF (HPNDIA.NE.0)   CALL FREHP(HPNDIA,IREAL8(DIMATM))
         IF (HPNPSE.NE.0)   CALL FREHP(HPNPSE,INTEG4(DIMATM))
         IF (HPNHGL.NE.0)   CALL FREHP(HPNHGL,INTEG4(DIMATM))
         IF (NNATOM .GT. 0) THEN
            HPNHGL = ALLHP(INTEG4(NNATOM))
            HPNPSE = ALLHP(INTEG4(NNATOM))
            HPNDIA = ALLHP(IREAL8(NNATOM))
            HPNCAL = ALLHP(INTEG4(NNATOM))
            HPIGR = ALLHP(INTEG4(NNATOM))
            HPNPR = ALLHP(INTEG4(NNATOM))
            HPPOPU = ALLHP(IREAL8(NNATOM))
         ELSE
            HPNHGL=0
            HPNPSE=0
            HPNDIA=0
            HPNCAL=0
            HPPOPU=0
            HPNPR=0
            HPIGR=0
         END IF
      END IF
C        
C     store the dimensions values in common block
      DIMATM=NNATOM
      DIMNHG=NNHG
      DIMMAT=NMATDIM
C     
      RETURN 
      END

