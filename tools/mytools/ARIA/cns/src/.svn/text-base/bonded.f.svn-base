      FUNCTION BONDED (IATOM)
C returns name of the last atom bonded to iatom
C Michael Nilges, EMBL
C
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER IATOM
      INTEGER BONDED
C local
      INTEGER I,BATOM,IBT,JBT
C
C begin
C
      BATOM = 0
C
      DO I=1,NBOND
        IBT=IB(I)
        JBT=JB(I)
        IF (IBT.GT.0.AND.JBT.GT.0) THEN
          IF (IBT.EQ.IATOM) THEN
            BATOM = JBT
          ELSEIF (JBT.EQ.IATOM) THEN
            BATOM = IBT
          END IF
        END IF
      END DO
      BONDED=BATOM
      END
