REMARKS  TOPH19.pep -MACRO for protein sequence
SET ECHO=FALSE END

! this is a macro to define standard protein peptide bonds
! and termini to generate a protein sequence.
! it should be added as @TOPH11.PEP in the SEGMent SEQUence
! level.
! 
! Author: Axel Brunger, 19-JAN-84



LINK PEPP    HEAD - *     TAIL + PRO     END  { LINK to PRO }
LINK PEPT    HEAD - *     TAIL + *       END

FIRSt PROP                TAIL + PRO     END { nter for PRO }
FIRSt NTER                TAIL + *       END

LAST  CTER   HEAD - *                    END


SET ECHO=TRUE END
