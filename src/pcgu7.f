      MODULE PCGUMODULE
C     ******************************************************************
C     PCGU LINEAR SOLVER CALLED FROM SMS
C     ******************************************************************
        USE PCGUSOLVERMODULE, ONLY: PCGUSOLVER
        PRIVATE
        PUBLIC :: PCGU7U1AR
        PUBLIC :: PCGU7U1AP
        PUBLIC :: PCGU7U1DA
        PUBLIC :: SET_PCGUINPUT
        TYPE(PCGUSOLVER), SAVE :: SOLVER

      CONTAINS

      SUBROUTINE PCGU7U1AR(IN, NJA, NEQS, MXITER, HICLOSE, ITER1,
     1   IPRSMS,IFDPARAM,IPCGUM)
C     ******************************************************************
C     READ PCGU DATA AND INITIALIZE THE SOLVER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT, IA, JA
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER, INTENT(IN) :: IN
      INTEGER, INTENT(IN) :: NJA
      INTEGER, INTENT(IN) :: NEQS
      INTEGER, INTENT(IN) :: MXITER
      DOUBLE PRECISION, INTENT(IN) :: HICLOSE
      INTEGER, INTENT(IN) :: ITER1
      INTEGER, INTENT(IN) :: IPRSMS
      INTEGER, INTENT(IN) :: IFDPARAM
      INTEGER, INTENT(INOUT) :: IPCGUM
C     + + + LOCAL DEFINITIONS + + +
      CHARACTER*200 LINE
      CHARACTER (LEN= 10) :: clin(1:2)
      CHARACTER (LEN= 20) :: cipc(0:3)
      CHARACTER (LEN= 20) :: cscale(0:2)
      CHARACTER (LEN= 25) :: corder(0:2)
      CHARACTER (LEN= 16) :: ccnvgopt(0:1)
      CHARACTER (LEN=  4) :: cval
      INTEGER :: LLOC, ISTART, ISTOP, I
      INTEGER :: ILINMETH, IPC, ISCL, IORD, IAORD, ICNVGOPT
      REAL :: R, HCLOSE, RCLOSE, RELAX
      REAL, PARAMETER :: RZERO = 0.0
C     + + + DATA + + +
      DATA clin  /'CG        ',
     2            'BCGS      '/
      DATA cipc  /'NONE                ',
     2            'JACOBI              ',
     3            'INCOMPLETE LU       ',
     4            'MOD. INCOMPLETE LU  '/
      DATA cscale/'NO SCALING          ',
     2            'SYMMETRIC SCALING   ',
     3            'L2 NORM SCALING     '/
      DATA corder/'ORIGINAL ORDERING        ',
     2            'RCM ORDERING             ',
     3            'MINIMUM DEGREE ORDERING  '/
      DATA ccnvgopt/'INFINITY NORM   ',
     2              'L2 NORM         '/
C     + + + FORMATS + + +
02000 FORMAT (1X,/1X,'PCGU -- UNSTRUCTURED CONJUGATE-GRADIENT SOLUTION',
     &        ' PACKAGE, VERSION 7.02, 08/13/2013')
02010 FORMAT (1X,/,14X,'SOLUTION BY THE CONJUGATE-GRADIENT METHOD',
     &        /,1X,66('-'),/,
     &        ' MAXIMUM OF ',I6,' CALLS OF SOLUTION ROUTINE',/,
     &        ' MAXIMUM OF ',I6,
     &        ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE',/,
     &        ' LINEAR ACCELERATION METHOD            =',1X,A,/,
     &        ' MATRIX PRECONDITIONING TYPE           =',1X,A,/,
     &        ' MATRIX SCALING APPROACH               =',1X,A,/,
     &        ' MATRIX REORDERING APPROACH            =',1X,A,/,
     &        ' SOLVE ONLY ACTIVE NODES               =',1X,L1,/,
     &        ' HEAD CHANGE CRITERION FOR CLOSURE     =',E15.5,/,
     &        ' RESIDUAL CHANGE CRITERION FOR CLOSURE =',E15.5,/,
     &        ' RESIDUAL CONVERGENCE NORM             =',1X,A,/,
     &        ' RELAXATION FACTOR                     =',E15.5,/,
     &        '  ONLY USED WITH MILU0 PRECONDITIONER',//)
02020 FORMAT (///,1X,'PCGU DATA INPUT ERROR:',
     &          /,2X,'SCALING MUST BE USED (ISCL.GT.0) IF USING',
     &          /,2X,'THE ILU0 OR MILU0 PRECONDITIONERS (IPC.EQ.2 OR',
     &          /,2X,'IPC.EQ.3) WITH MATRIX REORDERING (IORD.GT.0)')
C     ------------------------------------------------------------------
      ILINMETH = 0
      IPC = 0
      ISCL = 0
      IORD = 0
      IAORD = 0
      ICNVGOPT = 0
      RCLOSE = RZERO
      RELAX = RZERO
      HCLOSE = REAL(HICLOSE)
C
C-------PRINT A MESSAGE IDENTIFYING UPCG PACKAGE
      WRITE (IOUT,2000)
C
C-------READ AND PRINT COMMENTS
      CALL URDCOM(IN,IOUT,LINE)
      IF ( IFDPARAM.EQ.0 ) THEN
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IPC,R,IOUT,IN)
        cval = LINE(ISTART:ISTOP)
        SELECT CASE (cval)
          CASE ( 'CG' )
            ILINMETH = 1
          CASE ( 'BCGS' )
            ILINMETH = 2
          CASE DEFAULT
            ILINMETH = 1
            READ (CVAL,*) IPC
        END SELECT
        IF ( cval.EQ.'CG  ' .OR. cval.EQ.'BCGS' ) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPC,R,IOUT,IN)
        END IF
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISCL,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IORD,R,IOUT,IN)
C         IORD PLUS 10 SOLVES ONLY THE ACTIVE NODES
        IF ( IORD.GE.10 ) THEN
          IORD  = IORD - 10
          IAORD = 1
        END IF
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RCLOSE,IOUT,IN)
C         A NEGATIVE RCLOSE TESTS THE L2 NORM OF THE RESIDUAL
        IF ( RCLOSE.LT.RZERO ) THEN
          RCLOSE = ABS( RCLOSE )
          ICNVGOPT = 1
        END IF
        IF ( IPC.EQ.3 ) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RELAX,-IOUT,IN)
          IF ( ISTART.EQ.200 ) THEN
            IF ( ISTOP.EQ.200 ) RELAX = 0.97
          END IF
          IF ( RELAX.EQ.RZERO ) THEN
            IF ( LINE(200:200).EQ.'E' ) RELAX = 0.97
          END IF
        ELSE
          RELAX = 0.0
        END IF
      ELSE
        CALL SET_PCGUINPUT(IFDPARAM,ILINMETH,IPC,ISCL,IORD,
     2                     RCLOSE,RELAX)
      END IF
      IPCGUM = ILINMETH
C
C-------ERROR CHECKING FOR OPTIONS
      IF ( IPC.LT.0  ) IPC  = 0
      IF ( IPC.GT.3  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: IPC  MUST BE .LE. 3'
        CALL USTOP('PCGU7AR: IPC  MUST BE .LE. 3')
      END IF
      IF ( ISCL.LT.0 ) ISCL = 0
      IF ( ISCL.GT.2  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: ISCL MUST BE .LE. 2'
        CALL USTOP('PCGU7AR: ISCL MUST BE .LE. 2')
      END IF
      IF ( IORD.LT.0 ) IORD = 0
      IF ( IORD.GT.2  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: IORD MUST BE .LE. 2'
        CALL USTOP('PCGU7AR: IORD MUST BE .LE. 2')
      END IF
      IF ( RELAX.LT.0.0 ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: RELAXPCGU MUST BE .GE. 0.0'
        CALL USTOP('PCGU7AR: RELAXPCGU MUST BE .GE. 0.0')
      END IF
      IF ( RELAX.GT.1.0 ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: RELAXPCGU MUST BE .LE. 1.0'
        CALL USTOP('PCGU7AR: RELAXPCGU MUST BE .LE. 1.0')
      END IF
C
C-------PRINT THE SOLVER OPTIONS
      WRITE (IOUT,2010) MXITER, ITER1,
     2                  clin(ILINMETH), cipc(IPC),
     3                  cscale(ISCL), corder(IORD), (IAORD.GT.0),
     4                  HCLOSE, RCLOSE, ccnvgopt(ICNVGOPT),
     5                  RELAX
C
C-------ENSURE THAT SCALING IS USED WITH THE ILU0 AND MILU0
C       PRECONDITIONERS IF RCM OR MINIMUM DEGREE ORDERING IS USED
      IF ( IPC.EQ.2 .OR. IPC.EQ.3 ) THEN
        IF ( IORD.NE.0 ) THEN
          IF ( ISCL.EQ.0 ) THEN
            WRITE ( IOUT,2020 )
            CALL USTOP('SCALING MUST BE USED FOR ILU0 AND MILU0 '//
     2                 'WITH REORDERING')
          END IF
        END IF
      END IF
C
C-------INITIALIZE THE SOLVER
      CALL SOLVER%INIT(IOUT, IPRSMS, ITER1, ILINMETH, IPC,
     2                 ISCL, IORD, IAORD, ICNVGOPT,
     3                 HCLOSE, RCLOSE, RELAX,
     4                 NEQS, NJA, IA, JA)
C
C-------RETURN
      RETURN
      END SUBROUTINE PCGU7U1AR
C
      SUBROUTINE PCGU7U1AP(AC,RHS,HNEW,IBOUND,ICNVG,KSTP,KITER,IN_ITER)
C     ******************************************************************
C     SOLVE THE SYSTEM OF EQUATIONS FOR HNEW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      DOUBLEPRECISION, DIMENSION(*), INTENT(IN)    :: AC
      DOUBLEPRECISION, DIMENSION(*), INTENT(IN)    :: RHS
      DOUBLEPRECISION, DIMENSION(*), INTENT(INOUT) :: HNEW
      INTEGER, DIMENSION(*), INTENT(IN)            :: IBOUND
      INTEGER, INTENT(INOUT)                       :: ICNVG
      INTEGER, INTENT(IN)                          :: KSTP
      INTEGER, INTENT(IN)                          :: KITER
      INTEGER, INTENT(INOUT)                       :: IN_ITER
C     ------------------------------------------------------------------
      CALL SOLVER%SOLVE(KSTP, KITER, IBOUND, AC, HNEW, RHS, ICNVG,
     2                  IN_ITER)
C
C-------RETURN
      RETURN
      END SUBROUTINE PCGU7U1AP
C
      SUBROUTINE PCGU7U1DA()
C     ******************************************************************
C     DEALLOCATE PCGU DATA
C     ******************************************************************
      CALL SOLVER%DESTROY()
C
C-------RETURN
      RETURN
      END SUBROUTINE PCGU7U1DA
C
      SUBROUTINE SET_PCGUINPUT(IFDPARAM,ILINMETH,IPC,ISCL,IORD,
     2                         RCLOSE,RELAX)
C     ******************************************************************
C     SET THE PCGU OPTIONS FOR THE SIMPLE, MODERATE, AND COMPLEX OPTIONS
C     ******************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: IFDPARAM
      INTEGER, INTENT(INOUT) :: ILINMETH
      INTEGER, INTENT(INOUT) :: IPC
      INTEGER, INTENT(INOUT) :: ISCL
      INTEGER, INTENT(INOUT) :: IORD
      REAL, INTENT(INOUT)    :: RCLOSE
      REAL, INTENT(INOUT)    :: RELAX
C Simple option
      SELECT CASE ( IFDPARAM )
      CASE(1)
        ILINMETH=1
        IPC = 2
        ISCL = 0
        IORD = 0
        RCLOSE = 1.0e-1
        RELAX = 0.0
C Moderate
      CASE(2)
        ILINMETH=2
        IPC = 3
        ISCL = 0
        IORD = 0
        RCLOSE = 1.0e-1
        RELAX = 0.97
C Complex
      CASE(3)
        ILINMETH=2
        IPC = 3
        ISCL = 0
        IORD = 0
        RCLOSE = 1.0e-1
        RELAX = 0.97
      END SELECT
      RETURN
      END SUBROUTINE SET_PCGUINPUT

      END MODULE PCGUMODULE
