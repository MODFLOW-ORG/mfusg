      MODULE PCGULIBMODULE
C     ******************************************************************
C     SCALING, PRECONDITIONER, REORDERING, AND SPARSE MATRIX ROUTINES
C     USED BY THE PCGU SOLVER
C     ******************************************************************
      CONTAINS
C
C
C-------ROUTINE TO SCALE THE COEFFICIENT MATRIX (AC), 
C       THE RHS (B), AND THE ESTIMATE OF HNEW (X)
      SUBROUTINE SPCGU_SCALE(IOPT,ISCL,NIAC,NNZC,IA,JA,AC,X,B,
     2                       DSCALE,DSCALE2)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IOPT
        INTEGER, INTENT(IN) :: ISCL
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IA
        INTEGER, DIMENSION(NNZC),   INTENT(IN) :: JA
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(INOUT) :: AC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: X
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: B
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: DSCALE
        DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: DSCALE2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, n
        INTEGER :: id, jc
        INTEGER :: i0, i1
        DOUBLEPRECISION :: v, c1, c2
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C---------SCALE SCALE AC, X, AND B
        IF ( IOPT.EQ.0 ) THEN
C-----------SYMMETRIC SCALING
          SELECT CASE ( ISCL )
            CASE ( 1 )
              DO n = 1, NIAC
                id   = IA(n)
                v    = AC(id)
                c1   = 1.0D0 / SQRT( ABS( v ) )
                DSCALE(n)  = c1
                DSCALE2(n) = c1
              END DO
C               SCALE AC -- AC = DSCALE(row) * AC(i) * DSCALE2(col)
              DO n = 1, NIAC
                c1 = DSCALE(n)
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = DSCALE2(jc)
                  AC(i) = c1 * AC(i) * c2 
                END DO
              END DO
C-----------L-2 NORM SCALING
            CASE ( 2 )
C               SCALE EACH ROW SO THAT THE L-2 NORM IS 1
              DO n = 1, NIAC
                c1 = 0.0D0
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  c1 = c1 + AC(i) * AC(i)
                END DO
                c1 = SQRT( c1 )
                IF ( c1.EQ.0.0D0 ) THEN
                  c1 = 1.0D0
                ELSE
                  c1 = 1.0D0 / c1
                END IF
                DSCALE(n) = c1 
C                 INITIAL SCALING OF AC -- AC = DSCALE(row) * AC(i)              
                DO i = i0, i1
                  AC(i) = c1 * AC(i)
                END DO
              END DO
C               SCALE EACH COLUMN SO THAT THE L-2 NORM IS 1
              DO n = 1, NIAC
                DSCALE2(n) = 0.0D0
              END DO
              c2 = 0.0D0
              DO n = 1, NIAC
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = AC(i)
                  DSCALE2(jc) = DSCALE2(jc) + c2 * c2
                END DO
              END DO
              DO n = 1, NIAC
                c2 = DSCALE2(n)
                IF ( c2.EQ.0.0D0 ) THEN
                  c2 = 1.0D0
                ELSE
                  c2 = 1.0D0 / SQRT( c2 )
                END IF
                DSCALE2(n) = c2
              END DO
C               FINAL SCALING OF AC -- AC = DSCALE2(col) * AC(i)              
              DO n = 1, NIAC
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = DSCALE2(jc)
                  AC(i) = c2 * AC(i)
                END DO
              END DO
          END SELECT
C-----------SCALE X AND B
          DO n = 1, NIAC
            c1    = DSCALE(n)
            c2    = DSCALE2(n)
            X(n)  = X(n) / c2
            B(n)  = B(n) * c1
          END DO
C---------UNSCALE SCALE AC, X, AND B
        ELSE
          DO n = 1, NIAC
            c1 = DSCALE(n)
            i0 = IA(n)
            i1 = IA(n+1) - 1
C             UNSCALE AC
            DO i = i0, i1
              jc = JA(i)
              c2 = DSCALE2(jc)
              AC(i) = ( 1.0D0 / c1 ) * AC(i) * ( 1.0D0 / c2 ) 
            END DO
C             UNSCALE X AND B
            c2   = DSCALE2(n)
            X(n) = X(n) * c2
            B(n) = B(n) / c1
          END DO     
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_SCALE
C
C-------ROUTINE TO UPDATE THE PRECONDITIONER
      SUBROUTINE SPCGU_PCU(IOUT,NNZC,NIAC,NIAPC,NJAPC,NNZAPC,IPC,RELAX,
     2                     AC,IAC,JAC,APC,IAPC,JAPC,IWC,WC)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IOUT
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NIAPC
        INTEGER, INTENT(IN) :: NJAPC
        INTEGER, INTENT(IN) :: NNZAPC
        INTEGER, INTENT(IN) :: IPC
        REAL, INTENT(IN) :: RELAX
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)     :: AC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)    :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)      :: JAC
        DOUBLEPRECISION, DIMENSION(NNZAPC), INTENT(INOUT) :: APC
        INTEGER, DIMENSION(NIAPC+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NJAPC), INTENT(INOUT)   :: JAPC
        INTEGER, DIMENSION(NIAPC), INTENT(INOUT)   :: IWC
        DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT) :: WC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: izero
        DOUBLEPRECISION :: delta
C     + + + FUNCTIONS + + +
C     + + + FORMATS + + +
2000    FORMAT (/,' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.',
     &          /,' REPLACING ZERO AND SIGN-CHANGING PIVOTS WITH A',
     &          ' SMALL VALUE (SPCGU_PCU)')
C     + + + CODE + + +
        izero = 0
        delta = 0.0D0
        SELECT CASE(IPC)
C           NO PRE-CONDITIONER
          CASE (0)
C           JACOBI PRE-CONDITIONER
          CASE (1)
            CALL SPCGU_PCJ(NNZC,NIAC,AC,APC,IAC,JAC)
C           ILU0 AND MILU0
          CASE (2,3)
            LUPC: DO
              CALL SPCGU_PCILU0(NNZC,NIAC,AC,IAC,JAC,
     2                          APC,IAPC,JAPC,IWC,WC,
     3                          RELAX,izero,delta)
              IF ( izero.NE.1 ) THEN
                EXIT LUPC
              END IF
              delta = 1.5D0 * delta + 0.001
C               USE THE LARGEST SHIFT AND SMALL PIVOTS ON THE LAST PASS
              IF ( delta.GT.0.5D0 ) THEN
                WRITE(IOUT,2000)
                delta = 0.5D0
                izero = 2
              END IF
            END DO LUPC
C           ADDITIONAL PRECONDITIONERS - ILUT, etc.
        END SELECT
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCU
C
C-------JACOBI PRECONDITIONER - INVERSE OF DIAGONAL 
      SUBROUTINE SPCGU_PCJ(NNZC,NIAC,AC,APC,IAC,JAC)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)      :: AC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)   :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAC
        INTEGER, DIMENSION(NNZC),   INTENT(IN) :: JAC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, n
        INTEGER :: ic0, ic1
        INTEGER :: id
        DOUBLEPRECISION :: t
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NIAC
            ic0 = IAC(n)
            ic1 = IAC(n+1) - 1
            id = IAC(n)
            DO i = ic0, ic1
              IF ( JAC(i).EQ.n ) THEN
                id = i
                EXIT
              END IF
            END DO
            t  = AC(id)
            IF ( ABS( t ).GT.DZERO ) t = DONE / t
            APC(n) = t
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCJ

      SUBROUTINE SPCGU_JACA(NIAC,A,D1,D2)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: t
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NIAC
          t     = A(n) * D1(n)
          D2(n) = t
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_JACA

      SUBROUTINE SPCGU_PCILU0(NNZC,NIAC,AC,IAC,JAC,
     2                        APC,IAPC,JAPC,IWC,WC,
     3                        RELAX,IZERO,DELTA)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)     :: AC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)    :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)      :: JAC
        DOUBLEPRECISION, DIMENSION(NNZC), INTENT(INOUT)   :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NNZC), INTENT(INOUT)   :: JAPC
        INTEGER, DIMENSION(NIAC), INTENT(INOUT)   :: IWC
        DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)   :: WC
        REAL, INTENT(IN) :: RELAX
        INTEGER, INTENT(INOUT) :: IZERO
        DOUBLEPRECISION, INTENT(IN) :: DELTA
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1, id0
        INTEGER :: iic0, iic1, iid0
        INTEGER :: iu, iiu
        INTEGER :: j, n
        INTEGER :: jj
        INTEGER :: jcol, jw
        INTEGER :: jjcol
        DOUBLEPRECISION :: drelax
        DOUBLEPRECISION :: c1, c2
        DOUBLEPRECISION :: sd1
        DOUBLEPRECISION :: tl
        DOUBLEPRECISION :: rs
        DOUBLEPRECISION :: d
        LOGICAL :: lsmall
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
        DOUBLEPRECISION, PARAMETER :: DTINY = 1.0D-06
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        drelax = DBLE( RELAX )
C         IZERO OF 2 REPLACES BAD PIVOTS INSTEAD OF RETURNING
        lsmall = IZERO.EQ.2
        IZERO    = 0
        DO n = 1, NIAC
          IWC(n)  = 0
          WC(n)   = DZERO
        END DO
        MAIN: DO n = 1, NIAC
          ic0 = IAC(n)
          ic1 = IAC(n+1) - 1
          DO j = ic0, ic1
            jcol      = JAC(j)
            IWC(jcol) = 1
            WC(jcol) = WC(jcol) + AC(j)
          END DO
          ic0 = IAPC(n)
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n)
          rs   = DZERO
          LOWER: DO j = ic0, iu-1
            jcol     = JAPC(j)
            iic0     = IAPC(jcol) 
            iic1     = IAPC(jcol+1) - 1
            iiu      = JAPC(jcol)
            tl       = WC(jcol) * APC(jcol)
            WC(jcol) = tl
            DO jj = iiu, iic1
              jjcol = JAPC(jj)
              jw    = IWC(jjcol)
              IF ( jw.NE.0 ) THEN
                WC(jjcol) = WC(jjcol) - tl * APC(jj)
              ELSE
                rs = rs + tl * APC(jj)
              END IF
            END DO
          END DO LOWER
C           DIAGONAL - CALCULATE INVERSE OF DIAGONAL FOR SOLUTION
          d   = WC(n)
          tl  = ( DONE + DELTA ) * d - ( drelax * rs )
C-----------ENSURE THAT THE SIGN OF THE DIAGONAL HAS NOT CHANGED AND IS NOT ZERO
          sd1 = SIGN(d,tl)
          IF ( sd1.NE.d .OR. ABS(tl).EQ.DZERO ) THEN
            IF ( lsmall ) THEN
              tl = SIGN(DTINY,d)
            ELSE
              IZERO = 1
              EXIT MAIN
            END IF
          END IF
          APC(n) = DONE / tl
C           RESET POINTER FOR IW TO ZERO
          IWC(n) = 0
          WC(n)  = DZERO
          DO j = ic0, ic1
            jcol = JAPC(j)
            APC(j) = WC(jcol)
            IWC(jcol) = 0
            WC(jcol) = DZERO
          END DO
        END DO MAIN
C---------SCALE THE UPPER ENTRIES BY THE INVERSE OF THE DIAGONAL FOR
C         SPCGU_ILU0A, AFTER THE FACTORIZATION HAS USED THEM
        IF ( IZERO.EQ.0 ) THEN
          DO n = 1, NIAC
            iu  = JAPC(n)
            ic1 = IAPC(n+1) - 1
            DO j = iu, ic1
              APC(j) = APC(j) * APC(n)
            END DO
          END DO
        END IF
!        do n = 1, niac
!          ic0 = IAPC(n)
!          ic1 = IAPC(n+1) - 1
!          iu  = JAPC(n)
!          do j = ic0,iu-1
!            jcol = JAPC(j)
!            write (10991,'(i10,i10,g15.7)' ) n, jcol, 
!     2                                       APC(j)
!          end do
!          write (10991,'(i10,i10,g15.7)' ) n, n, 
!     2                                     APC(n)
!          do j = iu,ic1
!            jcol = JAPC(j)
!            write (10991,'(i10,i10,g15.7)' ) n, jcol, 
!     2                                       APC(j)
!          end do
!        end do
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCILU0

      SUBROUTINE SPCGU_ILU0A(NNZC,NIAC,APC,IAPC,JAPC,R,D)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(INOUT)  :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAPC
        INTEGER, DIMENSION(NNZC), INTENT(IN)   :: JAPC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)     :: R
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: D
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1, id0
        INTEGER :: iu
        INTEGER :: jcol
        INTEGER :: j, n
        DOUBLEPRECISION :: t
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C         FORWARD SOLVE - APC * D = R
        FORWARD: DO n = 1, NIAC
          t   = R(n)
          ic0 = IAPC(n)
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n) - 1
          LOWER: DO j = ic0, iu
            jcol = JAPC(j)
            t    = t - APC(j) * D(jcol)
          END DO LOWER
          D(n) = t
        END DO FORWARD
C         BACKWARD SOLVE - D = D / U, WITH THE UPPER ENTRIES SCALED BY
C         THE INVERSE OF THE DIAGONAL. THE NEAREST COLUMN, SOLVED JUST
C         BEFORE THIS ROW, IS SUBTRACTED LAST TO SHORTEN THE DEPENDENCY
        BACKWARD: DO n = NIAC, 1, -1
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n)
          t   = D(n) * APC(n)
          UPPER: DO j = ic1, iu, -1
            jcol = JAPC(j)
            t    = t - APC(j) * D(jcol)
          END DO UPPER
          D(n) = t
        END DO BACKWARD
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_ILU0A
C
C---------TEST FOR SOLVER CONVERGENCE
C         ICNVGOPT OF 1 ALSO ACCEPTS A RESIDUAL REDUCED BY EPFACT
        SUBROUTINE SPCGU_TESTCNVG( ICNVGOPT,ICNVG,
     2                             HMAX,RMAX,
     3                             RMAX0,EPFACT,HCLOSE,RCLOSE )
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)         :: ICNVGOPT
        INTEGER, INTENT(INOUT)      :: ICNVG
        DOUBLEPRECISION, INTENT(IN) :: HMAX
        DOUBLEPRECISION, INTENT(IN) :: RMAX
        DOUBLEPRECISION, INTENT(IN) :: RMAX0
        DOUBLEPRECISION, INTENT(IN) :: EPFACT
        DOUBLEPRECISION, INTENT(IN) :: HCLOSE
        DOUBLEPRECISION, INTENT(IN) :: RCLOSE
C     + + + CODE + + +
        IF ( ICNVGOPT.EQ.0 ) THEN
          IF ( ABS(HMAX).LE.HCLOSE .AND. ABS(RMAX).LE.RCLOSE ) THEN
            ICNVG = 1
          END IF
        ELSE
          IF ( RMAX.LE.RCLOSE ) THEN
            ICNVG = 1
          ELSE IF ( RMAX.LE.RMAX0*EPFACT ) THEN
            ICNVG = 1
          END IF
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_TESTCNVG
C
C---------GENERATE IAPC AND JAPC FROM IAC AND JAC
C         JAPC(1:NIAC) HAS THE POSITION OF THE UPPER ENTRY FOR A ROW
C         JAPC(NIAC+1:NNZC) IS THE COLUMN POSITION FOR ENTRY
C         APC(1:NIAC) PRECONDITIONED INVERSE OF THE DIAGONAL
C         APC(NIAC+1:NNZC) PRECONDITIONED ENTRIES FOR OFF DIAGONALS
        SUBROUTINE SPCGU_PCCRS(NIAC,NNZC,IAC,JAC,
     2                         IAPC,JAPC,IWC)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)         :: NIAC
        INTEGER, INTENT(IN)         :: NNZC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)    :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)      :: JAC
        INTEGER, DIMENSION(NIAC+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NNZC), INTENT(INOUT)   :: JAPC
        INTEGER, DIMENSION(NIAC+1), INTENT(INOUT) :: IWC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n, j
        INTEGER :: i0, i1
        INTEGER :: nlen
        INTEGER :: ic,ip
        INTEGER :: jcol
        INTEGER, DIMENSION(:), ALLOCATABLE :: iarr
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NIAC
          IWC(n) = 0
        END DO
        ip = NIAC + 1
        DO n = 1, NIAC
          i0 = IAC(n)
          i1 = IAC(n+1) - 1
          nlen = i1 - i0
          ALLOCATE( iarr(nlen) )
          ic = 0
          DO j = i0, i1
            jcol = JAC(j)
            IF ( jcol.EQ.n ) CYCLE
            IWC(jcol) = j
            ic = ic + 1
            iarr(ic) = jcol
          END DO
          CALL SPCGU_ISORT(nlen,iarr)
          IAPC(n) = ip
          DO j = 1, nlen
            jcol = iarr(j)
            JAPC(ip) = jcol
            ip = ip + 1
          END DO
          DEALLOCATE(iarr)
        END DO
        IAPC(NIAC+1) = NNZC + 1
C---------POSITION OF THE FIRST UPPER ENTRY FOR ROW         
        DO n = 1, NIAC
          i0 = IAPC(n)
          i1 = IAPC(n+1) - 1
          JAPC(n) = IAPC(n+1)
          DO j = i0, i1
            jcol = JAPC(j)
            IF ( jcol.GT.n ) THEN
              JAPC(n) = j
              EXIT
            END IF
          END DO
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCCRS
C
C-------SIMPLE INPLACE SORTING ROUTINE FOR AN INTEGER ARRAY      
      SUBROUTINE SPCGU_ISORT(NVAL,IARRAY)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER,INTENT(IN) :: NVAL
        INTEGER,DIMENSION(NVAL),INTENT(INOUT) :: IARRAY
C     + + + LOCAL DEFINITIONS + + +
        integer :: i, j, itemp
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO i = 1, NVAL-1
            DO j = i+1, NVAL
                if(IARRAY(i).GT.IARRAY(j)) then
                    itemp = IARRAY(j)
                    IARRAY(j) = IARRAY(i)
                    IARRAY(i) = itemp
                END IF
            END DO
        END DO
      END SUBROUTINE SPCGU_ISORT
      
      SUBROUTINE SPCGU_MV(NNZC,NIAC,A,D1,D2,IAC,JAC)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)   :: JAC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1
        INTEGER :: icol
        INTEGER :: m, n
        DOUBLEPRECISION :: t
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NIAC
C           ADD DIAGONAL AND OFF-DIAGONAL TERMS
          t     = DZERO
          ic0   = IAC(n)
          ic1   = IAC(n+1)-1
          DO m = ic0, ic1
            icol = JAC(m)
            t  = t + A(m) * D1(icol)
          END DO
          D2(n) = t
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_MV

      DOUBLEPRECISION FUNCTION SPCGU_DP(NIAC,A,B) RESULT(C)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: B
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        C = DZERO
        DO n = 1, NIAC
          C = C + A(n) * B(n)
        END DO
C---------RETURN
        RETURN
      END FUNCTION SPCGU_DP
C
C-------CALCULATE LORDER AND IORDER FOR RCM OR MINIMUM DEGREE ORDERING
C       USING THE GENRCM AND ODRV ROUTINES IN xmdlib_2.f
      SUBROUTINE SPCGU_CALC_ORDER(IORD,NIAC,NNZC,IA,JA,LORDER,IORDER)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IORD
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)  :: IA
        INTEGER, DIMENSION(NNZC),   INTENT(IN)  :: JA
        INTEGER, DIMENSION(NIAC), INTENT(INOUT) :: LORDER
        INTEGER, DIMENSION(NIAC), INTENT(INOUT) :: IORDER
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        INTEGER :: nsp
        INTEGER :: iflag
        INTEGER, DIMENSION(:), ALLOCATABLE :: iwork0, iwork1
C     + + + CODE + + +
        DO n = 1, NIAC
          LORDER(n) = 0
          IORDER(n) = 0
        END DO
        ALLOCATE ( iwork0(NIAC) )
        SELECT CASE ( IORD )
          CASE ( 1 )
            ALLOCATE ( iwork1(NIAC) )
            CALL GENRCM(NIAC, NNZC, IA, JA,
     2                  LORDER, iwork0, iwork1 )
          CASE ( 2 )
            nsp = 3 * NIAC + 4 * NNZC
            ALLOCATE ( iwork1(nsp) )
            CALL ODRV( IA, JA, LORDER, iwork0, iwork1,
     2                 NIAC, NNZC, nsp, iflag )
            IF ( iflag.NE.0 ) THEN
              CALL USTOP('ERROR CREATING MINIMUM DEGREE ORDER'//
     2                   'PERMUTATION ')
            END IF
        END SELECT
C         GENERATE INVERSE OF LORDER
        DO n = 1, NIAC
          IORDER( LORDER(n) ) = n
        END DO
        DEALLOCATE ( iwork0, iwork1 )
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_CALC_ORDER
C
C-------MOVE INACTIVE NODES TO THE END OF THE SYSTEM OF EQUATIONS
C       AND DROP THEIR COLUMNS FROM THE ACTIVE ROWS
      SUBROUTINE SPCGU_CALC_IAORDER(NIA_ALL,NNZ_ALL,IA_ALL,JA_ALL,
     2                              AC,IBOUND,NIAC,NNZC,
     3                              IA_ACT,JA_ACT,A_ACT,
     4                              LORDER_ACT,IORDER_ACT)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)    :: NIA_ALL
        INTEGER, INTENT(IN)    :: NNZ_ALL
        INTEGER, DIMENSION(NIA_ALL+1), INTENT(IN) :: IA_ALL
        INTEGER, DIMENSION(NNZ_ALL), INTENT(IN)   :: JA_ALL
        DOUBLEPRECISION, DIMENSION(NNZ_ALL), INTENT(IN) :: AC
        INTEGER, DIMENSION(NIA_ALL), INTENT(IN)   :: IBOUND
        INTEGER, INTENT(INOUT) :: NIAC
        INTEGER, INTENT(INOUT) :: NNZC
        INTEGER, DIMENSION(NIA_ALL+1), INTENT(INOUT) :: IA_ACT
        INTEGER, DIMENSION(NNZ_ALL), INTENT(INOUT)   :: JA_ACT
        DOUBLEPRECISION, DIMENSION(NNZ_ALL), INTENT(INOUT) :: A_ACT
        INTEGER, DIMENSION(NIA_ALL), INTENT(INOUT) :: LORDER_ACT
        INTEGER, DIMENSION(NIA_ALL), INTENT(INOUT) :: IORDER_ACT
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: j, n
        INTEGER :: i0, i1
        INTEGER :: jcol
        INTEGER :: nt
        INTEGER :: icnt
C     + + + CODE + + +
C---------ACTIVE NODES FIRST, THEN INACTIVE NODES
        icnt = 0
        DO n = 1, NIA_ALL
          IF ( IBOUND(n).NE.0 ) THEN
            icnt = icnt + 1
            IORDER_ACT(icnt) = n
          END IF
        END DO
        NIAC = icnt
        DO n = 1, NIA_ALL
          IF ( IBOUND(n).EQ.0 ) THEN
            icnt = icnt + 1
            IORDER_ACT(icnt) = n
          END IF
        END DO
        DO n = 1, NIA_ALL
          LORDER_ACT( IORDER_ACT(n) ) = n
        END DO
C---------FILL IA_ACT AND JA_ACT, WITH ONLY THE DIAGONAL IN INACTIVE ROWS
        IA_ACT(1) = 1
        icnt = 0
        DO n = 1, NIA_ALL
          nt = IORDER_ACT(n)
          i0 = IA_ALL(nt)
          i1 = IA_ALL(nt+1) - 1
          icnt = icnt + 1
          JA_ACT(icnt) = n
          IF ( IBOUND(nt).NE.0 ) THEN
            DO j = i0+1, i1
              jcol = JA_ALL(j)
              IF ( IBOUND(jcol).EQ.0 ) CYCLE
              icnt = icnt + 1
              JA_ACT(icnt) = LORDER_ACT(jcol)
            END DO
          END IF
          IA_ACT(n+1) = icnt + 1
        END DO
        NNZC = IA_ACT(NIAC+1) - 1
C---------FILL A_ACT
        CALL SPCGU_IAPERM(NIA_ALL,NNZ_ALL,IA_ALL,JA_ALL,AC,
     2                    IBOUND,NIAC,A_ACT,IORDER_ACT)
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_CALC_IAORDER
C
C-------FILL THE ACTIVE ROWS OF A_ACT FROM AC
      SUBROUTINE SPCGU_IAPERM(NIA_ALL,NNZ_ALL,IA_ALL,JA_ALL,AC,
     2                        IBOUND,NIAC,A_ACT,IORDER_ACT)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)    :: NIA_ALL
        INTEGER, INTENT(IN)    :: NNZ_ALL
        INTEGER, DIMENSION(NIA_ALL+1), INTENT(IN) :: IA_ALL
        INTEGER, DIMENSION(NNZ_ALL), INTENT(IN)   :: JA_ALL
        DOUBLEPRECISION, DIMENSION(NNZ_ALL), INTENT(IN) :: AC
        INTEGER, DIMENSION(NIA_ALL), INTENT(IN)   :: IBOUND
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZ_ALL), INTENT(INOUT) :: A_ACT
        INTEGER, DIMENSION(NIA_ALL), INTENT(IN) :: IORDER_ACT
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: j, n
        INTEGER :: i0, i1
        INTEGER :: nt
        INTEGER :: icnt
C     + + + CODE + + +
        icnt = 0
        DO n = 1, NIAC
          nt = IORDER_ACT(n)
          i0 = IA_ALL(nt)
          i1 = IA_ALL(nt+1) - 1
          DO j = i0, i1
            IF ( IBOUND(JA_ALL(j)).EQ.0 ) CYCLE
            icnt = icnt + 1
            A_ACT(icnt) = AC(j)
          END DO
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_IAPERM
C
C-------ROUTINES FROM SPARSKIT TO PERMUTATE A LINEAR SYSTEM OF EQUATIONS
C       IN ORDER TO REORDER THE MATRIX TO MINIMIZE THE BANDWIDTH USING
C       THE REVERSE CUTHILL MCKEE ALGORITHM
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR
! format. i.e., it computes P A Q, where P, Q are permutation matrices.
! P maps row i into row perm(i) and Q maps column j into column qperm(j):
!      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
! In the particular case where Q is the transpose of P (symmetric
! permutation of A) then qperm is not needed.
! note that qperm should be of length ncol (number of columns) but this
! is not checked.
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja,
!    ia = input matrix in a, ja, ia format
! perm       = integer array of length n containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix -- also the destination of column i in case
!         permutation is symmetric (job .le. 2)
!
! qperm      = same thing for the columns. This should be provided only
!         if job=3 or job=4, i.e., only in the case of a nonsymmetric
!        permutation of rows and columns. Otherwise qperm is a dummy
!
! job      = integer indicating the work to be done:
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
!             job = 1      permute a, ja, ia into ao, jao, iao
!             job = 2 permute matrix ignoring real values.
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
!             job = 3      permute a, ja, ia into ao, jao, iao
!             job = 4 permute matrix ignoring real values.
!
! on return:
!-----------
! ao, jao, iao = input matrix in a, ja, ia format
!
! in case job .eq. 2 or job .eq. 4, a and ao are never referred to
! and can be dummy arguments.
! Notes:
!-------
!  1) algorithm is in place
!  2) column indices may not be sorted on return even  though they may be
!     on entry.
!----------------------------------------------------------------------c
! local variables
      integer locjob, mod
!
!     locjob indicates whether or not real values must be copied.
!
      locjob = mod(job,2)
!
! permute rows first
!
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
!
! then permute columns
!
      locjob = 0
!
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob)
      else
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob)
      endif
!
      return
!-------end-of-dperm----------------------------------------------------
      end subroutine
!-----------------------------------------------------------------------
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm       = integer array of length nrow containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1)
!
!     determine pointers for output matix.
!
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
!
! get pointers from lengths
!
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
!
! copying
!
      do 100 ii=1,nrow
!
! old row = ii  -- new row = iperm(ii) -- ko = new pointer
!
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
!
      return
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine

!-----------------------------------------------------------------------
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.
! the result is written in the output matrix  ao, jao, iao.
! cperm computes B = A P, where  P is a permutation matrix
! that maps column j into column perm(j), i.e., on return
!      a(i,j) becomes a(i,perm(j)) in new matrix
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow       = row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm      = integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values ao and ignore iao.
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!----------------------------------------------------------------------c
! local parameters:
      integer k, i, nnz
!
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k))
 100  continue
!
!     done with ja array. return if no need to touch values.
!
      if (job .ne. 1) return
!
! else get new pointers -- and copy values too.
!
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
!
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
!
      return
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine
!----------------------------------------------------------------------- 
      subroutine vperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x 
! according to the permutation array perm(*), i.e., on return, 
! the vector x satisfies,
!
!	x(perm(j)) :== x(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! x	= input vector
!
! on return:
!---------- 
! x	= vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables 
      real*8 tmp, tmp1
!
      init      = 1
      tmp       = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
!     
! loop
! 
 6    k = k+1
!
! save the chased element --
! 
      tmp1      = x(ii) 
      x(ii)     = tmp
      next      = perm(ii) 
      if (next .lt. 0 ) goto 65
!     
! test for end 
!
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
!
! end loop 
!
      goto 6
!
! reinitilaize cycle --
!
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp       = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
!     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
!     
      return
!-------------------end-of-vperm--------------------------------------- 
!-----------------------------------------------------------------------
      end subroutine
      END MODULE PCGULIBMODULE
