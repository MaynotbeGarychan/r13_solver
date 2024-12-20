      
      SUBROUTINE LUDCMP (A, N, NP, INDX, D)
C-----  LU decomposition
C-----  Use single precision on cray
C
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NMAX=200, TINY=1.0E-20)
        DIMENSION A(NP,NP), INDX(N), VV(NMAX)
        D=1.
        DO I=1,N
            AAMAX=0.
            DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
            END DO
            IF (AAMAX.EQ.0.) THEN
            PRINT *,A(1,1),A(1,2),A(1,3)
            PRINT *,A(2,1),A(2,2),A(2,3)
            PRINT *,A(3,1),A(3,2),A(3,3)
c            PAUSE 'Singular matrix.'
            ENDIF
            VV(I)=1./AAMAX
        END DO
        DO J=1,N
            DO I=1,J-1
            SUM=A(I,J)
            DO K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
            END DO
            AAMAX=0.
            DO I=J,N
            SUM=A(I,J)
            DO K=1,J-1
                SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
                IMAX=I
                AAMAX=DUM
            END IF
            END DO
            IF (J.NE.IMAX) THEN
            DO K=1,N
                DUM=A(IMAX,K)
                A(IMAX,K)=A(J,K)
                A(J,K)=DUM
            END DO
            D=-D
            VV(IMAX)=VV(J)
            END IF
            INDX(J)=IMAX
            IF (A(J,J).EQ.0.) A(J,J)=TINY
            IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
                A(I,J)=A(I,J)*DUM
            END DO
            END IF
        END DO
        RETURN
        END

      SUBROUTINE LUBKSB (A, N, NP, INDX, B)
C-----  Linear equation solver based on LU decomposition
C-----  Use single precision on cray
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NP,NP), INDX(N), B(N)
        II=0
        DO I=1,N
            LL=INDX(I)
            SUM=B(LL)
            B(LL)=B(I)
            IF (II.NE.0) THEN
            DO J=II,I-1
                SUM=SUM-A(I,J)*B(J)
            END DO
            ELSE IF (SUM.NE.0.) THEN
            II=I
            END IF
            B(I)=SUM
        END DO
        DO I=N,1,-1
            SUM=B(I)
            IF (I.LT.N) THEN
            DO J=I+1,N
                SUM=SUM-A(I,J)*B(J)
            END DO
            END IF

            B(I)=SUM/A(I,I)
        END DO
        RETURN
        END