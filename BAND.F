C   ==============
      SUBROUTINE ABDGXY(ABD, bigG, delcout, Nin,NJin)
C   ==============
      implicit double precision(a-h,o-z)
      integer*4 Nin,NJin,N,NJ,ii,jj,kk
      REAL*8 ABD(Nin,3*Nin,NJin), bigG(Nin, NJin), 
     & delCout(Nin, NJin)
      COMMON A(12,12),B(12,12),delC(12,10003),D(12,25),G(12),
     &X(12,12),Y(12,12),N,NJ

      N=Nin
      NJ=NJin
    
      do j=1,NJ

        if(j.eq.1) then
            do 1 ii=1,n
                do 10 kk=1,n
                    B(ii,kk)=ABD(ii,kk,    j)
                    D(ii,kk)=ABD(ii,n+kk,  j)
                    X(ii,kk)=ABD(ii,2*n+kk,j)
10              continue
                G(ii)=bigG(ii,j) 
1           continue

        else if (j.eq.NJ) then 
            do 2 ii=1,n
                do 20 kk=1,n
                    Y(ii,kk)=ABD(ii,kk,    j)
                    A(ii,kk)=ABD(ii,n+kk,  j)
                    B(ii,kk)=ABD(ii,2*n+kk,j)
20              continue
                G(ii)=bigG(ii,j) 
2           continue

        else
            do 3 ii=1,n
                do 30 kk=1,n
                    A(ii,kk)=ABD(ii,kk,    j)
                    B(ii,kk)=ABD(ii,n+kk,  j)
                    D(ii,kk)=ABD(ii,2*n+kk,j)
30              continue
                G(ii)=bigG(ii,j) 
3           continue

        end if

      call BAND(J)

      end do

      do ii=1,n
        do jj=1,nj
            delcout(ii,jj)=delc(ii,jj)
        end do
      end do

      return
 999  format(1x,12(g12.4))
      end

C     ===================
      SUBROUTINE BAND(J)
C     ===================
      implicit double precision (A-H,O-Z)
      integer j
      character*40 determwarning
      DIMENSION E(12,13,10003)
      COMMON A(12,12),B(12,12),delC(12,10003),D(12,25),G(12),
     &X(12,12),Y(12,12),N,NJ

      SAVE E,NP1
      IF (J-2) 1,6,8

 1    NP1=N+1
      DO 2 I=1,N
          D(I,2*N+1)=G(I)
          DO 2 L=1,N
              LPN=L+N
              D(I,LPN)=X(I,L)
 2    CONTINUE

      CALL MATINV(N,2*N+1,DETERM)

      IF (DETERM) 4,3,4
 3    write(determwarning,101) J

 4    DO 5 K=1,N
          E(K,NP1,1)=D(K,2*N+1)
          DO 5 L=1,N
          E(K,L,1)=-D(K,L)
          LPN=L+N
 5        X(K,L)=-D(K,LPN)
      RETURN

 6    DO 7 I=1,N
          DO 7 K=1,N
              DO 7 L=1,N
                  D(I,K)= D(I,K) + A(I,L) * X(L,K)
 7    CONTINUE

 8    IF (J-NJ) 11,9,9

 9    DO 10 I=1,N
          DO 10 L=1,N
              G(I)= G(I) - Y(I,L) * E(L,NP1,J-2)
              DO 10 M=1,N
                  A(I,L)= A(I,L) + Y(I,M) * E(M,L,J-2)
 10   CONTINUE

 11   DO 12 I=1,N
          D(I,NP1)=-G(I)
          DO 12 L=1,N
              D(I,NP1)= D(I,NP1) + A(I,L) * E(L,NP1,J-1)
              DO 12 K=1,N
                B(I,K)= B(I,K) + A(I,L) * E(L,K,J-1)
 12   CONTINUE

      CALL MATINV(N,NP1,DETERM)

      IF (DETERM) 14,13,14

 13   write(determwarning,101) J

 14   DO 15 K=1,N
          DO 15 M=1,NP1
              E(K,M,J)=-D(K,M)
 15   CONTINUE

      IF (J-NJ) 20,16,16

 16   DO 17 K=1,N
          delC(K,J)=E(K,NP1,J)
 17   CONTINUE

      DO 18 JJ=2,NJ
          M=NJ-JJ+1
          DO 18 K=1,N
              delc(K,M)= E(K,NP1,M)
              DO 18 L=1,N
                  delC(K,M)= delC(K,M) + E(K,L,M) * delC(L,M+1)  
 18   CONTINUE

      DO 19 L=1,N
          DO 19 K=1,N
              delC(K,1)= delC(K,1) + X(K,L) * delC(L,3)
 19   CONTINUE
 
 20   RETURN

 101  FORMAT(15H DETERM=0 AT J=,I4)

      END

C     ===================
      SUBROUTINE MATINV(N,M,DETERM)
C     ===================
      implicit double precision (A-H,O-Z)
      double precision determ
      integer n, m
      DIMENSION ID(12)
      COMMON A(12,12),B(12,12),delC(12,10003),D(12,25)

      DETERM=1.0

      DO 1 I=1,N
 1        ID(I)=0

      DO 18 NN=1,N
          BMAX=1.1

          DO 6 I=1,N
              IF (ID(I).NE.0) GOTO 6
              BNEXT=0.0
              BTRY=0.0

              DO 5 J=1,N
                  IF (ID(J).NE.0) GOTO 5
                  IF (DABS(B(I,J)).LE.BNEXT) GOTO 5
                  BNEXT=DABS(B(I,J))
                  IF (BNEXT.LE.BTRY) GOTO 5
                  BNEXT=BTRY
                  BTRY=DABS(B(I,J))
                  JC=J
 5            CONTINUE 

              IF (BNEXT.GE.BMAX*BTRY) GOTO 6
              BMAX=BNEXT/BTRY
              IROW=I
              JCOL=JC 
 6        CONTINUE

          IF (ID(JC).EQ.0) GOTO 8
          DETERM=0.0
          RETURN

 8        ID(JCOL)=1
          IF (JCOL.EQ.IROW) GOTO 12

 9        DO 10 J=1,N
              SAVE=B(IROW,J)
              B(IROW,J)=B(JCOL,J)
 10           B(JCOL,J)=SAVE

          DO 11 K=1,M
              SAVE=D(IROW,K)
              D(IROW,K)=D(JCOL,K)
 11           D(JCOL,K)=SAVE

 12       F=1.0/B(JCOL,JCOL)
          DO 13 J=1,N
 13            B(JCOL,J)=B(JCOL,J)*F
          DO 14 K=1,M
 14           D(JCOL,K)=D(JCOL,K)*F

          DO 18 I=1,N
              IF (I.EQ.JCOL) GOTO 18
              F=B(I,JCOL)
              DO 16 J=1,N
 16               B(I,J)=B(I,J)-F*B(JCOL,J)
              DO 17 K=1,M
 17               D(I,K)=D(I,K)-F*D(JCOL,K)
 18   CONTINUE

      RETURN
      END
      
