C********************************************************************C
C                                                                    C
C                              'TRACER'                              C
C                                                                    C
C  A PROGRAM TO CALCULATE THE DISPERSION NO./NO. OF TANKS IN SERIES  C
C  EQUIVALENT TO A REAL SYSTEM FROM WHICH LITHIUM CONCENTRATION/TIME C
C                DATA HAS BEEN MEASURED                              C
C  DATA SHOULD BE INPUT VIA A DATA FILE WITH 2 COLUMNS:-             C
C  TIME AND CONCENTRATION RESPECTIVELY                               C
C                                                                    C
C********************************************************************C
      CHARACTER NAME *20,DATFILE *20,INPUTFILE *20,ANS1,ANS2
      CHARACTER RESULTSFILE*24
      INTEGER COUNT
      REAL   NDUL,MEANTIME
      COMMON N
      DIMENSION CONC(300),TIME(1:300),THETA(1:300),SEE(300)

C               TITLE PAGE
C               **********

      CHARACTER Y*4,YY*4,CLR*4
C     CALL SPACE(HEADER)
      Y=CHAR(27)//'[7m'
      YY=CHAR(27)//'[0m'
      CLR= CHAR(27)//'[2J'
      WRITE(*,*)CLR
      WRITE(*,'(1X,A4,A,A4)')Y,'  Copyright WRc (C) 1987  ',YY
      WRITE(*,'(1X,A4,A,A4)')Y,'Activated Sludge Processes',YY
      WRITE(*,'(/////)')
      WRITE(*,'(33X,A4,A18,A4)')Y,'TRACER VERSION 1.1',YY
      WRITE(*,*)
      WRITE(*,'(37X,A4,A10,A4)')Y,'(12.11.87)',YY
      DO 2 K = 1,5
    2 WRITE(*,*)
      CALL RETURN
C     HEADER='      SPDG...   LITHIUM TRACER TEST ANALYSIS   ...SPDG'
C     CALL SPACE(HEADER)

C               TERMINAL OR FILE INPUT?
C               ***********************
C $ in the format statement is a Microsoft extension to prevent
C a carriage return.
C
      WRITE(*,'(///1X,''  ENTER NAME OF PLANT  '',$)')
      READ(*,'(A)') NAME
      RESULTSFILE = NAME(1:)//'.DAT'
      OPEN(UNIT=9,FILE=RESULTSFILE,STATUS='NEW')
C     CALL SPACE(HEADER)
      WRITE(*,'(///1X,''  ENTER NO. OF DATA POINTS  '',\)')
      READ(*,'(I5)') N
      IF(N .GT. 300)THEN
        WRITE(*,'(1X,''        TOO MANY DATA POINTS !!!!'')')
        CALL RETURN
        GO TO 1000
      ENDIF
C     CALL SPACE(HEADER)
      WRITE(*,130)
  130 FORMAT(1X,' IS THE DATA ALREADY ON FILE? (Y OR N) ',\)
      READ(*,140)ANS1
  140 FORMAT(A1)
      IF(ANS1.EQ.'Y' .OR. ANS1 .EQ. 'y')THEN
        CALL READFILE(DATFILE,TIME,CONC)
      ELSEIF(ANS1.EQ.'N' .OR. ANS1 .EQ. 'n')then
        CALL TERMINAL(INPUTFILE,TIME,CONC)
      ENDIF
      CALL RETURN

C           REQUIRED RELATIONS FOR CALCULATION
C           **********************************

      LL        = 0
      SUM1      = 0.0
      SUM2      = 0.0
      DO 10 I = 1,N-1
        SUM1 = SUM1+(0.5*(CONC(I+1)+CONC(I))*(TIME(I+1)-TIME(I)))
        SUM2 = SUM2+(0.25*(CONC(I+1)+CONC(I))*(TIME(I+1)**2-TIME(I)**2))
   10 CONTINUE


C     CALCULATE MEAN RESIDENCE TIME, NORMALISED TIMES AND CONCENTRATIONS
C     ******************************************************************

      MEANTIME=SUM2/SUM1
      AREA    = SUM1
      CO = AREA/MEANTIME
      DO 20 I=1,N
         THETA(I)=TIME(I)/MEANTIME
         SEE(I)=CONC(I)/CO
   20 CONTINUE
      SIGMA = 2.
      IF(LL .EQ. 1)GO TO 100
  151 I=1
  150 IF(THETA(I).LT.SIGMA)THEN
         I=I+1
         IF(I.LE.N)THEN
            I=I+1
            GO TO 150
         ENDIF
      ELSEIF(THETA(I).GT.2.00)THEN
C        CALL SPACE(HEADER)
         WRITE(*,160) SIGMA
         READ(*,170)ANS2
      ENDIF
  160 FORMAT(//1X,'DO YOU WANT TO IGNORE DATA POINTS PAST THETA = ',
     +F3.1,' ? (Y/N)') 
  170 FORMAT(A1)
      IF (ANS2.EQ.'Y' .OR. ANS2 .EQ. 'y')THEN
      CA  LL DOCK(THETA,SEE,CO,MEANTIME,SIGMA)
      ENDIF
C     CALL SPACE(HEADER)


C            CALCULATE THE VARIANCE (SIGMA SQUARED)
C            **************************************

  100 SUM4 = 0.
      SUM3 = 0.      
      DO 40 I = 1,N-1
        SUM3 = SUM3+(0.5*(SEE(I+1)+SEE(I))*(THETA(I+1)-THETA(I)))
        TEMP2= THETA(I+1)**2-THETA(I)**2
        SUM4=SUM4+0.125*(SEE(I+1)+SEE(I))*(THETA(I+1)+THETA(I))*TEMP2
   40 CONTINUE
      AREA= SUM3*CO*MEANTIME
      VAR = (SUM4/SUM3) - 1.
      IF (VAR.LE.0.)THEN
         WRITE(6,250)
         WRITE(9,250)
         GO TO 55
      ELSE
      PRINT*,'       DISPERSION NUMBER WILL NOW BE CALCULATED'
      CALL RETURN
      ENDIF


C            CALCULATE D/uL AND NUMBER OF TANKS
C            **********************************

      DUL= VAR/2.
      KOUNT  = 0
      IF(VAR .GE. 1.0)THEN
        WRITE(6,'('' This plant has a high degree of backmixing.'')')
        WRITE(6,'('' There is also a lot of shortcircuiting and'')')
        WRITE(6,'('' deadspace. Because of these problems it is'')')
        WRITE(6,'('' to know the precise end of the tail of the'')')
        WRITE(6,'('' output response.'',/)')
        WRITE(6,'('' The variance will be recalculated cutting off'')')
        WRITE(6,'('' the tail lying beyond THETA = 2 STD. DEV.'')')
        WRITE(6,'('' This ensures that at least 95% of the data'')')
        WRITE(6,'('' available  is used.'')')
        CALL RETURN
        SIGMA = 1.96*(SQRT(VAR)) + 1.
        IF(SIGMA  .LT.  THETA(N))THEN
          GO TO 151
        ELSE
          DUL = 1./(3.*(1. - .1E-10))
        ENDIF
      ELSEIF(VAR .GT. 0.98  .AND.  VAR .LT. .1E1)THEN
        DUL  = 1./(3.*(1. - VAR))
        GO TO 54
      ELSE
   50   FUNC=2.*DUL-(2.*DUL**2)*(1.-EXP(-1./DUL))-VAR
        DERIV= 2.-4.*DUL+((4.*DUL+2.)*EXP(-1./DUL))
        NDUL = DUL-(FUNC/DERIV)
        DIFF = NDUL-DUL
        IF (ABS(DIFF) .GT. 0.0001)THEN
          IF(KOUNT .GT. 100)GO TO 53
          KOUNT=KOUNT+1
          DUL  = NDUL
          GO TO 50
        ENDIF
      ENDIF
      GO TO 54
   53 PRINT*,' CANNOT CONVERGE TO A D/uL VALUE AFTER 100 ITERATIONS'
      CALL RETURN
      DUL = 1./(3.*(1. - VAR))
   54 TANKS = (1./VAR)


C                   OUTPUT OF RESULTS
C                   *****************

   55 L=6
   60 WRITE(L,180)NAME
      WRITE(L,'(//5X,''TIME'',8X,''CONC'',15X,''THETA'',8X,''CNORM'')')
      WRITE(L,200)(TIME(I),CONC(I),THETA(I),SEE(I),I=1,N)
      WRITE(L,210)MEANTIME,AREA,CO
      WRITE(L,'(//1X,70(''*''))')
      WRITE(L,230)VAR,DUL,TANKS
      IF(L.EQ.6)THEN
         L=L+3
         GO TO 60
      ENDIF
      WRITE(6,240) RESULTSFILE
  180 FORMAT(//1X,70('*')/23X,'TRACER RESULTS'/23X,'FOR ',A20,
     +/1X,70('*'))
  200 FORMAT(/1X,F9.4,5X,F8.4,13X,F6.4,7X,F6.4)
  210 FORMAT(//1X,'MEAN RESIDENCE TIME=',F8.4,2X,'AREA=',
     +F8.2,2X,'Co=',F9.4)
  230 FORMAT(/1X,'VARIANCE  =  ',F6.4,10X,'D/UL=',G11.4,10X,
     +'N=',F8.4//1X,70('*'))
  240 FORMAT(//1X,'OUTPUT IS STORED IN FILE ',A24)
  250 FORMAT(//1X,'NEGATIVE VARIANCE... PLEASE CHECK YOUR DATA !!!!')
      CLOSE(9)
 1000 STOP
      END      


C               DATA ENTRY FROM ALREADY EXISTING FILE
C               *************************************

      SUBROUTINE READFILE(DATFILE,TIME,CONC)
      COMMON N
      CHARACTER DATFILE *20
      DIMENSION TIME(1:300),CONC(1:300)
      WRITE(*,'(/1X,''ENTER NAME OF INPUT FILE '',\)')
      READ(*,'(A20)')DATFILE
      OPEN(UNIT=8,FILE=DATFILE,STATUS='OLD')
      DO 5015 J=1,N
         READ(8,*)TIME(J),CONC(J)
 5015 CONTINUE
      CLOSE(8)
      WRITE(*,'(///1X,I3,'' DATA POINTS ENTERED'')') N
      RETURN
      END


C               DATA ENTRY VIA THE TERMINAL
C               ***************************

      SUBROUTINE TERMINAL(INPUTFILE,TIME,CONC)
      COMMON N
      CHARACTER INPUTFILE *20
      DIMENSION TIME(1:300),CONC(1:300)
      WRITE(*,'(//1X,''GIVE A FILE NAME FOR THE INPUT DATA '',\)')
      READ(*,'(A20)')INPUTFILE
      OPEN(UNIT=7,FILE=INPUTFILE,STATUS='NEW')
      DO 620 J=1,N
        WRITE(*,'(//1X,'' ENTER TIME '',\)')
        READ*,TIME(J)
        WRITE(*,'(//1X,'' ENTER CONCENTRATION OF Li '',\)')
        READ*,CONC(J)
        WRITE(7,'(/1X,F9.4,1X,F9.4)') TIME(J),CONC(J)
  620 CONTINUE
      CLOSE(7)
      WRITE(*,'(//1X,I3,1X,'' DATA POINTS ENTERED '')') N
      RETURN
      END


C               DOCK POINTS PAST THETA=2
C               ************************

      SUBROUTINE DOCK(TH,SE,C0,MIU,SIG)
      DIMENSION TH(1:300),SE(1:300),C(1:300),T(1:300)
      COMMON N
      REAL MIU
      DO  5 I=1,N
        IF(TH(I)  .EQ.  SIG)THEN
          CTHA=SE(I)
        ELSEIF(TH(I) .GT. SIG .AND. TH(I-1) .LT. SIG)THEN
          DELTA=(TH(I)-TH(I-1))
          CTHA=SE(I-1)+(SE(I)-SE(I-1))*(SIG-TH(I-1))/DELTA
        ENDIF
    5 CONTINUE
      DO 10 J=1,N
        IF (TH(J) .GT. SIG)THEN
          SE(J)=0.0
        ELSE
          SE(J) = SE(J) - CTHTA*TH(J)/SIG
        ENDIF
   10 CONTINUE

C     CALCULATE NEW MEAN RESIDENCE TIME AND AREA UNDER CURVE
C     ======================================================

      DO 15 K = 1,N
        C(K) = C0*SE(K)
   15   T(K) = TH(K)*MIU
      SUM1  = 0.
      SUM2  = 0.
      DO 20 L = 1,N-1
        SUM1  = SUM1+(0.5*(C(L+1)+C(L))*(T(L+1)-T(L)))
        SUM2  = SUM2+(0.25*(C(L+1)+C(L))*(T(L+1)**2-T(L)**2))
   20 CONTINUE
      MIU   = SUM2/SUM1
      C0    = SUM1/MIU
      DO 25 II = 1,N
        SE(II)= C(II)/C0
   25   TH(II)= T(II)/MIU
      RETURN
      END
 
*************************************************************************
*                                                                       *
*       Subroutine SPACE clears the screen an returns the cursor        *
*       to the top of the page.                                         *
*                                                                       *
*************************************************************************

C       SUBROUTINE SPACE(CONTROL)
C        CHARACTER *80 CONTROL
C        CHARACTER *8 Y,ZA
C        CHARACTER *3 Z/'[H'/,ZZ/'[7m'/
C        CHARACTER *3 X/'[2J'/,YY/'[0m'/,AZ/'[1m'/
C
C        ZA=CHAR(27)//AZ
C        PRINT *,ZA
C        Y=CHAR(27)//X
C        PRINT*,Y
C        Y=CHAR(27)//Z//CHAR(27)//ZZ
C        PRINT*,Y,CONTROL
C        Y=CHAR(27)//YY
C        PRINT*,Y
C        ZA=CHAR(27)//AZ
C        PRINT *,ZA
C
C       RETURN
C       END

 

*************************************************************************
*                                                                       *
*       Subroutine RETURN stops the screen form clearing so that        *
*       the page can be read. The return key has to be pressed to       *
*       continue.                                                       *
*                                                                       *
*************************************************************************

        SUBROUTINE RETURN
        CHARACTER*1 ANS
        CHARACTER*4 X,X1
        WRITE(*,*)
        X=CHAR(27)//'[7m'
        X1=CHAR(27)//'[0m'
        WRITE(*,'(1X,A4,A25,A4)')X,' Press Return to continue',X1
        READ(5,1)ANS
  1     FORMAT(A)
        RETURN
        END
