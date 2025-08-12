!
!     Mr. LUAN BIN
!
!     Euler Method FOR SEDIMENT TRANSPORT OVER IMPERVIOUS SURFACE: - (C) COPYRIGHT 2022.10.09
!
! ---------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------

		IMPLICIT NONE
		REAL(8), ALLOCATABLE ::Ct(:),Qt(:),Mt(:), TMDATA(:), VMDATA(:), PMDATA(:), Ct0(:), Cet(:), Cet1(:)
		REAL(kind=8):: h0, tc, R, a, W0, L, B, Ce, b0, kL, kU, Dk, Td, DT, t, fCe, DCT, EDCT, Ctc, PCt, CCt, TP, Ctt, Qtc, Mtc, NSE, k, NSE0, TV0, TC0, PC0, AV, k0, Ctc0, CF, CFL, CFU, DCF, CF0
		INTEGER:: I, J, JMAX, N, ITMP, M, MDATA, O
		CHARACTER:: FILNAM*200

! --------------------------------------------------------------------------------

!		READ IN INPUT FILE NAME
100		WRITE(*,*) 'ENTER INPUT FILE NAME [Test1.dat]. 0 TO EXIT:'
		READ(*,'(A)') FILNAM
		FILNAM = ADJUSTL(FILNAM)
!		IF(LEN_TRIM(FILNAM).EQ.1 .AND. FILNAM(1:1).EQ.'0')  STOP
!		IF(LEN_TRIM(FILNAM).EQ.0)   FILNAM="Test1.dat"

!		OPEN INPUT FILE
		OPEN(11, FILE=FILNAM, STATUS='OLD', IOSTAT=ITMP)
		IF(ITMP .NE. 0) THEN
			WRITE(*,*) 'FILE DOES NOT EXIST !'
			GOTO 100
		END IF

!		OPEN FILES FOR OUTPUT
		N = LEN_TRIM(FILNAM)
		IF(N .GT. 4) THEN
			N = N - 4
		END IF
		OPEN(UNIT = 21, FILE=FILNAM(1:N)//'_LST.DAT', STATUS = 'UNKNOWN')		!LIST FILE
        OPEN(UNIT = 31, FILE=FILNAM(1:N)//'_NSE.DAT', STATUS = 'UNKNOWN')       !calibration file


! ------------------------------------------------------------------------------
!     READ IN COMPUTATIONAL PARAMETERS, ALLOCATE ARRIES
! ------------------------------------------------------------------------------

		WRITE(*,*) 'READING CONTROL PARAMETERS ...'
		READ(11,'(A)')
		READ(11,'(15(10X,F10.8/),10X,F10.8)')   h0, R, a, W0, L, B, kL, kU, Dk, CFL, CFU, DCF, b0, Td, DT, TP
        READ(11,'(10X,I10)') MDATA
		ALLOCATE( TMDATA(MDATA), VMDATA(MDATA), PMDATA(MDATA) )
		READ(11,'(A)')
		DO O = 1, MDATA
			READ(11,*) TMDATA(O), VMDATA(O)
		END DO
        CONTINUE
		CLOSE(11)

		WRITE(*,'(1X,A,//)') 'COMPLETE READING CONTROL FILE'


		JMAX = INT(Td/DT)+1
		WRITE(*,*) JMAX
		tc = (L/(a*(R**(2.0/3.0))))**0.6
		WRITE(*,*) tc
		I = INT(tc/DT)+1
		WRITE(*,*) I
		M = INT(TP/DT)+1
		WRITE(*,*) M

		ALLOCATE( Ct(JMAX+1), Qt(JMAX+1), Mt(JMAX+1), Ct0(JMAX+1), Cet(JMAX+1), Cet1(JMAX+1) )
        WRITE(31,'(//3A13)') 'CF', 'k', 'NSE'

! ------------------------------------------------------------------------------
!     COMPUTATION PART
! ------------------------------------------------------------------------------

        Ct(1) = 0.0
        Qt(1) = 0.0
        NSE = 0.0

!       BEFORE EQUILIBRATION
        DO CF = CFL, CFU, DCF
        DO k = kL, kU, Dk
		Cet(1) = CF*W0*(1.0-exp(-k*h0))/(B*L*h0)
        DO J = 1, I
            t =(J-1)*DT
            Ce = CF*W0*k*exp(-k*h0)*exp(-k*R*t)/(B*L*(b0-k*h0))+(CF*W0*(1.0-exp(-k*h0))/(B*L*h0)-CF*W0*k*exp(-k*h0)/(B*L*(b0-k*h0)))*exp(-b0*t*R/h0)
            fCe = CF*W0*k*exp(-k*h0)*(exp(-k*R*t)-1.0)/(B*L*(b0-k*h0)*(-k*R))+(CF*W0*(1.0-exp(-k*h0))/(B*L*h0)-CF*W0*k*exp(-k*h0)/(B*L*(b0-k*h0)))*(-h0/(b0*R))*(exp(-b0*R*t/h0)-1.0)
            PCt = Ct(J)+DT*DCT(Ct(J),t,h0,R,W0,L,B,Ce,b0,fCe)
            CCt = Ct(J)+DT*DCT(PCt,t,h0,R,W0,L,B,Ce,b0,fCe)
            Cet(J+1) = Ce
            Ct(J+1) = 0.5*(PCt+CCt)
            Qt(J+1) = a*B*(R*t)**(5.0/3.0)
            WRITE(*,*) t ,  Ct(J+1)
        END DO

        PCt=Ct(I+1)+(tc-I*DT)*DCT(Ct(I+1),I*DT,h0,R,W0,L,B,Ce,b0,fCe)
        CCt=Ct(I+1)+(tc-I*DT)*DCT(PCt,I*DT,h0,R,W0,L,B,Ce,b0,fCe)
        Ctc=0.5*(PCt+CCt)
        Qtc=B*L*R

!       AFTER EQUILIBRATION
        DO J = I+1, JMAX
            t=(J-1)*DT
            Ce = CF*W0*k*exp(-k*h0)*exp(-k*R*t)/(B*L*(b0-k*h0))+(CF*W0*(1.0-exp(-k*h0))/(B*L*h0)-CF*W0*k*exp(-k*h0)/(B*L*(b0-k*h0)))*exp(-b0*t*R/h0)
            PCt = Ct(J)+DT*EDCT(Ct(J),t,h0,R,W0,L,B,Ce,b0,a)
            CCt = Ct(J)+DT*EDCT(PCt,t,h0,R,W0,L,B,Ce,b0,a)
            Cet(J+1) = Ce
            Ct(J+1) = 0.5*(PCt+CCt)
            Qt(J+1)=Qtc
            WRITE(*,*) t ,  Ct(J+1)
        END DO

!       calibration
        DO O = 1, MDATA
            J= INT((60.0*TMDATA(O))/DT)+1
            PMDATA(O) = Ct(J)+(Ct(J+1)-Ct(J))*(60.0*TMDATA(O)-DT*(J-1))/DT  !  deviation calculation
        END DO
        TV0 = 0.0
		TC0 = 0.0
		PC0 = 0.0
        DO O = 1, MDATA
            TV0=VMDATA(O)+TV0
			TC0=(VMDATA(O)-PMDATA(O))*(VMDATA(O)-PMDATA(O))+TC0
        END DO
        AV=TV0/MDATA   !average value
        DO O = 1, MDATA
			PC0=(VMDATA(O)-AV)*(VMDATA(O)-AV)+PC0
		END DO
		NSE0=1-TC0/PC0
        WRITE(*,*) CF, k ,  NSE0
		WRITE(31,'(10F14.8)') CF, k ,  NSE0
		IF ( NSE0.GE.NSE ) THEN
            NSE=NSE0
            CF0=CF
            k0=k
            Ctc0=Ctc
            Mtc=Qtc*Ctc0
            DO J=1,JMAX+1
                Ct0(J)=Ct(J)
                Mt(J) = Qt(J)*Ct0(J)
                Cet1(J)=Cet(J)
            END DO
        END IF
        END DO
        END DO
        
! ------------------------------------------------------------------------------
!     WRITE OUTPUT FILE
! ------------------------------------------------------------------------------

		WRITE(21,10)
		WRITE(21,20) JMAX
		WRITE(21,30) h0, R, a, W0, L, B, Td, DT, tc/60.0, Ctc0, Ct0(1), 1000000*Qtc, 1000*Mtc, k0, b0, NSE, CF0
 10		FORMAT(/1X,'INPUT DATA'/)
 20		FORMAT( 10X,'JMAX   =',I8 )
 30		FORMAT( 10X,'h0     =', F14.8/10X,'R      =', F14.8/10X,'a      =', F14.8/	&
	&			10X,'W0     =', F14.8/10X,'L      =', F14.8/10X,'B      =', F14.8/	&
	&			10X,'Td     =', F14.8/10X,'DT     =', F14.8/10X,'tc     =', F14.8/	&
	&			10X,'Ctc    =', F14.8/10X,'Ct0    =', F14.8/10X,'Qtc    =', F14.8/	&
	&			10X,'Mtc    =', F14.8/10X,'k      =', F14.8/10X,'b      =', F14.8/	&
	&			10X,'NSE    =', F14.8/10X,'CF0    =', F14.8 )

		WRITE(21,'(//5A16,1000(A8,I8))') 'TIME(MIN)', 'Ct(g/L)', 'Qt(mL/s)', 'Mt(g/s)', 'Ce(g/L)'

		DO J =1, INT(JMAX/M+1)
            t= (J-1)*DT*M
            WRITE(21,'(1800000E16.8)') t/60.0, Ct0(M*(J-1)+1), 1000000*Qt(M*(J-1)+1), 1000*Mt(M*(J-1)+1), Cet1(M*(J-1)+1)
        END DO
        CLOSE(21)
        DEALLOCATE( Ct, Qt, Mt, Ct0, Cet, Cet1 )
        DEALLOCATE( TMDATA, VMDATA, PMDATA )
    END

!       TSS CONCENTRATION BEFORE EQUILIBRATION
        REAL(8) FUNCTION DCT(Ctt,t,h0,R,W0,L,B,Ce,b0,fCe)
        IMPLICIT NONE
        REAL(kind=8):: h0, R, W0, L, B, t, Ctt, Ce, b0, fCe
        IF( t.EQ.0.0 ) THEN
            DCT = 0.0
        ELSE
            DCT = (t*R*b0*Ce+5.0/3.0*b0*R*fCe)/(h0*t+0.625*R*t*t)-Ctt*(5.0*h0+8.0*R*t)/(3.0*t*(h0+0.625*R*t))
        END IF
		END

!       TSS CONCENTRATION AFTER EQUILIBRATION
        REAL(8) FUNCTION EDCT(Ctt,t,h0,R,W0,L,B,Ce,b0,a)
        IMPLICIT NONE
        REAL(kind=8):: h0, R, W0, L, B, t, Ctt, Ce, b0, a
        IF( t.EQ.0.0 ) THEN
            EDCT = 0.0
        ELSE
            EDCT = R*b0*Ce/(h0+0.625*((R*L/a)**0.6))-Ctt*R/(h0+0.625*((R*L/a)**0.6))
        END IF
		END
