!##################################################################################################
! This module has the attributes and methods for the Ideal Donnan Osmotic Pressure Model
! -------------------------------------------
! Date: 2021/06
!
! Authors:  Rafael Geronimo
!           
!    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
    
module ModIdealDonnanOsmoticPressureModel

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModOsmoticPressureModel
    use ModLoadHistoryData

    implicit none

    
    !==================================================
    type ClassCF0History

        type(ClassLoadCase), pointer, dimension(:) :: LoadCase => null()
        integer :: nLoadCases
        character*100 :: TableName

        contains
            procedure :: ReadCF0TimeDiscretization
            procedure :: ReadCF0ValueDiscretization
            procedure :: CreateCF0ConstantLoadHistory

        end type
    !==================================================
        
    !==================================================
    type ClassCextHistory

        type(ClassLoadCase), pointer, dimension(:) :: LoadCase => null()
        integer :: nLoadCases
        character*100 :: TableName

        contains
            procedure :: ReadCextTimeDiscretization
            procedure :: ReadCextValueDiscretization
            procedure :: CreateCextConstantLoadHistory

        end type
    !==================================================
        
        
    
    type IdealDonnanOsmoticPressureModelProperties

        !Class atributes
        ! -----------------------------------------------------------------------------------
            real(8)                 :: phi_i                            !Internal osmotic coefficient
            real(8)                 :: phi_e                            !External osmotic coefficient
            real(8)                 :: c_e                              !External bath concentration
            real(8)                 :: c_F0                             !Initial fixed charge density multiplier
            real(8)                 :: n0                               !Initial fluid volume fraction
            real(8)                 :: Temp                             !Absolute temperature
            real(8)                 :: ConstR                           !Universal gas constant
            character(len=100)      :: CF0_LoadCurveTime                !Time discretization of CF0 load curve
            character(len=100)      :: CF0_LoadCurveValue               !Value of CF0 load curve      
            type (ClassCF0History)  :: CF0_CurveData                    !CF0 curve data class
            type (ClassCextHistory) :: Cext_CurveData                   !Cext curve data class
            
    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfOsmoticPressureModel": Attributes and methods of the Osmotic Pressure model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassOsmoticPressureModel) :: ClassIdealDonnanOsmoticPressureModel

    ! Class Attributes : Usually the state variables (instant and internal variables)
    !----------------------------------------------------------------------------------------
        type (IdealDonnanOsmoticPressureModelProperties), pointer :: Properties => null()
        real(8) :: c_F  
      
    contains

        ! Class Methods
        !----------------------------------------------------------------------------------
        procedure :: GetLambda                                   => GetLambda_IdealDonnan
        procedure :: GetFixedChargeDensityAndOsmoticGradient     => GetIdealDonnanOPFixedChargeDensityAndOsmoticGradient
        procedure :: CopyOsmoticPressureProperties               => CopyProperties_IdealDonnanOP
        procedure :: ReadOsmoticPressureParameters               => ReadOsmoticPressureParameters_IdealDonnan

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    

    contains
    
        !==========================================================================================
        ! Specific Procedures
        !==========================================================================================
        subroutine GetIdealDonnanOPFixedChargeDensityAndOsmoticGradient(this, LC, ST)
        
            use ModMathRoutines
            use ModLoadHistoryData
            
            class(ClassIdealDonnanOsmoticPressureModel) :: this
            integer :: LC, ST
            real(8) :: c_F0, c_e
            
            
            c_e = this%Properties%c_e*this%Properties%Cext_CurveData%LoadCase(LC)%Step(ST)%FinalVal
            
            c_F0 = this%Properties%c_F0*this%Properties%CF0_CurveData%LoadCase(LC)%Step(ST)%FinalVal
            
            this%c_F = c_F0*(this%Properties%n0/(this%Properties%n0 - 1 + this%J))
            
            !this%DeltaPi = this%Properties%phi_i*constR*this%Properties%Temp*(SQRT((this%c_F**(2.0d0)) + 4.0d0*(this%Properties%c_e**(2.0d0)))) - 2.0d0*this%Properties%phi_e*constR*this%Properties%Temp*this%Properties%c_e
            
            !this%DeltaPi = this%Properties%phi_i*constR*this%Properties%Temp * SQRT(this%c_F**2.0d0 + 4.0d0*(this%Properties%c_e**2.0d0)) - 2.0d0*this%Properties%phi_e*constR*this%Properties%Temp*this%Properties%c_e
            
            this%DeltaPi = this%Properties%phi_i*this%Properties%constR*this%Properties%Temp * SQRT((this%c_F*this%c_F) + (c_e*c_e)) - (this%Properties%phi_e*this%Properties%constR*this%Properties%Temp*c_e)
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine CopyProperties_IdealDonnanOP(this,Reference)

             class(ClassIdealDonnanOsmoticPressureModel)     :: this
             class(ClassOsmoticPressureModel)                :: Reference

             select type ( Reference )
            
                 class is ( ClassIdealDonnanOsmoticPressureModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_IdealDonnanOP"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine ReadOsmoticPressureParameters_IdealDonnan(this,DataFile)

            use ModParser
            use ModLoadHistoryData

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassIdealDonnanOsmoticPressureModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type (ClassParser) :: DataFile

		    !************************************************************************************
		    character(len=100),dimension(11) :: ListOfOptions, ListOfValues
		    logical,dimension(11)            :: FoundOption
		    integer                         :: i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)
            
            ListOfOptions=['phi_i', 'phi_e', 'c_e', 'c_F0', 'n0', 'Temp', 'gas constant R' ,'CF0 Curve time', 'CF0 Load Curve', 'Cext Curve time', 'Cext Load Curve']

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadOsmoticPressureParameters_IdealDonnanModel :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
                       
            this%Properties%phi_i  = ListOfValues(1)
            this%Properties%phi_e  = ListOfValues(2)
            this%Properties%c_e    = ListOfValues(3)
            this%Properties%c_F0   = ListOfValues(4)
            this%Properties%n0     = ListOfValues(5)
            this%Properties%Temp   = ListOfValues(6)
            this%Properties%constR = ListOfValues(7)
            
            call this%Properties%CF0_CurveData%ReadCF0TimeDiscretization (ListOfValues(8))
            
            call this%Properties%CF0_CurveData%ReadCF0ValueDiscretization (ListOfValues(9))
            
            call this%Properties%Cext_CurveData%ReadCextTimeDiscretization (ListOfValues(10))
            
            call this%Properties%Cext_CurveData%ReadCextValueDiscretization (ListOfValues(11))
            
            
            !************************************************************************************


        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetLambda_IdealDonnan(this, LC, ST)
        
        
            use ModMathRoutines
        
            class(ClassIdealDonnanOsmoticPressureModel)  :: this
            integer  :: LC, ST
            real (8) :: c_F0, c_e
            
            
            c_F0 = this%Properties%c_F0*this%Properties%CF0_CurveData%LoadCase(LC)%Step(ST)%FinalVal
            
            c_e = this%Properties%c_e*this%Properties%Cext_CurveData%LoadCase(LC)%Step(ST)%FinalVal
            
            !Lambda (Constant needed to obtain the Osmotic Pressure Gradient derivative in relation to displacement u)
            
            !this%Lambda = (-1.0d0)*this%J*this%Properties%phi_i*constR*this%Properties%Temp * (( (this%Properties%n0**2.0d0)*(this%Properties%c_F0**2.0d0) ) / ( (this%Properties%n0 - 1.0d0 + this%J )**3.0d0 )) / ( SQRT( ( (this%Properties%c_F0**2.0d0) *(this%Properties%n0 /  (this%Properties%n0 - 1.0d0 + this%J))**2.0d0 ) + 4.0d0*(this%Properties%c_e**2.0d0 ) ))

            this%Lambda = (-1.0d0)*this%J*this%Properties%phi_i*this%Properties%constR*this%Properties%Temp * (( (this%Properties%n0**2.0d0)*(c_F0**2.0d0) ) / ( (this%Properties%n0 - 1.0d0 + this%J )**3.0d0 )) / ( SQRT( ( (c_F0**2.0d0) *(this%Properties%n0 /  (this%Properties%n0 - 1.0d0 + this%J))**2.0d0 ) + (c_e*c_e) ))

            
        end subroutine
        !==========================================================================================
        
        
        
        !==========================================================================================
        subroutine ReadCF0TimeDiscretization(this,FileName)

            use ModParser

            implicit none

            class(ClassCF0History) :: this
            character (len=*) :: FileName

            type (ClassParser) :: TimeData
            character(len=255) :: String
            integer            :: nLC, i, n
            real(8)            :: DeltaTime
            real(8), allocatable, dimension(:) :: InitialTimes, FinalTimes
            integer, allocatable, dimension(:) :: Steps


            call TimeData%Setup(FileName,26)

            call TimeData%GetNextString(string) ; call TimeData%CheckError

            if ( TimeData%CompareStrings(String,'Number of Load Cases') ) then

                call TimeData%GetNextString(string) ; call TimeData%CheckError

                nLC=string
                allocate(InitialTimes(nLC), FinalTimes(nLC), Steps(nLC))

                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if ( TimeData%CompareStrings(String,'Load Case Initial Time Final Time Steps') ) then

                    do i = 1,nLC
                        read(26,*) n, InitialTimes(i), FinalTimes(i), Steps(i)
                    enddo

                else

                    call TimeData%RaiseError('Expected (Load Case	Initial Time	Final Time	Steps Time) in Time Discretization File')

                endif

            else

                call TimeData%RaiseError('Expected (Number of Load Cases) in Time Discretization File')

            endif

            call TimeData%CloseFile


            allocate(this%LoadCase(nLC))

            ! Guardando para usar no QuasiStaticAnalysis
            this%nLoadCases = nLC

            do i = 1,nLC

                allocate( this%LoadCase(i)%Step(Steps(i)) )

                ! Guardando para usar no QuasiStaticAnalysis
                this%LoadCase(i)%nSteps = Steps(i)

                DeltaTime = (FinalTimes(i)-InitialTimes(i))/ dble(Steps(i))
                do n = 1,Steps(i)

                    this%LoadCase(i)%Step(n)%InitTime  = InitialTimes(i) + dble(n-1)*DeltaTime

                    this%LoadCase(i)%Step(n)%FinalTime = this%LoadCase(i)%Step(n)%InitTime + DeltaTime

                enddo

            enddo


        end subroutine
        !=================================================================================================


        !=================================================================================================
        subroutine ReadCF0ValueDiscretization(this,TableName)

            use ModParser
            use ModMathRoutines        
            implicit none

            class(ClassCF0History) :: this
            character (len=*)       :: TableName
            type (ClassParser)      :: TimeData
           
            character(len=255)                    :: String, FileName
            integer                               :: cont, flagInitial, flagFinal, ST, LC
            real(8), allocatable, dimension(:,:)  :: TimeAndValue
            logical                               :: FileExist
            character(len=255), allocatable, dimension(:)  :: AuxString


            this%TableName = TableName

            call Split(TableName,AuxString,'.')

            if (size(AuxString) == 1) then
                FileName = trim(TableName)//'.tab'
            else
                FileName = trim(TableName)
            endif


            inquire(file=trim(FileName),exist=FileExist)

            if (.not.FileExist) then
                write(*,*) 'Table Could not be found: '//trim(TableName)
                stop
            endif


            ! Leitura da Tabela "contínua" informada
            !---------------------------------------------------------------------

            call TimeData%Setup(FileName,31)

            ! contar a quantidade de números da tabela informada
            cont = 0
            do while(.true.)
                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if (EOF(TimeData)) exit
                if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                    cont = cont + 1
                endif
            enddo

            call TimeData%CloseFile

            allocate( TimeAndValue(cont,2))

            ! Leitura da Tabela Informada
            call TimeData%Setup(FileName,31)
            cont = 0
            do while(.true.)
                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if(eof(TimeData)) exit
                if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                    cont = cont + 1
                    call TimeData%GetOriginalLine(String)
                    read(String,*) TimeAndValue(cont,1) , TimeAndValue(cont,2)
                endif
            enddo

            call TimeData%CloseFile
            !---------------------------------------------------------------------
            ! Interpolação dos valores segundo a discretização do tempo requerida
            ! no arquivo "time discretization"
            !---------------------------------------------------------------------

            do LC = 1, this%nLoadCases

                do ST = 1, this%LoadCase(LC)%nSteps

                    call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%InitTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%InitVal, flagInitial )

                    call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%FinalTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%FinalVal, flagFinal )

                    if (flagInitial == 1 .and. flagFinal == 1 ) then
                        this%LoadCase(LC)%Step(ST)%active = .true.
                    else
                        this%LoadCase(LC)%Step(ST)%active = .false.
                    endif

                enddo

            enddo

        end subroutine
       !=================================================================================================

       !=================================================================================================
        subroutine CreateCF0ConstantLoadHistory(this,value)


            implicit none

            class(ClassCF0History) :: this
            real(8) :: value

            integer                ::  ST, LC


            ! Interpolação dos valores segundo a discretização do tempo requerida
            ! no arquivo "time discretization"
            !---------------------------------------------------------------------

            do LC = 1, this%nLoadCases

                do ST = 1, this%LoadCase(LC)%nSteps

                    this%LoadCase(LC)%Step(ST)%InitVal = value

                    this%LoadCase(LC)%Step(ST)%FinalVal = value

                    this%LoadCase(LC)%Step(ST)%active = .true.

                enddo

            enddo
        
        end subroutine
        !=================================================================================================  
        
        !==========================================================================================
        subroutine ReadCextTimeDiscretization(this,FileName)

            use ModParser

            implicit none

            class(ClassCextHistory) :: this
            character (len=*) :: FileName

            type (ClassParser) :: TimeData
            character(len=255) :: String
            integer            :: nLC, i, n
            real(8)            :: DeltaTime
            real(8), allocatable, dimension(:) :: InitialTimes, FinalTimes
            integer, allocatable, dimension(:) :: Steps


            call TimeData%Setup(FileName,26)

            call TimeData%GetNextString(string) ; call TimeData%CheckError

            if ( TimeData%CompareStrings(String,'Number of Load Cases') ) then

                call TimeData%GetNextString(string) ; call TimeData%CheckError

                nLC=string
                allocate(InitialTimes(nLC), FinalTimes(nLC), Steps(nLC))

                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if ( TimeData%CompareStrings(String,'Load Case Initial Time Final Time Steps') ) then

                    do i = 1,nLC
                        read(26,*) n, InitialTimes(i), FinalTimes(i), Steps(i)
                    enddo

                else

                    call TimeData%RaiseError('Expected (Load Case	Initial Time	Final Time	Steps Time) in Time Discretization File')

                endif

            else

                call TimeData%RaiseError('Expected (Number of Load Cases) in Time Discretization File')

            endif

            call TimeData%CloseFile


            allocate(this%LoadCase(nLC))

            ! Guardando para usar no QuasiStaticAnalysis
            this%nLoadCases = nLC

            do i = 1,nLC

                allocate( this%LoadCase(i)%Step(Steps(i)) )

                ! Guardando para usar no QuasiStaticAnalysis
                this%LoadCase(i)%nSteps = Steps(i)

                DeltaTime = (FinalTimes(i)-InitialTimes(i))/ dble(Steps(i))
                do n = 1,Steps(i)

                    this%LoadCase(i)%Step(n)%InitTime  = InitialTimes(i) + dble(n-1)*DeltaTime

                    this%LoadCase(i)%Step(n)%FinalTime = this%LoadCase(i)%Step(n)%InitTime + DeltaTime

                enddo

            enddo


        end subroutine
        !=================================================================================================


        !=================================================================================================
        subroutine ReadCextValueDiscretization(this,TableName)

            use ModParser
            use ModMathRoutines        
            implicit none

            class(ClassCextHistory) :: this
            character (len=*)       :: TableName
            type (ClassParser)      :: TimeData
           
            character(len=255)                    :: String, FileName
            integer                               :: cont, flagInitial, flagFinal, ST, LC
            real(8), allocatable, dimension(:,:)  :: TimeAndValue
            logical                               :: FileExist
            character(len=255), allocatable, dimension(:)  :: AuxString


            this%TableName = TableName

            call Split(TableName,AuxString,'.')

            if (size(AuxString) == 1) then
                FileName = trim(TableName)//'.tab'
            else
                FileName = trim(TableName)
            endif


            inquire(file=trim(FileName),exist=FileExist)

            if (.not.FileExist) then
                write(*,*) 'Table Could not be found: '//trim(TableName)
                stop
            endif


            ! Leitura da Tabela "contínua" informada
            !---------------------------------------------------------------------

            call TimeData%Setup(FileName,31)

            ! contar a quantidade de números da tabela informada
            cont = 0
            do while(.true.)
                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if (EOF(TimeData)) exit
                if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                    cont = cont + 1
                endif
            enddo

            call TimeData%CloseFile

            allocate( TimeAndValue(cont,2))

            ! Leitura da Tabela Informada
            call TimeData%Setup(FileName,31)
            cont = 0
            do while(.true.)
                call TimeData%GetNextString(string) ; call TimeData%CheckError
                if(eof(TimeData)) exit
                if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                    cont = cont + 1
                    call TimeData%GetOriginalLine(String)
                    read(String,*) TimeAndValue(cont,1) , TimeAndValue(cont,2)
                endif
            enddo

            call TimeData%CloseFile
            !---------------------------------------------------------------------
            ! Interpolação dos valores segundo a discretização do tempo requerida
            ! no arquivo "time discretization"
            !---------------------------------------------------------------------

            do LC = 1, this%nLoadCases

                do ST = 1, this%LoadCase(LC)%nSteps

                    call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%InitTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%InitVal, flagInitial )

                    call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%FinalTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%FinalVal, flagFinal )

                    if (flagInitial == 1 .and. flagFinal == 1 ) then
                        this%LoadCase(LC)%Step(ST)%active = .true.
                    else
                        this%LoadCase(LC)%Step(ST)%active = .false.
                    endif

                enddo

            enddo

        end subroutine
       !=================================================================================================

       !=================================================================================================
        subroutine CreateCextConstantLoadHistory(this,value)


            implicit none

            class(ClassCextHistory) :: this
            real(8) :: value

            integer                ::  ST, LC


            ! Interpolação dos valores segundo a discretização do tempo requerida
            ! no arquivo "time discretization"
            !---------------------------------------------------------------------

            do LC = 1, this%nLoadCases

                do ST = 1, this%LoadCase(LC)%nSteps

                    this%LoadCase(LC)%Step(ST)%InitVal = value

                    this%LoadCase(LC)%Step(ST)%FinalVal = value

                    this%LoadCase(LC)%Step(ST)%active = .true.

                enddo

            enddo
        
        end subroutine
        !=================================================================================================  
    
    end module


