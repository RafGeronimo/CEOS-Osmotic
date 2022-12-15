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
    
module ModFixedOsmoticPressureModel

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModOsmoticPressureModel

    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheOsmoticPressureModel": Attributes and methods of the Osmotic Pressure model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type FixedOsmoticPressureModelProperties

        !Class atributes
        ! -----------------------------------------------------------------------------------
            real(8) :: phi_i                            !Internal osmotic coefficient
            real(8) :: phi_e                            !External osmotic coefficient
            real(8) :: c_e                              !External bath concentration
            real(8) :: c_F0                             !Initial fixed charge density
            real(8) :: Temp                             !Absolute temperature
            real(8) :: constR                           !Ideal gas constant
            real(8), dimension(3,3) :: F                !Deformation gradient
            real(8) :: n0                               !Initial fluid volume fraction

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheOsmoticPressureModel": Attributes and methods of the Osmotic Pressure model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassOsmoticPressureModel) :: ClassFixedOsmoticPressureModel

    ! Class Attributes : Usually the state variables (instant and internal variables)
    !----------------------------------------------------------------------------------------
        type (FixedOsmoticPressureModelProperties), pointer :: Properties => null()
        real(8) :: c_F                              !Fixed charge density
       !real (8) :: Lambda                          !Osmotic pressure gradient derivative constant
        
    contains

        ! Class Methods
        !----------------------------------------------------------------------------------
        procedure :: GetLambda                                   => GetLambda_Fixed
        procedure :: GetFixedChargeDensityAndOsmoticGradient     => GetFixedOsmoticPressureGradient
        procedure :: CopyOsmoticPressureProperties               => CopyProperties_FixedOP    
        procedure :: ReadOsmoticPressureParameters               => ReadOsmoticPressureParameters_Fixed
        
    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains
    
    !==========================================================================================
    ! Specific Procedures
    !==========================================================================================
    
        subroutine GetFixedOsmoticPressureGradient(this, LC, ST)
        
            class(ClassFixedOsmoticPressureModel) :: this
            integer :: LC, ST
            
            this%c_F = this%Properties%c_F0
            
            !this%DeltaPi = this%Properties%phi_i*this%Properties%constR*this%Properties%Temp*(((this%c_F**2.0d0) + 4.0d0*(this%Properties%c_e**2.0d0))**(1.0d0/2.0d0)) - 2.0d0*this%Properties%phi_e*this%Properties%constR*this%Properties%Temp*this%Properties%c_e
            
            !this%DeltaPi = this%Properties%phi_i*this%Properties%constR*this%Properties%Temp * SQRT((this%c_F**2.0d0) + (this%Properties%c_e**2.0d0)) - (this%Properties%phi_e*this%Properties%constR*this%Properties%Temp*this%Properties%c_e)
            
            !Perfect Osmometer (FEBio - model)
            
            this%DeltaPi = -this%Properties%constR*this%Properties%Temp*this%Properties%c_e
            
            
        end subroutine
        
        !==========================================================================================

        !==========================================================================================
        subroutine CopyProperties_FixedOP(this,Reference)

             class(ClassFixedOsmoticPressureModel)     :: this
             class(ClassOsmoticPressureModel)          :: Reference

             select type ( Reference )
            
                 class is ( ClassFixedOsmoticPressureModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_FixedOP"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine ReadOsmoticPressureParameters_Fixed(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassFixedOsmoticPressureModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser) :: DataFile

		    !************************************************************************************
		    character(len=100),dimension(6) :: ListOfOptions, ListOfValues
		    logical,dimension(6)            :: FoundOption
		    integer                         :: i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=['phi_i', 'phi_e', 'c_e', 'c_F0', 'n0', 'Temp']

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadOsmoticPressureParameters_Fixed :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
                       
            this%Properties%phi_i = ListOfValues(1)
            this%Properties%phi_e = ListOfValues(2)
            this%Properties%c_e   = ListOfValues(3)
            this%Properties%c_F0  = ListOfValues(4)
            this%Properties%n0    = ListOfValues(5)
            this%Properties%Temp  = ListOfValues(6)

            !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetLambda_Fixed(this, LC, ST)
        
            use ModMathRoutines
            
            class(ClassFixedOsmoticPressureModel)  :: this
            integer :: LC, ST
            
            !Lambda (Constant needed to obtain the Osmotic Pressure Gradient derivative in relation to displacement u)
            this%Lambda = 0.0d0
            
            
        end subroutine
        !==========================================================================================

    end module

