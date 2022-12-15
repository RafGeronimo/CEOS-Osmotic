!##################################################################################################
! This module is used to register a new Permeability Model.
!--------------------------------------------------------------------------------------------------
! Date: 2022/09
!
! Authors:  Rafael Geronimo
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModOsmoticPressureModelLibrary

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModAnalysis
    use ModOsmoticPressureModel                            !Main osmotic pressure module
    use ModIdealDonnanOsmoticPressureModel                 !Ideal Donnan osmotic pressure model
    use ModFixedOsmoticPressureModel                       !Fixed Osmotic pressure model
   
    ! Permeability Models ID registered:
    type ClassOsmoticPressureModels                                
        integer   :: IdealDonnanOsmoticPressureModel                 = 1
        integer   :: FixedOsmoticPressureModel                       = 2
              
    end type

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	type(ClassOsmoticPressureModels), parameter :: OsmoticPressureModels=ClassOsmoticPressureModels()

    contains

		!==========================================================================================
        ! Routine 
        ! PermeabilityModel: Routine that allocates the Osmotic Pressure
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:           Author:
        !==========================================================================================
        subroutine AllocateOsmoticPressureModel( OsmoticPressureModel , AnalysisSettings , nGP , GaussPoints_Osmotic)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            integer , intent(in) :: OsmoticPressureModel , nGP

            ! Output variables
            ! -----------------------------------------------------------------------------------
            class(ClassOsmoticPressureModel), pointer, dimension(:), intent(out) :: GaussPoints_Osmotic 
            
            ! Internal variables: Instance of each available Osmotic Pressure Model.
            ! -----------------------------------------------------------------------------------
            type(ClassIdealDonnanOsmoticPressureModel), pointer , dimension(:)   :: IDOPM
            type(ClassFixedOsmoticPressureModel),       pointer , dimension(:)   :: FOPM
		    !************************************************************************************

            !************************************************************************************
            ! CONSTRUCT THE OSMOTIC PRESSURE VARIABLES IN THE GAUSS POINTS
		    !************************************************************************************
             if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then
                
                select case (OsmoticPressureModel)
                    ! -------------------------------------------------------------------------
                    ! Ideal Donnan Osmotic Pressure
                    ! -------------------------------------------------------------------------------
                    case (OsmoticPressureModels % IdealDonnanOsmoticPressureModel)
                            allocate( IDOPM(nGP) )
                            GaussPoints_Osmotic => IDOPM
                    ! -------------------------------------------------------------------------------
                    ! Fixed Osmotic Pressure
                    ! -------------------------------------------------------------------------------
                    case (OsmoticPressureModels % FixedOsmoticPressureModel)
                            allocate( FOPM(nGP) )
                            GaussPoints_Osmotic => FOPM
                    ! -------------------------------------------------------------------------------
                    case default
                        call Error( "Error: Osmotic Pressure Model not registered.")
                end select
            else
                call Error("Error: Construct Osmotic Pressure Models: Analysis type must be 3D.")
            endif 
		    !************************************************************************************
        end subroutine
        !==========================================================================================

		!==========================================================================================
        ! Routine OsmoticPressureModelIdentifier: Routine that identifies the Osmotic Pressure Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine OsmoticPressureModelIdentifier( model, AnalysisSettings, modelID )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModParser
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            character(len=*) , intent(in)    :: model

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(out) :: modelID
            
            ! Internal Variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser)     :: Comp

            !************************************************************************************

            !************************************************************************************
            ! DECODE THE STRING SUPPLIED BY GiD
		    !************************************************************************************
            call Comp%Setup()

            if (Comp%CompareStrings('Ideal_Donnan_Osmotic_Pressure', model).and.(AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration)) then

                modelID = OsmoticPressureModels % IdealDonnanOsmoticPressureModel                

            ! -----------------------------------------------------------------------------------   
            elseif (Comp%CompareStrings('Fixed_Osmotic_Pressure', model).and.(AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration)) then

                modelID = OsmoticPressureModels % FixedOsmoticPressureModel              
            
            ! -----------------------------------------------------------------------------------     
            else

                call Error( "Error: Osmotic Pressure Model not identified: "//trim(model))
                
            endif

		    !************************************************************************************

        end subroutine
        !==========================================================================================

end module



