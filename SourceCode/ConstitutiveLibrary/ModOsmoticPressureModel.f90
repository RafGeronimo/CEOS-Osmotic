!##################################################################################################
! This module has the attributes and methods to all Osmotic Gradient Models.
!--------------------------------------------------------------------------------------------------
! Date: 2022/08
!
! Authors:  Rafael Geronimo
!!------------------------------------------------------------------------------------------------
! Modifications:
!
!##################################################################################################
module ModOsmoticPressureModel
        
        use ModContinuumMechanics
        use ModMathRoutines
        use ModStatus
        use ModLoadHistoryData

        !real(8), parameter :: constR = 8.3145d-6!8.3145d0         !Universal gas constant [J/mol*K]
        
        type ClassOsmoticPressureModel
            real(8) :: DeltaPi                          !Osmotic Pressure Gradient
            real(8) :: J                                !Volumetric jacobian J = det(F)
            real(8) :: Lambda                           !Osmotic Pressure constant of Stifness Matix Kuu

            contains
        
                ! Class Methods
                !------------------------------------------------------------------------------------
                procedure :: GetLambda                                   => GetLambdaBase
                procedure :: GetFixedChargeDensityAndOsmoticGradient     => GetFixedChargeDensityAndOsmoticGradientBase
                procedure :: CopyOsmoticPressureProperties               => CopyOsmoticPressurePropertiesBase
                procedure :: ReadOsmoticPressureParameters               => ReadOsmoticPressureParametersBase
            
            
        end type
            
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        type ClassOsmoticPressureModelWrapper

            class(ClassOsmoticPressureModel) , pointer , dimension(:) :: Mat => null()
            integer                                                   :: MaterialID = -999
            integer                                                   :: ModelEnumerator = -999

        end type
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    
            
        
    contains 
    
        !==========================================================================================
        ! Dummy Procedures: To be used by the superclasses
        !==========================================================================================
    
        subroutine GetFixedChargeDensityAndOsmoticGradientBase(this, LC, ST)
            class(ClassOsmoticPressureModel) :: this
            integer :: LC, ST
            stop "Error: OsmoticGradientConstructor: no fixed charge density and osmotic pressure obtained"
        end subroutine
        
        !==========================================================================================
        
        subroutine CopyOsmoticPressurePropertiesBase(this,Reference)
            class(ClassOsmoticPressureModel) :: this , Reference
            stop "Error: CopyOsmoticPressureProperties"
        end subroutine
        
        !==========================================================================================
        
        subroutine ReadOsmoticPressureParametersBase(this,DataFile)
            use ModParser
            class(ClassOsmoticPressureModel)::this
            type(ClassParser)::DataFile
            !integer , intent(in):: FileNum
            stop "Error: ReadMaterialParameters (OsmoticPressure)"
        end subroutine
        
        !==========================================================================================
        
        !==========================================================================================
        
        subroutine GetLambdaBase(this, LC, ST)
            class(ClassOsmoticPressureModel):: this
            integer :: LC, ST
            stop "Error: Lambda (constant of Osmotic Pressure Gradient derivative)"
        end subroutine
        
        !==========================================================================================
        
end module
