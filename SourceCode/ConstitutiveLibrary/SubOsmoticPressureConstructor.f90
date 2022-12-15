!##################################################################################################
! This routine constructs the element.
!--------------------------------------------------------------------------------------------------
! Date: 2022
!
! Authors:  Rafael Geronimo
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine OsmoticPressureConstructor( ElementBiphasic, ElementList, GlobalNodesList, Permeability, OsmoticPressure, AnalysisSettings )

	!************************************************************************************
	! DECLARATIONS OF VARIABLES
	!************************************************************************************
	! Modules and implicit declarations
	! -----------------------------------------------------------------------------------
	use ModAnalysis
    use ModElementBiphasic
	use ModElementLibrary
	use ModNodes
    use ModPermeabilityModelLibrary
	use ModOsmoticPressureModelLibrary
    
	implicit none

	! Input variables
	! -----------------------------------------------------------------------------------
    class(ClassElementBiphasic), pointer					 :: ElementBiphasic
    type (ClassElementsWrapper) , pointer , dimension(:)	 :: ElementList
	type(ClassAnalysis)										 :: AnalysisSettings
	type(ClassNodes) , dimension(:) , pointer                :: GlobalNodesList
	class(ClassPermeabilityModelWrapper)  , pointer          :: Permeability
    class(ClassOsmoticPressureModelWrapper)  , pointer       :: OsmoticPressure

	! Internal variables
	! -----------------------------------------------------------------------------------
	integer :: i, nNodes, nGP, gp

	!************************************************************************************

	!************************************************************************************
	! CONSTRUCT OF MATERIALS
	!************************************************************************************

    call ElementBiphasic%AllocateGaussPoints_Osmotic(nGP)

    if ( ( nGP <= 0 ) .or. (nGP > 1000) ) then
        call Error ("Error: Number of the Gauss Points <=0 or >1000")
    endif

    
    ! Allocate the contitutive model for all element's Osmotic Gauss point
    call AllocateOsmoticPressureModel( OsmoticPressure%ModelEnumerator , AnalysisSettings , nGP , ElementBiphasic%GaussPoints_Osmotic)

    ! Copy Osmotic Pressure properties from reference material (read in the settings file) to
    ! Osmotic Gauss points
    ! And Construct the Osmotic Pressure Model
    ! -----------------------------------------------------------------------------------
	do gp = 1,nGP
        call ElementBiphasic%GaussPoints_Osmotic(gp)%CopyOsmoticPressureProperties(OsmoticPressure%Mat(1))
        !call ElementBiphasic%GaussPoints_Osmotic(gp)%OsmoticPressureModelDestructor()
        !call ElementBiphasic%GaussPoints_Osmotic(gp)%OsmoticPressureModelConstructor(AnalysisSettings)
        
    enddo
    
	!************************************************************************************
end subroutine
