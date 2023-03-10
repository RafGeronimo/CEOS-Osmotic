!##################################################################################################
! This module has the system of equations of monolithic FEM (Biphasic Analysis with Osmotic Pressure)
!--------------------------------------------------------------------------------------------------
! Date: 2022/09
!
! Author:  Rafael Geronimo
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!                               
!##################################################################################################
module ModFEMSystemOfEquationsMonolithicBiphasicOP

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions
    use ModBoundaryConditionsFluid
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic
    use ModOsmoticPressureModel

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassFEMSystemOfEquationsMonolithicBiphasicOP
        
        ! Solid-Fluid variables
        real(8), dimension(:), allocatable                     :: Xbar
        real(8), dimension(:),allocatable                      :: XConverged          ! Global U/P vector
        integer, dimension(:) , pointer                        :: DirichletDOF        
        integer, dimension(:), allocatable                     :: PrescDirichletSparseMapONE
        integer, dimension(:), allocatable                     :: PrescDirichletSparseMapZERO
        
        ! Fluid variables
        real(8),dimension(:), allocatable                      :: Fint_fluid , Fext_fluid  !(Pbar ~~ Ubar)
        integer , dimension(:) , pointer                       :: PresDOF
        
        ! Solid variables
        real(8),dimension(:),allocatable                       :: Fint_solid , Fext_solid
        real(8),dimension(:),allocatable                       :: VSolid              ! Global solid velocity
        integer, dimension(:) , pointer                        :: DispDOF
        integer, dimension(:), allocatable                     :: FixedSupportSparseMapZERO
        integer, dimension(:), allocatable                     :: FixedSupportSparseMapONE
        
        ! FEM Objects
        type (ClassElementsWrapper)  , dimension(:) , pointer  :: ElementList
        type (ClassNodes)            , dimension(:) , pointer  :: GlobalNodesList
        type (ClassAnalysis)                                   :: AnalysisSettings
        class (ClassBoundaryConditions)             , pointer  :: BCSolid
        class (ClassBoundaryConditionsFluid)        , pointer  :: BCFluid
        type (ClassGlobalSparseMatrix)              , pointer  :: Kg
        
        ! Time discretization variables
        real (8)                                               :: Time
        real (8)                                               :: DeltaT

    contains

        procedure :: EvaluateSystem         => EvaluateROP
        procedure :: EvaluateGradientSparse => EvaluateKtOP
        procedure :: PostUpdate             => FEMUpdateMeshOP

    end type

    contains

    !=================================================================================================
    subroutine EvaluateROP(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsMonolithicBiphasicOP) :: this
        real(8),dimension(:) :: X,R
        integer :: NDOFsolid
        integer :: NDOFfluid

            ! X -> Global Incognite of the biphasic analysis (U / P)
        
            NDOFsolid = this%AnalysisSettings%NDOFsolid
            NDOFfluid = this%AnalysisSettings%NDOFfluid
            
            ! Update and internal variables (In case that the consitutive model depends of the pore pressure, the SolveConstitutiveModel must be uptaded)
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X(1:NDOFsolid), this%Status)

            ! Update the deformation gradient and permeability on fluid gauss points
            call SolvePermeabilityModel( this%ElementList , this%AnalysisSettings , X(1:NDOFsolid), this%Status)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif
            
            ! Internal Force solid
            call InternalForceSolid(this%ElementList , this%AnalysisSettings , X((NDOFsolid+1):(NDOFsolid+NDOFfluid)), &
                                    this%Fint_solid, this%Status)
            
            ! Updating the solid velocity with a time integration method
            call GetSolidVelocity (this%XConverged(1:NDOFsolid), X(1:NDOFsolid), this%DeltaT, this%VSolid)
            
            ! Internal Force fluid
            call InternalForceFluid(this%ElementList , this%AnalysisSettings , X((NDOFsolid+1):(NDOFsolid + NDOFfluid)), &
                                    this%VSolid , this%Fint_fluid , this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Residual
            R(1:NDOFsolid)                           = this%Fint_solid - this%Fext_solid
            R((NDOFsolid+1):(NDOFsolid + NDOFfluid)) = this%Fint_fluid - this%Fext_fluid 

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine EvaluateKtOP(this,X,R,G)

        use ModMathRoutines
        use ModSparseMatrixRoutines

        class(ClassFEMSystemOfEquationsMonolithicBiphasicOP)        :: this
        class (ClassGlobalSparseMatrix), pointer :: G
        real(8),dimension(:) :: X , R
        real(8) :: norma
        integer :: NDOFsolid
        integer :: NDOFfluid
        
        ! Difference Centrar Variables - For comparation reasons
        !real(8), dimension(68,68) :: KgDFCFull 
        !real(8), dimension(4624) :: normDFC
        !type (ClassGlobalSparseMatrix), pointer :: KgSparseDFC
        !type (ClassGlobalSparseMatrix), pointer :: KgSparseAnalytic
        !type(SparseMatrix) :: KgSparseDFC_Raw
        !integer :: i ,j, k
        
        !allocate(KgSparseDFC)
        
        ! X -> Global Incognite of the biphasic analysis (U / P)      
        
        NDOFsolid = this%AnalysisSettings%NDOFsolid
        NDOFfluid = this%AnalysisSettings%NDOFfluid

        call TangentStiffnessMatrixMonolithic( this%AnalysisSettings , this%ElementList , this%DeltaT, this%VSolid, &
                                                X((NDOFsolid+1):(NDOFsolid + NDOFfluid)), this%Kg)
       
        ! ---------------------------------------------------------------------
        !! Difference Centrar Comparation - Use only 1 element
        !KgSparseAnalytic => this%Kg
        !
        !call TangentStiffnessDFC(this, X, R, KgDFCFull)
        !
        !k=0
        !do i= 1, 68
        !    do j = 1, 68
        !        k=k+1
        !        write(*,*) i, "," , j, ",", KgDFCFull(i,j), ",", KgSparseAnalytic%Val(k)
        !    end do
        !end do
        ! ---------------------------------------------------------------------
        
        ! The dirichelet BC (Fluid -> pressure) are being applied in the system Kx=R and not in Kx = -R
        R = -R
        !****************************************************************************************
        call this%BCSolid%ApplyBoundaryConditions(  this%Kg , R , this%DirichletDOF, this%Xbar , X, this%PrescDirichletSparseMapZERO, &
                                                    this%PrescDirichletSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        !****************************************************************************************
        R = -R

        G => this%Kg
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine FEMUpdateMeshOP(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsMonolithicBiphasicOP) :: this
        real(8),dimension(:)::X
        integer :: NDOFsolid

        ! X -> Global Incognite of the biphasic analysis (U / P)
        NDOFsolid = this%AnalysisSettings%NDOFsolid

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings, X(1:NDOFsolid))
        endif

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    !subroutine  TangentStiffnessDFCOP(this, X, R, KgDFCFull)
    !    
    !    ! Input variables
    !    real(8),dimension(:), intent(in) :: X , R
    !    class(ClassFEMSystemOfEquationsMonolithicBiphasicOP) :: this
    !
    !    ! Output variables
    !    real(8), dimension(:,:), intent(inout) :: KgDFCFull
    !    
    !    ! Internal variables
    !    real(8),allocatable, dimension(:) :: Xpert , Rpert_frente,Rpert_tras
    !    integer :: i, j, k, FileID_CompareKgDFC
    !    real(8) :: h
    !    real(8), dimension(68,68) :: KgDFCFull2
    !    
    !    allocate(Rpert_frente(size(R,1)))
    !    allocate(Rpert_tras(size(R,1)))
    !    allocate(Xpert(size(X,1)))
    !    
    !    FileID_CompareKgDFC = 13
    !    open (FileID_CompareKgDFC,file='CompareKgDFC.result',status='unknown')
    !            
    !    Xpert = X
    !    
    !    do i = 1, size(X,1)
    !        h = 1.0e-8
    !        do j = 1, size(X,1)
    !            if (j.eq.61) then
    !                h = h*1.0d0
    !            endif
    !            
    !            Xpert(j) = X(j) + h
    !            
    !            call FEMUpdateMesh(this,Xpert)
    !            
    !            call EvaluateR(this,Xpert,Rpert_frente)
    !            
    !            Xpert(j) = X(j) - h
    !
    !            call FEMUpdateMesh(this,Xpert)
    !            
    !            call EvaluateR(this,Xpert,Rpert_tras)
    !            
    !            KgDFCFull(i,j) = (Rpert_frente(i) - Rpert_tras(i))/(2*h)
    !                            
    !            call FEMUpdateMesh(this,X)
    !            
    !        end do
    !    end do
    !    
    !    !do i = 1, size(X,1)
    !    !    do j = 1, size(X,1)
    !    !        write(FileID_CompareKgDFC,*) 'K', j,i, '=', KgDFCFull(j,i)        
    !    !    end do
    !    !end do     
    !    
    !end subroutine
    !!=================================================================================================
    
end module

