        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 10 17:34:24 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATERIALCONSTRUCTOR__genmod
          INTERFACE 
            SUBROUTINE MATERIALCONSTRUCTOR(ELEMENT,ELEMENTLIST,         &
     &GLOBALNODESLIST,MATERIAL,ANALYSISSETTINGS)
              USE MODANALYSIS
              USE MODNODES
              USE MODELEMENTLIBRARY
              CLASS (CLASSELEMENT) ,POINTER :: ELEMENT
              TYPE (CLASSELEMENTSWRAPPER) ,POINTER :: ELEMENTLIST(:)
              TYPE (CLASSNODES) ,POINTER :: GLOBALNODESLIST(:)
              CLASS (CLASSCONSTITUTIVEMODELWRAPPER) ,POINTER :: MATERIAL
              TYPE (CLASSANALYSIS) :: ANALYSISSETTINGS
            END SUBROUTINE MATERIALCONSTRUCTOR
          END INTERFACE 
        END MODULE MATERIALCONSTRUCTOR__genmod
