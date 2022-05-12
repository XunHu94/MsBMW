        !COMPILER-GENERATED INTERFACE MODULE: Thu May 12 23:40:16 2022
        MODULE CONVOLUTION__genmod
          INTERFACE 
            SUBROUTINE CONVOLUTION(NUM1,XTMP,YTMP,ZTMP,DEN,SPEED,RARRY)
              USE VARMOD
              INTEGER(KIND=4) :: NUM1
              INTEGER(KIND=4) :: XTMP
              INTEGER(KIND=4) :: YTMP
              INTEGER(KIND=4) :: ZTMP
              REAL(KIND=4) :: DEN(80)
              REAL(KIND=4) :: SPEED(80)
              REAL(KIND=4) :: RARRY(NUM1)
            END SUBROUTINE CONVOLUTION
          END INTERFACE 
        END MODULE CONVOLUTION__genmod
