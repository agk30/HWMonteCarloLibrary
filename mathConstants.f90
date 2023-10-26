module mathConstants

    ! General use mathematical constants for use in all subroutines
    ! The r14 parameter is very important and the program will not compile without it
    ! It defines the kind for the real variables as double precision (14 decimal places, max exponent of 30)
    ! This method is more portable than the double precision variable
    double precision, parameter :: pi = 3.141592653589793D0
    double precision, parameter :: boltzmannConstant =1.38064852D-23
    double precision, parameter :: Degree180 = 180.0
    double precision, parameter :: R_to_D = Degree180/pi
    double precision, parameter :: D_to_R = pi/Degree180

    public :: pi, boltzmannConstant, R_to_D, D_to_R

end module mathConstants