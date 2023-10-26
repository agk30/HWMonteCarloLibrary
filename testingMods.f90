module mod_tests
    use mathConstants
    use directions
    use speeds

    integer, dimension(-90:90) :: angleDist

    contains

        ! produces a file of binned angle distributions. 
        subroutine angleDistribution(outgoing)
            implicit none

            double precision, dimension(3), intent(in) :: outgoing
            double precision, dimension(3) :: ingoing
            integer :: i
            double precision :: angle

            ingoing(1) = 0 ; ingoing(2) = 0 ; ingoing(3) = 1

            angle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360/(2*pi))
            
            if (outgoing(1) .lt. 0) then

                angle = -angle
                
            end if

            do i = -90, 90

                if ((angle .gt. i) .and. (angle .lt. i+1)) then

                    angleDist(i) = angleDist(i) + 1

                    exit

                end if

            end do

        end subroutine angleDistribution

        subroutine writeAngleDistribution
            implicit none

            integer :: i

            open(unit=500,file='angledist.txt')

            do i = -90, 90

                write(500,*) angleDist(i), (i)

            end do

        end subroutine writeAngleDistribution

        subroutine angle_speed_distribution(vector, speed, distributionBin)
            implicit none

            double precision, intent(in) :: vector(3), speed
            double precision :: normal(3), normalxz(2), normalyz(2), vectorxz(2), vectoryz(2)
            double precision :: anglexz, angleyz, maxSpeed
            integer :: distributionBin(:,:)
            integer :: i, j, angleBinSize, speedBinSize, intSpeed

            angleBinSize = 2

            speedBinSize = 10

            maxSpeed = 2500

            normal(1) = 0 ; normal(2) = 0 ; normal(3) = 1.0
            normalxz(1) = 0 ; normalxz(2) = 1
            normalyz(1) = 0 ; normalyz(2) = 1

            vectorxz(1) = vector(1) ; vectorxz(2) = vector(3)
            vectoryz(1) = vector(2) ; vectoryz(2) = vector(3)

            if ((speed .lt. maxSpeed) .and. (vector(1) .gt. 0)) then

                anglexz = acos((dot_product(vectorxz,normalxz)/(norm2(vectorxz)*norm2(normalxz))))
                angleyz = acos((dot_product(vectoryz,normalyz)/(norm2(vectoryz)*norm2(normalyz))))
                anglexz = ceiling((anglexz*(360.0/(pi*2)))/angleBinSize)
                angleyz = ceiling((angleyz*(360.0/(pi*2)))/angleBinSize)

                if ((anglexz .ge. 1) .and. (anglexz .le. 45)) then

                    intSpeed = ceiling(speed/speedBinSize)

                    !print *, speed

                    distributionBin(int(anglexz),intSpeed) = distributionBin(int(anglexz),intSpeed) + 1
                end if

            end if

        end subroutine angle_speed_distribution

        subroutine write_angle_speed(bin)

            integer, dimension(:,:) :: bin
            integer :: i, j, angleBinSize, speedBinSize, maxSpeed
            character :: c

            angleBinSize = 2

            speedBinSize = 10

            maxSpeed = 2500

            open(600,file='angle_speed_distribution.csv')

            write(600,'(a)',advance='no') "Speed m/s,"

            do j = 1, (90/angleBinSize)
                do i = 1, (maxSpeed/speedBinSize)
                    write(600,'(I4,a)',advance='no') i*speedBinSize,","
                    write(600,'(I4,a)',advance='no') j*angleBinSize,","
                    write(600,'(I5,a)',advance='no') bin(j,i),","
                    write(600,'(a)',advance='no') new_line(c)
                end do
            end do

        end subroutine write_angle_speed

        subroutine point_source (vector, point, finalSpeed, time)
            implicit none

            double precision :: vector(3), point(3), speed, finalSpeed, time, mass, internalRatio, surfaceMass, initialSpeed
            double precision :: deflectionAngle, ingoing(3), rand, theta
            
            speed = 1800D0

            internalRatio = 0
            mass = 17D-3
            surfaceMass = 100D0

            call cosine_distribution(0, vector, theta)

            ingoing(1) = 0.707106781
            ingoing(2) = 0
            ingoing(3) = -0.707106781

            !call deflection_angle(ingoing, vector, deflectionAngle)

            !call soft_sphere_speed(mass, internalRatio, surfaceMass, speed, deflectionAngle, finalSpeed)

            finalSpeed = 100

            point = 0

            call random_number(rand)

            time = 1.2D-4 + ((10D-6)*rand)

        end subroutine point_source

        integer function find_bin_index(value, bin_range, bin_size, num_bins)

            implicit none
            double precision, dimension(2) :: bin_range
            double precision :: bin_size, value
            integer :: num_bins, i, index

            do i = 1, num_bins
                if (value .lt. (bin_range(1) + (i-1)*bin_size)) then
                    index = i
                    if (i .eq. 1) then
                        !print *, particleSpeed
                    end if
                    EXIT
                end if
            end do

            find_bin_index = index

        end function find_bin_index

        subroutine two_vector_angle(vector1, vector2, angle)
            implicit none

            double precision, dimension(3) :: vector1, vector2
            double precision :: angle

            angle = acosd(dot_product(vector1,vector2)/(norm2(vector1)*norm2(vector2)))

        end subroutine two_vector_angle

        subroutine cosine_distribution_fixed_phi(phi_in, scatteredDirection, theta)

            implicit none

            double precision :: rand1, phi
            double precision, intent(in) :: phi_in
            double precision, intent(out) :: theta
            double precision, dimension(3), intent(out) :: scatteredDirection
        
            phi = phi_in*D_to_R

            call random_number(rand1)

            !theta = asin(SQRT(rand1))
            theta = rand1*(pi/2)

            scatteredDirection(1) = sin(theta)*cos(phi)
            scatteredDirection(2) = sin(theta)*sin(phi)
            scatteredDirection(3) = cos(theta)
            
        end subroutine cosine_distribution_fixed_phi

        subroutine cosine_distribution_fixed_phi2(cosinePower, scatteredDirection, phi_in, theta)

            implicit none
        
            double precision :: rand1, rand2, theta, x, phi
            integer, intent(in) :: cosinePower, phi_in
            double precision, dimension(3), intent(out) :: scatteredDirection
        
            call random_number(rand2)

            phi = phi_in*D_to_R
        
            x = rand2**(1.0/dble(cosinePower+1.0))
        
            theta = dacos(x)

            scatteredDirection(1) = sin(theta)*cos(phi)
            scatteredDirection(2) = sin(theta)*sin(phi)
            scatteredDirection(3) = cos(theta)
            
        end subroutine cosine_distribution_fixed_phi2

end module mod_tests
