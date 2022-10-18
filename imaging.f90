module imaging    
    use mathConstants
    use mod_tests
    use m_config

    ! Array shared by entire class - REMEMBER ALWAYS TO ALLOCATE BEFORE USE
    ! image2 is used in the writing of images viewed from along z axis (see write_imageTest)
    integer*8, dimension(:,:,:), allocatable :: image2!, image
    logical :: anglestats

    contains

            ! Uses the entry time and exit time to find the corresponding timepoint for imaging
        subroutine start_end_timepoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, &
             tStep, startTimePoint, endTimePoint)
            implicit none

            integer, intent(in) :: NumberOfTimePoints
            double precision, intent(in) :: entryTime, exitTime, probeStart, probeEnd, tStep
            integer, intent(out) :: startTimePoint, endTimePoint
            
            ! If entry time is less than probe start time, then imaging for that particle starts from the beginning of the process
            if (entryTime .lt. probeStart) then
                startTimePoint = 1
            else
                ! Finds nearest timepoint above entry time, but + 1 is added due to the fact that logically,
                ! the image sequence starts from image 1 rather than image 0
                startTimePoint = ceiling((entryTime - probeStart) / tStep) + 1
            end if

            ! If exit time is greater than end of probe, then end timepoint then particle is imaged 
            ! all the way til the end of the probe time
            if (exitTime .gt. probeEnd) then
                endTimePoint = NumberOfTimePoints
            else
                ! Similarly to entry time, exit timepoint is found as the nearest timepoint above exit time + 1 for the same reasons
                endTimePoint = floor((exitTime - probeStart) / tStep) + 1
            end if
        end subroutine start_end_timepoints

        ! Finds the position a particle is in at any given timepoint, and finds its corresponding pixel position then writes it
        ! to the image array, adding intensity to that pixel region
        subroutine position_in_probe(image, startTimePoint, &
             endTimePoint, xPx, zPx, t0, probeStart, tStep, &
             particleSpeed, pxMmRatio, particleVector, particleStartPos, &
             sheetDimensions, testMods, scatterIntensity, fLifeTime, &
             captureGateOpen, captureGateClose, surface_z)

            implicit none

            double precision, intent(inout), dimension(:,:,:) :: image
            integer, intent(in) :: startTimePoint, endTimePoint, xPx, zPx, surface_z
            double precision, intent(in) :: probeStart, tStep, particleSpeed, pxMmRatio, t0, scatterIntensity, fLifeTime, &
             captureGateOpen, captureGateClose
            double precision, dimension(3), intent(in) :: particleVector, particleStartPos, sheetDimensions
            integer :: t, posInProbexPx, posinProbeyPx, posInProbezPx, yPx
            double precision :: currentTime, emissionTime
            double precision, dimension(3) :: posInProbe
            logical, intent(in) :: testMods
            logical :: zImage

            ! testing purpose: should be left as false for normal operation, other inputs would be required in normal 
            ! operation to enable this properly.
            zImage = .false.

            ! for sake of creating image viewed along z axis, y pixels are the same as z but can be changed if needed
            yPx = zPx

            ! Loops from entry timepont to exit timepoint to avoid wasting cycles when particle is not within sheet
            do t = startTimePoint, endTimePoint
                ! currentTime refers to the time it has taken the particle to travel from its starting point
                ! to the point in space at the given timepoint
                currentTime = probeStart + (t-1)*tStep - t0

                ! Finds the time the molecule will emit its photon after absorbing probe light
                call fluorescence_time(fLifeTime, emissionTime)

                ! Finds position in space that molcule will emit its photon
                currentTime = currentTime + emissionTime

                ! Real space position for particle
                posInProbe(1) = particleStartPos(1) + (particleVector(1)*particleSpeed*currentTime)
                posInProbe(2) = particleStartPos(2) + (particleVector(2)*particleSpeed*currentTime)
                posInProbe(3) = particleStartPos(3) + (particleVector(3)*particleSpeed*currentTime)

                ! Relative pixel position for particle
                ! Note: the subtraction at the end of each statement alters the position of the particle within the image array.
                ! Altering the x-postion by half the width of the image centres the beam
                ! Similarly, the z position can be altered however a factor of 1.3 was found to centre the sheet within the middle
                ! of the image quite well
                posInProbexPx = (ceiling(posInProbe(1)/pxMmRatio) + floor(real(xPx/2)))
                posInProbeyPx = (ceiling(posInProbe(2)/pxMmRatio) + floor(real(yPx/2)))
                posInProbezPx = abs(ceiling(posInProbe(3)/pxMmRatio) - (real(surface_z)))
                !TODO put calculation of this factor earlier somewhere

                ! Only writes to array if particle is within bounds of the image
                if ((posInProbexPx .lt. xPx) .and. (posInProbexPx .gt. 0) .and. (posInProbe(3) .ge. 0)) then
                    ! Mimics gating process. Emission only detected if it occurs between gate open and gate close
                    if ((emissionTime .gt. captureGateOpen) .and. (emissionTime .lt. captureGateClose)) then            
                        if (particleVector(3) .gt. 0) then
                            ! If particle scatters from surface (vector component in z direction > 0) use relative scattering intensity
                            image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + scatterIntensity
                        else
                            ! else just use unit intensity
                            image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + 1D0
                        end if
                    end if
                    
                    ! for testing purposes to view an image along the z axis
                    if (zImage) then                  
                        image2(posInProbeyPx,posInProbexPx,t) = image2(posInProbeyPx,posInProbexPx,t) + 1D0
                    end if
                end if
            end do

        end subroutine position_in_probe

        subroutine directory_setup(path, date_time, linux, full_path, raw_path, blur_path, if_path, input_string, parent_path)
            implicit none

            character(200), intent(in) :: path
            character(200) :: string1, string2
            character(:), allocatable :: trim_path, delim, minus_extension
            character(:), allocatable, intent(out) :: full_path, raw_path, blur_path, if_path, parent_path
            character(:), allocatable, intent(in) :: input_string
            character(17), intent(in) :: date_time
            integer :: length, i
            logical :: linux

            trim_path = trim(path)
            length = len(trim_path)
            delim = "."
            call split_string(input_string, string1, string2, delim)
            delim = "/"
            minus_extension = string1
            call split_string(minus_extension, string1, string1, delim)

            ! Uses the lengtt of the string to find last character in string
            ! If last character is a "/" then subdirectory mdkir command does not need a new one
            ! in the run number directory
            if (trim_path(length:length) .eq. "/") then
                full_path = trim_path//trim(date_time)
            else
                full_path = trim_path//'/'//trim(date_time)
            end if

            parent_path = full_path//"_"//trim(string1)
            raw_path = full_path//"_"//trim(string1)//'/Raw Images'
            blur_path = full_path//"_"//trim(string1)//'/Blurred Images'
            if_path = full_path//"_"//trim(string1)//'/IF Adjusted Images'

            if (linux .eqv. .TRUE.) then

                call execute_command_line('mkdir -p "'//raw_path//'"')
                call execute_command_line('mkdir -p "'//blur_path//'"')
                call execute_command_line('mkdir -p "'//if_path//'"')
            else
                call execute_command_line('mkdir "'//full_path//'/Raw Images'//'"')
                call execute_command_line('mkdir "'//full_path//'/Blurred Images'//'"')
                call execute_command_line('mkdir "'//full_path//'/IF Adjusted Images'//'"')
            end if

        end subroutine directory_setup

        ! Writes out image array into a sequence of images
        subroutine write_image(image, xPx, zPx, startDelay, stopDelay, tstep, NumberOfTimePoints, date_time, raw_path, blur_path, if_path)
            implicit none

            double precision, intent(inout), dimension(:,:,:,:) :: image
            double precision :: startDelay, stopDelay, tstep
            integer :: t, i, j, k, xPx, zPx, NumberOfTimePoints, start_int, stop_int, tstep_int
            character(200) :: fileName, runID
            character(:), allocatable, intent(in) :: raw_path, blur_path, if_path
            character(17), intent(in) :: date_time
            character(3) :: imageNumber

            start_int = nint(startDelay*1E6)
            stop_int = nint(stopDelay*1E6)
            tstep_int = nint(tstep*1E6)

            print "(a)", adjustr(trim("Writing image files"))

            do k = 1, 2     
                do t = 1, NumberOfTimePoints
                    write(imageNumber, '(I0.3)') ((t*tstep_int)-(1*tstep_int)+start_int)
                    if (k == 1) then                
                        fileName = trim(raw_path)//"/image_"//imageNumber//".txt"
                    else if (k == 2) then
                        fileName = blur_path//"/image_"//imageNumber//".txt"
                    else
                        fileName = if_path//"/image_"//imageNumber//".txt"
                    end if

                    open(unit=20+t,file=filename)

                    do i = 1, zPx
                        do j = 1, xPx
                            write(20+t,'(ES12.5,a)',advance='no') image(i,j,t,k)," "
                        end do

                        write(20+t,*)
                    end do
                    close(unit=20+t)
                end do
            end do
        end subroutine write_image

        subroutine convim(imin,nx,ny,gaussdev,imout)
            !Convolutes input image imin with a gaussian of st. dev. gaussdev (in pixels), to produce imout.
                implicit none
                double precision, dimension(:,:), intent(in) :: imin(nx,ny)
                double precision, dimension(:,:), intent(out) :: imout(nx,ny)
                integer, intent(in) :: nx,ny
                double precision, intent(in) :: gaussdev
                
                double precision, dimension(:), allocatable :: gauss
                double precision, dimension(:,:) :: immid(nx,ny)
                
                double precision, parameter :: sqrt2 = 1.414213562d0
                
                integer :: gsize
                double precision :: cent
                integer :: icent
                integer :: j,k,l
                double precision :: lim1,lim2
                double precision :: dblej
                        
            !Calculate size of gaussian array required. This can just be a 1D gaussian since 2D Gaussian is separable
            ! and hence convolution can be applied stepwise for each dimension.
                if ((gaussdev .eq. 0.0d0) .and. (sum(imin) .gt. 0)) then
                        imout = imin
                else
                        gsize = floor(6.0d0*gaussdev)+1
                        if (mod(gsize,2) .eq. 0) gsize = gsize + 1 !Gaussian array will always have an odd-numbered size
                        allocate(gauss(gsize))
            !Calculate gaussian:
                        cent  = (dble(gsize)+1.0d0)*0.50d0
                        icent = int(cent)
                        
                        do j = 1,gsize
                            if (j .eq. icent) then
                                    lim1 = 0.5d0/(sqrt2*gaussdev)
                                    gauss(j) = erf(lim1)
                            else
                                    dblej = dble(j)
                                    lim1 = (abs(dblej - cent) +0.50d0)/(sqrt2*gaussdev)
                                    lim2 = lim1 - 0.5d0/(sqrt2*gaussdev)
                                    gauss(j) = erf(lim1)-erf(lim2)
                            end if
                        end do
            
            !Convolute along x.
                        do k = 1, ny
                            do j = 1, nx
                                    immid(j,k) = 0.0d0
                                    do l = 1, gsize
                                        if (((j+l-icent) .gt. 0) .and. ((j+l-icent) .le. nx)) then
                                                immid(j,k) = immid(j,k) + imin(j+l-icent,k)*gauss(l)
                                        end if
                                    end do
                            end do
                        end do
            !Convolute along y.
                        do j = 1, nx
                            do k = 1, ny
                                    imout(j,k) = 0.0d0
                                    do l = 1, gsize
                                        if (((k+l-icent) .gt. 0) .and. ((k+l-icent) .le. ny)) then
                                                imout(j,k) = imout(j,k) + immid(j,k+l-icent)*gauss(l)
                                        end if
                                    end do
                            end do
                        end do
            !Normalize output image to have same intensity as input image:
                        if ((sum(imin) .gt. 0) .or. (sum(imin) .lt. 0)) then
                            imout = imout*sum(imin)/sum(imout)
                        else
                            imout = 0  
                        end if     
                end if
        end subroutine convim

        subroutine fluorescence_time(fLifeTime, emissionTime)
            implicit none

            double precision :: rand
            double precision, intent(out) :: emissionTime
            double precision, intent(in) :: fLifeTime

            call random_number(rand)

            emissionTime = -fLifeTime*log(rand)

        end subroutine fluorescence_time

        subroutine date_time_string(string)
            character(17), intent(out) :: string
            character(4) :: year, month, day, hour, min, sec
            integer, dimension(8) :: values

            ! Values array holds: (1) year, (2) month, (3) day, (4) UTC difference in mins, (5) hour, (6) minute, (7) seconds, (8) milliseconds
            call date_and_time(VALUES=values)

            write(year, '(i4)') values(1)
            write(month, '(i2.2)') values(2)
            write(day, '(i2.2)') values(3)
            write(hour, '(i2.2)') values(5)
            write(min, '(i2.2)') values(6)
            write(sec, '(i2.2)') values(7)

            string = trim(year)//"-"//trim(month)//"-"//trim(day)//"_"//trim(hour)//trim(min)//trim(sec)
        end subroutine date_time_string

        SUBROUTINE split_string(instring, string1, string2, delim)
            CHARACTER(:), allocatable :: instring,delim
            CHARACTER(200),INTENT(OUT):: string1,string2
            INTEGER :: index
        
            instring = TRIM(instring)
        
            index = SCAN(instring,delim)
            string1 = instring(1:index-1)
            string2 = instring(index+1:)
        
        END SUBROUTINE split_string

           
        
end module imaging