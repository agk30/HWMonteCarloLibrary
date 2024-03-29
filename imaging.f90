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
             captureGateOpen, captureGateClose)

            implicit none

            double precision, intent(inout), dimension(:,:,:) :: image
            integer, intent(in) :: startTimePoint, endTimePoint, xPx, zPx
            double precision, intent(in) :: probeStart, tStep, particleSpeed, pxMmRatio, t0, scatterIntensity, fLifeTime, &
             captureGateOpen, captureGateClose
            double precision, dimension(3), intent(in) :: particleVector, particleStartPos, sheetDimensions
            integer :: t, posInProbexPx, posinProbeyPx, posInProbezPx, sheetCentrePx, yPx, i
            double precision :: currentTime, angle, emissionTime
            double precision, dimension(3) :: posInProbe
            logical, intent(in) :: testMods
            logical :: zImage

            ! testing purpose: should be left as false for normal operation, other inputs would be required in normal 
            ! operation to enable this properly.
            zImage = .false.

            ! for sake of creating image viewed along z axis, y pixels are the same as z but can be changed if needed
            yPx = zPx

            ! Loops from entry timepont to exit timepoint to avoic wasting cycles when particle is not within sheet
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
                posInProbezPx = abs(ceiling(posInProbe(3)/pxMmRatio) - (real(294)))
                !TODO put calculation of this factor earlier somewhere

                ! Only writes to array if particle is within bounds of the image
                if ((posInProbexPx .lt. xPx) .and. (posInProbexPx .gt. 0) .and. (posInProbe(3) .ge. 0)) then
                    ! Mimics gating process. Emission only detected if it occurs between gate open and gate close
                    if ((emissionTime .gt. captureGateOpen) .and. (emissionTime .lt. captureGateClose)) then            
                        if (particleVector(3) .gt. 0) then
                            image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + scatterIntensity
                        else
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

        subroutine directory_setup(path, runID, input_param, linux)
            implicit none

            character(200), intent(in) :: path, runID
            character(:), allocatable :: trim_path, trim_runID, full_path
            integer :: length, i
            logical :: linux
            type(CFG_t) :: input_param

            trim_runID = trim(runID)

            print '(a)', "Creating output image directories..."

            trim_path = trim(path)
            length = len(trim_path)

            ! Uses the lenght of the string to find last character in string
            ! If last character is a "/" then subdirectory mdkir command does not need a new one
            ! in the run number directory
            if (trim_path(length:length) .eq. "/") then
                full_path = trim_path//'Run '//trim_runID
            else
                full_path = trim_path//'/Run '//trim_runID
            end if

            if (linux .eqv. .TRUE.) then
                call execute_command_line('mkdir -p "'//full_path//'/Raw Images'//'"')
                call execute_command_line('mkdir -p "'//full_path//'/Blurred Images'//'"')
                call execute_command_line('mkdir -p "'//full_path//'/IF Adjusted Images'//'"')
            else
                call execute_command_line('mkdir "'//full_path//'/Raw Images'//'"')
                call execute_command_line('mkdir "'//full_path//'/Blurred Images'//'"')
                call execute_command_line('mkdir "'//full_path//'/IF Adjusted Images'//'"')
            end if

            call CFG_write(input_param, full_path//"/input_values.cfg", .FALSE., .FALSE.)
        end subroutine directory_setup

        ! Writes out image array into a sequence of images
        subroutine write_image(image, xPx, zPx, startDelay, stopDelay, tstep, NumberOfTimePoints, runNumber, imagePath)
            implicit none

            double precision, intent(inout), dimension(:,:,:,:) :: image
            double precision :: startDelay, stopDelay, tstep
            integer :: t, i, j, k, xPx, zPx, runNumber, NumberOfTimePoints, start_int, stop_int, tstep_int
            character(200) :: fileName, runID
            character(200), intent(in) :: imagePath
            character(3) :: imageNumber

            start_int = int(startDelay*1E6)
            stop_int = int(stopDelay*1E6)
            tstep_int = int(tstep*1E6)

            print "(a)", 'Entering write'

            write(runID, '(i0)') runNumber

            do k = 1, 3     
                do t = 1, NumberOfTimePoints, tstep_int
                    write(imageNumber, '(I0.3)') ((t*tstep_int)-(1*tstep_int)+start_int)
                    if (k == 1) then                
                        fileName = trim(imagePath)//"Run "//trim(runID)//"/Raw Images/Image_"//imageNumber//".txt"
                    else if (k == 2) then
                        fileName = trim(imagePath)//"Run "//trim(runID)//"/Blurred Images/Image_"//imageNumber//".txt"
                    else
                        fileName = trim(imagePath)//"Run "//trim(runID)//"/IF Adjusted Images/Image_"//imageNumber//".txt"
                    end if

                    open(unit=20+t,file=filename)

                    do i = 1, zPx
                        do j = 1, xPx
                            write(20+t,'(ES12.5,a)',advance='no') image(i,j,t,k)," "
                        end do

                        write(20+t,*)
                    end do
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
        
end module imaging