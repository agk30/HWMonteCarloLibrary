!TODO fix this damn thing. Bad hardcoded variables - get it done pls

! Important to note that the image coming into this subroutine should be zeroed
! already or weird stuff happens
module sgconv
    contains
        ! builds the kernel used in 3D SG smoothing subroutine 
        ! returns 1D array in wraparound order from the 2D padded input array   
        subroutine construct_kernel(i, j, ksize, nd, columnKernel, padInput)
            implicit none
            
            integer :: k
            integer, intent(in) :: ksize, nd, i, j
            double precision, intent(inout), dimension(:,:) :: padInput
            double precision, intent(inout), dimension(:) :: columnKernel
            
            do k = 1, ksize
                columnKernel(((k-1)*ksize)+1 : k*ksize) = padInput((i+nd+nd-k+1), j:(j+ksize-1))
            end do

        end subroutine construct_kernel

        subroutine sg_array(xPx, zPx, ksize, matrixPath, input, output)
            implicit none
        
            integer :: xPx, zPx, i, j, ksize, nd, nk
            double precision, intent(in), dimension(:,:) :: input
            double precision, allocatable, dimension(:,:) :: diffinput
            double precision, dimension(:,:), intent (out) :: output
            double precision, allocatable, dimension(:,:) :: padInput
            double precision, allocatable, dimension(:) :: sgmatrix
            double precision, allocatable, dimension(:,:) :: kernel
            double precision, allocatable, dimension(:) :: columnKernel
            double precision :: dotprod
            character(100) :: matrixPath
            
            nd = (ksize-1)/2
            nk = ((ksize**2) - 1)/2
        
            allocate(diffinput(xPx,zPx))
            allocate(padInput(xPx+ksize-1,zPx+ksize-1))
            allocate(sgmatrix(ksize**2))
            allocate(kernel(-nd:nd,-nd:nd))
            allocate(columnKernel(ksize**2))

            ! open(12,file='SG Matrices/CC_027x027_00'//char(polyOrder)//'x00'//char(polyOrder)//'.dat')
            open(1000,file=matrixPath)
        
            read(1000,*) sgmatrix(1:(((ksize**2)-1)/2)+1)
        
            ! input SG matrices only contain half + 1 of the required array, symmetric about the last entry. this loop builds
            ! the rest of the array
            do i = 1, nk     
                sgmatrix(i + nk + 1) = sgmatrix(nk+1-i)    
            end do
    
            padInput = 0D0
        
            ! builds a padded array. the image array requires padding outside of its bounds with the addition of nd fields on 
            ! either side of the image and also above and below it
            do i = 1, 420     
                do j = 1, 420
        
                    padInput(i+nd,j+nd) = input(i,j)
        
                end do       
            end do
      
            do i = 1, 420        
                do j = 1, 420      
                    call construct_kernel(i, j, ksize, nd, columnKernel, padInput)
       
                    ! convolutes image with SG cooefficients
                    dotprod = dot_product(columnKernel, sgmatrix)
        
                    output(i,j) = dotprod

                    !ensures no negative numbers in output (not really needed, should be fixed)
                    if (output(i,j) .lt. 0) then
                        !output(i,j) = 0
                    end if  
                end do  
            end do

        end subroutine sg_array

end module sgconv



