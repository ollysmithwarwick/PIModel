module PIInit
! **************************************************************** !
! Module containing tools to get initial conditions from a file or !
! initialise with a type of initial condition                      !
! **************************************************************** !
  use parameters
  use pi_io

  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

contains
  subroutine getVec(size, initInfoFile, vec)
    ! Converts filie into PI state as vector
    ! initMode
    !     0 From file
    !     1 Gaussian at center
    !     2 Random white noise
    !     3 Initialize with Fourier values

    integer, intent(in) :: size
    character(len = 400), intent(in) :: initInfoFile
    character(len = 400) :: initFile
    type(C_PTR) :: planf, planb
    complex(dp), dimension(:), allocatable :: in, out
    real(dp) :: initAmp, param1, phaseP=0d0,globalN=0d0, globalE=0d0
    integer  :: lineN, LBound=-1, UBound=-1, sizef=-1
    integer :: seed, initMode
    integer, dimension(1) :: seed_
    real(dp) :: sigma, fphase
    real(dp) :: translateL=0.0
    real(dp) :: S_, shift_, phase_
    real(dp),    dimension(6*size), intent(out) :: vec
    real(dp),    dimension(6*size) :: vecTmp
    complex(dp), dimension(size,2) :: fvec
    real(dp),    dimension(size, 4)             :: matR, matI
    complex(dp), dimension(size, 4)             :: mat
    integer :: j, j1, j2
    
    namelist /initInfo/ initMode, initAmp, seed, lineN, initFile, param1, LBound, UBound, translateL, phaseP, globalN, globalE, sizef
    initFile = ""
    initAmp = 1.0_dp
    param1  = 1.0_dp
    lineN   = 0
    seed    = 0
    vec     = 0.0_dp
    matR    = 0.0_dp
    matI    = 0.0_dp
    mat     = 0.0_dp


    write(*,*) 'PIInit: InitInfoFile= '//trim(adjustl(initInfoFile))
    open(21, file = initInfoFile)
    read(21, nml=initInfo)
    close(21)
    mat = 0.0_dp

    if (initMode ==0) then ! Read from file
       if (initFile /= "" .AND. lineN /= 0) then
          if (UBound > 0 .and. LBound > 0) then
             call file2vec(size, initFile, lineN, vec, LBound = LBound, UBound = UBound)
             write(*,*) 'vec:', vec
          else
             call file2vec(size, initFile, lineN, vec)
          end if
        else
          write(*,*) 'Error: initMode = 0 requires filename and line number'
          stop 0
       end if
       vec=vec*initAmp
    else if (initMode == 1) then ! Gaussian Blob
       sigma = param1
       do j = 1, size
          mat(j, 3:4) = exp(-((j-1-(size/2))*(j-1-(size/2)))/(2.0*sigma*sigma))
       end do
       mat = initAmp * mat
       call mat2vec(size, mat, vec)
    else if (initMode == 2) then ! Randomised
       write(*,*) 'init: initMode=2. Random noise init'
       if (seed /= 0) then
          seed_(1) = seed
          call random_seed(put = seed_)
       else
          call random_seed()
       end if
       do j = LBound, UBound
          call random_number(matR(j, 3))
          call random_number(matR(j, 4))
          call random_number(matI(j, 3))
          call random_number(matI(j, 4))
       end do
       matR(LBound:UBound, 3:4) = matR(LBound:UBound,3:4) * 2.0_dp - 1.0_dp
       matI(LBound:UBound, 3:4) = matI(LBound:UBound,3:4) * 2.0_dp - 1.0_dp
       mat(LBound:UBound, 3:4)  = matR(LBound:UBound, 3:4) + im*matI(LBound:UBound, 3:4)
       mat = initAmp * mat
       call mat2vec(size, mat, vec)
       write(*,*) vec
    else if (initMode==3) then
       ! Fourier Initialisation Mode 

       ! Initialise random seed
       if (sizef == -1) then
          error stop
       end if

       if (seed /= 0) then
          seed_(1) = seed
          call random_seed(put=seed_)
       else
          call random_seed()
       end if
       mat(:,:) = 0d0
       ! Read in mode amplitudes.
       ! Format phi: A(0), A(1), A(-1), A(2) etc
       !          n: A(0), A(1), A(-1), A(2) etc

       call file2fourier(size, sizef, initfile, fvec)

       ! Apply random phases
       do j1=1,size
          do j2=1,2
             call random_number(fphase)
             fvec(j1,j2) = fvec(j1,j2)*exp(im*fphase*2.0d0*pi)
          end do
       end do

       write(*,*) 'Fourier init mode: fvec: ', fvec

       ! Setup FFTW3
       allocate(in(size))
       allocate(out(size))
       planb=fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)

       ! Execute FFTW3
       do j1=1,2
          in(:)=fvec(:,j1)
          write(*,*) 'in: ', in(:)
          call fftw_execute_dft(planb,in,out)
          write(*,*) 'out: ', out(:)
          fvec(:,j1) = out(:)
       end do
       write(*,*) 'Fourier init mode: fvec: ', fvec
       call fftw_destroy_plan(planb)
       ! Convert to single vector
       mat(LBound:UBound, 3:4)=initAmp*fvec(:,:)
       call mat2vec(size,mat,vec)
    end if

    call translate(vec,vecTmp,translateL, phaseP, globalN, globalE)
    vec = vecTmp

  end subroutine getVec

  function norm(N, vec)
    real(dp) :: norm
    integer, intent(in) :: N
    real(dp), dimension(N) :: vec
    integer :: j
    
    norm = 0.0_dp
    do j = 1, N
       norm = norm + vec(j)*vec(j)
    end do
  end function norm
  
  subroutine translate(vecIn, vecOut, translateL,phaseP,globalN, globalE)
    real(dp), dimension(:), intent(in)  :: vecIn
    real(dp), dimension(:), intent(out) :: vecOut
    type(C_PTR) :: planf, planb
    complex(dp), dimension(:), allocatable :: in, out
    complex(dp), dimension(:,:), allocatable :: mat,tmp
    integer :: m,n,j, i1, ind
    real(dp) :: k
    real(dp), intent(in) :: translateL, phaseP, globalN, globalE

    n=size(vecIn)/6

    allocate(in(n))
    allocate(out(n))
    allocate(mat(n,4))
    allocate(tmp(n,4))
    tmp(:,:)=0d0
    call vec2mat(n, vecIn, mat)
    
    if (translateL /= 0) then
       planf=fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
       planb=fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
       do m = 1,4
          in(:)=mat(:,m)
          call fftw_execute_dft(planf, in, out)
          do j=1,n/2
             k=-2.0*pi*j/n
             out(n+1-j)=out(n+1-j)*exp(im*k*translateL)
          end do
          do j=1,n/2+1
             k=2.0*pi*(j-1)/n
             out(j)=out(j)*exp(im*k*translateL)
          end do
          in(:)=out(:)
          call fftw_execute_dft(planb,in,out)
          mat(:,m)=out(:)/n
       end do
    end if

    !write(*,*) 'mat0 = ', mat(:,1)
    
    do i1 =1,n
       tmp(i1, 1) = mat(i1,1)+globalN
       tmp(i1, 2) = mat(i1,2)+globalE

       tmp(i1, 3) = mat(i1,3)*exp(im*phaseP)
       tmp(i1, 4) = mat(i1,4)*exp(im*phaseP)
    end do

    !write(*,*) 'mat1 = ',tmp(:,1)

    call mat2vec(n, tmp, vecOut)

  end subroutine
end module PIInit
