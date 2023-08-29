module PI_io
! Tools for i/o relating to PIModel !

  use parameters
  use PIRunModule
  implicit none
  character(len = 20) :: mainfmt, mainfmtin
!  integer :: size, ierr
  
contains 

  subroutine setupIo(n, ierr)
    integer, intent(in)  :: n
    integer, intent(out) :: ierr
!    size = n
    ierr = 0
  end subroutine setupIo

  subroutine file2vec(sizex, fileName, lineN, vec, LBound, UBound)
    ! Read a vector from lineN of fileName. !

    integer, intent(in) :: sizex
    character(len = 100), intent(in) :: fileName
    real(dp), dimension(:) :: vec

    integer, intent(in) :: lineN
    integer, intent(in), optional :: UBound, LBound ! Can be used to take a range of values that is not the whole line in filename


    real(dp), dimension(:), allocatable :: vect
    integer :: U, L
    integer :: i1, i2
    allocate(vect(size(vec)))
    vec(:) = 0.0_dp
    vect(:) = 0.0_dp
    if(present(UBound) .and. present(LBound)) then
       U = UBound
       L = LBound
    else
       U = size(vec)
       L = 1
    end if
    write(mainfmt, '(A1,I8,A11)') '(', 6*sizex,'(E30.15E3))'
    open(12, file = fileName)
    do i1 = 1, lineN
       read(12, fmt = mainfmt) vect(:)
    end do
    close(12)
    vec(L:U) = vect(L:U)
    vec(L+sizex:U+sizex) =     vect(L+sizex:U+sizex)
    vec(2*L-1+2*sizex:2*U+2*sizex) =     vect(2*L-1+2*sizex:2*U+2*sizex)
    vec(2*L-1+4*sizex:2*U+4*sizex) =     vect(2*L-1+4*sizex:2*U+4*sizex)
    
  end subroutine file2vec

  subroutine file2fourier(sizex,sizef,filename,fvec)
    ! Reads in fourier coefficients from a file !
    integer, intent(in) :: sizex, sizef
    character(len = 100), intent(in) :: fileName
    complex(dp), dimension(:,:), intent(out) :: fvec
    real(dp), dimension(:,:), allocatable :: vect
    integer :: i1
    integer :: half
    allocate(vect(sizef,2))
    open(12,file=fileName)

    do i1=1,2
       read(12,*) vect(:,i1)
    end do
    close(12)
    half=(sizef+1)/2
    fvec(1:half,:)=vect(1:half,:)
    fvec(sizex-half:sizex,:)=vect(sizef-half:sizef,:)
  end subroutine file2fourier

  subroutine vec2file(size, fileName, vec, app)
    ! Print contents of vec to filename !
    ! If app is true the line will be appended !
    character, intent(in) :: fileName(:)
    integer, intent(in) :: app, size
    integer :: i1, i2
    real(dp), dimension(6*size) :: vec

    write(mainfmt, '(A1,I8,A11)') '(', 6*size,'(E30.15E3))'

    if (app == 1) then
       open(12, file = fileName, position = 'APPEND')
    else
       open(12, file = fileName)
    end if
    write(12, fmt = mainfmt) vec
    close(12)
  end subroutine vec2file

  subroutine writeFile(fileName, mat)
    ! Write contents of a matrix to filename !
    character(*), intent(in) :: fileName
    character(len = 21) :: outfmt
    real(dp), intent(in), dimension(1:, 1:) :: mat
    integer, dimension(2) :: shp
    integer :: cols, rows
    integer :: j

    shp = shape(mat)
    cols = shp(1)
    rows = shp(2)
    open(12, file = trim(adjustl(fileName)))
    write(outfmt, '(A1, I9, A11)') '(', cols, '(E30.15E3))'

    do j = 1, rows
       write(12, fmt=outfmt) mat(:, j)
    end do

    close(12)

  end subroutine writeFile

  subroutine writeFile1D(filename, vec, append)
    character(*), intent(in) :: filename
    character(len=21) :: outfmt
    logical, optional :: append
    real(dp), intent(in), dimension(1:) :: vec
    integer, dimension(1) :: shp
    integer :: cols

    shp=shape(vec)
    cols = shp(1)
    write(*,*) 'cols:', cols
    write(outfmt, '(A1, I9, A11)') '(', cols, '(E30.15E3))'
    write(*,*) FILENAME
    write(*,*) adjustl(trim(filename))
    if (present(append)) then
       if (append == .true.) then
          open(12, file=trim(adjustl(filename)), action='write',position='append')
       end if
    else
       open(12, file=trim(adjustl(filename)))
    end if
    write(12, fmt=outfmt) vec(:)
    close(12)

  end subroutine writeFile1D

  subroutine writeFileComplex(fileName, mat)
    character(*), intent(in) :: fileName
    character(len = 21) :: outfmt
    complex(dp), intent(in), dimension(1:, 1:) :: mat
    integer, dimension(2) :: shp
    integer :: cols, rows
    integer :: j

    shp = shape(mat)
    cols = shp(1)
    rows = shp(2)
    open(12, file = fileName)
    write(outfmt, '(A1, I9, A11)') '(', cols*2, '(E30.15E3))'

    do j = 1, rows
       write(12, fmt=outfmt) mat(:, j)
    end do

    close(12)

  end subroutine writeFileComplex

  subroutine vec2mat(size, vec, mat)
    ! Convert from 1D vec to 2D mat format
    integer, intent(in) :: size
    real(dp), dimension(6*size)    :: vec
    complex(dp), dimension(size, 4) :: mat
    integer :: i1, i2

    do i1 = 1, size
       mat(i1, 1) = vec(i1) ! NBAR
       mat(i1, 2) = vec(i1+size) !E
       mat(i1, 3) = vec(2*i1 + 2*size - 1) + im*vec(2*i1 + 2*size) !PHITILDE
       mat(i1, 4) = vec(2*i1 + 4*size - 1) + im*vec(2*i1 + 4*size) ! NTILDE
    end do
  end subroutine vec2mat

  subroutine mat2vec(size, mat, vec)
    integer, intent(in) :: size
    real(dp), dimension(6*size) :: vec
    complex(dp), dimension(size, 4), intent(in) :: mat
    integer :: i1, i2
!    write(*,*) 'mat2vec: mat: ', mat
    do i1 = 1, size
       vec(i1         ) = dreal(mat(i1, 1))
       vec(i1 +   size) = dreal(mat(i1, 2))
       vec(2*i1 + 2*size - 1) = dreal(mat(i1, 3))
       vec(2*i1 + 2*size    ) = dimag(mat(i1, 3))
       vec(2*i1 + 4*size - 1) = dreal(mat(i1, 4))
       vec(2*i1 + 4*size    ) = dimag(mat(i1, 4))
    end do
  end subroutine mat2vec

  subroutine writePlotInfo(runIn, statesFile, plotDir, fileName)
    type(PIRun) :: runIn
    character(len = 100) :: statesFile, plotDir, fileName
    real(dp) :: S, dt, L
    integer :: N_x, N_t, spf

    namelist /plotinfo/ S, dt, L, N_x, N_t, spf, statesFile, plotDir

    S = runIn%S
    dt = runIn%dt
    L = runIn%L
    N_x  = runIn%N_x
    N_t = runIn%N_t
    spf = runIn%spf
     
    open(12, file = fileName)
    write(12, nml=plotInfo)
    close(12)   
  end subroutine writePlotInfo

  subroutine GetStats(states, amps, fluxs, avgAmps, avgFluxs, pkAmps, pkFluxs)
    real(dp), dimension(1:,1:), intent(in) :: states

    real(dp), dimension(1:,1:) :: amps, fluxs
    real(dp), dimension(1:) :: avgAmps, avgFluxs, pkAmps, pkFluxs

    complex(dp), dimension(:), allocatable :: phi, n
    integer :: Nx,Nt, i1, i2
    integer, dimension(2) :: shp

    shp=shape(states)
    Nx=shp(1)/6
    Nt=shp(2)
    allocate(phi(Nx))
    allocate(n(Nx))

    do i1=1,Nt
       n(:)=0d0
       phi(:)=0d0

       do i2=1,Nx
          phi(i2) = states(2*i2 + 2*Nx -1, i1) + im * states(2*i2 +2*Nx, i1)
          n(i2)   = states(2*i2 + 4*Nx -1, i1) + im * states(2*i2 +4*Nx, i1)
       end do
       amps(:,i1) = real((phi(:)*conjg(phi(:)) + n(:)*conjg(n(:))))
       fluxs(:,i1) = real(-im*(phi(:)*conjg(n(:)) - n(:)*conjg(phi(:))))
       avgamps(i1) = sum(amps(:,i1))/Nx
       avgfluxs(i1) = sum(fluxs(:,i1))/Nx
       pkamps(i1) = maxval(amps(:,i1))
       pkfluxs(i1) = maxval(fluxs(:,i1))
    end do
  end subroutine GetStats

end module PI_io
