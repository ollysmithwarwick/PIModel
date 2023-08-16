program Bisect
! *********************************************** !
! Reads in an initialisation file, and parameters !
! and computes the time evolution of the PI model !
! system. The time evolution is stored in an      ! 
! output file                                     !
! *********************** Oliver Smith 31/07/2023 !

  use PIInit     
  use PIRunModule
  use PIModule

  implicit none
  type(PIRun) :: PIRunIn
  character (len=500) :: setupFile = 'BisectSetup.in', OutDir, PIRunFile, PIInitFile, StatesFile, TFile, ampfile, bisfile
  real(dp), dimension(:,:), allocatable :: states
  real(dp), dimension(:),   allocatable :: tVals

  ! Iterators
  integer :: jt, jf

  ! Bisection Params
  real(dp) :: lam, turb, initAmp, step, tmin
  integer :: maxits

  ! Bisection amplitude
  real(dp), dimension(:),   allocatable :: initamps
  real(dp), dimension(:), allocatable :: init
  real(dp), dimension(:,:), allocatable :: amps
  real(dp) :: R, upper, lower
  integer :: jb = 0

!  namelist /setup/ OutDir, PIRunFile, PIInitFile, StatesFile, TFile
  namelist /setup/  OutDir, PIRunFile, PIInitFile, StatesFile, TFile, AmpFile, BisFile
  namelist /bisectParams/ lam, turb, tmin, maxits, initAmp, step
  
  write(*,*) 'Bisect: Started'

  ! Read setupfile for locations of initialisation files
  write(*,*) 'Bisect: Reading setup file'
  open(12, file = setupFile)
  read(12, nml = setup)
  close(12)
  write(*,*) 'Bisect: Setup file read'

  ! Read PIRunFile into PIRunIn
  write(*,*) 'PI: Reading into PIRunIn from '//trim(adjustl(PIRunfile))
  call file2run(PIRunFile, PIRunIn)
  write(*,*) 'PI: PIRunIn read complete'
  write(*,*) 'PI: PIRunIn = '
  write(*,*) PIRunIn

  ! Read BisFile
  write(*,*) 'Bisect: Reading bisfile'
  open(12, file = BisFile)
  read(12, nml = bisectParams)
  close(12)
  write(*,*) 'Bisect: bisfile read'
  
  ! Read PIInitFile and uses it to determine initial conditions
  write(*,*) 'PI: Determining initial conditions from '//trim(adjustl(PIInitFile))
  call getVec(PIRunIn%N_x, PIInitFile, PIRunIn%init)
  write(*,*) 'PI: Initial conditions determined'

  ! Setup states and tVals with correct size
  write(*,*) 'PI: Allocating memory'
  allocate(states(6*PIRunIn%N_x, (PIRunIn%N_t/PIRunIn%spf) + 2))
  allocate(amps((PIRunIn%N_t/PIRunIn%spf) + 2, maxits))
  allocate(initamps(maxits))
  allocate(init(6*PIRunIn%N_x))
  allocate(tVals (               (PIRunIn%N_t/PIRunIn%spf) + 2))
  write(*,*) 'PI: Memory allocation complete'

  ! Bisection Setup
  R = initamp
  init = PIRunIn%init/amplitude(PIRunIn%init)

  states = 0d0
  amps = 0d0
  initamps = 0d0
  tVals = 0d0


  upper = -1
  lower = -1
  ! Main Loop
  do jb = 1, maxits
     write(*,*) 'Bisect: Iteration ', jb
     write(*,*) 'Bisect: Trying R = ', R
     
     PIRunIn%init = R * init
     call PI_Main(PIRunIn, states)
     do jf = 1, (PIRunIn%N_t/PIRunIn%spf) + 1
        amps(jf, jb) = amplitude(states(:,jf))
        if ((jf-1)*spf*dt < tmin) then
           cycle
        end if
        if (amplitude(states(:, jf)) < lam) then
           lower = R
           write(*,*) 'Bisect: Laminar'
           exit
        else if ((amplitude(states(:, jf)) > turb) .or. (isnan(states(0,jf)))) then
           upper = R
           write(*,*) 'Bisect: Turbulent'
           exit
        end if
     end do
     write(*,*) jf,(PIRunIn%N_t/PIRunIn%spf + 2) 
     if (jf == (PIRunIn%N_t/PIRunIn%spf)+2) then
        ! Not running long enough to decide and/or lam/turb need to be adjusted
        write(*,*) 'Bisect: Unable to determine if lam or turb'
        exit
     end if

     if (lower == -1) then
        R = R - step
     else if (upper == -1) then
        R = R + step
     else
        R = (upper + lower)/2d0
     end if
  end do

  ! Fill tVals with time values
  write(*,*) 'PI: Calculating tVals'
  jf = 0
  tVals(0) = 0d0
  do jt = 1, PIRunIn%N_t
    if (mod(jt,PIRunIn%spf) == 0) then
       jf = jf + 1
       tVals(jf + 1) = jt * PIRunIn%dt
    end if
  end do
  tVals(jf + 2) = PIRunIn%dt * PIRunIn%N_t
  write(*,*) 'PI: tVals calculation complete'
!  write(*,*) states(:,1000)
!  write(*,*) shape(states)
  call system('mkdir -p '//trim(adjustl(OutDir)))
  write(*,*) 'PI: Storing output.', TFile
  call writeFile( trim(adjustl(OutDir))//trim(adjustl(StatesFile)), states)
  call writeFile( trim(adjustl(OutDir))//trim(adjustl(AmpFile)), amps)
  call writeFile1D(trim(adjustl(OutDir))//trim(adjustl(TFile)),  tvals)
  write(*,*) 'PI: Copying input files'
  call system('cp '//trim(adjustl(PIRunFile))//' '//trim(adjustl(OutDir)))
  call system('cp '//trim(adjustl(PIInitFile))//' '//trim(adjustl(OutDir)))
  call system('cp '//trim(adjustl(bisFile))//' '//trim(adjustl(OutDir)))
  write(*,*) 'PI: Done'

  !deallocate(states)
  !deallocate(tvals)
  !deallocate(init)
  !deallocate(initamps)
  !deallocate(amps)

  contains

    function amplitude(vec)
      real(dp), dimension(:) :: vec
      integer :: j
      real(dp) :: tot
      real(dp) :: amplitude

      tot = 0d0
      do j = size(vec)/3, size(vec)
         tot = tot + vec(j)*vec(j)
      end do
      amplitude = sqrt(tot)/size(vec)
    end function amplitude

end program Bisect
