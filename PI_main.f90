program PIExample
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
  character (len=500) :: setupFile = 'PISetup.in', OutDir, PIRunFile, PIInitFile, StatesFile, TFile
  real(dp), dimension(:,:), allocatable :: states
  real(dp), dimension(:),   allocatable :: tVals

  ! Iterators
  integer :: jt, jf

  namelist /setup/ OutDir, PIRunFile, PIInitFile, StatesFile, TFile
  
  write(*,*) 'PI: Started'

  ! Read setupfile for locations of initialisation files
  write(*,*) 'PI: Reading setup file'
  open(12, file = setupFile)
  read(12, nml = setup)
  close(12)
  write(*,*) 'PI: Setup file read'

  ! Read PIRunFile into PIRunIn
  write(*,*) 'PI: Reading into PIRunIn from '//trim(adjustl(PIRunfile))
  call file2run(PIRunFile, PIRunIn)
  write(*,*) 'PI: PIRunIn read complete'
  write(*,*) 'PI: PIRunIn = '
  write(*,*) PIRunIn

  ! Read PIInitFile and uses it to determine initial conditions
  write(*,*) 'PI: Determining initial conditions from '//trim(adjustl(PIInitFile))
  call getVec(PIRunIn%N_x, PIInitFile, PIRunIn%init)
  write(*,*) 'PI: Initial conditions determined'

  ! Setup states and tVals with correct size
  write(*,*) 'PI: Allocating memory'
  allocate(states(6*PIRunIn%N_x, (PIRunIn%N_t/PIRunIn%spf) + 2))
  allocate(tVals (               (PIRunIn%N_t/PIRunIn%spf) + 2))
  write(*,*) 'PI: Memory allocation complete'

  ! Call main PI function
  write(*,*) 'PI: Calling PI_main'
  call PI_Main(PIRunIn, states)
  write(*,*) 'PI: PI_main complete'


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
  write(*,*) states(:,1000)
  write(*,*) shape(states)
  call system('mkdir -p '//trim(adjustl(OutDir)))
  write(*,*) 'PI: Storing output'
  call writeFile( trim(adjustl(OutDir))//trim(adjustl(StatesFile)), states)
  call writeFile1D( trim(adjustl(OutDir))//trim(adjustl(TFile))     ,  tvals)
  write(*,*) 'PI: Copying input files'
  call system('cp '//trim(adjustl(PIRunFile))//' '//trim(adjustl(OutDir)))
  call system('cp '//trim(adjustl(PIInitFile))//' '//trim(adjustl(OutDir)))
  write(*,*) 'PI: Done'
!  deallocate(states)
!  deallocate(tvals)
  
end program PIExample
