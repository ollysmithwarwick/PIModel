module PIRunModule
! *********************************************************** !
! Contains the derived type 'PIRun' which allows for the easy !
! storage and passing of model parameters. Also conatins tool !
! routines related to input and output of associated namelist !
! files                                                       !
! ********************************* Oliver Smith 31/07/2023   ! 

  use parameters
  implicit none

  type PIRun
     real(dp), dimension(:), allocatable :: init          ! The initial conditions of the PI model in a 1D vector
     real(dp) :: S              ! Background shear
     real(dp) :: dt             ! Timestep size
     real(dp) :: L              ! Domain size
     real(dp) :: shift = 0.0_dp ! Translation distance of final state (allows one to compare initial and translated final states - e.g RPOSolve)
     real(dp) :: phase = 0.0_dp ! Phase rotation of final state
     
     integer  :: N_x            ! Number of x grid points
     integer  :: N_t            ! Number of timesteps to perform
     integer  :: spf            ! Number of timesteps per output frame (spf = steps-per-frame)
     integer  :: nonlin=1       ! Flag for nonlinear terms. 0 = off (not tested - avoid). 1 = on
  end type PIRun
  
  interface write(formatted)
     module procedure writePIRun ! Overwrite default write for PIRun
  end interface write(formatted)

  contains

  subroutine assignRun(in, out)
  ! Acts like assignment operator out = in !
    type(PIRun) :: in, out
    out%S      = in%S
    out%dt     = in%dt
    out%L      = in%L
    out%shift  = in%shift
    out%phase  = in%phase
    out%N_x    = in%N_x
    out%N_t    = in%N_t
    out%spf    = in%spf
    out%nonlin = in%nonlin
    if (allocated(out%init)) then
       deallocate(out%init)
    end if
    allocate(out%init(size(in%init)))
    out%init   = in%init
  end subroutine assignRun

  subroutine run2file(runIn, fileName)
  ! Prints PIRun to a file via namelist !
    type(PIRun), intent(in) :: runIn
    character(len = 100), intent(in) :: fileName
    real(dp) :: S, dt, L, shift, phase
    integer :: N_x, N_t, spf, nonlin

    namelist /PIRun/ S, dt, L, shift, phase, N_x, N_t, spf, nonlin

    S          = runIn%S
    dt         = runIn%dt
    L          = runIn%L
    shift      = runIn%shift
    phase      = runIn%phase
    N_x        = runIn%N_x
    N_t        = runIn%N_t
    spf        = runIn%spf
    nonlin     = runIn%nonlin

    open(12, file = fileName)
    write(12, nml=PIRun)
    close(12)
  end subroutine run2file   
  
  subroutine file2run(fileName, runOut)
  ! Reads file into PIRun via namelist
    type(PIRun), intent(out) :: runOut
    character(len = 100) :: fileName
    
    real(dp) :: S, dt, L, shift, phase
    integer :: N_x, N_t, spf, nonlin, initMode
    
    namelist /runParams/ S, dt, L, N_x, N_t, spf, nonlin, shift, phase
    open(12, file = fileName)
    read(12, nml = runParams)
    close(12)
    
    runOut%S=S
    runOut%dt=dt
    runOut%L=L
    runOut%N_x = N_x
    runOut%N_t = N_t
    runOut%spf = spf
    runOut%nonlin = nonlin
    runOut%shift = shift
    runOut%phase = phase

    allocate(runOut%init(6*N_x))
    runOut%init(:)=0d0
  end subroutine file2run

  subroutine writePIRun(dtv, unit,iotype, v_list,iostat, iomsg)
    !! User_Defined Derived-Type I/O Standard
    class(PIRun), intent(in) :: dtv
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'shape(init) = ', shape(dtv%init),new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'S           = ', dtv%S,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'L           = ', dtv%L,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'shift       = ', dtv%shift,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'phase       = ', dtv%phase,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'N_x         = ', dtv%N_x,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'N_t         = ', dtv%N_t,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'spf         = ', dtv%spf,new_line('a')
    write(unit, fmt=*, iostat=iostat, iomsg=iomsg) 'nonlin      = ', dtv%nonlin    ,new_line('a')
  end subroutine writePIRun
end module PIRunModule
