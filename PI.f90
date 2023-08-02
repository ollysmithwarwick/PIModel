! PI.f90  - 2.0.0 - 01/08/2022                                                      !
! Olly Smith (olly.smith@warwick.ac.uk)                                             !

module PIModule
  ! A module containing subroutines which implement the Plasma Interchange (PI) model !
  ! using pseudospectral methods.                                                     !
  use, intrinsic :: iso_c_binding
  use PI_io
  use parameters
  use PIRunModule
  implicit none
  include 'fftw3.f03'

  ! PLANS FOR FFTS !
  type(C_PTR) :: planf, planb, plan32f, plan32b

  ! STATE DATA AND CALCULATION VARIABLES !
  complex(C_DOUBLE_COMPLEX), dimension(:),   allocatable :: in,   out,  in32, out32           ! Allocated memory for FFTs. Done this way to allow easier alignment of data                   
  complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: rval, fval, dval, tempval, newval ! Arrays containing current state and data used in calculations

  ! ITERATORS / COUNTERS !
  integer :: m, njump                                ! Integers of do loops and iterators

  ! GLOBAL INFO ABOUT SYSTEM/RUN !
  integer :: N_x, N_t, nonlin, N_frames                       ! Number of gridpoints (must be even for now), timesteps. Nonlinear flag (1 includes nonlin terms, 0 only linear - this may need to be checked) 
  real(dp) :: L, S, dt

  ! DERIVED QUANTITIES !
  integer :: N32, spf                                        ! Number of grid points in the expanded 3/2 space
  real(dp), dimension(:), allocatable :: k                   ! Vector containing k (wavenumber) values for each point in Fourier space
  real(dp) :: f, f2, f32, L_inv, t, t_total, dx              ! Helpful divisor terms

  ! CONSTANTS !
  integer, parameter :: NBAR = 1, E = 2, PHITILDE = 3, NTILDE = 4, DXNBAR = 5, DXPHITILDE = 6, B_PLUS = 5, B_MINUS = 6 ! These are labels to make identifying rval entries easier

contains

  subroutine PI_Main(PIRun0, states)
    ! The code works using the general sequence of steps: !
    !       1. Read all input information from setup.in, init.dat                       !
    !       2. Convert the initial conditions to Fourier space                          !
    !       3. Calculate nonlinear terms at the current time                            ! 
    !       4. Perorm a half-timestep dt/2 using the nonlinear terms found in 3         !
    !       5. Calculate the nonlinear terms at t  + dt/2                               !
    !       6. Use the terms found in step 5. to perform a full timestep t              !
    !       7. Remap k values if the shear causes them to go be closer to another       !
    !           grid point.                                                             !
    !       8. Go back to 3 until the full time is done, then print results.            !

    type(PIRun) :: PIRun0
    real(dp), dimension(1:,1:), intent(inout) :: states
    real(dp), dimension(:), allocatable :: init, vec
    real(dp) :: shift, phase
    integer :: nf, j

    write(*,*) 'PI_main: Called'
    
    S = PIRun0%S
    dt = PIRun0%dt
    L = PIRun0%L
    N_x = PIRun0%N_x
    N_t = PIRun0%N_t
    spf = PIRun0%spf
    nonlin = PIRun0%nonlin
    shift = PIRun0%shift
    phase = PIRun0%phase
    states = 0.0_dp

    call PI_setup()
!    write(*,*) 'PI_main: PI_setup complete' 

!    allocate(states(6*N_x, N_frames))
    allocate(init(6*N_x))
    allocate(vec(6*N_x))
!    write(*,*) 'PI_main: Allocation complete'
    vec = 0.0_dp
    init = PIRun0%init
    call vec2mat(N_x, init, rval)
!    write(*,*) 'vec2mat complete'
    states(:, 1) = init(:)

    t = 0
    t_total = 0
    nf = 0
    call fft(1)
!    write(*,*) 'shp = ', shape(states)

    do j = 1, N_t
       call time_step(dt)
       if (mod(j, spf) == 0) then
          nf = nf + 1
          call fft(-1)
          call mat2vec(N_x, rval, vec)
          states(:,nf + 1) = vec(:)
       end if
    end do

!    write(*,*) 'PI_main: Shifting final state'
    do j = 1, N_x
       fval(j,     NBAR) = fval(j,     NBAR) * exp(im*(k(j)*shift))
       fval(j,        E) = fval(j,        E) * exp(im*(k(j)*shift))
       fval(j, PHITILDE) = fval(j, PHITILDE) * exp(im*(k(j)*shift))
       fval(j,   NTILDE) = fval(j,   NTILDE) * exp(im*(k(j)*shift))
    end do

    call fft(-1)
    rval(:,  NTILDE) = rval(:,  NTILDE)*exp(im*phase)
    rval(:,PHITILDE) = rval(:,PHITILDE)*exp(im*phase)
    call mat2vec(N_x, rval, vec)
    states(:, nf + 2) = vec(:)

!    write(*,*) 'PI_main: Done'

!    if (PIRun0%statesFile /= '') then
!       call writeFile(PIRun0%statesFile, states)
!    end if

!    if (PIRun0%plotInfoFile /='') then
!       call run2file(PIRun0, PIRun0%plotInfoFile)
!    end if

    call end_pi
!    write(*,*) 'PI_main: Done'

  end subroutine PI_Main

  subroutine PI_setup()
    ! Sets up variables for use in the code.                  !
    !   ext: a flag determining whether the file 'init.dat'   !
    !        will be used (1) for the initial conditions      !
    !        or if these will be determined otherwise (0)     !

    ! integer :: ext
    ! character(len = 50) :: ConfigFile
    ! type(run) :: run
    ! ! Read Setup namelist for the variables: Name, S, N_x, N_t, nframes, dt, L, nonlin!
    ! open(21, file = ConfigFile)
    ! read(21, nml = PI_nml)
    ! close(21)
    
    ! Determine N32 = 3/2 * N_x, rounded up.
    N32 = CEILING(N_x*(1.0_dp+nonlin*1.0_dp/2.0_dp)) 

    ! Calculate commonly used fractions and store them as variables to save time
    f = (1.0_dp/(L*L))                               
    f2 = (1.0_dp/N_x)                                
    f32 = (1.0_dp/N32)
    L_inv = (1.0_dp/L)

    dx = L/(N_x)                                     ! dx is the length per grid point
!    steps_per_frame = int(N_t/N_frames)              ! Number of timesteps between each frame captured

    ! Allocate sizes to all complex arrays
    allocate(in(N_x))
    allocate(out(N_x))
    allocate(in32(N32))
    allocate(out32(N32))
    allocate(rval(N_x, 4))
    allocate(fval(N_x, 4))
    allocate(dval(N_x,4))
    allocate(tempval(N_x,4))
    allocate(newval(N_x,4))
    allocate(k(N_x))
    ! ---------- CREATE AN ARRAY CONTAINING k VALUES ------------ !
    ! Array index     : 1 2 3 4 5 6 7 8 ... N_x-3 N_x-2 N_x-1 N_x !
    ! k value (*2pi/L): 0 1 2 3 4 5 6 7 ...  -4     -3   -2    -1 !

    ! Negative k's first to allow k(N_x/2 + 1) to be overwritten !
    do m = 1, N_x/2
       k(N_x+1-m) = -2*pi*L_inv*m
    end do

    ! Positive k's !
    do m = 1, N_x/2 + 1            
       k(m) = 2*pi*L_inv*(m-1)
    end do
    
    
    ! -------------------------------------------------- !
    ! ------- READ INITIAL CONDITIONS FROM FILE ---------!
    ! ------ WRITE INFO FILE FOR PLOTTING SCRIPT --------!
    ! ------------------ PLAN FFTW ----------------------!
    rval(:,:) = 0.0_dp

!    call writeplotinfo(plotInfoFile)
    call plan
  end subroutine PI_setup

  
  subroutine readinit(inFileName, vec)
    ! Reads initial conditions from inFileName !
    character(len = 100), intent(in)  :: inFileName
    real(dp), dimension(6*N_x) :: vec
    call file2vec(N_x, inFileName, 1, vec)
    call vec2mat(N_x, vec, rval)
  end subroutine readinit
  
  ! subroutine writeplotinfo(plotInfoFile) ! Formerly writeinfo
  !   ! ------- WRITES GENERAL METADATA TO A FILE FOR PLOTTING SCRIPT PLOT.PY TO READ ------- !
  !   character(len = 50), intent(in) :: plotInfoFile
  !   open(31, file = plotInfoFile)
  !   write(31, '(A50)') NAME
  !   write(31, '(F23.15)') S
  !   write(31, '(I5)') N_frames
  !   write(31, '(I6)') N_x
  !   write(31, '(F23.15)') L
  !   write(31, '(F23.15)') (dt*N_t)/N_frames
  !   write(31, '(F23.15)')  dt
  !   write(31, '(I8)'    )  N_t
  !   close(31)
  ! end subroutine writeplotinfo

  subroutine plan
    ! DECLARE PLANS FOR FFTW !
    planf   = fftw_plan_dft_1d(N_x, in, out, FFTW_FORWARD,  FFTW_ESTIMATE)
    planb   = fftw_plan_dft_1d(N_x, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
    plan32f = fftw_plan_dft_1d(N32, in32, out32, FFTW_FORWARD, FFTW_ESTIMATE)
    plan32b = fftw_plan_dft_1d(N32, in32, out32, FFTW_BACKWARD, FFTW_ESTIMATE)
  end subroutine plan

  subroutine writeout(app, outFile)
    !  WRITES OUTPUT INTO OUT.DAT !
    !  flag: DETERMINES WHETHER TO OPEN FILE AS NEW !
    !        (REMOVING ALL PREVIOUS DATA, USED AT   !
    !        START OF A RUN) OR APPEND THE CURRENT  !
    !        STATE ONTO THE END                     !
    integer              :: app
    character(len = 100) :: outFile
    if (app == 0) then
       open(11, file = outFile,      form = 'formatted')
    else if (app == 1) then
       open(11, file = outFile,      form = 'formatted', position = 'append')
    end if
    
    ! Writes each pair of rows row as real part then imag part of a quantity !
    ! Each column is a grid point !
    do m = 1, 4
       write(11,1) real(rval(:,m))
       write(11,1) DIMAG(rval(:,m))
    end do
    close(11)
1   format(100000f30.15)
  end subroutine writeout

  subroutine shift(d, t_)
    ! Applies a forward or backward phase shift !
    ! exp{+-iStx) across all real space. See    !
    ! guide for motivation.                     !
    !   d: 1 for forward, -1 for backward       !
    !  t_: Time since last repeat kshift        !
    implicit none
    integer :: d 
    real(dp) :: t_
    complex(dp) :: p, q
    p = 1.0_dp ! Shift to apply
    q = zexp(d*im*S*t_*dx) ! Relative shift between neigbouring grid points
    do m = 1, N_x
       rval(m, NTILDE  ) = rval(m,NTILDE  )*p
       rval(m, PHITILDE) = rval(m,PHITILDE)*p
       p = p * q
    end do
  end subroutine shift

  subroutine fft(d)
    ! FFT Forward (d = 1) or Backward (d = -1) including relevant phase shifts
    implicit none
    integer :: d
    if (d == 1) then
       call shift(d, t)

       do m = 1, 4
          in(:) = rval(:, m)
          call fftw_execute_dft(planf, in, out)
          fval(:,m) = out(:)
       end do

       fval(:,:) = fval(:,:)*f2
!       fval(N_x/2+1,:) = 0.0_dp ! THIS MAY NOT BE NEEDED BUT REMOVES SINGLE UNSTABLE MODE

    else if (d == -1) then
       do m = 1, 4
          in(:) = fval(:,m)
          call fftw_execute_dft(planb, in, out)
          rval(:,m) = out(:)
       end do
       call shift(d, t)
    end if
  end subroutine fft

  subroutine zeropadding(fval_, fval32_)
    ! TAKES A VECTOR fval_ AND ADDS ZEROS IN ADDED k-VALUE COMPONENTS TO MAKE IT N32 LONG
    implicit none
    complex(C_DOUBLE_COMPLEX), dimension(N_x, 6), intent(in) :: fval_
    complex(C_DOUBLE_COMPLEX), dimension(N32, 6), intent(out) :: fval32_

    fval32_(:,:) = 0.0_dp
    do m = 1, N_x/2+1
       fval32_(m,2:6) = fval_(m,2:6) 
    end do

    do m = N_x/2+2, N_x       
       fval32_(N32-N_x+m,2:6) = fval_(m,2:6)
    end do
  end subroutine zeropadding

  subroutine unpadding(a,b)
    ! PERFORMS INVERSE OF zeropadding
    complex(C_DOUBLE_COMPLEX), dimension(N32) :: a
    complex(C_DOUBLE_COMPLEX), dimension(N_x) :: b
    do m = 1, N_x/2+1
       b(m) = a(m)
    end do
    
    do m = N_x/2+2, N_x
       b(m) = a(N32-N_x+m)
    end do
  end subroutine unpadding

  subroutine init_time
    ! SETS UP TIME
    ! t IS THE TIME SINCE LAST SHIFT
    ! t_total IS THE TOTAL SIMULATED TIME
    ! njump IS THE TOTAL NUMBER OF k JUMPS
    t = 0
    t_total = 0
    njump = 0
!    nf = 0
  end subroutine init_time

  subroutine time_step(dt_)
    double precision :: dt_
    ! Perfroms a single time step !
    dval(:,:) = 0.0_dp
    
    call step_forward(dt_)
    t = t + dt_
    ! If k values shift (due to shear, k(t) = k(0) - St) to be closer to nearby k values perform shift !
    if (abs(S*t) > L_inv * pi) then
       call kshift
    end if
    t_total = t_total + dt_

  end subroutine time_step

  subroutine korder
    ! RESHAPES fval VALUES INTO ORDER    !
    ! THIS ALLOWS A KSHIFT TO BE SIMPLE  ! 
    do m = 1, N_x/2+1
       tempval(N_x-N_x/2+m-1,3:4)=fval(m,3:4)
    end do
    do m = N_x/2+2,N_x
       tempval(m-N_x/2-1,3:4)=fval(m,3:4)
    end do
  end subroutine korder

  subroutine norder
    ! PERFORMS THE INVERSE OF KORDER !
    do m = 1, N_x/2+1
       fval(m,3:4) = tempval(N_x-N_x/2+m-1,3:4)
    end do
    do m = N_x/2+2,N_x
       fval(m,3:4) = tempval(m-N_x/2-1,3:4)
    end do
  end subroutine norder

  subroutine kshift
    ! Shift every fval along one space !
    implicit none
    integer :: nshift
    call korder

    if (S > 0) then
       nshift = int((t+L_inv*pi/S)/(2*L_inv*pi/S))
       t = t - nshift * 2 * L_inv * pi / S
       do m = 1, N_x-nshift
          tempval(m,3:4)=tempval(m+nshift,3:4)
       end do
       do m = 1, nshift
          tempval(N_x+1-m,3:4) = 0.0_dp
       end do
    else
       nshift = int(-(t+L_inv*pi/S)/2*L_inv*pi*S)
       t = t + nshift * 2 * L_inv * pi / S
       do m = N_x, 1 + nshift
          tempval(m,3:4) = tempval(m-nshift ,3:4)
       end do
       do m = 1, nshift
          tempval(1,3:4) = 0.0_dp
       end do
    end if
    njump = njump+nshift

    call norder
  end subroutine kshift

  subroutine nonlinear(fval, t_, outval)
    ! DEALS WITH THE NONLINEAR TERMS OF PI MODEL        !
    ! outval STORES THE NONLINEAR TERMS' FOURIER MODES  !
    ! outval_NBAR(x) = iD_x(n\phi*-\phin*)              !
    ! outval_E(x)    = iD_x(\phiD_x\phi*-\phi*D_x\phi)  ! 
    ! outval_phi(x)  = -i\phi(E)                        !
    ! outval_n(x)    = -inE + i\phiD_x(nbar)            !

    implicit none
    complex(C_DOUBLE_COMPLEX), dimension(N_x,4)  :: fval
    complex(C_DOUBLE_COMPLEX), dimension(N_x,6)  :: fval_
    complex(C_DOUBLE_COMPLEX), dimension(N32,6) :: rval32, fval32
    complex(C_DOUBLE_COMPLEX), dimension(N32)   :: tempval32, dval32
    complex(C_DOUBLE_COMPLEX), dimension(N_x,4), intent(out) :: outval
    real(dp), intent(in) :: t_
    
    fval_(:,1:4) = fval(:,1:4)
    do m = 1, N_x
       fval_(m,DXNBAR)     = im *    k(m)     * fval_(m, NBAR    )
       fval_(m,DXPHITILDE) = im * (k(m)-S*t_) * fval_(m, PHITILDE)
    end do
    
    call zeropadding(fval_, fval32)
    do m = 2, 6
       in32(:) = fval32(:,m)
       call fftw_execute_dft(plan32b, in32, out32)
       rval32(:,m) = out32(:)
    end do

    ! --------- NBAR CALCULATION -------------!
    tempval32(:) = conjg(rval32(:,NTILDE))*rval32(:,PHITILDE) - conjg(rval32(:,PHITILDE))*rval32(:,NTILDE)
    in32(:) = tempval32(:)
    call fftw_execute_dft(plan32f, in32, out32)
    dval32(:) = out32(:)
    call unpadding(dval32(:),outval(:,NBAR))
    do m = 1, N_x
       outval(m,NBAR) = -k(m)*outval(m,NBAR)
    end do

    ! --------- E CALCULATION ------------------!
    tempval32(:) = rval32(:, PHITILDE)*conjg(rval32(:,DXPHITILDE)) - conjg(rval32(:,PHITILDE))*rval32(:,DXPHITILDE)
    in32(:) = tempval32(:)
    call fftw_execute_dft(plan32f, in32, out32)
    dval32(:) = out32(:)
    call unpadding(dval32(:),outval(:,E))
    do m = 1, N_x
       outval(m, E) = -k(m) * outval(m,E)
    end do

    ! --------- PHITILDE CALCULATION -----------!
    tempval32(:) = -im * rval32(:, PHITILDE)*rval32(:,E)
    in32(:) = tempval32(:)
    call fftw_execute_dft(plan32f, in32, out32)
    dval32(:) = out32(:)
    call unpadding(dval32(:),outval(:,PHITILDE))

    ! --------- NTILDE CALCULATION -------------!
    tempval32(:) = -im * rval32(:, NTILDE) * rval32(:,E) + im * rval32(:,PHITILDE) * rval32(:,DXNBAR)
    in32(:) = tempval32(:)
    call fftw_execute_dft(plan32f, in32, out32)
    dval32(:) = out32(:)
    call unpadding(dval32(:),outval(:,NTILDE))

    outval(:,:) = f32*outval(:,:) ! Rescale due to transform back and forth

    do m = 1, N_x/2 -1
!       outval(N_x + 1 - m, NBAR) = conjg(outval(m + 1, NBAR))
!       outval(N_x + 1 - m,    E) = conjg(outval(m + 1,    E))
    end do
!    outval(N_x/2 + 1, NBAR) = 0.0_dp
!    outval(N_x/2 + 1,    E) = 0.0_dp
  end subroutine nonlinear

 subroutine linear_trap(fval_, dval_, t_, dt_, newval_)
   ! 'LINEAR TRAPEZIUM RULE - USES THE TRAPEZIUM RULE TO DEAL WITH LINEAR TERMS !
   ! INCLUDES THE NONLINEAR TERMS GENERATED IN NONLIN AS dval_                  !
   real(dp) :: a_minus, a_plus, frac, phicoeff, t_, dt_
   complex(dp) :: ncoeff
   complex(C_DOUBLE_COMPLEX), dimension(N_x, 4) :: fval_, dval_, newval_
   
   do m = 1, N_x
      a_plus  = 1.0_dp + dt_*k(m)*k(m)*0.5_dp
      a_minus = 1.0_dp - dt_*k(m)*k(m)*0.5_dp
      frac    = a_minus/a_plus
      newval_(m,     NBAR)  = fval_(m,NBAR)        + dt_*dval_(m, NBAR)
      newval_(m,        E)  = fval_(m,   E)*(frac) + dt_*dval_(m,E)/(a_plus)
      
      a_plus  = 1.0_dp + dt_*(k(m)-S*(t_+dt_))*(k(m)-S*(t_+dt_))*0.5_dp
      a_minus = 1.0_dp - dt_*(k(m)-S*t_)*(k(m)-S*t_)*0.5_dp
      
      frac = a_minus/a_plus
      
      phicoeff = frac + dt_*dt_/(4*a_plus*a_plus)
      ncoeff   = (im*dt_/(2.0_dp*a_plus))*(1.0_dp + frac)
      
      newval_(m, PHITILDE) = (4.0_dp*a_plus*a_plus/(4.0_dp*a_plus*a_plus-dt_*dt_))*&
           (phicoeff*fval_(m,PHITILDE) + ncoeff*fval_(m, NTILDE) + (im*dt_*dt_/(2*a_plus*a_plus))*dval_(m,NTILDE)&
           + (1.0_dp/a_plus)*dt_*dval_(m,PHITILDE))
      
      newval_(m, NTILDE) = fval_(m,NTILDE)*frac - im*dt_*0.5_dp*(1.0_dp/a_plus)*&
           (fval_(m, PHITILDE) + newval_(m, PHITILDE)) + dt_*dval_(m,NTILDE)/a_plus
   end do
 end subroutine linear_trap
  
 subroutine step_forward(dt_)
   double precision :: dt_
   dval = 0d0
   if (nonlin==1) then
      ! COMPUTES A HALF TIMESTEP TO GET FVAL AT t+dt/2 !
      ! USES THIS FVAL TO GET dval FOR FULL TIMESTEP   !
      call nonlinear(fval,t,dval)
      call linear_trap(fval, dval, t, dt_*0.5_dp, tempval)
      call nonlinear(tempval, t+dt_*0.5_dp, dval)
   end if

   call linear_trap(fval, dval, t, dt_, newval) ! Full step !
   
   ! Shift temporary data into permanent
   fval(:,1:4) = newval(:,:)
   fval(N_x/2+1, :) = 0.0_dp ! This mode can get messy - easiest to set to 0 !

 end subroutine step_forward
 
 subroutine write_frame(v)
   integer :: v ! Verbose mode flag
   ! FT BACK TO REAL VALUES AND PRINTS OUTPUT USING WRITEOUT !
   if (v == 1) then
      write(*,*) t, t_total, njump , fval(1,PHITILDE)
   end if
   call fft(-1)
!   call writeout(1)
1  format(100000f30.15) ! NOTE THIS FORMAT MUST BE LARGE ENOUGH TO CONTAIN ALL POINTS
 end subroutine write_frame
 
 subroutine end_pi

   ! CLEAN UP OF MEMORY TO END PROGRAM !

   call fftw_destroy_plan(planf)
   call fftw_destroy_plan(planb)
   call fftw_destroy_plan(plan32f)
   call fftw_destroy_plan(plan32b)
   call fftw_cleanup()
   deallocate(in)
   deallocate(out)
   deallocate(in32)
   deallocate(out32)
   deallocate(rval)
   deallocate(fval)
   deallocate(dval)



   deallocate(tempval)
   deallocate(newval)
   deallocate(k)
 end subroutine end_pi
end module PIMODULE


