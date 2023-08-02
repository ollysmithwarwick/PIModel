module parameters
  implicit none
  integer    , parameter :: dp = selected_real_kind(15,307)
  real(dp)   , parameter :: pi = 4.0_dp * atan(1.0_dp)
  complex(dp), parameter :: im = (0.0_dp, 1.0_dp)  
end module parameters
