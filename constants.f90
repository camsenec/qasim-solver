module constants_m
  use iso_c_binding
  implicit none

  ! << f90 constants >>
  integer, parameter :: SI = selected_int_kind(8)
  integer, parameter :: DI = selected_int_kind(16)
  integer, parameter :: SR = selected_real_kind(6)
  integer, parameter :: DR = selected_real_kind(12)

  ! << Mathematical constants >
  real(DR), parameter :: EPS = 1e-5
  real(DR), parameter :: INF = 1e5

  ! <<IO constants >
  ! in : 読み込みファイルのための装置番号
  integer(SI),parameter :: in = 17
  ! out : 読み込みファイルのための装置番号
  integer(SI),parameter :: out = 19

  integer(SI),parameter :: DIV = 10
  integer(SI), parameter :: DIV_LIGHT = 5

  integer(SI), parameter :: COLOR_NUM = 4
  integer(SI), parameter :: A = 0

end module constants_m
