module ut_m
  use constants_m
  implicit none

  public :: del_spaces

contains
  subroutine del_spaces(s)
    implicit none
    character(*), intent(inout) :: s
    character(len=len(s)) :: tmp
    integer(SI) :: i, j
    j = 1
    do i = 1, len(s)
      if (s(i:i)==' ') cycle
      tmp(j:j) = s(i:i)
      j = j + 1
    end do
    s = tmp(1:j-1)
  end subroutine del_spaces

end module ut_m
