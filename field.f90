!-------------------------------------------------------------------
!  spinglass_qa.f: quantum annealing by the path-integral monte carlo method
!
!       by tanaka tomoya, Kobe University.
!-------------------------------------------------------------------
module field_m
  use constants_m
  use calc_energ_m
  implicit none

  public :: rnd_seed, &
             init_sg, &
             init_coupling, &
             choose, &
             reverse_spin

contains
  ! set random seed
  subroutine rnd_seed
    implicit none
    integer(SI) ::  i , seedize
    integer(SI), dimension(:), allocatable :: seed

    call random_seed(size=seedize)
    allocate(seed(seedize))
    call random_seed(get = seed)
    do i = 1, seedize
      call system_clock(count = seed(i))
    end do

    call random_seed(put = seed(:))
    deallocate(seed)
  end subroutine rnd_seed


  ! initialize spinglass
  subroutine init_sg(spin, m, n)
    implicit none
    integer(SI) :: x, y, c
    integer(SI), intent(in) :: n
    integer(DI) :: k
    integer(DI), intent(in) :: m
    integer(SI),dimension(COLOR_NUM,n,n,m), intent(inout) :: spin
    real(SR) tmp

    do k = 1, m
      do y = 1, n
        do x = 1, n
          do c = 1,COLOR_NUM
            call random_number(tmp)
            spin(c, x, y, k) = nint(tmp)
            if (abs(spin(c, x, y, k)) == 0) then
              spin(c, x, y, k) = -1
            end if
          end do
        end do
      end do
    end do
  end subroutine init_sg


  ! initialize coupling
  subroutine init_coupling(j_couple,n,in)
    implicit none
    real(DR), dimension(n,n,n,n), intent(inout) :: j_couple
    integer(SI) :: ix, iy, jx, jy
    integer(SI),intent(in) :: n
    real(DR) :: tmpj
    integer(SI) :: in
    integer(DI) :: count

    count = 0
    j_couple = 0

    do
      read(in,*,end=100) ix, iy, jx, jy, tmpj
      j_couple(ix,iy,jx,jy) = tmpj
      j_couple(jx,jy,ix,iy) = tmpj
      !print *, j_couple(ix,iy,jx,jy)
      count = count + 1
    end do
    100 close(in)
    print * ,count

  end subroutine init_coupling


  ! choose update site
  subroutine choose(site_x, site_y, n)
    implicit none
    integer(SI), intent(inout) :: site_x, site_y
    integer(SI), intent(in) :: n
    real(DR) :: tmp

    call random_number(tmp)
    site_x = ceiling(tmp * n)
    call random_number(tmp)
    site_y = ceiling(tmp * n)

  end subroutine choose


  ! reverse spin based on spin_old
  subroutine reverse_spin(site_x, site_y, spin_old, spin_new, c, k, m, n)
    implicit none
    integer(SI), intent(in) :: n, site_x, site_y, c
    integer(DI), intent(in) :: k, m
    integer(SI), dimension(COLOR_NUM,n,n,m), intent(inout) :: spin_old, spin_new

    spin_new = spin_old
    if (spin_old(c, site_x, site_y, k) == 1) then
      spin_new(c, site_x, site_y, k) = -1
    else
      spin_new(c, site_x, site_y, k) = 1
    end if
  end subroutine reverse_spin

  subroutine output_spin(filename, spin, k, m, n)
    implicit none
    integer(SI), intent(in) :: n
    integer(SI), dimension(COLOR_NUM,n,n,m), intent(in) :: spin
    integer(SI) :: ix,iy,c
    integer(DI) :: k, m
    character(len=256)::linebuf
    integer(SI), parameter :: iw = 5000
    character(len=128), intent(in) :: filename

    open(iw,file=filename, status="replace")
    do iy = 1, n
      do ix = 1, n
        do c = 1, COLOR_NUM
          write(linebuf, *) ix, ',', iy, ',' ,c , ',', spin(c,ix,iy,k)
          call del_spaces(linebuf)
          write (iw, '(a)') trim(linebuf)
        enddo
      enddo
    enddo
    close(iw)

  end subroutine output_spin

  subroutine del_spaces(s)
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


end module field_m
