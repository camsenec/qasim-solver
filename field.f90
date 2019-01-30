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
    integer(SI) :: x, y
    integer(SI), intent(in) :: n
    integer(DI) :: k
    integer(DI), intent(in) :: m
    integer(SI),dimension(n,n,m), intent(inout) :: spin
    real(SR) tmp

    do k = 1, m
      do y = 1, n
        do x = 1, n
          call random_number(tmp)
          spin(x,y,k) = nint(tmp)
          if (abs(spin(x, y, k)) == 0) then
            spin(x, y, k) = -1.0
          end if
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
      ! print *, j_couple(ix,iy,jx,jy)
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
  subroutine reverse_spin(site_x, site_y, spin_old, spin_new, k, m, n)
    implicit none
    integer(SI), intent(in) :: n, site_x, site_y
    integer(DI), intent(in) :: k, m
    integer(SI), dimension(n,n,m), intent(inout) :: spin_old, spin_new

    spin_new = spin_old
    if (spin_old(site_x, site_y, k) == -1) then
      spin_new(site_x, site_y, k) = 1
    else
      spin_new(site_x, site_y, k) = -1
    end if
  end subroutine reverse_spin

  subroutine output_spin(filename, spin, k, m, n)
    implicit none
    integer(SI), intent(in) :: n
    integer(SI), dimension(n,n,m), intent(in) :: spin
    integer(SI) :: ix,iy
    integer(DI) :: k, m
    integer(SI), parameter :: iw = 5000
    character(len=128), intent(in) :: filename

    open(iw,file=filename, status="replace")
    do iy = 1, n
      do ix = 1, n
        write(iw,fmt='(i4, a ,i4, a, i4)') ix, ',', iy, ',', spin(ix,iy,k)
      enddo
      ! 改行
      write(iw,*)
    enddo
    close(iw)

  end subroutine output_spin


end module field_m
