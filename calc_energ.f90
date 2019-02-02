module calc_energ_m
  use constants_m
  implicit none

contains
  ! calculate energy of a slice on sa
  function energ_sa(j_couple, spin, c, k, m, n)
    implicit none
    real(DR), dimension(n,n,n,n), intent(in) :: j_couple
    integer(SI), dimension(COLOR_NUM,n,n,m) :: spin
    integer(SI) :: ix, iy, jx, jy, c
    integer(SI), intent(in) :: n
    integer(DI), intent(in) :: k, m
    real(DR) :: energ_sa

    energ_sa = 0.0d0
    do jy = 1, n
      do jx = 1, n
        do iy = 1, n
          do ix = 1, n
            energ_sa = energ_sa - 1/2.0 * j_couple(ix ,iy ,jx, jy) &
            * spin(c, ix, iy, k) * spin(c, jx, jy, k)
          end do
        end do
      end do
    end do

  end function energ_sa

  ! calculate average of energy of all replica in qa
  function energ_qa(j_couple, spin, j_tilda, m, n)
    implicit none
    real(DR), dimension(n,n,n,n), intent(in) :: j_couple
    real(DR), intent(in) :: j_tilda
    real(DR), dimension(:), allocatable :: energ_0, energ_1, energ_2
    integer(SI), dimension(COLOR_NUM,n,n,m), intent(in) :: spin
    integer(SI), intent(in) :: n
    integer(SI) :: ix, iy, jx, jy, c, tmp
    integer(DI), intent(in) :: m
    integer(DI) :: k, sum_color
    real(DR) :: energ_qa

    allocate(energ_0(m))
    allocate(energ_1(m))
    allocate(energ_2(m))

    do k = 1, m
      energ_0(k) = 0.0d0
      energ_1(k) = 0.0d0
      energ_2(k) = 0.0d0
    end do


    do k = 1, m
      do jy = 1, n
        do jx = 1, n
          do iy = 1, n
            do ix = 1, n
              do c = 1, COLOR_NUM
                energ_0(k) = energ_0(k) - (1/2.0) * j_couple(ix, iy, jx, jy) &
                * spin(c, ix, iy, k) * spin(c, jx, jy, k)
              end do
            end do
          end do
        end do
      end do
    end do

    !カラーごとに相互作用を計算する(カラーごとにスピンが揃っていく)
    do k = 2, m - 1
      do iy = 1, n
        do ix = 1, n
          do c = 1, COLOR_NUM
            energ_1(k) = energ_1(k) - spin(c, ix, iy, k - 1) * spin(c, ix, iy, k)
            energ_1(k) = energ_1(k) - spin(c, ix, iy, k) * spin(c, ix, iy, k + 1)
          end do
        end do
      end do
    end do

    !k = m
    do iy = 1, n
      do ix = 1, n
        do c = 1, COLOR_NUM
          energ_1(m) = energ_1(m) - spin(c, ix, iy, m-1) * spin(c, ix ,iy, m)
          energ_1(m) = energ_1(m) - spin(c, ix, iy, m) * spin(c, ix, iy, 1)
        end do
      end do
    end do

    !k = 1
    do iy = 1, n
      do ix = 1, n
        do c = 1,COLOR_NUM
          energ_1(1) = energ_1(m) - spin(c, ix, iy, m) * spin(c, ix ,iy, 1)
          energ_1(1) = energ_1(m) - spin(c, ix, iy, 1) * spin(c, ix ,iy, 2)
        end do
      end do
    end do
    !  print *, "energ_1", sum(energ_0) / size(energ_0)
    !  print *, "energ_2", j_tilda * sum(energ_1)

    !for penalty
    do k = 1,m
      do iy = 1,n
        do ix = 1,n
          sum_color = 0
          do c = 1, COLOR_NUM
            tmp = spin(c, ix, iy, k)
            if(tmp == -1) then
              tmp = 0
            end if
            sum_color = sum_color + tmp
          end do
          energ_2(k) = energ_2(k) + A * (1 - sum_color) * (1 - sum_color)
        end do
      end do
    end do

    energ_qa = sum(energ_0) / m + j_tilda * sum(energ_1) + sum(energ_2)
  end function energ_qa

  function delta_qa(j_couple, spin, j_tilda, site_x, site_y, k ,m, n)
    implicit none
    real(DR), intent(in) ::  j_tilda
    real(DR), dimension(n,n,n,n), intent(in) :: j_couple
    integer(SI), dimension(COLOR_NUM,n,n,m), intent(in) :: spin
    integer(SI) :: x, y, c, tmp
    integer(SI), intent(in) :: site_x, site_y, n
    integer(DI), intent(in) :: m, k
    integer(DI) :: sum_color
    real(DR) :: delta_energ_0, delta_energ_1, delta_energ_2
    real(DR) :: delta_qa

    delta_energ_0 = 0.0d0
    delta_energ_1 = 0.0d0
    delta_energ_2 = 0.0d0

    !一項目
    do c = 1, COLOR_NUM
      do y = 1, n
        do x = 1, n
          delta_energ_0 = delta_energ_0 - (1/2.0) * &
          j_couple(x, y, site_x, site_y) * spin(c, x, y, k) * spin(c, site_x, site_y, k)

          delta_energ_0 = delta_energ_0 - (1/2.0) * &
          j_couple(site_x, site_y, x, y) * spin(c, x, y, k) * spin(c, site_x, site_y, k)
        end do
      end do
    end do

    !二項目
    if(k == m) then
      do c = 1, COLOR_NUM
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, m-1) * &
        spin(c, site_x, site_y, m)
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, m) * &
        spin(c, site_x, site_y, 1)
      end do
    else if(k == 1) then
      do c = 1, COLOR_NUM
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, m) * &
        spin(c, site_x, site_y, 1)
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, 1) * &
        spin(c, site_x, site_y, 2)
      end do
    else
      do c = 1, COLOR_NUM
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, k - 1) * &
        spin(c, site_x, site_y, k)
        delta_energ_1 = delta_energ_1 - spin(c, site_x, site_y, k) * &
        spin(c, site_x, site_y, k + 1)
      end do
    end if

    !for penalty
    sum_color = 0
    do c = 1, COLOR_NUM
      tmp = spin(c, site_x, site_y, k)
      if(tmp == -1) then
        tmp = 0
      end if
      sum_color = sum_color + tmp
    end do
    delta_energ_2 = delta_energ_2 + A * (1 - sum_color) * (1 - sum_color)

    delta_qa = 2 * ((delta_energ_0 / m) + (j_tilda * delta_energ_1)) + delta_energ_2

  end function delta_qa
end module calc_energ_m
