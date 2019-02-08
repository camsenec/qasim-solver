function qa(n, m)
  use ut_m
  use field_m
  use constants_m
  use calc_energ_m
  use iso_c_binding
  implicit none
  !=========parameter　declarelation========
  !---------parameter for spinglass and general---------
  ! n : 1スライスにおけるサイト数, m : トロッター数
  integer(SI),intent(in) :: n, m
  ! i : 汎用イテレーター, tmp : 汎用一時変数, count : 汎用カウンタ
  integer(DI) :: i, tmp, count, count_color
  ! tau : モンテカルロステップ数
  integer(DI) :: tau
  !x:サイトのx座標, y:サイトのy座標
  integer(SI) :: x, y, c
  ! spin_old(i,:,:): i番目のトロッタースライスの遷移前の状態
  ! spin_new(i,:,:): i番目のトロッタースライスの遷移後の状態
  integer(SI), allocatable, dimension(:,:,:,:) :: spin_old, spin_new
  ! energ(i,c) : 色c, i番目のスライスの遷移前のエネルギー
  real(DR), allocatable, dimension(:,:) :: energ
  ! energ_old_qa : 遷移前の合計エネルギー
  ! energ_new_qa : 遷移後の合計エネルギ-
  ! energ_delta : 遷移前から遷移後のエネルギー差
  real(DR) :: energ_old_qa, energ_new_qa, energ_delta
  ! j_couple : カップリング
  real(DR), dimension(:,:,:,:), allocatable :: j_couple
  ! p: 反転させる確率, p_base : p>p_baseであるときに反転させる
  real(DR) :: prob, prob_base
  ! site_x, site_y : 反転させるサイトのx座標, y座標
  integer(SI) :: site_x, site_y
  ! qa : 返り値
  real(C_float) :: qa

  !--------parameter for qa--------
  ! qa_step : 量子アニーリングステップ数, sa_step : 古典アニーリングステップ数
  integer(DI) :: qa_step
  ! j_tilda : トロッタースライスごとの相互作用におけるカップリング
  real(DR) :: j_tilda
  ! gamma : アニーリング係数
  real(DR) :: gamma, gamma_init
  ! beta : 逆温度
  real(DR) :: beta
  ! k : スライス
  integer(SI):: k

  !======== initialize ========

  !-------- parameter for qa------

  ! set beta(becaues mt = m  / beta))
  beta = 10_DR
  ! set gamma_init
  gamma_init = 3_DR

  !-------- parameter set for qa------

  !set initial gamma
  gamma = gamma_init
  ! set qa_step
  qa_step = 100000

  !check parameter
  print *, 'beta:', beta
  print *, 'initial_gamma:', gamma_init
  print *, 'qa_step:', qa_step
  print *, 'n', n
  print *, 'm', m

  !-------- initialize --------
  ! set random seed
  call rnd_seed

  ! allocat memorye
  allocate(j_couple(n,n,n,n))
  allocate(spin_old(COLOR_NUM,n,n,m))
  allocate(spin_new(COLOR_NUM,n,n,m))
  allocate(energ(m, COLOR_NUM))

  ! initialize spin of all slice
  call init_sg(spin_old, m, n)

  ! initialize coupling
  call init_coupling(j_couple, n)

  ! initialize energ
  j_tilda = -1 / (2 * beta) * log(tanh(beta * gamma / m))
  energ_old_qa = energ_qa(j_couple, spin_old, j_tilda, m, n)

  !======== Quantumn Annealing ========
  do tau = 1, qa_step

    if(tau > 1) then
      ! remove TMF_term
      energ_old_qa = energ_old_qa - TMF_term(spin_old, j_tilda, m, n)
      ! calculate j_tilda
      j_tilda = -1 / (2 * beta) * log(tanh(beta * gamma / m))
      ! add TMF_term
      energ_old_qa = energ_old_qa + TMF_term(spin_old, j_tilda, m, n)
    end if

    do k = 1, m
      do i = 1,n*n
        ! select reversed site
        call choose(site_x, site_y, n)

        do c = 1, COLOR_NUM
          ! reverse spin
          call reverse_spin(site_x, site_y, spin_old, spin_new, c, k, m, n)

          ! calculate energ_new
          !energ_old_qa = energ_qa(j_couple, spin_old, j_tilda, m, n)
          !energ_new_qa = energ_qa(j_couple, spin_new, j_tilda, m, n)
          energ_delta = delta_qa(j_couple, spin_old, spin_new, j_tilda, site_x, site_y, c, k, m, n)
          energ_new_qa = energ_old_qa + energ_delta
          !print *,  k,i,c,tau, energ_new_qa - energ_old_qa, energ_delta

          ! calculate p
          if (energ_delta <= 0) then
            prob = 1
          else
            prob = exp(-beta * energ_delta)
          end if

          ! update sg based on probability p
          call random_number(prob_base)

          if (prob >= prob_base) then
            spin_old = spin_new
            energ_old_qa = energ_new_qa
          end if
        end do
      end do
    end do


    !end judge and outout
    count_color = 0
    if(mod(tau, DIV) == 0) then
      do c = 1, COLOR_NUM

        !calculate energy
        do k =  1, m
          energ(k, c) = energ_sa(j_couple, spin_old, c, k, m, n)
        end do

        print *, "j_tilda : ", j_tilda
        print *, "qa_step : ", tau
        count = 0
        do k = 1, m
          print *, beta, gamma, energ(k, c)

          !end judge among slices
          if (k < m .and. abs(energ(k, c) - energ(k + 1, c)) .le. EPS*1e-4) then
            count = count + 1
          end if
        end do

        !end judge among color
        if(count .ge. m - 1) then
          count_color = count_color + 1
        end if

      end do

    end if

    if(count_color .ge. COLOR_NUM) then
      exit
    end if

    !update gamma
    gamma = gamma*0.99

  end do

  !call output_spin(spin_old, 1_DI, m, n)
  qa = minval(energ(:,1))
  deallocate(j_couple)
  deallocate(spin_old, spin_new)
  deallocate(energ)

end function qa
