program fem2d
  use linarg
  implicit none

  integer, parameter:: iunit = 10

  !要素内の節点は右上の角から反時計回りで要素座標xiとeta
  double precision, parameter:: xi(8) = (/1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0/)
  double precision, parameter:: eta(8) = (/1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0/)

  !ガウス-ルジャンドル4点積分の積分座標と重み
  double precision, parameter:: uip = sqrt(3.0 + 2.0*sqrt(6.0/5.0))/7.0
  double precision, parameter:: uim = sqrt(3.0 - 2.0*sqrt(6.0/5.0))/7.0
  double precision, parameter:: gu(4) = (/-uip, -uim, uim, uip/)
  double precision, parameter:: wip = (18.0 + sqrt(30.0))/36.0
  double precision, parameter:: wim = (18.0 - sqrt(30.0))/36.0
  double precision, parameter:: gw(4) = (/wim, wip, wip, wim/)

  !応力型　sigma_x, sigma_y, tau_xy
  type:: stress
    double precision:: xx, yy, xy
  end type stress
  !要素型 構成する節点, 縦弾数係数, ポアソン比, 厚さ
  type:: element
    integer:: node(8)
    double precision:: ym, po, thick
  end type element
  !節点型 節点座標, 境界条件値, 変位or荷重(0なら変位, 1なら荷重)
  type:: node
    double precision:: x, y, bcx, bcy
    integer:: xbt, ybt
  end type node

  !全節点数と全要素数と計算する応力点数(要素ごと)
  integer:: i, j
  integer:: el, no, st
  integer:: iargc
  integer:: count
  character(256):: filename
  double precision:: pivot, factor

  double precision, allocatable:: dvec(:) !節点変位ベクトル
  double precision, allocatable:: kmat(:, :) !剛性マトリックス
  double precision, allocatable:: ddvec(:)
  integer, allocatable:: ipiv(:)

  type(element), allocatable:: elem(:)
  type(node), allocatable:: nd(:)

!init
  !オペランド(ファイル名)を読み取る
  !iargc: オペランドの数を取得するfunction, getarg(i, arg): オペランドのi番目をargに代入するsubroutine
  count = iargc()
  if(count /= 1) then
    write(*, *) "input data file name"
    stop
  else
    call getarg(1, filename)
  end if

  !全節点数と全要素数をファイルから読み込み
  call read_param(filename, el, no)

  !allocate
  allocate(dvec(2*no), kmat(2*no, 2*no), ddvec(2*no), ipiv(2*no))
  allocate(elem(el), nd(no))

  !節点座標と境界条件を読み込み
  call read_data(filename, el, no, elem, nd)
  !initされたデータのサマリー
  write(*, *) el, no
  do i=1, el
    write(*, *) (elem(i)%node(j), j=1, 8), elem(i)%ym, elem(i)%po, elem(i)%thick
  end do
  do i=1, no
    write(*, *) nd(i)%x, nd(i)%y, nd(i)%bcx, nd(i)%bcy, nd(i)%xbt, nd(i)%ybt
  end do

  !剛性マトリックスの作成
  call make_stiffness(no, el, nd, elem, kmat)

  !境界条件の設定
  call set_bc(no, el, nd, elem, kmat, dvec)

  !連立方程式を解く
  call solve(2*no, kmat, dvec, ddvec, ipiv)

  !do i=1, 2*no
    !pivot = dvec(ipiv(2*no - i + 1))
    !dvec(ipiv(2*no - i + 1)) = dvec(2*no - i + 1)
    !dvec(2*no - i + 1) = pivot
  !end do

  factor = 0.0
  call print_result(no, el, nd, elem, ddvec, factor)
  call print_stress(no, el, nd, elem, ddvec)
  call internal_info(no, el, nd, elem, ddvec, factor)
  stop

contains

  subroutine read_param(filename, el, no)
    implicit none
    character(256), intent(in):: filename
    integer, intent(inout):: el, no
    integer:: i, ios

    open(iunit, iostat=ios, file=filename, status="old")
    if(ios /= 0) then
      write(*, *) 'file open error'
      close(iunit)
      stop
    else
      read(iunit, *) no
      do i=1, no
        read(iunit, '()') !節点座標を読み飛ばし
      end do
      read(iunit, *) el
      rewind(iunit)
      close(iunit)
    end if
  end subroutine read_param

  subroutine read_data(filename, el, no, elem, nd)
    implicit none
    character(256), intent(in):: filename
    integer, intent(in):: el, no
    type(element), intent(inout):: elem(el)
    type(node), intent(inout):: nd(no)
    integer:: i, ios, j, k
    double precision:: value

    open(iunit, iostat=ios, file=filename, status='old')
    if(ios /= 0) then
      write(*, *) 'file open error'
      close(iunit)
      stop
    else
      read(iunit, '()') !noを読み飛ばし
      do i=1, no
        read(iunit, *) nd(i)%x, nd(i)%y
      end do
      read(iunit, '()') !el読み飛ばし
      do i=1, el
        read(iunit, *) elem(i)%node, elem(i)%ym, elem(i)%po, elem(i)%thick
      end do

      !境界条件デフォルト指定
      do i=1, no
        nd(i)%bcx = 0.0
        nd(i)%xbt = 1 !デフォルトでは荷重指定
        nd(i)%bcy = 0.0
        nd(i)%ybt = 1 !デフォルトでは荷重指定
      end do
      !境界条件の読み込み
      read(iunit, *) j !x方向への変位が指定された節点の数を読み込み
      do i=1, j
        read(iunit, *) k, value
        nd(k)%bcx = value
        nd(k)%xbt = 0 !節点kのx方向に変位境界条件が指定されている
      end do
      read(iunit, *) j !y方向への変位が指定された節点の数を読み込み
      do i=1, j
        read(iunit, *) k, value
        nd(k)%bcy = value
        nd(k)%ybt = 0 !節点kのy方向に変位境界条件が指定されている
      end do
      read(iunit, *) j !x方向への0でない荷重が指定された節点を読み込み
      do i=1, j
        read(iunit, *) k, value
        nd(k)%bcx = value
      end do
      read(iunit, *) j !y方向への0でない荷重が指定された節点を読み込み
      do i=1, j
        read(iunit, *) k, value
        nd(k)%bcy = value
      end do
      close(iunit)
    end if
  end subroutine read_data

  subroutine print_result(no, el, nd, elem, dvec, factor)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(in):: dvec(2*no)
    double precision, intent(inout):: factor

    integer:: i, j, ios, p
    double precision:: xmax = 0.0, xmin = 0.0, ymax = 0.0, ymin = 0.0
    double precision:: d = 0.0, maxd = 0.0

    open(iunit, iostat=ios, file="displ.txt", status="replace")
    write(*, *) "------displacement------"
    write(iunit, *) "x_value(1)", "y_value(2)", "u_value(3)", "v_value(4)"
    do i=1, no
      write(*, *) i, "x=", nd(i)%x, "y=", nd(i)%y, "u=", dvec(2*i-1), "v=", dvec(2*i)
      if (nd(i)%x < xmin) xmin = nd(i)%x
      if (nd(i)%x > xmax) xmax = nd(i)%x
      if (nd(i)%y < ymin) ymin = nd(i)%y
      if (nd(i)%y > ymax) ymax = nd(i)%y
      d = sqrt(dvec(2*i-1)*dvec(2*i-1) + dvec(2*i)*dvec(2*i))
      !maxd: 変位の最大値
      if(d > maxd) maxd = d
    end do
    !代表長さ
    d = sqrt((xmax - xmin)*(ymax - ymin))
    factor = d/maxd*0.05
    do i=1, el
      do j=1, 8
        p = elem(i)%node(j)
        write(iunit, *) nd(p)%x, nd(p)%y, dvec(2*p-1), dvec(2*p)
      end do
      p = elem(i)%node(1)
      write(iunit, *) nd(p)%x, nd(p)%y, dvec(2*p-1), dvec(2*p)
      write(iunit, *)
    end do
    close(iunit)
    open(iunit, iostat=ios, file="gpd.plot", status="replace")
    write(iunit, *) "f =", factor
    write(iunit, *) "plot 'displ.txt' using 1:2 w l, 'displ.txt' using ($1+f*$3):($2+f*$4) w l, &
    'displ.txt' using 1:2 w p pt 6 ps 1"
    write(iunit, *) "pause -1"
    close(iunit)
  end subroutine print_result

  !剛性マトリックスを作る
  subroutine make_stiffness(no, el, nd, elem, kmat)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(inout):: kmat(2*no, 2*no)

    integer:: i, j, k, l, m, ie, p, q
    double precision:: f, xi, eta, det, sum
    double precision:: dmat(3, 3), bmat(3, 16), db(3, 16), emat(16, 16), tmat(16, 16)

    !Kmat = \int{B^t D B}
    !init kmax
    kmat(:, :) = 0.0
    !要素のループ
    do ie=1, el
      !dmatを作る
      call make_dmat(elem(ie), dmat)

      emat(:, :) = 0.0; tmat(:, :) = 0.0
      !積分ループ
      do i=1, 4
        xi = gu(i)
        do j=1, 4
          eta = gu(j)
          !bmatを作る
          call make_bmat(elem(ie), no, nd, xi, eta, bmat, det)
          !B^t D Bを計算する
          emat = matmul(transpose(bmat), matmul(dmat, bmat))
          !重みとjacobianのdetをかけて積分
          tmat = tmat + det*emat*gw(i)*gw(j)*elem(ie)%thick
        end do
      end do

      !tmat: 要素ごとの剛性マトリックス, kmat:　全体の要素マトリックス
      do i=1, 8
        p = elem(ie)%node(i)
        do j=1, 8
          q = elem(ie)%node(j)
          kmat(2*q-1, 2*p-1) = kmat(2*q-1, 2*p-1) + tmat(2*j-1, 2*i-1)
          kmat(2*q, 2*p-1) = kmat(2*q, 2*p-1) + tmat(2*j, 2*i-1)
          kmat(2*q-1, 2*p) = kmat(2*q-1, 2*p) + tmat(2*j-1, 2*i)
          kmat(2*q, 2*p) = kmat(2*q, 2*p) + tmat(2*j, 2*i)
        end do
      end do

    end do
  end subroutine make_stiffness

  subroutine make_dmat(nelem, dmat)
    implicit none
    type(element), intent(in):: nelem
    double precision, intent(inout):: dmat(3, 3)

    double precision:: f

    f = nelem%ym/(1.0 - nelem%po*nelem%po)
    dmat(1, 1) = f; dmat(1, 2) = nelem%po*f; dmat(1, 3) = 0.0
    dmat(2, 1) = nelem%po*f; dmat(2, 2) = f; dmat(2, 3) = 0.0
    dmat(3, 1) = 0.0; dmat(3, 2) = 0.0; dmat(3, 3) = ((1.0 - nelem%po)/2.0)*f
  end subroutine make_dmat

  !bマトリックスを作る
  subroutine make_bmat(nelem, no, nd, xi, eta, bmat, det)
    implicit none
    type(element), intent(in):: nelem
    integer, intent(in):: no
    type(node), intent(in):: nd(no)
    double precision, intent(in):: xi, eta
    double precision, intent(inout):: bmat(3, 16)
    double precision, intent(inout):: det

    integer:: i, j, p
    double precision:: jacobian(2, 2)
    double precision:: dndx, dndy

    !dndx, dndyを計算する
    !jacobianを作る
    jacobian(:, :) = 0.0

    do i=1, 8
      p = nelem%node(i)
      jacobian(1, 1) = jacobian(1, 1) + nd(p)%x * dndxi(i, xi, eta)
      jacobian(1, 2) = jacobian(1, 2) + nd(p)%y * dndxi(i, xi, eta)
      jacobian(2, 1) = jacobian(2, 1) + nd(p)%x * dndet(i, xi, eta)
      jacobian(2, 2) = jacobian(2, 2) + nd(p)%y * dndet(i, xi, eta)
    end do
    !(dndx, dndy)^t = J(dndxi, dndet)^t
    det = jacobian(1, 1)*jacobian(2, 2) - jacobian(1, 2)*jacobian(2, 1)

    bmat(:, :) = 0.0
    !bmatを作る
    do i=1, 8
      dndx = (jacobian(2, 2)*dndxi(i, xi, eta) - jacobian(1, 2)*dndet(i, xi, eta))/det
      dndy = (-jacobian(2, 1)*dndxi(i, xi, eta) + jacobian(1, 1)*dndet(i, xi, eta))/det
      bmat(1, 2*i-1) = dndx
      bmat(2, 2*i) = dndy
      bmat(3, 2*i-1) = dndy
      bmat(3, 2*i) = dndx
    end do

    det = abs(det)
  end subroutine make_bmat

  !bmat内のNiとdNi/dxi, dNi/detaを計算するfunctions
  !Ni(xi, eta)
  function ni(i, nxi, neta) result(z)
    implicit none
    integer, intent(in):: i
    double precision, intent(in):: nxi, neta
    double precision:: z

    double precision:: x2, y2, xxi, yyi, xi2, yi2
    x2 = nxi*nxi; y2 = neta*neta; xxi = nxi*xi(i); yyi = neta*eta(i); xi2 = xi(i)*xi(i); yi2 = eta(i)*eta(i)

    z = (1.0 + xxi)*(1.0 + yyi)*(xxi + yyi - 1.0)*xi2*yi2/4.0 &
    + (1.0 - x2)*(1.0 + yyi)*(1.0 - xi2)/2.0 + (1.0 - y2)*(1.0 + xxi)*(1.0 - yi2)/2.0
    return
  end function ni

  !dNi(xi, eta)/dxi
  function dndxi(i, nxi, neta) result(z)
    implicit none
    integer, intent(in):: i
    double precision, intent(in):: nxi, neta
    double precision:: z

    double precision:: x2, y2, xxi, yyi, xi2, yi2
    y2 = neta*neta; xxi = nxi*xi(i); yyi = neta*eta(i); xi2 = xi(i)*xi(i); yi2 = eta(i)*eta(i)

    z = (1.0 + yyi)*xi2*xi(i)*yi2*(2.0*xxi + yyi)/4.0 &
    - nxi*(1.0 + yyi)*(1.0 - xi2) + (1.0 - y2)*xi(i)*(1.0 - yi2)/2.0
    return
  end function dndxi

  !dNi(xi, eta)/deta
  function dndet(i, nxi, neta) result(z)
    implicit none
    integer, intent(in):: i
    double precision, intent(in):: nxi, neta
    double precision:: z

    double precision:: x2, y2, xxi, yyi, xi2, yi2
    x2 = nxi*nxi; xxi = nxi*xi(i); yyi = neta*eta(i); xi2 = xi(i)*xi(i); yi2 = eta(i)*eta(i)

    z = (1.0 + xxi)*xi2*yi2*eta(i)*(xxi + 2.0*yyi)/4.0 &
    + (1.0 - x2)*eta(i)*(1.0 - xi2)/2.0 - neta*(1.0 + xxi)*(1.0 - yi2)
    return
  end function dndet

  subroutine make_stress(s, no, nd, nelem, uvec, xi, eta)
    implicit none
    type(stress), intent(inout):: s
    integer, intent(in):: no
    type(node), intent(in):: nd(no)
    type(element), intent(in):: nelem
    double precision, intent(in):: uvec(16)
    double precision, intent(in):: xi, eta

    double precision:: jacobian(2, 2), bmat(3, 16), dmat(3, 3), db(3, 16)
    double precision:: det
    integer:: i, j

    jacobian(:, :) = 0.0

    call make_dmat(nelem, dmat)
    call make_bmat(nelem, no, nd, xi, eta, bmat, det)

    s%xx = 0.0; s%yy = 0.0; s%xy = 0.0

    db = matmul(dmat, bmat)

    do i=1, 8
      j = nelem%node(i)
      s%xx = s%xx + db(1, 2*i-1)*uvec(2*j-1) + db(1, 2*i)*uvec(2*j)
      s%yy = s%yy + db(2, 2*i-1)*uvec(2*j-1) + db(2, 2*i)*uvec(2*j)
      s%xy = s%xy + db(3, 2*i-1)*uvec(2*j-1) + db(3, 2*i)*uvec(2*j)
    end do
  end subroutine make_stress

  subroutine print_stress(no, el, nd, elem, uvec)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(in):: uvec(16)

    integer:: i, j, ie
    double precision:: nxi, neta, x, y
    type(stress):: s

    write(*, *) "stress"
    do ie=1, el
      do i=1, 8
        nxi = xi(i)
        neta = eta(i)
        call make_stress(s, no, nd, elem(ie), uvec, nxi, neta)
        x = 0.0; y = 0.0
        do j=1, 8
          x = x + nd(elem(ie)%node(j))%x*ni(j, nxi, neta)
          y = y + nd(elem(ie)%node(j))%y*ni(j, nxi, neta)
        end do
        write(*, *) ie, "x=", x, "y=", y, "sx=", s%xx, "sy=", s%yy, "txy=", s%xy
      end do
    end do
  end subroutine print_stress

  subroutine internal_info(no, el, nd, elem, uvec, factor)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(in):: uvec(16)
    double precision, intent(in):: factor

    integer:: i, j, k, ie, divi, p
    double precision:: xi, eta, x, y, u, v
    type(stress):: s

    divi = 5

    open(iunit, file="internal.txt", status="replace")
    do ie=1, el
      do i=-divi, divi
        eta = dble(i)/dble(divi)
        do j=-divi, divi
          xi = dble(j)/dble(divi)
          x = 0.0; y = 0.0; u = 0.0; v = 0.0
          do k=1, 8
            p = elem(ie)%node(k)
            x = x + nd(p)%x*ni(k, xi, eta)
            y = y + nd(p)%y*ni(k, xi, eta)
            u = u + uvec(2*p-1)*ni(k, xi, eta)
            v = v + uvec(2*p)*ni(k, xi, eta)
          end do
          call make_stress(s, no, nd, elem(ie), uvec, xi, eta)
          write(iunit, *) "x=", x, "y=", y, "u=", u, "v=", v, "sx=", s%xx, "sy=", s%yy, "txy=", s%xy
        end do
        write(iunit, *) !splotのための空行
      end do
      write(iunit, *) !要素ごとの空行
    end do
    close(iunit)
    open(iunit, file="gps.plot", status="replace")
    write(iunit, *) "set grid"
    write(iunit, *) "set pm3d"
    write(iunit, *) "set pm3d map"
    write(iunit, *) "set pm3d interpolate 3,3"
    !write(iunit, *) "set format cb '%%4.1le{%%L}'"
    write(iunit, *) "set size ratio -1"
    write(iunit, *) "set view 0,0,1"
    write(iunit, *) "set contour"
    write(iunit, *) "set palette defined ( -1 '#000030', 0 '#000090', 1 '#000fff', 2 &
    '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 &
    '#7f0000')"
    write(iunit, *) "unset surface"
    write(iunit, *) "unset key"
    write(iunit, *) "f =", factor
    write(iunit, *) "set title 'Sxx plot(deform. factor =", factor, ")'"
    write(iunit, *) "splot 'internal.txt' u ($2+f*$6):($4+f*$8):($10) w l"
    write(iunit, *) "pause -1"
    write(iunit, *) "set title 'Syy plot (deform. factor =", factor, ")'"
    write(iunit, *) "splot 'internal.txt' u ($2+f*$6):($4+f*$8):($12) w l"
    write(iunit, *) "pause -1"
    write(iunit, *) "set title 'Sxy plot (deform. factor =", factor, ")'"
    write(iunit, *) "splot 'internal.txt' u ($2+f*$6):($4+f*$8):($14) w l"
    write(iunit, *) "pause -1"
    close(iunit)
  end subroutine internal_info

  !境界条件の設定
  subroutine set_bc(no, el, nd, elem, kmat, dvec)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(inout):: kmat(2*no, 2*no)
    double precision, intent(inout):: dvec(2*no)

    integer:: i, j

    !init dvec
    do i=1, no
      dvec(2*i-1) = nd(i)%bcx
      dvec(2*i) = nd(i)%bcy
      !bcxが変位として指定されている場合,　kmatのi行がi列目だけ1であと0
      if(nd(i)%xbt == 0) then
        do j=1, 2*no
          kmat(2*i-1, j) = 0.0
        end do
        kmat(2*i-1, 2*i-1) = 1.0
      end if
      !bcyが変位として指定されている場合,　kmatのi行がi列目だけ1であと0
      if(nd(i)%ybt == 0) then
        do j=1, 2*no
          kmat(2*i, j) = 0.0
        end do
        kmat(2*i, 2*i) = 1.0
      end if
    end do
  end subroutine set_bc

end program fem2d
