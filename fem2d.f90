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

  !全節点数と全要素数
  integer:: i, j
  integer:: el, no
  integer:: iargc
  integer:: count
  character(256):: filename
  double precision:: pivot

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
  !call gauss(2*no, kmat, dvec)

  !do i=1, 2*no
    !pivot = dvec(ipiv(2*no - i + 1))
    !dvec(ipiv(2*no - i + 1)) = dvec(2*no - i + 1)
    !dvec(2*no - i + 1) = pivot
  !end do

  call print_result(no, el, nd, elem, ddvec)
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

  subroutine print_result(no, el, nd, elem, dvec)
    implicit none
    integer, intent(in):: no, el
    type(node), intent(in):: nd(no)
    type(element), intent(in):: elem(el)
    double precision, intent(in):: dvec(2*no)

    integer:: i, j, ios, p
    double precision:: xmax = 0.0, xmin = 0.0, ymax = 0.0, ymin = 0.0
    double precision:: factor = 0.0, d = 0.0, maxd = 0.0

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
    kmat = 0.0
    !要素のループ
    do ie=1, el
      !dmatを作る
      f = elem(ie)%ym/(1.0 - elem(ie)%po*elem(ie)%po)
      dmat(1, 1) = f; dmat(1, 2) = elem(ie)%po*f; dmat(1, 3) = 0.0
      dmat(2, 1) = elem(ie)%po*f; dmat(2, 2) = f; dmat(2, 3) = 0.0
      dmat(3, 1) = 0.0; dmat(3, 2) = 0.0; dmat(3, 3) = ((1.0 - elem(ie)%po)/2.0)*f

      emat = 0.0; tmat = 0.0
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
          kmat(2*p-1, 2*q-1) = kmat(2*p-1, 2*q-1) + tmat(2*i-1, 2*j-1)
          kmat(2*p, 2*q-1) = kmat(2*p, 2*q-1) + tmat(2*i, 2*j-1)
          kmat(2*p-1, 2*q) = kmat(2*p-1, 2*q) + tmat(2*i-1, 2*j)
          kmat(2*p, 2*q) = kmat(2*p, 2*q) + tmat(2*i, 2*j)
        end do
      end do

      open(iunit, file="test.txt", status="replace")
      do i=1, 2*no
        write(iunit, *) (kmat(i, j), j=1, 2*no)
      end do
      close(iunit)

    end do
  end subroutine make_stiffness

  !bマトリックスを作る
  subroutine make_bmat(nelem, no, nd, xi, eta, bmat, det)
    implicit none
    type(element), intent(in):: nelem
    integer, intent(in):: no
    type(node), intent(in):: nd(no)
    double precision, intent(in):: xi, eta
    double precision, intent(inout):: bmat(3, 16)
    double precision, intent(inout):: det

    integer:: i, p
    double precision:: jacobian(2, 2) = 0.0
    double precision:: dndx, dndy

    !dndx, dndyを計算する
    !jacobianを作る
    do i=1, 8
      p = nelem%node(i)
      jacobian(1, 1) = jacobian(1, 1) + nd(p)%x * dndxi(i, xi, eta)
      jacobian(1, 2) = jacobian(1, 2) + nd(p)%y * dndxi(i, xi, eta)
      jacobian(2, 1) = jacobian(2, 1) + nd(p)%x * dndet(i, xi, eta)
      jacobian(2, 2) = jacobian(2, 2) + nd(p)%y * dndet(i, xi, eta)
    end do
    !(dndx, dndy)^t = J(dndxi, dndet)^t
    det = jacobian(1, 1)*jacobian(2, 2) - jacobian(1, 2)*jacobian(2, 1)

    bmat = 0.0
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
