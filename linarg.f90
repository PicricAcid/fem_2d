module linarg
  implicit none
  private

  public:: solve

contains

!linear_solver
  subroutine solve(n, a, b, x, ipiv)
!intent(in) parameter
    integer, intent(in):: n
!intent(in) arrays
    double precision, intent(inout):: a(:, :), b(:), x(:)
    integer, intent(inout):: ipiv(:)
!variables
    integer:: i, j, k, unit, imax, jmax, t
    double precision:: sum, cir1, cir3, pivot
!arrays
    double precision, allocatable:: l(:, :), u(:, :), y(:), jp(:)
!init variables
    allocate(l(n, n), u(n, n), y(n), jp(n))
    l = 0.0d0
    u = 0.0d0

    do i=1, n
      jp(i) = i
    end do

!solve the equations Ax = b for x
    do k=1, n
!pivotting
      pivot = abs(a(k, k))
      imax = k
      jmax = k
      do i=k, n
        do j=k, n
          if(abs(a(i, k)) > pivot) then
            pivot = abs(a(i, k))
            imax = i
            jmax = j
          end if
        end do
      end do

      if(jmax /= k) then
        t = jp(k)
        jp(k) = jmax
        jp(jmax) = t
      end if

      do i=1, n
        pivot = a(i, k)
        a(i, k) = a(i, jmax)
        a(i, jmax) = pivot
      end do

      do j=1, n
        pivot = a(k, j)
        a(k, j) = a(imax, j)
        a(imax, j) = pivot
      end do

      pivot = b(k)
      b(k) = b(imax)
      b(imax) = pivot

      ipiv(k) = imax

!LU decomposition
      cir1 = 1.0d0/a(k, k)
      u(k, k) = a(k, k) !cir1
      do i=k+1, n
        a(i, k) = a(i, k)*cir1
      end do
      do j=k+1, n
        cir3 = a(k, j)
        u(k, j) = cir3
        do i=k+1, n
          a(i, j) = a(i, j) - a(i, k)*cir3
        end do
      end do
    end do

  do i=1, n
    do j=1, n
      if(i == j) then
        l(i, j) = 1.0d0
      else if(i > j) then
        l(i, j) = a(i, j)
      end if
    end do

  end do

! forward substitution
    do i=1, n
      sum = 0
      do k=1, i-1
        sum = sum + l(i, k)*y(k)
      end do
      y(i) = b(i) - sum
    end do

!back substitution
    do i=n, 1, -1
      sum = 0
      do k=i+1, n
        sum = sum + u(i, k)*x(k)
      end do
      x(i) = (y(i) - sum)/u(i, i)
    end do

  end subroutine solve

end module linarg
