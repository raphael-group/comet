! Find the kth entry of one in the vth row or column of the matrix A.

function kth_one(arr, v, direction, k, m, n) result(i)

  implicit none

  integer, intent(in) :: direction, m, n, k, v
  integer, intent(in) :: arr(m, n)
  integer :: total, i

  ! Keep track of the number of ones we've seen so far.
  total = 0

  ! If direction==1, then search along column k.
  if (direction .eq. 1) then

     do i=1, m
        if (arr(i, v) .eq. 1) then
           total = total + 1
        end if
        if (total .eq. k) then
           exit
        end if
     end do

  ! Otherwise, if direction==2, then search along row k.
  else

     do i=1, n
        if (arr(v, i) .eq. 1) then
           total = total + 1
        end if
        if (total .eq. k) then
           exit
        end if
     end do

  end if

end function kth_one

! Seed the PRNG with the integer n.

! The random_seed routine requires m seeds, but we just pass the same
! seed n a total of m times for convenience.

subroutine seed_prng(n)

  implicit none

  integer, intent(in) :: n
  integer :: m
  integer, allocatable :: x(:)

  call random_seed(size = m)
  allocate(x(m))
  x = n
  call random_seed(put = x)
  deallocate(x)

end subroutine seed_prng

! Find the smallest index i such that dist(i)>p.

! We use a binary search.

function choice(dist, n, p) result(i)

  implicit none

  integer, intent(in) ::  n
  double precision, intent(in) :: dist(n)
  double precision, intent(in) :: p
  integer :: i, j

  i = n/2
  j = n/4

  do while (j>0)
    if (dist(i)<p) then
      i = i+j
    else
      i = i-j
    end if
    j = j/2
  end do

  do while (i<n .and. dist(i+1)<=p)
    i = i+1
  end do

  do while (i>1 .and. dist(i-1)>p)
    i = i-1
  end do

end function choice

! Construct the normalized cumulative sum of arr.

subroutine cdf(output, arr, n)

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: arr(n)
  integer :: i
  double precision :: output(n)

  output(1) = arr(1)
  do i=2,n
     output(i) = arr(i) + output(i-1)
  end do
  output = output / sum(arr)

end subroutine cdf

subroutine bipartite_edge_swap(B, A, nswap, max_tries, seed, verbose, m, n)

  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: m, n, nswap, max_tries, seed
  integer, intent(in) :: A(m,n)
  integer, intent(out) :: B(m,n)
  integer :: u, v, x1, y1, x2, y2, swapcount, iter, i
  integer :: kth_one, choice
  integer :: x_degrees(m), y_degrees(n)
  double precision :: r(4), x_dist(m), y_dist(n)

  ! Initialize the permuted matrix.
  B = A
  swapcount = 0
  iter = 0

  ! Find the degrees of the nodes in each partition.
  do i=1,m
    x_degrees(i) = sum(B(i,:))
  end do
  do i=1,n
    y_degrees(i) = sum(B(:,i))
  end do

  ! Convert the degree into degree distributions.
  call cdf(x_dist, x_degrees, m)
  call cdf(y_dist, y_degrees, n)

  ! Use the given random seed.
  call seed_prng(seed)

  ! Try to swap until reaching the maximum number of swaps or tries.
  do while (swapcount < nswap .and. iter < max_tries)

    ! Select four random numbers uniformly in [0,1] for each iteration.
    call random_number(r)

    ! Select random nodes in each partition.
    x1 = choice(x_dist, m, r(1))
    y1 = choice(y_dist, n, r(2))

    ! Choose a random neighbor of each node.
    u  = int(r(3) * x_degrees(x1)) + 1
    v  = int(r(4) * y_degrees(y1)) + 1

    y2 = kth_one(B, x1, 2, u, m, n)
    x2 = kth_one(B, y1, 1, v, m, n)

    ! If x1 is not already connected to y1 and x2 is not already
    ! connected to y2, then we can perform the edge swap.

    if ((B(x1, y1) .eq. 0) .and. (B(x2, y2) .eq. 0)) then
      swapcount = swapcount + 1

      ! Update the adjacency matrix.
      B(x1, y1) = 1
      B(x1, y2) = 0
      B(x2, y2) = 1
      B(x2, y1) = 0

    end if

    iter = iter + 1

  end do

  ! Warn if the maximum number of tries were reached.
  if (iter >= max_tries) then
     print *, '[Warning] Bipartite edge swap: Maximum tries exceeded.'
  end if

  if (verbose) then
    print *, 'Number of swaps:        ', swapcount
    print *, 'Maximum number of swaps:', nswap
    print *, 'Number of tries:        ', iter
    print *, 'Maximum number of tries:', max_tries
  end if

end subroutine bipartite_edge_swap
