function nth_one(arr, index, t, num, m, n)
  
  integer, intent(in) :: num, m, n, index, t
  integer, intent(in) :: arr(m, n)
  integer :: total, i, nth_one

  total = 0 ! number of ones we've seen so far

  ! If we are searching in a row
  if (t .eq. 0) then
     do i=1, n
        if (arr(index, i) .eq. 1) then
           total = total + 1
        end if
        if (total .eq. num) then
           nth_one = i
           exit
        end if
     end do
  ! Otherwise we're searching in a column
  else
     do i=1, m
        if (arr(i, index) .eq. 1) then
           total = total + 1
        end if
        if (total .eq. num) then
           nth_one = i
           exit
        end if
     end do
  end if
end function nth_one

! Seed the PRNG with the given number.
! random_seed requires a collection of m seeds, but
! we just pass in the same seed m times
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

function choice(dist, n, p)
  implicit none

  integer, intent(in) ::  n
  double precision, intent(in) :: dist(n)
  double precision, intent(in) :: p
  integer :: choice, i

  do i=1, n
     if (dist(i) < p) then
        choice = i
        exit
     end if
  end do
end function choice


subroutine cdf(output, arr, n)
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: arr(n)
  integer :: i
  double precision :: output(n), mass(n)
  double precision :: prev, total

  total = dble(sum(arr))
  mass = dble(arr) / total
  prev = 0
  do i=1, n
     output(i) = mass(i) + prev
     prev = output(i)
  end do
end subroutine cdf

subroutine bipartite_edge_swap(B, A, x_degrees, y_degrees, nswap, max_tries, seed, m, n)

    implicit none

    ! Declarations
    integer, intent(in) :: m, n, nswap, max_tries, seed
    integer, intent(in) :: A(m,n), x_degrees(m), y_degrees(n)
    integer, intent(out) :: B(m,n)
    integer :: u, v, x1, y1, x2, y2, swapcount, iter, nth_one, choice
    double precision :: r(4), x_dist(m), y_dist(n)

    ! Initialize the permuted matrix
    B = A
    swapcount = 0
    iter = 0

    ! Convert the degree distributions into distributions
    call cdf(x_dist, x_degrees, m)
    call cdf(y_dist, y_degrees, n)

    ! Use the given random seed
    call seed_prng(seed)

    do while (swapcount < nswap .and. iter < max_tries)
       ! select six random numbers, enough for this iteration
       call random_number(r)

       ! select random nodes in each partition
       x1 = choice(x_dist, m, r(1))
       y1 = choice(y_dist, n, r(2))

       ! choose a random neighbor of each node
       u  = int(r(3) * x_degrees(x1))
       v  = int(r(4) * y_degrees(y1))
       
       y2 = nth_one(B, x1, 0, u, m, n)
       x2 = nth_one(B, y1, 1, v, m, n)

       ! If x1 is not already connected to y1 and
       ! x2 is not already connected to y2 then we will
       ! perform the edge swap
       if (B(x1, y1) .eq. 0 .and. B(x2, y2) .eq. 0) then
          swapcount = swapcount + 1 

          ! Update the adjacency matrix
          B(x1, y1) = 1
          B(x1, y2) = 0
          B(x2, y2) = 1
          B(x2, y1) = 0

       end if

       iter = iter + 1
       
    end do

    if (iter > max_tries) then
       print *, '[Warning] Bipartite edge swap: Maximum tries exceeded.'
    end if
end subroutine bipartite_edge_swap
