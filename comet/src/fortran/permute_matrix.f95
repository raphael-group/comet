subroutine nth_one(loc, arr, dim, t, num, m, n)
  implicit none

  integer, intent(in) :: num, m, n, dim, t
  integer, intent(in) :: arr(m, n)
  integer, intent(out) :: loc
  integer :: count, i, length

  count = 0
  i = 1

  if (t .eq. 0) then
     length = n
  else
     length = m
  end if
  
  do while (i <= length)
     ! If we are searching in a row
     if (t .eq. 0) then
        if (arr(dim, i) .eq. 1) then
           count = count + 1
        end if
     ! Otherwise we're searching in a column
     else
        if (arr(i, dim) .eq. 1) then
           count = count + 1
        end if
     end if
     ! Exit when we've found the location of the num-th one
     if (count .eq. num) then
        loc = i
        exit
     end if
     i = i + 1
  end do
end subroutine nth_one

subroutine bipartite_edge_swap(B, A, degrees, nswap, max_tries, m, n, num_nodes)

    implicit none

    integer, intent(in) :: m, n, num_nodes, nswap, max_tries
    integer, intent(in) :: A(m,n), degrees(num_nodes)
    integer, intent(out) :: B(m,n)
    integer :: u, v, x1, y1, x2, y2, swapcount, iter
    real :: r(4)

    ! Initialize the permuted
    B = A

    ! Use the given random seed
    ! TO-DO: fix this later. See http://goo.gl/mleKan
    ! call random_seed()

    swapcount = 0
    iter = 0

    do while (swapcount < nswap .and. iter < max_tries)
       ! select six random numbers, enough for this iteration
       call random_number(r)

       ! select random nodes in each partition
       x1 = int(1 + r(1) * m)
       y1 = int(1 + r(2) * n)

       ! choose a random neighbor of each node
       u  = int(1 + r(3) * degrees(x1))
       v  = int(1 + r(4) * degrees(m + y1))
       
       call nth_one(y2, B, x1, 0, u, m, n)
       call nth_one(x2, B, y1, 1, v, m, n)

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
