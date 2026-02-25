program DIRSF
    implicit none
    integer :: i, j, it, n, m
    integer, allocatable :: adj_matrix(:,:), indegree(:), outdegree(:)
    integer, allocatable :: iter_indeg(:), iter_outdeg(:)
    real(8) :: alpha, beta, gamma, delta_in, delta_out
    real(8) :: c1, c2, exponent_in, exponent_out
    real(8), allocatable :: pk_in(:), pk_out(:)
    real(8) :: time1, time2
    integer :: N_final, kmax_in, kmax_out, niter

    open(10, file="output/directed_scale_free_graph.txt", status="replace")
    open(11, file="output/degree.txt", status="replace")
    !open(12, file="output/indegree_distribution.txt", status="replace")
    !open(13, file="output/outdegree_distribution.txt", status="replace")

    niter = 20 ! number of iterations to average over

    ! Initialize parameters for directed scale-free graph
    ! Method A: terminate when n == N_final (vertex count is the clock).
    ! Edge count m is a random variable; expected value ~ N_final / (alpha+gamma).
    N_final   = 30000   ! target number of vertices
    alpha     = 0.1d0
    beta      = 0.8d0
    gamma     = 0.1d0
    delta_in  = 1.0d0
    delta_out = 1.0d0

    c1 = (alpha + beta) / (1d0 + delta_in*(alpha + gamma))
    c2 = (beta + gamma) / (1d0 + delta_out*(alpha + gamma))
    exponent_in  = 1d0 + 1d0/c1
    exponent_out = 1d0 + 1d0/c2

    print *, "N_final=", N_final
    print *, "Expected edges ~", nint(real(N_final) / (alpha + gamma))
    print *, "Network density ~ 1/(alpha+gamma)", "= ", real(N_final) / (alpha + gamma) / real(N_final*(N_final-1))
    print *, "alpha=", alpha
    print *, "beta=", beta
    print *, "gamma=", gamma
    print *, "delta_in=", delta_in
    print *, "delta_out=", delta_out
    print *, "Expected in-degree  exponent=", exponent_in
    print *, "Expected out-degree exponent=", exponent_out
    print *, repeat("-", 50)

    allocate(adj_matrix(N_final, N_final))
    allocate(indegree(N_final), outdegree(N_final))

    ! iter arrays: each iteration contributes exactly N_final vertices
    allocate(iter_indeg(N_final*niter))
    allocate(iter_outdeg(N_final*niter))

    call random_seed() ! initialize RNG

    do it = 1, niter

        print *, "Iteration ", it
        call cpu_time(time1)

        call grow_directed_scale_free(N_final, alpha, beta, gamma, delta_in, delta_out, &
            & n, m, indegree, outdegree, adj_matrix)

        call cpu_time(time2)
        print *, '  Vertices: ', n, '  Edges: ', m
        print *, '  Generating time (seconds): ', time2 - time1

        ! n is guaranteed == N_final (method A terminates on vertex count)
        iter_indeg((it-1)*N_final+1 : it*N_final) = indegree(1:N_final)
        iter_outdeg((it-1)*N_final+1 : it*N_final) = outdegree(1:N_final)

    end do

    kmax_in  = maxval(iter_indeg)
    kmax_out = maxval(iter_outdeg)

    print *, "Max in-degree: ",  kmax_in
    print *, "Max out-degree: ", kmax_out

    allocate(pk_in(0:kmax_in), pk_out(0:kmax_out))

    call degree_distribution(N_final*niter, iter_indeg, iter_outdeg, pk_in, pk_out, kmax_in, kmax_out)    

    close(10)
    close(11)

end program DIRSF

!===========================================================================

subroutine degree_distribution(N, indegree, outdegree, pk_in, pk_out, kmax_in, kmax_out)

    implicit none
    integer, intent(in)  :: N                           ! total samples = N_final * niter
    integer, intent(in)  :: indegree(N), outdegree(N)
    integer, intent(in) :: kmax_in, kmax_out
    real(8), intent(inout) :: pk_in(0:kmax_in), pk_out(0:kmax_out)    

    integer :: i
    integer, allocatable :: hist_in(:), hist_out(:)

    open(12, file="output/indegree_distribution.txt",  status="replace")
    open(13, file="output/outdegree_distribution.txt", status="replace")

    allocate(hist_in (0:kmax_in))
    allocate(hist_out(0:kmax_out))
    hist_in  = 0
    hist_out = 0

    !--------------------------------------
    ! 2. Count
    !--------------------------------------
    do i = 1, N
        hist_in(indegree(i))  = hist_in(indegree(i))  + 1
        hist_out(outdegree(i)) = hist_out(outdegree(i)) + 1
    end do

    !--------------------------------------
    ! 3. Normalize by total number of samples (N_final * niter)
    !--------------------------------------

    do i = 0, kmax_in
        pk_in(i) = dble(hist_in(i)) / dble(N)
    end do
    do i = 0, kmax_out
        pk_out(i) = dble(hist_out(i)) / dble(N)
    end do

    do i = 1, kmax_in
        write(12, *) i, pk_in(i)
    end do
    do i = 1, kmax_out
        write(13, *) i, pk_out(i)
    end do

    deallocate(hist_in, hist_out)
    close(12)
    close(13)

end subroutine degree_distribution

!===========================================================================

subroutine grow_directed_scale_free(N_final, alpha, beta, gamma, delta_in, delta_out, &
    & n, t, indeg, outdeg, adj_matrix)
!-----------------------------------------------------------------------
! Directed scale-free graph growth (no self-loops, no multiple edges).
!
! Method A: terminate when n == N_final (vertex count is the clock).
! Edge count t is a random variable, expected ~ N_final / (alpha+gamma).
!
! At each discrete step we add exactly ONE directed edge.
! With probability alpha:  (A) add new vertex v, add edge v -> w,
!                              where w chosen ~ (d_in + delta_in)
! With probability beta:   (B) add edge v -> w between existing vertices,
!                              v chosen ~ (d_out + delta_out),
!                              w chosen ~ (d_in  + delta_in), independently.
! With probability gamma:  (C) add new vertex w, add edge v -> w,
!                              where v chosen ~ (d_out + delta_out)
! alpha + beta + gamma = 1, and delta_in, delta_out >= 0.
!
! Barabasi-Albert model (m=1) is a special case with
!   beta=gamma=delta_out=0, alpha=delta_in=1  (per Bollobas et al. 2003).
!
! Convention: adj_matrix(i,j) = 1 means directed edge j -> i.
!
! INPUT:
!   N_final              : target number of vertices (termination condition)
!   alpha,beta,gamma     : probabilities, must satisfy alpha+beta+gamma=1
!   delta_in,delta_out   : attractiveness parameters >= 0
!
! OUTPUT:
!   adj_matrix(i,j)      : adjacency matrix; (i,j)=1 means edge j -> i
!   outdeg(i), indeg(i)  : degree arrays consistent with adj_matrix
!   n                    : final vertex count (== N_final)
!   t                    : final edge count (random)
!
! NOTE: RNG must be initialized by caller via random_seed().
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in)    :: N_final
    real(8), intent(in)    :: alpha, beta, gamma, delta_in, delta_out
    integer, intent(inout) :: n, t
    integer, intent(inout) :: adj_matrix(N_final, N_final)
    integer, intent(inout) :: outdeg(N_final), indeg(N_final)

    integer :: v, w, newv
    real(8) :: r
    integer :: attempts, max_attempts

! ---------- sanity checks ----------
    if (N_final < 2)                                    stop "N_final must be >= 2"
    if (alpha < 0d0 .or. beta < 0d0 .or. gamma < 0d0) stop "alpha,beta,gamma must be >= 0"
    if (abs((alpha+beta+gamma) - 1d0) > 1d-10)         stop "alpha+beta+gamma must equal 1"
    if (delta_in < 0d0 .or. delta_out < 0d0)           stop "delta_in and delta_out must both be >= 0"
    if (alpha + gamma <= 0d0)                           stop "alpha+gamma must be > 0 (otherwise no new vertices ever)"

! ---------- initialize: G0 = two vertices, one edge 1->2 ----------
    n = 2
    t = 1
    adj_matrix       = 0
    indeg            = 0
    outdeg           = 0
    adj_matrix(2, 1) = 1   ! edge 1->2: stored at A(2,1) per convention
    outdeg(1)        = 1
    indeg(2)         = 1
    max_attempts     = 1000

! ---------- growth loop: run until n == N_final ----------
    do while (n < N_final)
        attempts = 0
        do
            attempts = attempts + 1
            if (attempts > max_attempts) then
                print *, "Warning: too many resampling attempts in step (B)."
                exit
            end if

            call random_number(r)

            if (r < alpha) then
            ! -------- (A): new vertex v, edge v -> w --------
                newv = n + 1
                w    = sample_vertex_by_indeg(n, indeg, delta_in)
                n    = newv
                call add_edge_update(newv, w)
                exit

            else if (r < alpha + beta) then
            ! -------- (B): edge between existing vertices --------
                v = sample_vertex_by_outdeg(n, outdeg, delta_out)
                w = sample_vertex_by_indeg (n, indeg,  delta_in)
                if (v == w)                cycle   ! no self-loops
                if (adj_matrix(w, v) == 1) cycle   ! no multiple edges
                call add_edge_update(v, w)
                exit

            else
            ! -------- (C): new vertex w, edge v -> w --------
                newv = n + 1
                v    = sample_vertex_by_outdeg(n, outdeg, delta_out)
                n    = newv
                call add_edge_update(v, newv)
                exit
            end if

        end do
    end do

contains

    subroutine add_edge_update(v, w)
        ! Add directed edge v->w; update adj_matrix, degrees, and edge count t.
        ! Convention: adj_matrix(i,j) = 1 means edge j -> i.
        integer, intent(in) :: v, w
        adj_matrix(w, v) = 1
        outdeg(v) = outdeg(v) + 1
        indeg(w)  = indeg(w)  + 1
        t = t + 1
    end subroutine add_edge_update

    integer function sample_vertex_by_indeg(n, indeg, delta_in) result(idx)
        ! Sample a vertex with probability proportional to (d_in + delta_in).
        integer, intent(in) :: n
        integer, intent(in) :: indeg(n)
        real(8), intent(in) :: delta_in
        integer :: i
        real(8) :: total, r, cumulative

        total = 0d0
        do i = 1, n
            total = total + indeg(i) + delta_in
        end do

        call random_number(r)
        r = r * total

        cumulative = 0d0
        do i = 1, n
            cumulative = cumulative + indeg(i) + delta_in
            if (r < cumulative) then
                idx = i
                return
            end if
        end do
        idx = n   ! fallback (floating-point rounding guard)
    end function sample_vertex_by_indeg

    integer function sample_vertex_by_outdeg(n, outdeg, delta_out) result(idx)
        ! Sample a vertex with probability proportional to (d_out + delta_out).
        integer, intent(in) :: n
        integer, intent(in) :: outdeg(n)
        real(8), intent(in) :: delta_out
        integer :: i
        real(8) :: total, r, cumulative

        total = 0d0
        do i = 1, n
            total = total + outdeg(i) + delta_out
        end do

        call random_number(r)
        r = r * total

        cumulative = 0d0
        do i = 1, n
            cumulative = cumulative + outdeg(i) + delta_out
            if (r < cumulative) then
                idx = i
                return
            end if
        end do
        idx = n   ! fallback (floating-point rounding guard)
    end function sample_vertex_by_outdeg

end subroutine grow_directed_scale_free