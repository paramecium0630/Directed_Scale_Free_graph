program DIRSF
    implicit none
    integer :: it
    integer, allocatable :: adj_matrix(:,:), indegree(:), outdegree(:)
    integer, allocatable :: iter_indeg(:), iter_outdeg(:)
    real(8) :: alpha, beta, gamma, delta_in, delta_out
    real(8) :: c1, c2, exponent_in, exponent_out
    real(8), allocatable :: pk_in(:), pk_out(:)
    real(8) :: time1, time2
    integer :: N_final, kmax_in, kmax_out, niter, m_links

    open(10, file="output/directed_scale_free_graph.txt", status="replace")
    open(11, file="output/degree.txt", status="replace")
    !open(12, file="output/indegree_distribution.txt", status="replace")
    !open(13, file="output/outdegree_distribution.txt", status="replace")

    niter = 20 ! number of iterations to average over

    ! Initialize parameters for directed scale-free graph
    ! Method A: terminate when n == N_final (vertex count is the clock).
    ! Edge count is a random variable; expected value ~ m_links*N_final/(alpha+gamma).
    N_final   = 20000   ! target number of vertices
    alpha     = 0.15d0
    beta      = 0.80d0
    gamma     = 0.05d0
    delta_in  = 1.0d0
    delta_out = 1.0d0
    m_links   = 10
    !!! test 
    c1 = real(m_links,8) * (alpha + beta) / (real(m_links,8) + delta_in*(alpha + gamma))
    c2 = real(m_links,8) * (beta + gamma) / (real(m_links,8) + delta_out*(alpha + gamma))
    exponent_in  = 1d0 + 1d0/c1
    exponent_out = 1d0 + 1d0/c2

    print *, "N_final=", N_final
    print *, "m_links=", m_links
    print *, "Expected edges ~", nint(real(m_links,8) * real(N_final,8) / (alpha + gamma))
    print *, "Network density ~ m/(alpha+gamma)", "= ", &
        real(m_links,8) * real(N_final,8) / (alpha + gamma) / real(N_final*(N_final-1),8)
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

        call grow_directed_scale_free(N_final, m_links, alpha, beta, gamma, delta_in, delta_out, &
            & indegree, outdegree, adj_matrix)

        call cpu_time(time2)
        print *, '  Vertices: ', N_final, '  Edges: ', sum(outdegree(1:N_final))
        print *, '  Generating time (seconds): ', time2 - time1

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

subroutine grow_directed_scale_free(N_final, m_links, alpha, beta, gamma, delta_in, delta_out, &
    & indeg, outdeg, adj_matrix)
!-----------------------------------------------------------------------
! Directed scale-free graph growth (no self-loops, no multiple edges).
!
! Method A: terminate when n == N_final (vertex count is the clock).
! Edge count is a random variable, expected ~ m_links*N_final/(alpha+gamma).
!
! At each discrete step we add m_links directed edges.
! With probability alpha:  (A) add new vertex v, add m_links edges v -> w,
!                              where w chosen ~ (d_in + delta_in)
! With probability beta:   (B) add m_links edges v -> w between existing vertices,
!                              v chosen ~ (d_out + delta_out),
!                              w chosen ~ (d_in  + delta_in), independently.
! With probability gamma:  (C) add new vertex w, add m_links edges v -> w,
!                              where v chosen ~ (d_out + delta_out)
! alpha + beta + gamma = 1, and delta_in, delta_out >= 0.
!
! Barabasi-Albert-like growth uses m_links>=1 where each new node connects
! to multiple existing nodes via preferential attachment.
!
! Convention: adj_matrix(i,j) = 1 means directed edge j -> i.
!
! INPUT:
!   N_final              : target number of vertices (termination condition)
!   m_links              : number of edges added per growth step
!   alpha,beta,gamma     : probabilities, must satisfy alpha+beta+gamma=1
!   delta_in,delta_out   : attractiveness parameters >= 0
!
! OUTPUT:
!   adj_matrix(i,j)      : adjacency matrix; (i,j)=1 means edge j -> i
!   outdeg(i), indeg(i)  : degree arrays consistent with adj_matrix
!
! NOTE: RNG must be initialized by caller via random_seed().
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in)    :: N_final, m_links
    real(8), intent(in)    :: alpha, beta, gamma, delta_in, delta_out
    integer, intent(inout) :: adj_matrix(N_final, N_final)
    integer, intent(inout) :: outdeg(N_final), indeg(N_final)

    integer :: n, v, w, newv
    real(8) :: r
    integer :: attempts, max_attempts

! ---------- sanity checks ----------
    if (N_final < 2)                                    stop "N_final must be >= 2"
    if (m_links < 1)                                    stop "m_links must be >= 1"
    if (alpha < 0d0 .or. beta < 0d0 .or. gamma < 0d0) stop "alpha,beta,gamma must be >= 0"
    if (abs((alpha+beta+gamma) - 1d0) > 1d-10)         stop "alpha+beta+gamma must equal 1"
    if (delta_in < 0d0 .or. delta_out < 0d0)           stop "delta_in and delta_out must both be >= 0"
    if (alpha + gamma <= 0d0)                           stop "alpha+gamma must be > 0 (otherwise no new vertices ever)"

! ---------- initialize: G0 = two vertices, one edge 1->2 ----------
    n = 2
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
            ! -------- (A): new vertex v, m_links edges newv -> w --------
                newv = n + 1
                n    = newv
                call add_edges_from_new_vertex(newv)
                exit

            else if (r < alpha + beta) then
            ! -------- (B): m_links edges between existing vertices --------
                call add_edges_between_existing(n)
                exit

            else
            ! -------- (C): new vertex w, m_links edges v -> newv --------
                newv = n + 1
                n    = newv
                call add_edges_to_new_vertex(newv)
                exit
            end if

        end do
    end do

contains

    subroutine add_edges_from_new_vertex(newv)
        integer, intent(in) :: newv
        integer :: added, w_try, tries, max_tries_local, m_eff

        m_eff = min(m_links, newv-1)
        added = 0
        tries = 0
        max_tries_local = max_attempts * max(1, m_eff)

        do while (added < m_eff .and. tries < max_tries_local)
            tries  = tries + 1
            w_try  = sample_vertex_by_indeg(newv-1, indeg, delta_in)
            if (adj_matrix(w_try, newv) == 1) cycle
            call add_edge_update(newv, w_try)
            added = added + 1
        end do
    end subroutine add_edges_from_new_vertex

    subroutine add_edges_to_new_vertex(newv)
        integer, intent(in) :: newv
        integer :: added, v_try, tries, max_tries_local, m_eff

        m_eff = min(m_links, newv-1)
        added = 0
        tries = 0
        max_tries_local = max_attempts * max(1, m_eff)

        do while (added < m_eff .and. tries < max_tries_local)
            tries  = tries + 1
            v_try  = sample_vertex_by_outdeg(newv-1, outdeg, delta_out)
            if (adj_matrix(newv, v_try) == 1) cycle
            call add_edge_update(v_try, newv)
            added = added + 1
        end do
    end subroutine add_edges_to_new_vertex

    subroutine add_edges_between_existing(n_now)
        integer, intent(in) :: n_now
        integer :: added, tries, max_tries_local

        added = 0
        tries = 0
        max_tries_local = max_attempts * max(1, m_links)

        do while (added < m_links .and. tries < max_tries_local)
            tries = tries + 1
            v = sample_vertex_by_outdeg(n_now, outdeg, delta_out)
            w = sample_vertex_by_indeg (n_now, indeg,  delta_in)
            if (v == w) cycle
            if (adj_matrix(w, v) == 1) cycle
            call add_edge_update(v, w)
            added = added + 1
        end do
    end subroutine add_edges_between_existing

    subroutine add_edge_update(v, w)
        ! Add directed edge v->w; update adj_matrix and degrees.
        ! Convention: adj_matrix(i,j) = 1 means edge j -> i.
        integer, intent(in) :: v, w
        adj_matrix(w, v) = 1
        outdeg(v) = outdeg(v) + 1
        indeg(w)  = indeg(w)  + 1
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
