program DIRSF
    implicit none
    integer :: i, j, it, n, m
    integer, allocatable :: adj_matrix(:,:), indegree(:), outdegree(:)
    integer, allocatable :: iter_indeg(:), iter_outdeg(:)
    real(8) :: alpha, beta, gamma, delta_in, delta_out
    real(8) :: density, c1, c2, exponent_in, exponent_out
    real(8), allocatable :: pk_in(:), pk_out(:)
    real(8) :: time1, time2
    integer :: N_final, E_final, kmax_in, kmax_out, niter

    open(10, file="output/directed_scale_free_graph.txt", status="replace")
    open(11, file="output/degree.txt", status="replace")
    !open(12, file="output/indegree_distribution.txt", status="replace")
    !open(13, file="output/outdegree_distribution.txt", status="replace")

    niter = 50 ! number of iterations to average over

    ! Initialize parameters for directed scale-free graph
    N_final = 10000 ! target number of vertices
    E_final = 200000 ! target number of edges
    alpha = 1d0/3
    beta = 1d0/3
    gamma = 1d0/3
    delta_in = 1.0d0
    delta_out = 1.0d0

    c1 = (alpha + beta) / (1 + delta_in*(alpha + gamma))
    c2 = (beta + gamma) / (1 + delta_out*(alpha + gamma))
    exponent_in = 1 + 1d0/c1
    exponent_out = 1 + 1d0/c2
    density = real(E_final) / (real(N_final) * real(N_final - 1))

    print *, "N_final=", N_final
    print *, "E_final=", E_final
    print *, "Density=", density
    print *, "alpha=", alpha
    print *, "beta=", beta
    print *, "gamma=", gamma
    print *, "delta_in=", delta_in
    print *, "delta_out=", delta_out
    print *, "Expected in-degree exponent=", exponent_in
    print *, "Expected out-degree exponent=", exponent_out
    print *, repeat("-", 50)

    allocate(adj_matrix(N_final, N_final))    
    allocate(indegree(N_final), outdegree(N_final))

    allocate(pk_in(N_final), pk_out(N_final))
    allocate(iter_indeg(N_final*niter), iter_outdeg(N_final*niter)) ! for collecting degree data across iterations

    call random_seed() ! initialize RNG

    do it = 1, niter

    print *, "Iteration ", it
    call cpu_time(time1)
    
    call grow_directed_scale_free(N_final, E_final, alpha, beta, gamma, delta_in, delta_out, &
        & n, m, indegree, outdegree, adj_matrix)
    call cpu_time(time2)
    print*, 'Generating time (seconds): ', time2 - time1

    iter_indeg((it-1)*N_final+1:it*N_final) = indegree(1:N_final)
    iter_outdeg((it-1)*N_final+1:it*N_final) = outdegree(1:N_final)

    enddo

    call degree_distribution(N_final, iter_indeg, iter_outdeg, pk_in, pk_out, kmax_in, kmax_out)

    print*, "Max in-degree: ", kmax_in
    print*, "Max out-degree: ", kmax_out

    ! call cpu_time(time1)
    ! ! Write adjacency matrix to file
    ! do i = 1, n
    !     do j = 1, n
    !         if (adj_matrix(i, j) == 1) then
    !             write(10, *) i, j
    !         end if
    !     end do
    ! end do
    ! call cpu_time(time2)
    ! print*, 'write time (seconds): ', time2 - time1

    ! ! Write degree information to file
    ! do i = 1, n
    !     write(11, '(I0, 1X, I0)') indegree(i), outdegree(i)
    ! end do

    close(10)
    close(11)
    close(12)
    close(13)

end program DIRSF

subroutine degree_distribution(N, indegree, outdegree, pk_in, pk_out, kmax_in, kmax_out)

    implicit none
    integer, intent(in) :: N
    integer, intent(in) :: indegree(N), outdegree(N)

    real(8), allocatable, intent(out) :: pk_in(:), pk_out(:)
    integer, intent(out) :: kmax_in, kmax_out

    integer :: i
    integer, allocatable :: hist_in(:), hist_out(:)

    open(12, file="output/indegree_distribution.txt", status="replace")
    open(13, file="output/outdegree_distribution.txt", status="replace")

    !--------------------------------------
    ! 1. Max degree
    !--------------------------------------
    kmax_in  = maxval(indegree)
    kmax_out = maxval(outdegree)

    allocate(hist_in(0:kmax_in))
    allocate(hist_out(0:kmax_out))
    hist_in  = 0
    hist_out = 0

    !--------------------------------------
    ! 2. Count
    !--------------------------------------
    do i = 1, N
        hist_in(indegree(i))   = hist_in(indegree(i))   + 1
        hist_out(outdegree(i)) = hist_out(outdegree(i)) + 1
    end do

    !--------------------------------------
    ! 3. Normalize
    !--------------------------------------
    allocate(pk_in(0:kmax_in))
    allocate(pk_out(0:kmax_out))

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

end subroutine degree_distribution

subroutine grow_directed_scale_free(N_final, E_final, alpha, beta, gamma, delta_in, delta_out, &
    & n, t, indeg, outdeg, adj_matrix)
!-----------------------------------------------------------------------
! Directed scale-free graph growth (no self-loops, no multiple edges).
!
! At each discrete step we add exactly ONE directed edge.
! With probability alpha:  (A) add new vertex v, add edge v -> w,
!                           where w chosen ~ (d_in + delta_in)
! With probability beta:   (B) add edge v -> w between existing vertices,
!                           v chosen ~ (d_out + delta_out),
!                           w chosen ~ (d_in  + delta_in), independently.
! With probability gamma:  (C) add new vertex w, add edge v -> w,
!                           where v chosen ~ (d_out + delta_out)
! alpha + beta + gamma = 1, and delta_in, delta_out >= 0 are tunable parameters.
!
! Barabási–Albert model is a special case with alpha=gamma=0, beta=1, and delta_in=delta_out=1.
!
! INPUT:
!   N_final, E_final     : target numbers of vertices and edges
!   alpha,beta,gamma     : probabilities, must satisfy alpha+beta+gamma=1
!   delta_in,delta_out   : attractiveness parameters >=0
!
! OUTPUT data structures:
!   A(i,j)   : adjacency matrix (0/1) for directed edges j -> i
!   outdeg(i), indeg(i): degrees consistent with A
!   n, t     : current numbers of vertices and edges
!
! IMPORTANT:
!   - This routine assumes RNG is initialized elsewhere (random_seed).
!   - Complexity: sampling is O(n) per draw; (B) may resample.
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: N_final, E_final ! target number of vertices and edges
    real(8), intent(in) :: alpha, beta, gamma, delta_in, delta_out ! parameters for preferential attachment
    integer, intent(inout) :: n, t ! current number of vertices and edges
    integer, intent(inout) :: adj_matrix(N_final, N_final)     ! 0/1 adjacency matrix
    integer, intent(inout) :: outdeg(N_final), indeg(N_final) ! degree arrays

    integer :: i, j, v, w, newv ! vertex indices
    real(8) :: r ! random number for choosing operation
    integer :: attempts, max_attempts ! for resampling in (B)

! ---------- sanity checks ----------
    if (N_final < 2) stop "N_final must be >= 2"
    if (E_final < 1) stop "E_final must be >= 1"
    if (alpha < 0d0 .or. beta < 0d0 .or. gamma < 0d0) stop "alpha,beta,gamma must be >=0"
    if (abs((alpha+beta+gamma)-1d0) > 1d-10) stop "alpha+beta+gamma must be 1"
    if (delta_in < 0d0 .and. delta_out < 0d0) stop "either delta_in or delta_out must be >=0"

! We start from G0: one vertex, no edges (N0=1, t0=0).
! Since each new vertex comes with one new edge, need enough edges to create vertices.
    if (E_final < (N_final - 1)) then
        stop "Impossible targets: need E_final >= N_final-1 (starting from 1 vertex)."
    end if

! ---------- initialize ----------
!   1 -> 2
    n = 2 ! start with one vertex
    t = 1 ! no edges yet
    adj_matrix = 0 ! clear adjacency matrix
    outdeg = 0
    indeg = 0
    adj_matrix(2, 1) = 1 ! add edge 1 -> 2
    outdeg(1) = 1
    indeg(2) = 1
    max_attempts = 1000

! ---------- growth loop ----------
    do while (n < N_final .or. t < E_final)
        attempts = 0
        do 
            attempts = attempts + 1
            if (attempts > max_attempts) then
                print *, "Warning: too many attempts to sample edge in (B). Consider increasing max_attempts."
                exit
            end if
            call random_number(r)

            ! If reached max vertices → force type (B)
            if (n >= N_final) then

                v = sample_vertex_by_outdeg(n, t, outdeg, delta_out) ! sample source vertex by out-degree
                w = sample_vertex_by_indeg(n, t, indeg,  delta_in) ! sample target vertex by in-degree

                if (v == w) cycle ! no self-loops
                if (adj_matrix(w, v) == 1) cycle ! no multiple edges

                call add_edge_update(v, w) ! add edge v->w and update degrees and edge count
                exit
            end if

            if (r < alpha) then
            ! -------- (A) --------
            newv = n + 1 ! new vertex index
            w = sample_vertex_by_indeg(n, t, indeg, delta_in)

            n = newv
            call add_edge_update(newv, w)
            exit

            else if (r < alpha + beta) then
            ! -------- (B) --------
            v = sample_vertex_by_outdeg(n, t, outdeg, delta_out)
            w = sample_vertex_by_indeg(n, t, indeg,  delta_in)

            if (v == w) cycle
            if (adj_matrix(w, v) == 1) cycle

            call add_edge_update(v, w)
            exit

            else
            ! -------- (C) --------
            newv = n + 1
            v = sample_vertex_by_outdeg(n, t, outdeg, delta_out)

            n = newv
            call add_edge_update(v, newv)
            exit
            end if

        end do
    end do

contains
    
    subroutine add_edge_update(v, w)
        integer, intent(in) :: v, w
        ! Add directed edge v->w and update adjacency matrix, degree counts, and edge count
        adj_matrix(w, v) = 1 ! note: A(i,j)=1 means edge j->i
        outdeg(v) = outdeg(v) + 1
        indeg(w) = indeg(w) + 1
        t = t + 1
    end subroutine add_edge_update

    integer function sample_vertex_by_indeg(n, t, indeg, delta_in) result(idx)
        integer, intent(in) :: n, t
        integer, intent(in) :: indeg(n)
        real(8), intent(in) :: delta_in
        integer :: i
        real(8) :: total_attractiveness, r, cumulative

        ! Sample a vertex with probability proportional to (d_in + delta_in).
        total_attractiveness = 0d0
        do i = 1, n
            total_attractiveness = total_attractiveness + indeg(i) + delta_in
        end do

        call random_number(r)
        r = r * total_attractiveness
        
        ! Walk through vertices to find which one corresponds to the random number
        cumulative = 0d0
        do i = 1, n
            cumulative = cumulative + indeg(i) + delta_in
            if (r < cumulative) then
                idx = i
                return
            end if
        end do

        idx = n ! fallback (should not happen if total_attractiveness > 0)
    end function sample_vertex_by_indeg

    integer function sample_vertex_by_outdeg(n, t, outdeg, delta_out) result(idx)
        integer, intent(in) :: n, t
        integer, intent(in) :: outdeg(n)
        real(8), intent(in) :: delta_out
        integer :: i
        real(8) :: total_attractiveness, r, cumulative

        total_attractiveness = 0d0
        do i = 1, n
            total_attractiveness = total_attractiveness + outdeg(i) + delta_out
        end do

        call random_number(r)
        r = r * total_attractiveness

        cumulative = 0d0
        do i = 1, n
            cumulative = cumulative + outdeg(i) + delta_out
            if (r < cumulative) then
                idx = i
                return
            end if
        end do

        idx = n ! fallback (should not happen if total_attractiveness > 0)
    end function sample_vertex_by_outdeg

end subroutine grow_directed_scale_free
