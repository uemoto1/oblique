subroutine calc_oblique
    use input_parameter
    use em_field
    use phys_constants
    use math_constants
    implicit none
    real(8), allocatable :: Ac_old(:, :)
    real(8), allocatable :: Ac_cur(:, :)
    real(8), allocatable :: Ac_new(:, :)
    real(8), allocatable :: J_old(:, :)
    real(8), allocatable :: J_cur(:, :)
    real(8), allocatable :: J_new(:, :)
    real(8), allocatable :: P_old(:, :)
    real(8), allocatable :: P_cur(:, :)
    real(8), allocatable :: P_new(:, :)
    real(8) :: theta, rmatrix(3, 3)

    real(8) :: cp1 
    character(256) :: file_ac
    integer :: it, iz

    real(8) :: c1, c2, c3

    theta = (pi / 180.0d0) * theta_deg

    allocate(Ac_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(P_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(P_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(P_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))

    ! Initial electromagnetic field
    call calc_Ac_ext_t(-dt, -hz_m / cspeed_au * cos(theta), -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_old)
    call calc_Ac_ext_t(0.0d0, -hz_m / cspeed_au * cos(theta), -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_cur)

    rmatrix(1, 1:3) = (/ cos(theta), 0.0d0, sin(theta)   /)
    rmatrix(2, 1:3) = (/ 0.0d0, 1.0d0, 0.0d0 /)
    rmatrix(3, 1:3) = (/ -sin(theta), 0.0d0, cos(theta)   /)

    Ac_old(:, :) = matmul(rmatrix, Ac_old(:, :))
    Ac_cur(:, :) = matmul(rmatrix, Ac_cur(:, :))

    Ac_new = 0.0d0
    J_new = 0.0d0
    J_cur = 0.0d0
    J_old = 0.0d0
    P_new = 0.0d0
    P_cur = 0.0d0
    P_old = 0.0d0

    c1 = (2.0d0 - (omega0 * dt) ** 2) / (1.0d0 + gamma * dt * 0.5d0)
    c2 = (1 - gamma * dt * 0.5d0) / (1.0d0 + gamma * dt * 0.5d0)
    c3 = alpha / (1.0d0 + gamma * dt * 0.5d0)

    do it = 0, nt
        Ac_new = 0.0d0

        do iz = -nzvacl_m, nz_m+nzvacr_m
             Ac_new(2, iz) = 2.0d0 * Ac_cur(2, iz) - Ac_old(2, iz) &
                & + (cspeed_au * dt / cos(theta)) ** 2 * (Ac_cur(2, iz + 1) - 2 * Ac_cur(2, iz) + Ac_cur(2, iz - 1)) / hz_m ** 2 &
                & + 4.0 * pi / (cos(theta) ** 2) * J_cur(2, iz) * dt ** 2
        end do

        do iz = 1, nz_m
            J_new(1:3, iz) = c1 * J_cur(1:3, iz) - c2 * J_old(1:3, iz) &
                & - c3 * (Ac_new(1:3, iz) - 2.0d0 * Ac_cur(1:3, iz) + Ac_old(1:3, iz))
            P_new(1:3, iz) = P_old(1:3, iz) + 2.0d0 * dt * J_cur(1:3, iz)
        end do


        Ac_old = Ac_cur
        Ac_cur = Ac_new
        J_old = J_cur
        J_cur = J_new
        P_old = P_cur
        P_cur = P_new

        if (mod(it, 100) == 0) then
            write(file_ac, "('Ac_', i6.6, '.txt')") it
            write(*, "(a)") trim(file_ac)
            open(10, file=trim(file_ac), action="write")
            do iz = -nzvacl_m, nz_m+nzvacr_m
                write(10, "(4es20.10e4)") iz * hz_m, ac_cur(:, iz)
            end do
            close(10)
        end if
    end do

    stop
end subroutine calc_oblique





subroutine calc_normal
    use input_parameter
    use em_field
    use phys_constants
    use math_constants
    implicit none
    real(8), allocatable :: Ac_old(:, :)
    real(8), allocatable :: Ac_cur(:, :)
    real(8), allocatable :: Ac_new(:, :)
    real(8), allocatable :: J_old(:, :)
    real(8), allocatable :: J_cur(:, :)
    real(8), allocatable :: J_new(:, :)

    real(8) :: cp1 
    character(256) :: file_ac
    integer :: it, iz

    real(8) :: c1, c2, c3

    allocate(Ac_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))

    ! Initial electromagnetic field
    call calc_Ac_ext_t(-dt, -hz_m / cspeed_au, -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_old)
    call calc_Ac_ext_t(0.0d0, -hz_m / cspeed_au, -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_cur)
    Ac_new = 0.0d0
    J_new = 0.0d0
    J_cur = 0.0d0
    J_old = 0.0d0

    c1 = (2.0d0 - (omega0 * dt) ** 2) / (1.0d0 + gamma * dt * 0.5d0)
    c2 = (1 - gamma * dt * 0.5d0) / (1.0d0 + gamma * dt * 0.5d0)
    c3 = alpha / (1.0d0 + gamma * dt * 0.5d0)

    do it = 0, nt
        do iz = -nzvacl_m, nz_m+nzvacr_m
             Ac_new(1, iz) = 2.0d0 * Ac_cur(1, iz) - Ac_old(1, iz) &
                & + dt ** 2 * cspeed_au ** 2 * (Ac_cur(1, iz+1) - 2.0d0 * Ac_cur(1, iz) + Ac_cur(1, iz-1)) / hz_m ** 2 &
                & + 4.0 * pi * J_cur(1, iz) * dt ** 2
             Ac_new(2, iz) = 2.0d0 * Ac_cur(2, iz) - Ac_old(2, iz) &
                & + dt ** 2 * cspeed_au ** 2 * (Ac_cur(2, iz+1) - 2.0d0 * Ac_cur(2, iz) + Ac_cur(2, iz-1)) / hz_m ** 2 &
                & + 4.0 * pi * J_cur(2, iz) * dt ** 2
        end do

        do iz = 1, nz_m
            J_new(1:3, iz) = c1 * J_cur(1:3, iz) - c2 * J_old(1:3, iz) &
                & - c3 * (Ac_new(1:3, iz) - 2.0d0 * Ac_cur(1:3, iz) + Ac_old(1:3, iz))
        end do


        Ac_old = Ac_cur
        Ac_cur = Ac_new
        J_old = J_cur
        J_cur = J_new

        if (mod(it, 100) == 0) then
            write(file_ac, "('Ac_', i6.6, '.txt')") it
            write(*, "(a)") trim(file_ac)
            open(10, file=trim(file_ac), action="write")
            do iz = -nzvacl_m, nz_m+nzvacr_m
                write(10, "(4es20.10e4)") iz * hz_m, ac_cur(:, iz)
            end do
            close(10)
        end if
    end do

    stop
end subroutine calc_normal










program main
    use input_parameter
    implicit none
    call read_input()

    call calc_oblique
    stop
end program main
