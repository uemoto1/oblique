subroutine calc_oblique
    use input_parameter
    use em_field
    use phys_constants
    use math_constants
    implicit none
    real(8), allocatable :: Ac_old(:, :)
    real(8), allocatable :: Ac_cur(:, :)
    real(8), allocatable :: Ac_new(:, :)
    real(8), allocatable :: Jld_old(:, :)
    real(8), allocatable :: Jld_cur(:, :)
    real(8), allocatable :: Jld_new(:, :)
    real(8), allocatable :: Pld_old(:, :)
    real(8), allocatable :: Pld_cur(:, :)
    real(8), allocatable :: Pld_new(:, :)
    real(8), allocatable :: J_cur(:, :)
    real(8), allocatable :: P_cur(:, :)
    real(8), allocatable :: w(:)
    real(8) :: theta, rmatrix(3, 3)

    character(256) :: file_ac
    integer :: it, iz

    real(8) :: c1, c2, c3

    integer :: ns

    theta = (pi / 180.0d0) * theta_oblique_deg
    ns = n_smooth_oblique

    allocate(Ac_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Pld_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Pld_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Pld_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(J_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(P_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(w(-nzvacl_m-1: nz_m+nzvacr_m+1))

    ! Initial electromagnetic field
    call calc_Ac_ext_t(-dt, -hz_m / cspeed_au * cos(theta), -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_old)
    call calc_Ac_ext_t(0.0d0, -hz_m / cspeed_au * cos(theta), -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_cur)

    ! rmatrix(1, 1:3) = (/ cos(theta), 0.0d0, sin(theta)   /)
    ! rmatrix(2, 1:3) = (/ 0.0d0, 1.0d0, 0.0d0 /)
    ! rmatrix(3, 1:3) = (/ -sin(theta), 0.0d0, cos(theta)   /)

    ! Ac_old(:, :) = matmul(rmatrix, Ac_old(:, :))
    ! Ac_cur(:, :) = matmul(rmatrix, Ac_cur(:, :))

    w = 0.0d0
    do iz = 1, nz_m
        if (iz < ns) then
            w(iz) = weight((dble(iz) - dble(0.5)) / dble(ns))
        else if (nz_m - ns - 1 < iz) then
            w(iz) = weight((dble(nz_m) - dble(iz) - dble(0.5)) / dble(ns))
        else
            w(iz) = 1.0d0
        end if
    end do

    do iz = lbound(w,1), ubound(w,1)
        write(10,*) w(iz)
    end do
    

    Ac_new = 0.0d0
    Jld_new = 0.0d0
    Jld_cur = 0.0d0
    Jld_old = 0.0d0
    Pld_new = 0.0d0
    Pld_cur = 0.0d0
    Pld_old = 0.0d0
    P_cur = 0.0
    J_cur = 0.0

    c1 = (2.0d0 - (omega0 * dt) ** 2) / (1.0d0 + gamma * dt * 0.5d0)
    c2 = (1 - gamma * dt * 0.5d0) / (1.0d0 + gamma * dt * 0.5d0)
    c3 = alpha / (1.0d0 + gamma * dt * 0.5d0)

    do it = 0, nt

        do iz = 1, nz_m
            P_cur(1:3, iz) = w(iz) * Pld_cur(1:3, iz)
            J_cur(1:3, iz) = w(iz) * Jld_cur(1:3, iz)
        end do

        Ac_new = 0.0d0

        do iz = -nzvacl_m, nz_m+nzvacr_m
            Ac_new(1, iz) = 2.0d0 * Ac_cur(1, iz) - Ac_old(1, iz) &
            & + (cspeed_au * dt / cos(theta)) ** 2 * (Ac_cur(1, iz + 1) - 2 * Ac_cur(1, iz) + Ac_cur(1, iz - 1)) / hz_m ** 2 &
            & + 4.0d0 * pi * dt ** 2 * J_cur(1, iz) &
            & + 4.0d0 * pi * sin(theta) / cos(theta) ** 2 * cspeed_au * dt ** 2 &
            & * (P_cur(3, iz+1) - P_cur(3, iz-1)) / (2.0d0 * hz_m)
            
            Ac_new(2, iz) = 2.0d0 * Ac_cur(2, iz) - Ac_old(2, iz) &
            & + (cspeed_au * dt / cos(theta)) ** 2 * (Ac_cur(2, iz + 1) - 2 * Ac_cur(2, iz) + Ac_cur(2, iz - 1)) / hz_m ** 2 &
            & + 4.0 * pi / (cos(theta) ** 2) * J_cur(2, iz) * dt ** 2
            
            Ac_new(3, iz) = Ac_old(3, iz) &
            & + 2.0d0 * cspeed_au * dt * sin(theta) / (cos(theta) ** 2) * (ac_cur(1, iz+1) - ac_cur(1, iz-1)) / (2.0d0 * hz_m) &
            & + 4.0d0 * pi * 2.0d0 * dt / (cos(theta) ** 2) * P_cur(3, iz) 
        end do


    !     Ac_new_ms(3,ix_m,iy_m,iz_m) = 2*Ac_ms(3,ix_m,iy_m,iz_m) - Ac_old_ms(3,ix_m,iy_m,iz_m) &
    !     & + (c_light*dt/cos_oblique)**2 &
    !     & * (Ac_ms(3,ix_m+1,iy_m,iz_m)-2*Ac_ms(3,ix_m,iy_m,iz_m)+Ac_ms(3,ix_m-1,iy_m,iz_m))/Hx_m**2 &
    !     & - 4*pi*dt**2*Jm_ms(3,ix_m,iy_m,iz_m) &
    !     & - 4*pi*sin_oblique/cos_oblique**2*c_light*dt**2 &
    !     & * (Pm_ms(1,ix_m+1,iy_m,iz_m)-Pm_ms(1,ix_m-1,iy_m,iz_m))/(2*Hx_m)
    !   Ac_new_ms(1,ix_m,iy_m,iz_m) = Ac_old_ms(1,ix_m,iy_m,iz_m) &
    !     & + 2*dt*c_light*sin_oblique/cos_oblique**2*(Ac_ms(3,ix_m+1,iy_m,iz_m)-Ac_ms(3,ix_m-1,iy_m,iz_m))/(2*Hx_m) &
    !     & - 4*pi*(2*dt)/cos_oblique**2*Pm_ms(1,ix_m,iy_m,iz_m)
  

        do iz = 1, nz_m
            Jld_new(1:3, iz) = c1 * Jld_cur(1:3, iz) - c2 * Jld_old(1:3, iz) &
                & - c3 * (Ac_new(1:3, iz) - 2.0d0 * Ac_cur(1:3, iz) + Ac_old(1:3, iz))
            Pld_new(1:3, iz) = Pld_old(1:3, iz) + 2.0d0 * dt * Jld_cur(1:3, iz)
        end do

        Ac_old = Ac_cur
        Ac_cur = Ac_new
        Jld_old = Jld_cur
        Jld_cur = Jld_new
        Pld_old = Pld_cur
        Pld_cur = Pld_new

        if (mod(it, 100) == 0) then
            write(file_ac, "('Ac_', i6.6, '.txt')") it
            open(10, file=trim(file_ac), action="write")
            do iz = -nzvacl_m, nz_m+nzvacr_m
                write(10, "(4es20.10e3)") iz * hz_m, ac_cur(:, iz)
            end do
            close(10)
            ! Log
            write(*, "(a,2es12.3e3)") trim(file_ac), minval(ac_cur), maxval(ac_cur)
        end if

    end do

    stop
contains

    real(8) function weight(tt)
        implicit none
        real(8) :: tt
        weight = -2.0d0 * tt ** 3 + 3.0d0 * tt ** 2
        return
    end function weight
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
    real(8), allocatable :: Jld_old(:, :)
    real(8), allocatable :: Jld_cur(:, :)
    real(8), allocatable :: Jld_new(:, :)

    real(8) :: cp1 
    character(256) :: file_ac
    integer :: it, iz

    real(8) :: c1, c2, c3

    allocate(Ac_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Ac_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_old(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_cur(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))
    allocate(Jld_new(1:3, -nzvacl_m-1: nz_m+nzvacr_m+1))

    ! Initial electromagnetic field
    call calc_Ac_ext_t(-dt, -hz_m / cspeed_au, -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_old)
    call calc_Ac_ext_t(0.0d0, -hz_m / cspeed_au, -nzvacl_m-1, nz_m+nzvacr_m+1, Ac_cur)
    Ac_new = 0.0d0
    Jld_new = 0.0d0
    Jld_cur = 0.0d0
    Jld_old = 0.0d0

    c1 = (2.0d0 - (omega0 * dt) ** 2) / (1.0d0 + gamma * dt * 0.5d0)
    c2 = (1 - gamma * dt * 0.5d0) / (1.0d0 + gamma * dt * 0.5d0)
    c3 = alpha / (1.0d0 + gamma * dt * 0.5d0)

    do it = 0, nt
        do iz = -nzvacl_m, nz_m+nzvacr_m
             Ac_new(1, iz) = 2.0d0 * Ac_cur(1, iz) - Ac_old(1, iz) &
                & + dt ** 2 * cspeed_au ** 2 * (Ac_cur(1, iz+1) - 2.0d0 * Ac_cur(1, iz) + Ac_cur(1, iz-1)) / hz_m ** 2 &
                & + 4.0 * pi * Jld_cur(1, iz) * dt ** 2
             Ac_new(2, iz) = 2.0d0 * Ac_cur(2, iz) - Ac_old(2, iz) &
                & + dt ** 2 * cspeed_au ** 2 * (Ac_cur(2, iz+1) - 2.0d0 * Ac_cur(2, iz) + Ac_cur(2, iz-1)) / hz_m ** 2 &
                & + 4.0 * pi * Jld_cur(2, iz) * dt ** 2
        end do

        do iz = 1, nz_m
            Jld_new(1:3, iz) = c1 * Jld_cur(1:3, iz) - c2 * Jld_old(1:3, iz) &
                & - c3 * (Ac_new(1:3, iz) - 2.0d0 * Ac_cur(1:3, iz) + Ac_old(1:3, iz))
        end do


        Ac_old = Ac_cur
        Ac_cur = Ac_new
        Jld_old = Jld_cur
        Jld_cur = Jld_new

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
