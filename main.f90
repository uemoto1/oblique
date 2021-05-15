program main
    use input_parameter
    use em_field
    implicit none
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: i
    call read_input()


    allocate(Ac_ext_t(1:3, 0:nt))

    call calc_Ac_ext_t(0.0d0, dt, 0, nt, Ac_ext_t)

    do i = 0, nt
        write(*, *) i*dt, Ac_ext_t(1:3, i)
    end do

    stop
end program main
