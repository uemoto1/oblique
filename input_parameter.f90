
module input_parameter
    implicit none
    ! calculation
    character(256) :: theory
    ! control
    character(256) :: sysname
    character(256) :: base_directory
    ! multiscale
    character(256) :: fdtddim
    integer :: nx_m
    integer :: ny_m
    integer :: nz_m
    real(8) :: hx_m
    real(8) :: hy_m
    real(8) :: hz_m
    integer :: nxvacl_m
    integer :: nxvacr_m
contains
subroutine read_input()
    implicit none
    character(256) :: tmp
    integer :: iret
    integer, parameter :: ifp = 9

    namelist/calculation/ &
    & theory
    namelist/control/ &
    & sysname, &
    & base_directory
    namelist/multiscale/ &
    & fdtddim, &
    & nx_m, &
    & ny_m, &
    & nz_m, &
    & hx_m, &
    & hy_m, &
    & hz_m, &
    & nxvacl_m, &
    & nxvacr_m

    theory = '1d'
    sysname = ''
    base_directory = ''
    fdtddim = ''
    nx_m = 1
    ny_m = 1
    nz_m = 1
    hx_m = 1e2
    hy_m = 1e2
    hz_m = 1e2
    nxvacl_m = 1000
    nxvacr_m = 1000

    open(ifp, file='.namelist.tmp', action='write', status='replace')
    do while (.true.)
        read(*, '(a)', iostat=iret) tmp
        if (iret < 0) exit
        tmp = adjustl(tmp)
        if (tmp(1:1) .eq. '#') cycle
        write(ifp, '(a)') trim(tmp)
    end do
    close(ifp)

    open(ifp, file='.namelist.tmp', action='read', status='old')
    read(ifp, nml=calculation, iostat=iret)
    rewind(ifp)
    read(ifp, nml=control, iostat=iret)
    rewind(ifp)
    read(ifp, nml=multiscale, iostat=iret)
    rewind(ifp)
    close(ifp)

    write(*, '(a)') '# calculation:'
    write(*, '(a,a)') '#     theory:', trim(theory)
    write(*, '(a)') '# control:'
    write(*, '(a,a)') '#     sysname:', trim(sysname)
    write(*, '(a,a)') '#     base_directory:', trim(base_directory)
    write(*, '(a)') '# multiscale:'
    write(*, '(a,a)') '#     fdtddim:', trim(fdtddim)
    write(*, '(a,100i10)') '#     nx_m:', nx_m
    write(*, '(a,100i10)') '#     ny_m:', ny_m
    write(*, '(a,100i10)') '#     nz_m:', nz_m
    write(*, '(a,100es25.15e3)') '#     hx_m:', hx_m
    write(*, '(a,100es25.15e3)') '#     hy_m:', hy_m
    write(*, '(a,100es25.15e3)') '#     hz_m:', hz_m
    write(*, '(a,100i10)') '#     nxvacl_m:', nxvacl_m
    write(*, '(a,100i10)') '#     nxvacr_m:', nxvacr_m

    return
end subroutine read_input
end module input_parameter

