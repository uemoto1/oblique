
module input_parameter
    implicit none
character(256) :: theory
integer :: nt
real(8) :: dt
real(8) :: e_impulse
character(256) :: ae_shape1
real(8) :: E_amplitude1
real(8) :: I_wcm2_1
real(8) :: tw1
real(8) :: omega1
real(8), dimension(3) :: epdir_re1
real(8), dimension(3) :: epdir_im1
real(8) :: phi_cep1
character(256) :: ae_shape2
real(8) :: E_amplitude2
real(8) :: I_wcm2_2
real(8) :: tw2
real(8) :: omega2
real(8), dimension(3) :: epdir_re2
real(8), dimension(3) :: epdir_im2
real(8) :: phi_cep2
real(8) :: t1_t2
real(8) :: t1_start
character(256) :: fdtddim
integer :: nx_m
integer :: ny_m
integer :: nz_m
real(8) :: hx_m
real(8) :: hy_m
real(8) :: hz_m
integer :: nzvacl_m
integer :: nzvacr_m
real(8) :: alpha
real(8) :: gamma
real(8) :: omega0
real(8) :: theta_deg
contains
subroutine read_input()
    implicit none
    character(256) :: tmp
    integer :: iret
    integer, parameter :: ifp = 9

namelist/calculation/theory
namelist/tgrid/nt, &
dt
namelist/emfield/e_impulse, &
ae_shape1, &
E_amplitude1, &
I_wcm2_1, &
tw1, &
omega1, &
epdir_re1, &
epdir_im1, &
phi_cep1, &
ae_shape2, &
E_amplitude2, &
I_wcm2_2, &
tw2, &
omega2, &
epdir_re2, &
epdir_im2, &
phi_cep2, &
t1_t2, &
t1_start
namelist/multiscale/fdtddim, &
nx_m, &
ny_m, &
nz_m, &
hx_m, &
hy_m, &
hz_m, &
nzvacl_m, &
nzvacr_m, &
alpha, &
gamma, &
omega0, &
theta_deg

theory = '1d'
nt = 0
dt = 0.0d0
e_impulse = 0.0d0
ae_shape1 = 'none'
E_amplitude1 = 0.0d0
I_wcm2_1 = -1.0d0
tw1 = 0.0d0
omega1 = 0.0d0
epdir_re1 = 0.0d0
epdir_im1 = 0.0d0
phi_cep1 = 0.0d0
ae_shape2 = 'none'
E_amplitude2 = 0.0d0
I_wcm2_2 = -1.0d0
tw2 = 0.0d0
omega2 = 0.0d0
epdir_re2 = 0.0d0
epdir_im2 = 0.0d0
phi_cep2 = 0.0d0
t1_t2 = 0.0d0
t1_start = 0.0d0
fdtddim = ''
nx_m = 1
ny_m = 1
nz_m = 1
hx_m = 1.0d2
hy_m = 1.0d2
hz_m = 1.0d2
nzvacl_m = 1000
nzvacr_m = 1000
alpha = 1.2d0
gamma = 1.0d-3
omega0 = 1.0d0
theta_deg = 0.0d0

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
read(ifp, nml=calculation, iostat=iret); rewind(ifp)
read(ifp, nml=tgrid, iostat=iret); rewind(ifp)
read(ifp, nml=emfield, iostat=iret); rewind(ifp)
read(ifp, nml=multiscale, iostat=iret); rewind(ifp)
    close(ifp)

write(*, '(a, 99a)') '# calculation.theory:', trim(theory)
write(*, '(a, 99i9)') '# tgrid.nt:', nt
write(*, '(a, 99es25.15)') '# tgrid.dt:', dt
write(*, '(a, 99es25.15)') '# emfield.e_impulse:', e_impulse
write(*, '(a, 99a)') '# emfield.ae_shape1:', trim(ae_shape1)
write(*, '(a, 99es25.15)') '# emfield.E_amplitude1:', E_amplitude1
write(*, '(a, 99es25.15)') '# emfield.I_wcm2_1:', I_wcm2_1
write(*, '(a, 99es25.15)') '# emfield.tw1:', tw1
write(*, '(a, 99es25.15)') '# emfield.omega1:', omega1
write(*, '(a, 99es25.15)') '# emfield.epdir_re1:', epdir_re1
write(*, '(a, 99es25.15)') '# emfield.epdir_im1:', epdir_im1
write(*, '(a, 99es25.15)') '# emfield.phi_cep1:', phi_cep1
write(*, '(a, 99a)') '# emfield.ae_shape2:', trim(ae_shape2)
write(*, '(a, 99es25.15)') '# emfield.E_amplitude2:', E_amplitude2
write(*, '(a, 99es25.15)') '# emfield.I_wcm2_2:', I_wcm2_2
write(*, '(a, 99es25.15)') '# emfield.tw2:', tw2
write(*, '(a, 99es25.15)') '# emfield.omega2:', omega2
write(*, '(a, 99es25.15)') '# emfield.epdir_re2:', epdir_re2
write(*, '(a, 99es25.15)') '# emfield.epdir_im2:', epdir_im2
write(*, '(a, 99es25.15)') '# emfield.phi_cep2:', phi_cep2
write(*, '(a, 99es25.15)') '# emfield.t1_t2:', t1_t2
write(*, '(a, 99es25.15)') '# emfield.t1_start:', t1_start
write(*, '(a, 99a)') '# multiscale.fdtddim:', trim(fdtddim)
write(*, '(a, 99i9)') '# multiscale.nx_m:', nx_m
write(*, '(a, 99i9)') '# multiscale.ny_m:', ny_m
write(*, '(a, 99i9)') '# multiscale.nz_m:', nz_m
write(*, '(a, 99es25.15)') '# multiscale.hx_m:', hx_m
write(*, '(a, 99es25.15)') '# multiscale.hy_m:', hy_m
write(*, '(a, 99es25.15)') '# multiscale.hz_m:', hz_m
write(*, '(a, 99i9)') '# multiscale.nzvacl_m:', nzvacl_m
write(*, '(a, 99i9)') '# multiscale.nzvacr_m:', nzvacr_m
write(*, '(a, 99es25.15)') '# multiscale.alpha:', alpha
write(*, '(a, 99es25.15)') '# multiscale.gamma:', gamma
write(*, '(a, 99es25.15)') '# multiscale.omega0:', omega0
write(*, '(a, 99es25.15)') '# multiscale.theta_deg:', theta_deg

    return
end subroutine read_input
end module input_parameter

