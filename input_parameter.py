#!/usr/bin/env python3

def construct():
    add_kw("calculation", "theory", "character(256)", default="'1d'")
    add_kw("tgrid", "nt", "integer", default="0")
    add_kw("tgrid", "dt", "real(8)", default="0.0d0")
    add_kw("emfield", "e_impulse", "real(8)", default="0.0d0")
    add_kw("emfield", "ae_shape1", "character(256)", default="'none'")
    add_kw("emfield", "E_amplitude1", "real(8)", default="0.0d0")
    add_kw("emfield", "I_wcm2_1", "real(8)", default="-1.0d0")
    add_kw("emfield", "tw1", "real(8)", default="0.0d0")
    add_kw("emfield", "omega1", "real(8)", default="0.0d0")
    add_kw("emfield", "epdir_re1", "real(8)", default="0.0d0", dimension=3)
    add_kw("emfield", "epdir_im1", "real(8)", default="0.0d0", dimension=3)
    add_kw("emfield", "phi_cep1", "real(8)", default="0.0d0")
    add_kw("emfield", "ae_shape2", "character(256)", default="'none'")
    add_kw("emfield", "E_amplitude2", "real(8)", default="0.0d0")
    add_kw("emfield", "I_wcm2_2", "real(8)", default="-1.0d0")
    add_kw("emfield", "tw2", "real(8)", default="0.0d0")
    add_kw("emfield", "omega2", "real(8)", default="0.0d0")
    add_kw("emfield", "epdir_re2", "real(8)", default="0.0d0", dimension=3)
    add_kw("emfield", "epdir_im2", "real(8)", default="0.0d0", dimension=3)
    add_kw("emfield", "phi_cep2", "real(8)", default="0.0d0")
    add_kw("emfield", "t1_t2", "real(8)", default="0.0d0")
    add_kw("emfield", "t1_start", "real(8)", default="0.0d0")
    add_kw("multiscale", "fdtddim", "character(256)", default="''")
    add_kw("multiscale", "nx_m", "integer", default="1")
    add_kw("multiscale", "ny_m", "integer", default="1")
    add_kw("multiscale", "nz_m", "integer", default="1")
    add_kw("multiscale", "hx_m", "real(8)", default="1.0d2")
    add_kw("multiscale", "hy_m", "real(8)", default="1.0d2")
    add_kw("multiscale", "hz_m", "real(8)", default="1.0d2")
    add_kw("multiscale", "nzvacl_m", "integer", default="1000")
    add_kw("multiscale", "nzvacr_m", "integer", default="1000")
    add_kw("multiscale", "alpha", "real(8)", default="1.2d0")
    add_kw("multiscale", "gamma", "real(8)", default="1.0d-3")
    add_kw("multiscale", "omega0", "real(8)", default="1.0d0")

















import collections

data = collections.defaultdict(dict)
def add_kw(group, name, f90type, default=None, dimension=0):
    global data
    data[group][name] = (f90type, default, dimension)

template = """
module input_parameter
    implicit none
{DEFINE}
contains
subroutine read_input()
    implicit none
    character(256) :: tmp
    integer :: iret
    integer, parameter :: ifp = 9

{NAMELIST}

{DEFAULT}

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
{READ}
    close(ifp)

{VARDUMP}

    return
end subroutine read_input
end module input_parameter
"""


def create_define():
    tmp = []
    for group in data:
        for name in data[group]:
            f90type, default, dimension = data[group][name]
            if dimension:
                tmp += ["{F90TYPE}, dimension({DIMENSION}) :: {NAME}".format(
                    F90TYPE=f90type,
                    DIMENSION=dimension,
                    NAME=name,
                )]
            else:
                tmp += ["{F90TYPE} :: {NAME}".format(
                    F90TYPE=f90type,
                    NAME=name,
                )]   
    return "\n".join(tmp)

def create_read():
    tmp = []
    for group in data:
        tmp += ["read(ifp, nml={GROUP}, iostat=iret); rewind(ifp)".format(GROUP=group)]
    return "\n".join(tmp)

def create_namelist():
    tmp = []
    for group in data:
        tmp += ["namelist/{GROUP}/{LIST}".format(
            GROUP=group,
            LIST=", &\n".join(data[group].keys())
        )]
    return "\n".join(tmp)


def create_default():
    tmp = []
    for group in data:
        for name in data[group]:
            f90type, default, dimension = data[group][name]
            if default:
                tmp += ["{NAME} = {DEFAULT}".format(
                    NAME=name,
                    DEFAULT=default
                )]
    return "\n".join(tmp)

def create_vardump():
    tmp = []
    for group in data:
        for name in data[group]:
            f90type, default, dimension = data[group][name]
            if "int" in f90type.lower():
                tmp += ["write(*, '(a, 99i9)') '# {GROUP}.{NAME}:', {NAME}".format(
                    GROUP=group,
                    NAME=name)
                ]
            elif "real" in f90type.lower():
                tmp += ["write(*, '(a, 99es25.15)') '# {GROUP}.{NAME}:', {NAME}".format(
                    GROUP=group,
                    NAME=name)
                ]
            elif "char" in f90type.lower():
                tmp += ["write(*, '(a, 99a)') '# {GROUP}.{NAME}:', trim({NAME})".format(
                    GROUP=group,
                    NAME=name)
                ]
    return "\n".join(tmp)

construct()

print(template.format(
    DEFINE=create_define(),
    DEFAULT=create_default(),
    NAMELIST=create_namelist(),
    READ=create_read(),
    VARDUMP=create_vardump(),
))

