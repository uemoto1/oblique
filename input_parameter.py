#!/usr/bin/env python3

data = [
    {
        "name": "calculation", "variable": [
            {"name": "theory", "type": "character(256)",  "default": "'1d'", "dimension": ""},
        ]
    },
    {
        "name": "control", "variable": [
            {"name": "sysname", "type": "character(256)",  "default": "''", "dimension": ""},
            {"name": "base_directory", "type": "character(256)",  "default": "''", "dimension": ""},
        ]
    },
    {
        "name": "multiscale", "variable": [
            {"name": "fdtddim", "type": "character(256)",  "default": "''", "dimension": ""},
            {"name": "nx_m", "type": "integer",  "default": "1", "dimension": ""},
            {"name": "ny_m", "type": "integer",  "default": "1", "dimension": ""},
            {"name": "nz_m", "type": "integer",  "default": "1", "dimension": ""},
            {"name": "hx_m", "type": "real(8)",  "default": "1e2", "dimension": ""},
            {"name": "hy_m", "type": "real(8)",  "default": "1e2", "dimension": ""},
            {"name": "hz_m", "type": "real(8)",  "default": "1e2", "dimension": ""},
            {"name": "nxvacl_m", "type": "integer",  "default": "1000", "dimension": ""},
            {"name": "nxvacr_m", "type": "integer",  "default": "1000", "dimension": ""},
        ]
    },
]




























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
        group_name = group["name"]
        tmp += ["    ! {GROUP}".format(GROUP=group_name)]
        for variable in group["variable"]:
            variable_name = variable["name"]
            variable_type = variable["type"]
            variable_dimension = variable.get("dimension", "")
            if variable_dimension:
                tmp += ["    {TYPE}, dimension({DIMENSION}) :: {NAME}".format(
                    TYPE=variable_type,
                    DIMENSION=variable_dimension,
                    NAME=variable_name,
                )]
            else:
                tmp += ["    {TYPE} :: {NAME}".format(
                    TYPE=variable_type,
                    NAME=variable_name,
                )]   
    return "\n".join(tmp)

def create_read():
    tmp = []
    for group in data:
        group_name = group["name"]
        tmp += ["    read(ifp, nml={GROUP}, iostat=iret)".format(
            GROUP=group_name
        )]
        tmp += ["    rewind(ifp)"]
    return "\n".join(tmp)

def create_namelist():
    tmp = []
    for group in data:
        group_name = group["name"]
        tmp2 = []
        for variable in group["variable"]:
            variable_name = variable["name"]
            tmp2 += ["    & {NAME}".format(NAME=variable_name)]
        tmp += ["    namelist/{GROUP}/ &".format(GROUP=group_name)]
        tmp += [", &\n".join(tmp2)]
    return "\n".join(tmp)

def create_default():
    tmp = []
    for group in data:
        group_name = group["name"]
        for variable in group["variable"]:
            variable_name = variable["name"]
            variable_default = variable.get("default", "")
            if variable_default:
                tmp += ["    {NAME} = {DEFAULT}".format(
                    NAME=variable_name,
                    DEFAULT=variable_default
                )]
    return "\n".join(tmp)

def create_vardump():
    tmp = []
    for group in data:
        group_name = group["name"]
        tmp += ["    write(*, '(a)') '# {GROUP}:'".format(
            GROUP=group_name
        )]
        for variable in group["variable"]:
            variable_name = variable["name"]
            variable_type = variable["type"]
            if "real" in variable_type:
                tmp += ["    write(*, '(a,100es25.15e3)') '#     {NAME}:', {NAME}".format(
                    NAME=variable_name
                )]
            elif "integer" in variable_type:
                tmp += ["    write(*, '(a,100i10)') '#     {NAME}:', {NAME}".format(
                    NAME=variable_name
                )]
            elif "character" in variable_type:
                tmp += ["    write(*, '(a,a)') '#     {NAME}:', trim({NAME})".format(
                    NAME=variable_name
                )]
    return "\n".join(tmp)


print(template.format(
    DEFINE=create_define(),
    DEFAULT=create_default(),
    NAMELIST=create_namelist(),
    READ=create_read(),
    VARDUMP=create_vardump(),
))

