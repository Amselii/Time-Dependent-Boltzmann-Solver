program Forward_Euler_steady_state
    use Foward_Euler

    implicit none
    INTEGER :: nu, l_idx
    integer :: file, ii
    DOUBLE PRECISION :: mass, q, n_E, Gamma, K_constant, sigma0
    DOUBLE PRECISION :: E_min, E_max, dE, PI
    DOUBLE PRECISION, ALLOCATABLE :: E(:), z_init(:), z_new(:)
    DOUBLE PRECISION, DIMENSION(1) :: n_i, I, E_avg

    ! Defining constants
    mass = 9.11e-28
    q = 4.8e-10
    n_E = 1e10 !1e3-1e10 cm^3
    Gamma = 20.0 !15-30
    PI = 3.1415
    K_constant = 4*PI*q**4/mass *n_E*Gamma
    sigma0 = 6e-7

    E_min = 1.6e-13 !low
    E_max = 1.6e-6!1.6e-6 !1Mev
    dE = 1.6e-12*5 !(E_max - E_min) 1eV
    nu = int((E_max - E_min)/dE)


    ALLOCATE(E(nu))
    E = [(E_min + (l_idx-1)*dE, l_idx=1, nu)]

    n_i = [1e10]
    I = [2.17e-11]
    E_avg = 0.6*I

    ALLOCATE(z_init(size(E)))
    z_init = 0.0

    ALLOCATE(z_new(size(E)))
    z_new = z_init
    z_new(size(z_new)) = 4e14

    call Numerical(E, E_avg, I, n_i, K_constant, z_new, sigma0, mass, E_max)

    open(newunit=file, file="data_test_steady_energy_res5eV_start4e14.txt", status="replace", action="write")
    do ii=1, nu
        write(file, *) E(ii), z_new(ii)
    end do
    close(file)


    !print *, z_new

end program