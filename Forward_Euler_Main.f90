program Forward_Euler_steady_state
    use Foward_Euler

    IMPLICIT NONE
    INTEGER ::  l_idx, nu
    INTEGER :: file, ii
    DOUBLE PRECISION :: m, q, n_E, Gamma, K_constant, sigma0, steady
    DOUBLE PRECISION :: E_min, E_max, dE, PI, dt, time_step
    DOUBLE PRECISION, ALLOCATABLE :: E(:), z_init(:), z_new(:)
    DOUBLE PRECISION, DIMENSION(1) :: n_i, I, E_avg

    ! Defining constants
    m = 9.11e-28
    q = 4.8e-10
    Gamma = 20.0 !15-30
    PI = 3.1415
    sigma0 = 6e-7

    E_min = 1.6e-13 !low
    E_max = 1.6e-6!1.6e-6 !1Mev
    dE = 1.6e-12*5 !(E_max - E_min) 1eV
    nu = int((E_max - E_min)/dE)

    n_E = 1e10
    K_constant = 4.0*PI*q**4/m *n_E*Gamma
    n_i = [1e10]
    I = [2.17e-11]
    E_avg = 0.6*I
    
    ALLOCATE(E(nu))
    E = [(E_min + (l_idx-1)*dE, l_idx=1, nu)]

    dt = 1.0

    ALLOCATE(z_init(size(E)))
    z_init = 0.0

    ALLOCATE(z_new(size(E)))
    z_new = z_init
    z_new(size(z_new)) = 0.0

    time_step = 0.0
    steady = 0.0

    open(newunit=file, file="data_testing_res_25eV.txt", status="replace", action="write")
    do while (time_step <= 0.0)

        !n_E = 1e10*(time_step/86400.0)**(-3)
        !K_constant = 4.0*PI*q**4/m *n_E*Gamma
        !n_i = [1e10*(time_step/86400.0)**(-3)]

        call Numerical(E, E_avg, I, n_i, K_constant, z_new, z_init, sigma0, m, dt, E_max, steady)

        z_init = z_new
        time_step = time_step + dt

        do ii=1, nu
            write(file, *) E(ii), z_new(ii)
        end do
    end do

    close(file)


    !print *, z_new

end program