program Simulation
    use time_dependence_relativistic

    IMPLICIT NONE
    INTEGER ::  l_idx, nu
    INTEGER :: file, ii
    DOUBLE PRECISION :: m, q, n_E, Lambda, K_constant, sigma0, steady, c_speed
    DOUBLE PRECISION :: E_min, E_max, dE, PI, start_time, time_total, dt
    DOUBLE PRECISION, ALLOCATABLE :: E(:), z_init(:), z_new(:)
    DOUBLE PRECISION, DIMENSION(1) :: n_i, I, E_avg

    !Defining constants in cgs units
    c_speed = 2.9e10 !speed of light
    m = 9.11e-28 !electron mass
    q = 4.8e-10 !electron charge
    Lambda = 20.0 !between 15-30
    PI = 3.1415
    sigma0 = 1e-19 !cross section normalisation constant

    I = [2.17e-11] !ionisation potential atom
    E_avg = 3.0*I !average energy

    E_min = 1.6e-13 !minimum energy
    E_max = 1.6e-6 !maximum energy
    dE = 1.6e-12*10.0 !energy resolution
    nu = int((E_max - E_min)/dE)

    !n_E = 1e8
    !K_constant = 4.0*PI*q**4/m*n_E*Gamma
    !n_i = [1e8]
    
    ALLOCATE(E(nu))
    E = [(E_min + (l_idx-1)*dE, l_idx=1, nu)]

    start_time = 2.0*86400.0

    ALLOCATE(z_init(size(E)))
    z_init = 0.0

    ALLOCATE(z_new(size(E)))
    z_new = z_init

    time_total = start_time

    !time step
    dt = 2.0*86400.0

    !steady = 1.0 means NO steady state
    !steady = 0.0 means steady state is on
    steady = 0.0
    
    open(newunit=file, file="data_steady_KN_10eV_3I_t2_rel_sigma1e-19.txt", status="replace", action="write")
    do while (time_total <= 2.0*86400.0)

        !time dependent number densities
        n_E = 1e8*(time_total/86400.0)**(-3)
        K_constant = 4.0*PI*q**4/m *n_E*Lambda
        n_i = [1e8*(time_total/86400.0)**(-3)]

        !call z(E,t)
        call Numerical(E, E_avg, I, n_i, K_constant, z_new, z_init, sigma0, m, c_speed, dt, E_max, steady, time_total)

        z_init = z_new
        dt = 2.0*86400.0!0.2*time_total
        time_total = time_total + dt

        do ii=1, nu
            write(file, *) E(ii), z_new(ii)
        end do

        print *, 'Finished an iteration'
        print *, time_total
    end do

    close(file)

end program