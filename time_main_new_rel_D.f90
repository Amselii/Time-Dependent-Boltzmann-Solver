program Simulation
    use time_dependence_relativistic

    IMPLICIT NONE
    INTEGER ::  l_idx, nu
    INTEGER :: file, file2, ii
    DOUBLE PRECISION :: m, q, n_E, Gamma, K_constant, sigma0, steady, c_speed, D, S_int
    DOUBLE PRECISION :: E_min, E_max, dE, PI, start_time, time_total, dt
    DOUBLE PRECISION, ALLOCATABLE :: E(:), z_init(:), z_new(:), P_heat(:), P_ion(:)
    DOUBLE PRECISION, DIMENSION(1) :: n_i, I, E_avg

    ! Defining constants
    c_speed = 2.9e10
    m = 9.11e-28
    q = 4.8e-10
    Gamma = 20.0 !15-30
    PI = 3.1415
    sigma0 = 1e-18!6e-7

    E_min = 1.6e-13 !low
    E_max = 1.6e-12*3000 !1.6e-6!1.6e-6 !1Mev
    dE = 1.6e-12*2!*10.0 !(E_max - E_min) 1eV
    nu = int((E_max - E_min)/dE)

    n_E = 1e8
    K_constant = 4.0*PI*q**4/m*n_E*Gamma
    n_i = [1e8]
    I = [2.17e-11]
    E_avg = 0.6*I
    
    ALLOCATE(E(nu))
    E = [(E_min + (l_idx-1)*dE, l_idx=1, nu)]

    start_time = 0.0!1.0*86400.0

    ALLOCATE(z_init(size(E)))
    z_init = 0.0

    ALLOCATE(z_new(size(E)))
    z_new = z_init
    z_new(size(z_new)) = 0.0

    ALLOCATE(P_heat(size(E)))
    P_heat = 0.0

    ALLOCATE(P_ion(size(E)))
    P_ion = 0.0

    time_total = start_time
    steady = 1.0

    dt =0.000001 !1.0*86400.0
    

    !open(newunit=file, file="data_testing_KN_10eV_3I_2d_rel_adaptivedt02_sigma1e-19.txt", status="replace", action="write")
    open(newunit=file2, file="Power_heating_ion_D_keV_sigma01e-18_2eV_res_3k_new.txt", status="replace", action="write")
    open(newunit=file, file="keV_test_sigma1e-18_2eV_res_3k_new.txt", status="replace", action="write")
    !open(newunit=file2, file="trash2.txt", status="replace", action="write")
    do while (time_total <= 10.0)

        !n_E = 1e8*(time_total/86400.0)**(-3)
        !K_constant = 4.0*PI*q**4/m *n_E*Gamma
        !n_i = [1e8*(time_total/86400.0)**(-3)]
        D = 0.0
        S_int = 0.0
        !P_heat = 0.0
        !P_ion = 0.0

        call Numerical(E, E_avg, I, n_i, K_constant, z_new, z_init, sigma0, m, &
             c_speed, dt, E_max, steady, time_total, D, S_int, P_heat, P_ion)

        z_init = z_new
        dt = 0.1 !0.2*time_total
        time_total = time_total + dt

        do ii=1, nu
            write(file, *) E(ii), z_new(ii), time_total
            write(file2, *) P_heat(ii), P_ion(ii), D, S_int, E(ii), time_total
            !print *, "L_heat", L_heat(ii)
            !print *, "L_heat", L_ion(ii)
        end do

         !P_heat, P_ion, D, S_int, time_total

        print *, 'Finished an iteration'
        print *, time_total
    end do

    close(file)
    close(file2)



end program