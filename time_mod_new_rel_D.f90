MODULE time_dependence_relativistic
    IMPLICIT NONE

    PRIVATE
    PUBLIC S_E, f, g, sigma, beta, gamma_func, Numerical

    contains

    !time dependent source function
    !E_j energy of the electron at bin j
    FUNCTION S_E(E_j, E_max, time) result(s_result)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, E_max, time
        DOUBLE PRECISION :: s_result
        
        if (E_max*0.98< E_j .and. E_j < E_max) then
            s_result = 1.0 / (E_max*0.1)!*(time/86400.0)**(-1.5)
        else
            s_result = 0.0
        end if
    END FUNCTION

    !cross section primaries
    !E_prime primary electron energy
    !E_j energy of the electron at bin j
    !I_i ionisation potential
    !E_avg average energy
    FUNCTION f(E_prime, E_j, I_i, E_avg) result(f_res)
        IMPLICIT NONE 
        DOUBLE PRECISION, intent(in) :: E_prime, E_j, I_i, E_avg
        DOUBLE PRECISION :: f_res, num, den

        num = 1.0
        den = 1.0 + ((E_prime - E_j - I_i)/E_avg)**2

        f_res = num/den
    END FUNCTION

    !cross section secondaries
    !E_j energy of the electron at bin j
    !E_avg average energy
    FUNCTION g(E_j, E_avg) result(g_res)
        IMPLICIT NONE
        DOUBLE PRECISION, intent(in) :: E_j, E_avg
        DOUBLE PRECISION :: g_res, num, den
        
        num = 1.0
        den = (1.0 + ((E_j)/E_avg)**2)
        g_res = num/den
    END FUNCTION

    FUNCTION sigma(E_j, I_i) result(sigma_res)
        IMPLICIT NONE 
        DOUBLE PRECISION, INTENT(IN) :: E_j, I_i
        DOUBLE PRECISION :: sigma_res
        sigma_res = 4.5e-14*(1.6e-12)**2*log(E_j/I_i)/(E_j*I_i)
    END FUNCTION

    !beta coefficient for relativistic speeds
    !E_j energy of the electron at bin j
    !mass electron mass
    !c_speed speed of light
    FUNCTION beta(E_j, mass, c_speed) result(beta_res)
        IMPLICIT NONE
        DOUBLE PRECISION, intent(in) :: E_j, mass, c_speed
        DOUBLE PRECISION :: beta_res

        beta_res = sqrt(1 - (1 /(1+E_j/(mass*c_speed**2))**2))
    END FUNCTION

    !gamma coefficient for relativity
    !E_j energy of the electron at bin j
    !mass electron mass
    !c_speed speed of light
    FUNCTION gamma_func(E_j, mass, c_speed) result(gamma_res)
        IMPLICIT NONE
        DOUBLE PRECISION, intent(in) :: E_j, mass, c_speed
        DOUBLE PRECISION :: gamma_res

        gamma_res = 1 + E_j/(mass*c_speed**2)
    END FUNCTION

    !The numerical calculation of z(E,t) using implicit Euler
    !number_dens atom number density n_i
    !K_constant from heat loss
    !dtime the time step
    !steady defines if one calculates steady state or not 
    !time total time
    SUBROUTINE Numerical(Energy, E_average, Ionization, number_dens, K_constant, z_new, z_initial, sigma0, &
                         mass, c_speed, dtime, E_max, steady, time, D, S_int, P_heat, P_ion)
        IMPLICIT NONE
        INTEGER :: j_idx, NE, k_idx, i_idx, i_idx2
        DOUBLE PRECISION :: Ii, E_avgi
        DOUBLE PRECISION :: c1, c2, c3, c4, LHS, Ej, Ep, integral1, integral2, dEnergy, n_idx
        DOUBLE PRECISION, INTENT(IN) :: Energy(:), E_average(:), Ionization(:), number_dens(:),&
                                        K_constant, E_max, sigma0, steady, mass, dtime, time, c_speed
        DOUBLE PRECISION, INTENT(OUT) :: D, S_int
        DOUBLE PRECISION, INTENT(INOUT):: P_heat(:), P_ion(:)
        DOUBLE PRECISION, INTENT(INOUT) :: z_new(:), z_initial(:)

        !Delta E
        dEnergy = Energy(2) - Energy(1)

        NE = size(Energy)

        !Back propagation from highest to lowest energy
        do j_idx = NE-1, 1, -1
            Ej = Energy(j_idx)
            c1 = 1/(beta(Ej, mass, c_speed)*c_speed)

            !Integral of the downscattered primaries
            integral1 = 0.0

            !Integral of the creation of secondaries
            integral2 = 0.0

            !calculation of the integrals using trapezoidal rule
            do k_idx = j_idx+1, size(Energy)-1 
                !Energy of the integration variable E'
                Ep = Energy(k_idx)

                !summing over the ith atoms and atom number densities
                do i_idx = 1, size(number_dens)
                    n_idx = number_dens(i_idx)

                    Ii = Ionization(i_idx)
                    E_avgi = E_average(i_idx)

                    !Integral of the downscattered primaries
                    !f(E_k, E_j, I_i, E_avg) is the Opal cross section without normalisation
                    if  (Ep > (Ii + Ej)) then 
                        integral1 = integral1 + n_idx*sigma0 &
                        * ( z_new(k_idx)*f(Ep, Ej, Ii, E_avgi)/(E_avgi*atan((Ep-Ii)/(E_avgi))) & !sigma(Ep, Ii)
                             + z_new(k_idx+1)*f(Energy(k_idx+1), Ej, Ii, E_avgi)& !sigma(Energy(k_idx+1), Ii)
                             /(E_avgi*atan((Energy(k_idx+1)-Ii)/(E_avgi))) )/2.0 &
                        * dEnergy
                    end if

                    !Integral of the creation of secondaries
                    !g(E_j, E_avg) is the Opal cross section without normalisation
                    if (Ep > (Ii + Ej)) then
                        integral2 = integral2 + n_idx*sigma0 &
                        * ( z_new(k_idx)/(E_avgi*atan((Ep-Ii)/(E_avgi))) & !sigma(Ep, Ii)
                            + z_new(k_idx+1)/(E_avgi*atan((Energy(k_idx+1)-Ii)/(E_avgi))) )/2.0 & !sigma(Energy(k_idx+1), Ii)
                            * g(Ej, E_avgi) * dEnergy
                    end if

                end do
            
            end do

            !ionisation loss term
            c2 = 0.0
            !summing over the ith atoms and atom number densities
            do i_idx2 = 1, size(number_dens)
                Ii = Ionization(i_idx2)
                n_idx = number_dens(i_idx2)
                E_avgi = E_average(i_idx2)
                
                !if (Ej > Ii) then
                c2 = c2 + n_idx * sigma0/(E_avgi*atan((Ej-Ii)/(E_avgi))) * E_avgi * atan((Ej-Ii)/(E_avgi)) !sigma(Ej, Ii)
                !end if
            end do

            !One part of the heat loss
            c3 = K_constant/(beta(Ej, mass, c_speed)**2*c_speed**2)&
                *(1.0 / dEnergy + 2.0 /(beta(Ej, mass, c_speed)**2 *mass*c_speed**2 *gamma_func(Ej, mass, c_speed)**3)) 
            
            !LHS of the equation
            LHS = c1/dtime*steady + c2 + c3

            !Other part of the heat loss
            c4 = K_constant/(beta(Ej, mass, c_speed)**2*c_speed**2)* z_new(j_idx+1)/ dEnergy 

            !Calculation of z(E,t)
            z_new(j_idx) = (S_E(Ej,E_max, time) + integral1 + integral2 + c4 + z_initial(j_idx)*c1/dtime*steady) / LHS

            P_heat(j_idx) = z_new(j_idx)*K_constant/(beta(Ej, mass, c_speed)**2*c_speed**2)
            !print *, "L_heat", L_heat(j_idx)
                
            Ii = 2.17e-11
            !if (Ej > Ii) then
            P_ion(j_idx) = z_new(j_idx)*n_idx*sigma0*Ii !sigma(Ej, Ii)
            !print *, "L_ion", L_ion(j_idx)
            !end if
            !print *, "Lheat", K_constant/(beta(Ej, mass, c_speed)**2*c_speed**2)&
            !    *((z_new(j_idx+1) - z_new(j_idx))/dEnergy&
            !     - (2*z_new(j_idx))/(beta(Ej, mass, c_speed)**2*gamma_func(Ej, mass, c_speed)**3*mass*c_speed**2) )
            !print *, "Lion", z_new(j_idx)*c2
            !print *, "c2", c2
            !print *, "c3", c3
            !print *, "c4", c4 !K_constant/(beta(Ej, mass, c_speed)**2*c_speed**2)/ dEnergy
            !print *, "I1", integral1
            !print *, "I2", integral2

            !Loss(j_idx) = (c2*z_new(j_idx) + c3*z_new(j_idx) - c4) !*Energy(j_idx)

        end do

        !P_heat = 0.0
        !P_ion = 0.0
        D = 0.0
        S_int = 0.0

        do j_idx = 1, NE, 1

           ! P_heat =  P_heat + L_heat(j_idx)!* dEnergy
           ! P_ion = P_ion + L_ion(j_idx)!* dEnergy

            D = D + (P_heat(j_idx) + P_ion(j_idx))* dEnergy

            S_int = S_int + S_E(Energy(j_idx), E_max, time)*Energy(j_idx)*dEnergy


        end do

        print *, "D", D
        print *, "S_int", S_int

        return

    end subroutine Numerical

end module

