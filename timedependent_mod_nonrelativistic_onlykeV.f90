MODULE time_dependence_nonrelativistic
    IMPLICIT NONE

    PRIVATE
    PUBLIC S_E, f, g, vel, Numerical

    contains

    !source function
    !E_j energy of the electron at bin j
    FUNCTION S_E(E_j, E_max) result(s_result)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, E_max
        DOUBLE PRECISION :: s_result
        
        if (4.6e-9< E_j .and. E_j < E_max) then
            s_result = 1.0 / (E_max-4.6e-9)
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

    !non-relativistic velocity
    !E_j energy of the electron at bin j
    !mass electron mass
    FUNCTION vel(E_j, mass) result(vel_res)
        IMPLICIT NONE
        DOUBLE PRECISION, intent(in) :: E_j, mass
        DOUBLE PRECISION :: vel_res
        
        vel_res = sqrt(mass/(2*E_j))
    END FUNCTION

    !The numerical calculation of z(E,t) using implicit Euler
    !number_dens atom number density n_i
    !K_constant from heat loss
    !dtime the time step
    !steady defines if one calculates steady state or not
    SUBROUTINE Numerical(Energy, E_average, Ionization, number_dens, K_constant, z_new, z_initial, &
                        sigma0, mass, dtime, E_max, steady)
        IMPLICIT NONE
        INTEGER :: j_idx, NE, k_idx, i_idx, i_idx2
        DOUBLE PRECISION :: Ii, E_avgi
        DOUBLE PRECISION :: c1, c2, c3, c4, LHS, Ej, Ep, integral1, integral2, dEnergy, n_idx          
        DOUBLE PRECISION, INTENT(IN) :: Energy(:), E_average(:), Ionization(:), number_dens(:),&
                                        K_constant, E_max, sigma0, mass, dtime, steady
        DOUBLE PRECISION, INTENT(INOUT) :: z_new(:), z_initial(:)

        !Delta E
        dEnergy = Energy(2) - Energy(1)

        NE = size(Energy)

        !Back propagation from highest to lowest energy
        do j_idx = NE-1, 1, -1
            Ej = Energy(j_idx)
            c1 = vel(Ej, mass)

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
                        * ( z_new(k_idx)*f(Ep, Ej, Ii, E_avgi)/(E_avgi*atan((Ep-Ii)/(E_avgi))) &
                             + z_new(k_idx+1)*f(Energy(k_idx+1), Ej, Ii, E_avgi)/(E_avgi*atan((Energy(k_idx+1)-Ii)/(E_avgi))) )/2.0 &
                        * dEnergy
                    end if

                    !Integral of the creation of secondaries
                    !g(E_j, E_avg) is the Opal cross section without normalisation
                    if (Ep > (Ii + Ej)) then
                        integral2 = integral2 + n_idx*sigma0 &
                        * ( z_new(k_idx)/(E_avgi*atan((Ep-Ii)/(E_avgi))) &
                            + z_new(k_idx+1)/(E_avgi*atan((Energy(k_idx+1)-Ii)/(E_avgi))) )/2.0 &
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
                
                c2 = c2 + n_idx * sigma0/(E_avgi*atan((Ej-Ii)/(E_avgi))) * E_avgi * atan((Ej-Ii)/(E_avgi))
            end do

            !One part of the heat loss
            c3 = ( mass / (2.0 * Ej) ) * K_constant *(1.0 / dEnergy + 1.0 / Ej)

            !LHS of the equation
            LHS = c1/dtime*steady + c2 + c3

            !Other part of the heat loss
            c4 = ( mass / (2.0 * Ej) ) * K_constant * z_new(j_idx+1)/ dEnergy

            !Calculation of z(E,t)
            z_new(j_idx) = (S_E(Ej,E_max) + integral1 + integral2 + c4 + z_initial(j_idx)*c1/dtime*steady) / LHS
            !print *, "s", S_E(Ej,E_max, time)
            !print *, "I1", integral1
            !print *, "I2", integral2
            !print *, "c2", c2
            !print *, "c3", c3
            !print *, "c4", c4 
        end do

        return

    END SUBROUTINE Numerical

END MODULE
