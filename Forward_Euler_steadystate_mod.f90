module Foward_Euler
    implicit none

    private
    public S_E, f, g, Numerical

    contains

    function S_E(E_j, E_max) result(s_result)
        implicit none
        DOUBLE PRECISION, INTENT(IN) :: E_j, E_max
        Double precision :: s_result
        
        if (E_max*0.98< E_j .and. E_j < E_max) then
            s_result = 1.0 / (E_max*0.1)
        else
            s_result = 0.0
        end if
    end function

    !cross section primaries
    function f(E_prime, E_j, E_avg) result(f_res)
        implicit none 
        Double precision, intent(in) :: E_prime, E_j, E_avg
        Double precision :: f_res, num, den

        num = 1.0
        den = 1.0 + ((E_prime - E_j)/E_avg)**2

        f_res = num/den
    end function

    !cross section secondaries
    function g(E_j, I_i, E_avg) result(g_res)
        implicit none 
        Double precision, intent(in) :: E_j, I_i, E_avg
        Double precision :: g_res, num, den
        
        num = 1.0
        den = (1.0 + ((I_i + E_j)/E_avg)**2)
        g_res = num/den
    end function

    subroutine Numerical(Energy, E_average, Ionization, number_dens, K_constant, z_new, sigma0, m, E_max)
        implicit none
        integer :: j_idx, NE, k_idx, i_idx, i_idx2
        double precision :: c2, c3, c4, LHS, sigma0, Ej, Ii, integral1, integral2, Ep, E_avgi, dEnergy, m, n_idx
        double precision, intent(in) :: Energy(:), E_average(:), Ionization(:), number_dens(:), K_constant, E_max
        double precision, allocatable :: z_new(:)

    

        dEnergy = Energy(2) - Energy(1)

        NE = size(Energy)

        do j_idx = NE-1, 1, -1
            Ej = Energy(j_idx)

            integral1 = 0.0
            integral2 = 0.0

            do k_idx = j_idx+1, size(Energy)

                Ep = Energy(k_idx)

                do i_idx = 1, size(number_dens)
                    n_idx = number_dens(i_idx)

                    Ii = Ionization(i_idx)
                    E_avgi = E_average(i_idx)

                    
                    if (Ep > Ii) then
                        integral1 =  integral1 + n_idx*sigma0*z_new(k_idx)*f(Ep, Ej, E_avgi)*dEnergy
                    end if

                    if (Ep > (Ii + Ej)) then

                        integral2 = integral2 + n_idx*sigma0*z_new(k_idx)*g(Ej, Ii, E_avgi)*dEnergy
                    end if

                end do
            
            end do

            c2 = 0.0
            do i_idx2 = 1, size(number_dens)
                
                n_idx = number_dens(i_idx2)
                E_avgi = E_average(i_idx2)
                
                c2 = c2 + n_idx * sigma0 * E_avgi * atan(Ej/E_avgi)
            end do

            c3 = ( m / (2.0 * Ej) ) * K_constant *(1.0 / dEnergy + 1.0 / Ej)

            LHS = c2 + c3

            c4 = ( m / (2.0 * Ej) ) * K_constant * z_new(j_idx+1)/ dEnergy

          
            z_new(j_idx) = (S_E(Ej, E_max) + integral1 + integral2 + c4) / LHS
            !print *, "S", S_E(Ej, E_max)
            !print *, "I1", integral1
            !print *, "I2", integral2
            !print *, "c4", c4
            !print *, "c3", c3
            !print *, "c2", c2
            !print *, "Znew", z_new

        end do

        return

    end subroutine Numerical

end module

