program test_toroidal_momentum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval, interval_t
    use milca_toroidal_momentum, only: toroidal_momentum_step, &
        toroidal_momentum_step_jvp, toroidal_momentum_step_vjp
    implicit none

    call exchange()
    call common_flux()
    call manufactured_history()
    call derivatives()
    call interval_extrema()
contains
    subroutine exchange()
        real(dp), parameter :: old(3) = [2.0_dp, 1.0_dp, -1.0_dp]
        real(dp), parameter :: face(4) = [0.0_dp, 4.0_dp, -2.0_dp, 0.0_dp]
        real(dp), parameter :: torque(3) = 0.0_dp
        real(dp) :: next(3), reversed(3)

        call toroidal_momentum_step(old, 0.25_dp, face, torque, next)
        call close(next, [1.0_dp, 2.5_dp, -1.5_dp], 'signed shell exchange')
        call close([sum(next)], [sum(old)], 'closed angular momentum ledger')
        call toroidal_momentum_step(-old, 0.25_dp, -face, -torque, reversed)
        call close(reversed, -next, 'toroidal convention reversal')
    end subroutine exchange

    subroutine common_flux()
        real(dp), parameter :: flux = 2.0_dp**60
        real(dp) :: next(1)

        call toroidal_momentum_step([0.0_dp], 0.5_dp, [flux, flux], [-2.0_dp], next)
        call close(next, [-1.0_dp], 'small source with a large common face flux')
    end subroutine common_flux

    subroutine manufactured_history()
        real(dp), parameter :: edge(4) = [0.0_dp, 0.2_dp, 0.7_dp, 1.0_dp]
        real(dp), parameter :: steps(3) = [0.125_dp, 0.375_dp, 0.25_dp]
        real(dp), parameter :: initial(3) = [1.0_dp, -2.0_dp, 4.0_dp]
        real(dp) :: width(3), squared_width(3), face(4), torque(3), old(3), next(3)
        real(dp) :: t, midpoint, dt, expected(3)
        type(interval_t) :: old_box(3), next_box(3), face_box(4), torque_box(3)
        integer :: k

        width = edge(2:) - edge(:3)
        squared_width = edge(2:)**2 - edge(:3)**2
        old = initial
        old_box = interval(initial)
        t = 0.0_dp
        do k = 1, size(steps)
            dt = steps(k)
            midpoint = t + dt/2
            face = (1 + 2*midpoint)*edge**2 - midpoint*edge
            torque = (4*midpoint - 3)*width
            face_box = interval(face - 1.0e-14_dp, face + 1.0e-14_dp)
            torque_box = interval(torque - 1.0e-14_dp, torque + 1.0e-14_dp)
            call toroidal_momentum_step(old, dt, face, torque, next)
            call toroidal_momentum_step(old_box, interval(dt), face_box, torque_box, next_box)
            t = t + dt
            expected = initial + (2.5_dp*t**2 - 3*t)*width - (t + t**2)*squared_width
            call close(next, expected, 'independent integrated history')
            call close([sum(next)], [sum(initial) - 4*t + 1.5_dp*t**2], 'open momentum ledger')
            if (.not. all(next_box%lo <= expected .and. expected <= next_box%hi)) &
                error stop 'interval enclosure of integrated history'
            old = next
            old_box = next_box
        end do
    end subroutine manufactured_history

    subroutine derivatives()
        real(dp), parameter :: old(3) = [0.5_dp, -1.25_dp, 4.0_dp], dt = 0.5_dp
        real(dp), parameter :: face(4) = [2.0_dp, -1.0_dp, 3.5_dp, 1.25_dp]
        real(dp), parameter :: torque(3) = [-1.0_dp, 0.25_dp, -0.5_dp]
        real(dp), parameter :: direction(3) = [1.0_dp, -0.5_dp, 0.25_dp], d_dt = 0.5_dp
        real(dp), parameter :: d_face(4) = [0.5_dp, 0.25_dp, -0.75_dp, 1.5_dp]
        real(dp), parameter :: d_torque(3) = [0.75_dp, -0.75_dp, 1.0_dp]
        real(dp), parameter :: weight(3) = [2.0_dp, -1.5_dp, 0.25_dp], h = 2.0_dp**(-12)
        real(dp) :: tangent(3), hi(3), lo(3), b_old(3), b_dt, b_face(4), b_torque(3)

        call toroidal_momentum_step(old + h*direction, dt + h*d_dt, &
                                    face + h*d_face, torque + h*d_torque, hi)
        call toroidal_momentum_step(old - h*direction, dt - h*d_dt, &
                                    face - h*d_face, torque - h*d_torque, lo)
        call toroidal_momentum_step_jvp(old, dt, face, torque, direction, d_dt, d_face, &
                                        d_torque, tangent)
        call close(tangent, (hi - lo)/(2*h), 'centred directional derivative')
        call toroidal_momentum_step_vjp(old, dt, face, torque, weight, b_old, b_dt, &
                                        b_face, b_torque)
        call close([dot_product(weight, tangent)], &
            [dot_product(b_old, direction) + b_dt*d_dt + dot_product(b_face, d_face) &
            + dot_product(b_torque, d_torque)], 'adjoint pairing over shared faces')
    end subroutine derivatives

    subroutine interval_extrema()
        type(interval_t) :: next(1)

        call toroidal_momentum_step([interval(12.0_dp, 14.0_dp)], &
            interval(0.25_dp, 0.5_dp), [interval(-1.0_dp, 2.0_dp), interval(4.0_dp, 6.0_dp)], &
            [interval(2.0_dp, 3.0_dp)], next)
        if (next(1)%lo > 9.5_dp .or. next(1)%hi < 14.5_dp) &
            error stop 'interval enclosure of both extrema'
    end subroutine interval_extrema

    subroutine close(value, expected, label)
        real(dp), intent(in) :: value(:), expected(:)
        character(*), intent(in) :: label
        if (all(abs(value - expected) <= 32*epsilon(1.0_dp)*max(1.0_dp, abs(expected)))) return
        print *, label
        print *, 'computed ', value
        print *, 'expected ', expected
        error stop 1
    end subroutine close
end program test_toroidal_momentum
