module milca_toroidal_momentum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t
    use milca_toroidal_momentum_leaves, only: momentum_shell, momentum_shell_interval, &
        momentum_shell_jvp, momentum_shell_vjp
    implicit none
    private
    public :: toroidal_momentum_step, toroidal_momentum_step_jvp, toroidal_momentum_step_vjp

    interface toroidal_momentum_step
        module procedure step, step_interval
    end interface
contains
    !> Cell-integrated angular momentum [kg m**2/s] over dt [s]. Face and
    !> shell torques [N m] are time averages. Faces are ordered inside out,
    !> positive flux transports positive angular momentum toward the exterior.
    pure subroutine step(old, dt, face, torque, next)
        real(dp), intent(in) :: old(:), dt, face(size(old) + 1), torque(size(old))
        real(dp), intent(out) :: next(size(old))
        integer :: n

        n = size(old)
        call momentum_shell(old, dt, face(:n), face(2:), torque, next)
    end subroutine step

    pure subroutine step_interval(old, dt, face, torque, next)
        type(interval_t), intent(in) :: old(:), dt, face(size(old) + 1), torque(size(old))
        type(interval_t), intent(out) :: next(size(old))
        integer :: n

        n = size(old)
        call momentum_shell_interval(old, dt, face(:n), face(2:), torque, next)
    end subroutine step_interval

    pure subroutine toroidal_momentum_step_jvp(old, dt, face, torque, &
                                              d_old, d_dt, d_face, d_torque, d_next)
        real(dp), intent(in) :: old(:), dt, face(size(old) + 1), torque(size(old))
        real(dp), intent(in) :: d_old(size(old)), d_dt, d_face(size(old) + 1), d_torque(size(old))
        real(dp), intent(out) :: d_next(size(old))
        integer :: n

        n = size(old)
        call momentum_shell_jvp(old, dt, face(:n), face(2:), torque, &
                                d_old, d_dt, d_face(:n), d_face(2:), d_torque, d_next)
    end subroutine toroidal_momentum_step_jvp

    pure subroutine toroidal_momentum_step_vjp(old, dt, face, torque, &
                                              b_next, b_old, b_dt, b_face, b_torque)
        real(dp), intent(in) :: old(:), dt, face(size(old) + 1), torque(size(old))
        real(dp), intent(in) :: b_next(size(old))
        real(dp), intent(out) :: b_old(size(old)), b_dt, b_face(size(old) + 1), b_torque(size(old))
        real(dp) :: b_dt_shell, b_inner
        integer :: i

        b_dt = 0.0_dp
        b_face(1) = 0.0_dp
        do i = 1, size(old)
            call momentum_shell_vjp(old(i), dt, face(i), face(i + 1), torque(i), b_next(i), &
                                    b_old(i), b_dt_shell, b_inner, b_face(i + 1), b_torque(i))
            b_face(i) = b_face(i) + b_inner
            b_dt = b_dt + b_dt_shell
        end do
    end subroutine toroidal_momentum_step_vjp
end module milca_toroidal_momentum
