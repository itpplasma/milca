module milca_ideal_fb_coil_field
    ! Coil field for the 3D vacuum on the plasma boundary
    ! (theory/derivations/free-boundary-plasma-boundary-sheet-3d.md): B and its
    ! first and second derivatives in (R, Z) at fixed toroidal angle, inside
    ! the plasma and on its boundary, and F, the toroidal circulation of B
    ! over 2 pi (R B_phi averaged over the angle).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, abstract, public :: fb_field_jets_t
        real(dp) :: f_vac = 0
    contains
        procedure(field_jets_eval), deferred :: eval
    end type

    abstract interface
        subroutine field_jets_eval(self, r, z, v, b, db, d2b)
            !! (B_R, B_phi, B_Z) at (R, Z, phi = v), derivatives in (R, Z)
            !! (db(:, 1:2)) and second derivatives (RR, RZ, ZZ).
            import :: fb_field_jets_t, dp
            class(fb_field_jets_t), intent(in) :: self
            real(dp), intent(in) :: r, z, v
            real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        end subroutine field_jets_eval
    end interface

end module milca_ideal_fb_coil_field
