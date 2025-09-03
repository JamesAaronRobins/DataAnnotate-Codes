program dihedral_force
! Script modified for correct Answer - Not the origionally produced script for the promot
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), dimension(3) :: bead1, bead2, bead3, bead4
    real(dp), dimension(3) :: b1, b2, b3, m, n
    real(dp), dimension(3) :: fi, fj, fk, fl
    real(dp) :: k_phi, phi0, phi, dV_dphi
    real(dp) :: m_norm, n_norm, b2_norm, b2_sq
    real(dp) :: b1_dot_b2, b3_dot_b2
    real(dp) :: x, y
    real(dp) :: f1_norm, f2_norm, f3_norm, f4_norm

    ! Input coordinates
    bead1 = [199.059_dp, 53.939_dp, 456.220_dp]
    bead2 = [193.109_dp, 53.492_dp, 457.693_dp]
    bead3 = [188.280_dp, 53.617_dp, 461.206_dp]
    bead4 = [184.850_dp, 54.600_dp, 466.799_dp]

    ! Force constants
    k_phi = 0.081_dp
    phi0  = 2.302_dp

    ! Bond vectors
    b1 = bead1 - bead2
    b2 = bead3 - bead2
    b3 = bead4 - bead3

    ! Cross products for planes
    m = cross_product(b1, b2)
    n = cross_product(b2, b3)

    ! Magnitudes
    m_norm = norm2(m)
    n_norm = norm2(n)
    b2_norm = norm2(b2)

    ! Dihedral angle
    x = dot_product(m, n)
    y = dot_product(cross_product(m, n), b2 / b2_norm)
    phi = atan2(y, x)

    ! Derivative of potential
    dV_dphi = -k_phi * sin(phi + phi0) * b2_norm

    ! Dot products for projection
    b2_sq = dot_product(b2, b2)
    b1_dot_b2 = dot_product(b1, b2)
    b3_dot_b2 = dot_product(b3, b2)

    ! Forces from Bekker's paper
    fi = (-dV_dphi / (m_norm**2)) * m
    fl = ( dV_dphi / (n_norm**2)) * n

    fj = (-1.0_dp + b1_dot_b2 / b2_sq) * fi - (b3_dot_b2 / b2_sq) * fl
    fk = (-b1_dot_b2 / b2_sq) * fi + (-1.0_dp + b3_dot_b2 / b2_sq) * fl


    ! Magnitudes
    f1_norm = norm2(fi)
    f2_norm = norm2(fj)
    f3_norm = norm2(fk)
    f4_norm = norm2(fl)

    ! Output in JSON format
    print *, '{'
    print '(A,F8.3)', '  "Bead 1": ', f1_norm, ','
    print '(A,F8.3)', '  "Bead 2": ', f2_norm, ','
    print '(A,F8.3)', '  "Bead 3": ', f3_norm, ','
    print '(A,F8.3)', '  "Bead 4": ', f4_norm
    print *, '}'

contains

    pure function cross_product(a, b) result(c)
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

end program dihedral_force
