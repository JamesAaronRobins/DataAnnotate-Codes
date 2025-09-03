program bbDIHforce
    !Use the intrinsic Fortran module for accurate floating points
    use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

    implicit none
    !Define all of the variables used in the code and calculations
    integer, parameter :: PREC = REAL64 
    real(PREC) :: imp1(3), imp2(3), imp3(3), imp4(3)
    real(PREC) :: f, delta, cosine, anglsign, d1, phi
    !real(PREC) :: c1(3), c2(3), c3(3), t1(3), t2(3)
    real(PREC) :: v12(3), v32(3),v34(3), fi(3), fj(3), fk(3), fl(3)
    !real(PREC) :: cos1
    real(PREC) :: m(3), n(3)
    real(PREC) :: dih, akj, akj2, pre, dvijvkj_akj2, dvklvkj_akj2
    real(PREC) :: dih_k = 0.081, dih_p0 = 2.302
    real(PREC) :: Forces1(3), Forces2(3), Forces3(3), Forces4(3)
    real(PREC) :: normF1, normF2, normF3, normF4

    ! Bead positions, given in prompt
    imp1 = [199.059,  53.939, 456.220]
    imp2 = [193.109,  53.492, 457.693]
    imp3 = [188.280,  53.617, 461.206]
    imp4 = [184.850,  54.600, 466.799]

    !Calculate the vectors between the beads
    v12(:) = vectors(imp1,imp2)
    v32(:) = vectors(imp3,imp2)
    v34(:) = vectors(imp3,imp4)

    !Calculate the 2 planes that the dihedral angle is calculated between 
    ! m = rij x rkj
    m(1) = v12(2)*v32(3) - v12(3)*v32(2)
    m(2) = v12(3)*v32(1) - v12(1)*v32(3)
    m(3) = v12(1)*v32(2) - v12(2)*v32(1)

    ! n = rkj x rkl
    n(1) = v32(2)*v34(3) - v32(3)*v34(2)
    n(2) = v32(3)*v34(1) - v32(1)*v34(3)
    n(3) = v32(1)*v34(2) - v32(2)*v34(1)

    ! This is the scalar rkj - Both the values are used later
    akj2 = dot_product(v32, v32)
    akj = sqrt(akj2)
    print*, akj
    ! Calculation of the dihedral angle
    dih = atan2(akj * dot_product(v12, n), dot_product(m, n))

    ! Here pre is the calculation of the value for dV/d(phi) normalised to akj 
    pre = - dih_k * sin(dih + dih_p0) * akj

    ! Integration of the force into each bead
    fi(:) = -pre / dot_product(m,m) * m(:)
    fl(:) =  pre / dot_product(n,n) * n(:)
    
    dvijvkj_akj2 = dot_product(v12, v32) / akj2
    dvklvkj_akj2 = dot_product(v34, v32) / akj2
    fj(:) = (-1.0_PREC + dvijvkj_akj2) * fi(:) + (          - dvklvkj_akj2) * fl(:)
    fk(:) = (          - dvijvkj_akj2) * fi(:) + (-1.0_PREC + dvklvkj_akj2) * fl(:)

    ! Sum the final force as the normal of the dimensions - Not just 1+2+3
    normF1 = norm2(fi)
    normF2 = norm2(fj)
    normF3 = norm2(fk)
    normF4 = norm2(fl)

    ! Output result as JSON-style
    write(*,'(A)') '{'
    write(*,'(A,F8.5)') '  "Bead 1": ', normF1, ','
    write(*,'(A,F8.5)') '  "Bead 2": ', normF2, ','
    write(*,'(A,F8.5)') '  "Bead 3": ', normF3, ','
    write(*,'(A,F8.5)') '  "Bead 4": ', normF4
    write(*,'(A)') '}'

contains

pure function vectors(v1, v2) result (vector)
    use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

    integer, parameter :: PREC = REAL64 

    real(PREC) :: vector(3)
    real(PREC), intent(in) :: v1(3), v2(3)
    integer :: i
    vector(:) = v1(:) - v2(:)

end function vectors

pure function cross(v1, v2) result(crossproduct)
    use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

    integer, parameter :: PREC = REAL64 

    real(PREC) :: crossproduct(3)
    real(PREC), intent(in) :: v1(3), v2(3)
    integer :: i

    crossproduct(1) = v1(2)*v2(3) - v1(3)*v2(2)
    crossproduct(2) = v1(3)*v2(1) - v1(1)*v2(3)
    crossproduct(3) = v1(1)*v2(2) - v1(2)*v2(1)
    
end function cross

end program
