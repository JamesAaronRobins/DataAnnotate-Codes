! dihedral_force.f90
! This program computes the forces on a 4‑atom dihedral using a periodic
! dihedral potential of the form V = k * (1 + cos(phi + phi0)).  The
! dihedral angle phi is computed from the four bead coordinates using
! cross products.  Forces are estimated via central finite differences on
! the potential energy, which avoids the need to code the rather
! complicated analytic derivatives.  The code prints the magnitude of
! the force acting on each bead to five decimal places.  The formulas
! for the periodic dihedral energy follow the GROMACS manual, where
! V(φ) = k (1 + cos(n φ − φ_s))【31514368016376†L7207-L7223】; here the phase
! shift φ_s has been replaced by −φ0 so that V = k (1 + cos(φ + φ0)).  The
! LAMMPS manual states a similar harmonic dihedral potential
! E = K [1 + d cos(n φ)]【458590365153931†L94-L122】, and reversing the sign of
! the phase shift corresponds to replacing φ − φ0 with φ + φ0.

program dihedral_force
  implicit none
  ! use double precision throughout
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: k = 0.081_dp        ! force constant (kcal/mol)
  real(dp), parameter :: phi0 = 2.302_dp     ! ideal angle (radians)
  real(dp), parameter :: h = 1.0e-6_dp       ! finite difference step size (Å)
  integer :: i, j
  real(dp) :: positions(4,3), pos_plus(4,3), pos_minus(4,3)
  real(dp) :: forces(4,3), mag(4)

  ! Initialise the positions (Å)
  !positions = reshape([ &
  !    199.059_dp, 53.939_dp, 456.220_dp, &
  !    193.109_dp, 53.492_dp, 457.693_dp, &
  !    188.280_dp, 53.617_dp, 461.206_dp, &
  !    184.850_dp, 54.600_dp, 466.799_dp  &
  !  ], shape(positions))
  positions(1,:) = [199.059_dp, 53.939_dp, 456.220_dp]
  positions(2,:) = [193.109_dp, 53.492_dp, 457.693_dp]
  positions(3,:) = [188.280_dp, 53.617_dp, 461.206_dp]
  positions(4,:) = [184.850_dp, 54.600_dp, 466.799_dp]

  ! Compute force on each bead by finite differences
  do i = 1,4
    do j = 1,3
      pos_plus  = positions
      pos_minus = positions
      pos_plus(i,j)  = pos_plus(i,j)  + h
      pos_minus(i,j) = pos_minus(i,j) - h
      forces(i,j) = -(compute_energy(pos_plus) - compute_energy(pos_minus)) / (2.0_dp*h)
    end do
  end do

  ! Compute magnitude of force on each bead
  do i = 1,4
    mag(i) = sqrt( forces(i,1)**2 + forces(i,2)**2 + forces(i,3)**2 )
  end do

  ! Print results
  print '(a)', '{'
  do i = 1,4
     ! Format bead index and force magnitude to five decimals
     write(*,'(a,i0,a,f0.5,a)', advance='no') '  "Bead ', i, '": ', mag(i), ''
     if (i < 4) then
        print '(a)', ','
     else
        print*, ' '
     end if
  end do
  print '(a)', '}'

contains

  function compute_energy(pos) result(energy)
    ! Compute the dihedral energy V = k * (1 + cos(phi + phi0))
    implicit none
    real(dp), intent(in) :: pos(4,3)
    real(dp) :: energy
    real(dp) :: p1(3), p2(3), p3(3), p4(3)
    real(dp) :: b1(3), b2(3), b3(3)
    real(dp) :: n1(3), n2(3), n1_u(3), n2_u(3), b2_u(3), m1(3)
    real(dp) :: x, y, phi

    ! Extract bead positions
    p1 = pos(1,:)
    p2 = pos(2,:)
    p3 = pos(3,:)
    p4 = pos(4,:)
    ! Bond vectors
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    ! Normals to the planes
    n1 = cross_product(b1, b2)
    n2 = cross_product(b2, b3)
    ! Normalise vectors
    if (norm(n1) /= 0.0_dp) then
       n1_u = n1 / norm(n1)
    else
       n1_u = 0.0_dp
    end if
    if (norm(n2) /= 0.0_dp) then
       n2_u = n2 / norm(n2)
    else
       n2_u = 0.0_dp
    end if
    if (norm(b2) /= 0.0_dp) then
       b2_u = b2 / norm(b2)
    else
       b2_u = 0.0_dp
    end if
    ! Auxiliary vector
    m1 = cross_product(n1_u, b2_u)
    ! Compute x and y components for atan2
    x = dot_product(n1_u, n2_u)
    y = dot_product(m1, n2_u)
    phi = atan2(y, x)
    ! Energy
    energy = k * (1.0_dp + cos(phi - phi0))
  end function compute_energy

  function cross_product(a,b) result(c)
    ! Return cross product of two 3‑vectors
    implicit none
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  function norm(a) result(res)
    ! Euclidean norm of a 3‑vector
    implicit none
    real(dp), intent(in) :: a(3)
    real(dp) :: res
    res = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
  end function norm

end program dihedral_force