! Fortran program to calculate dihedral forces based on the GROMACS manual
! and the principles from the paper "Force and Virial of Torsional-Angle-Dependent
! Potentials" by H. Bekker, H. J. C. Berendsen, W. F. van Gunsteren (1995).

program DihedralForceCalculator
  implicit none

  ! Bead coordinates
  real, dimension(3) :: r1 = [199.059, 53.939, 456.220]
  real, dimension(3) :: r2 = [193.109, 53.492, 457.693]
  real, dimension(3) :: r3 = [188.280, 53.617, 461.206]
  real, dimension(3) :: r4 = [184.850, 54.600, 466.799]

  ! Force constant and ideal angle
  real, parameter :: k_phi = 0.081 ! kcal/mol/rad^2
  real, parameter :: phi0 = 2.302    ! radians

  ! Bond vectors
  real, dimension(3) :: r21, r32, r43

  ! Normal vectors to the planes
  real, dimension(3) :: m, n

  ! Dihedral angle
  real :: phi

  ! Derivative of the potential
  real :: dV_dphi

  ! Forces on each bead
  real, dimension(3) :: f1, f2, f3, f4

  ! Scalar sum of force components
  real :: f1_scalar, f2_scalar, f3_scalar, f4_scalar

  ! Intermediate variables
  real :: dot_mn, norm_m_sq, norm_n_sq, norm_r32
  real :: cos_phi, sin_phi
  real, dimension(3) :: cross_mn_r32

  ! Pi
  real, parameter :: pi = 4.0 * atan(1.0)

  ! Calculate bond vectors
  r21 = r1 - r2
  r32 = r2 - r3
  r43 = r3 - r4

  ! Calculate normal vectors
  m = cross_product(r21, r32)
  n = cross_product(r32, r43)

  ! Calculate dihedral angle
  dot_mn = dot_product(m, n)
  norm_m_sq = dot_product(m, m)
  norm_n_sq = dot_product(n, n)
  cos_phi = dot_mn / sqrt(norm_m_sq * norm_n_sq)

  ! The sign of the angle is determined by the sign of (m x n) . r32
  cross_mn_r32 = cross_product(m, n)
  if (dot_product(cross_mn_r32, r32) < 0.0) then
     phi = -acos(cos_phi)
  else
     phi = acos(cos_phi)
  end if

  ! GROMACS potential (modified)
  ! V = 0.5 * k * (phi + phi0)^2
  ! The derivative with respect to phi is:
  dV_dphi = k_phi * (phi + phi0)

  ! Calculate forces based on the paper by Bekker et al.
  norm_r32 = sqrt(dot_product(r32, r32))
  f1 = -dV_dphi * (norm_r32 / norm_m_sq) * m
  f4 = dV_dphi * (norm_r32 / norm_n_sq) * n

  ! Forces on the central atoms
  f2 = -f1 + (dot_product(r21, r32) / dot_product(r32, r32)) * f1 - &
       (dot_product(r43, r32) / dot_product(r32, r32)) * f4
  f3 = -f4 - (dot_product(r21, r32) / dot_product(r32, r32)) * f1 + &
       (dot_product(r43, r32) / dot_product(r32, r32)) * f4


  ! Calculate the scalar sum of the force components
  f1_scalar = sum(f1)
  f2_scalar = sum(f2)
  f3_scalar = sum(f3)
  f4_scalar = sum(f4)

  ! Print the result in JSON format
  print '(A)', '{'
  print '("  ""Bead 1"": { ""Force"": """, F10.5, """ }")', f1_scalar
  print '("  ""Bead 2"": { ""Force"": """, F10.5, """ }")', f2_scalar
  print '("  ""Bead 3"": { ""Force"": """, F10.5, """ }")', f3_scalar
  print '("  ""Bead 4"": { ""Force"": """, F10.5, """ }")', f4_scalar
  print '(A)', '}'


contains

  function cross_product(v1, v2) result(res)
    real, dimension(3), intent(in) :: v1, v2
    real, dimension(3) :: res
    res(1) = v1(2) * v2(3) - v1(3) * v2(2)
    res(2) = v1(3) * v2(1) - v1(1) * v2(3)
    res(3) = v1(1) * v2(2) - v1(2) * v2(1)
  end function cross_product

end program DihedralForceCalculator