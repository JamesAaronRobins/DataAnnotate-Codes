
!  dihedral_force.f90  –  Modern Fortran (Fortran-2008)
!
!  Implements the Bekker-Berendsen–van Gunsteren torsion‐force algorithm
!  with a harmonic potential V = k (φ + φ0)**2    (note the sign change).
!-------------------------------------------------------------------------------
module vect
   implicit none
contains
   pure function dot(a,b)   result(d)
      real(8), intent(in) :: a(3), b(3)
      real(8)             :: d
      d = sum(a*b)
   end function
   pure function norm(a) result(n)
      real(8), intent(in) :: a(3)
      real(8)             :: n
      n = sqrt(dot(a,a))
   end function
   pure function cross(a,b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8)             :: c(3)
      c = [ a(2)*b(3) - a(3)*b(2),  &
            a(3)*b(1) - a(1)*b(3),  &
            a(1)*b(2) - a(2)*b(1) ]
   end function
end module

program dihedral_force
   use, intrinsic :: iso_fortran_env, only : wp=>real64
   use vect
   implicit none
   real(wp), parameter :: k = 0.081_wp        ! kcal mol-1 rad-2
   real(wp), parameter :: phi0 = 2.302_wp     ! rad
   !--- Cartesian coordinates of the 4 beads (Å) -------------------------------
   real(wp), dimension(3) :: r1 = [199.059_wp, 53.939_wp, 456.220_wp]
   real(wp), dimension(3) :: r2 = [193.109_wp, 53.492_wp, 457.693_wp]
   real(wp), dimension(3) :: r3 = [188.280_wp, 53.617_wp, 461.206_wp]
   real(wp), dimension(3) :: r4 = [184.850_wp, 54.600_wp, 466.799_wp]
   !--- Internal vectors -------------------------------------------------------
   real(wp) :: b1(3), b2(3), b3(3), a(3), b(3)
   real(wp) :: ra2, rb2, rcb, rcb2, sb, cosphi, sinphi, phi
   real(wp) :: dVdphi, factor
   real(wp) :: f1(3), f2(3), f3(3), f4(3), s12, s32
   real(wp) :: scalar(4)
   !---------------------------------------------------------------------------
   b1 = r1 - r2
   b2 = r3 - r2
   b3 = r4 - r3
   a  = cross(b1, b2)
   b  = cross(b2, b3)

   ra2 = dot(a,a)
   rb2 = dot(b,b)
   sb  = norm(b2)

   cosphi = dot(a,b)/(sqrt(ra2)*sqrt(rb2))
   sinphi = sb*dot(b1,b)/sqrt(ra2*rb2)
   phi    = atan2(sinphi, cosphi)            ! Bekker definition (−π,π]

   dVdphi = 2.0_wp*k*(phi + phi0)            ! derivative of modified potential

   !--- Forces on atom 1 and 4 -------------------------------------------------
   f1 = -dVdphi *  ( sb / ra2 ) * a
   f4 = -dVdphi * (-sb / rb2 ) * b           ! note extra minus from Bekker

   !--- Auxiliary dot products -------------------------------------------------
   s12 = dot(b1,b2)/ra2
   s32 = dot(b3,b2)/rb2

   !--- Forces on atom 2 and 3 (Bekker 1995, Eqs. 8 & 9) -----------------------
   f2 = -f1 + s12 * f1 - s32 * f4
   f3 = -f1 - f2 - f4                        ! momentum conservation

   !--- Scalar forces (requested: sum of components) ---------------------------
   scalar(1) = sum(f1)
   scalar(2) = sum(f2)
   scalar(3) = sum(f3)
   scalar(4) = sum(f4)

   print*, scalar
   !print(*,'(a,4(3f12.6))') "  F1,F2,F3,F4 (x,y,z) / kcal mol-1 Å-1:", f1, f2, f3, f4
   !print(*,'(a,3f12.6)')      "  ∑F (should be ≈0):", f1+f2+f3+f4
   !print(*,'(a)')              " "
   !print(*,'(a)')              "  Scalar (ΣFx+Fy+Fz) forces per bead:"
   !print(*,'(f10.5,2x,f10.5,2x,f10.5,2x,f10.5)') scalar
   !--- JSON block -------------------------------------------------------------
   !write(*,'(a)') 'JSON = { "Bead 1": '//trim(adjustl(tof(scalar(1))))// &
   !     ', "Bead 2": '//trim(adjustl(tof(scalar(2))))// &
   !     ', "Bead 3": '//trim(adjustl(tof(scalar(3))))// &
   !     ', "Bead 4": '//trim(adjustl(tof(scalar(4))))//' }'
!contains
!   pure function tof(x) result(s)
!!!      real(wp), intent(in) :: x
!      character(len=:), allocatable :: s
!      write(unit='(f10.5)',fmt='(f10.5)') x
!      read(unit='(f10.5)',fmt='(a)')  s
!   end function tof
end program dihedral_force