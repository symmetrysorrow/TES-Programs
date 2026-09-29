! Compiled membrane conductivity k(T) = C * max(T, Tfloor)**p.
!
! Replaces the per-integration-point MATC expression generated for the
! Membrane material.  The coefficient, exponent and floor are written into
! the Material section by build_cases.py (keyword names below), so the
! physics still comes from the project's material expression; this routine
! only evaluates it.  The floor reproduces the MATC guard
! 0.5*((T+1e-12)+abs(T-1e-12)) == max(T, 1e-12).
FUNCTION MembraneConductivity(Model, n, T) RESULT(k)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: T, k
  TYPE(ValueList_t), POINTER :: Material
  LOGICAL :: Found
  INTEGER :: i
  REAL(KIND=dp), SAVE :: Coeff = 0.0_dp, Expo = 0.0_dp, Tfloor = 1.0e-12_dp
  LOGICAL, SAVE :: Loaded = .FALSE.

  IF (.NOT. Loaded) THEN
    ! Search the materials instead of GetMaterial(): the Phase24 assembly
    ! fast path evaluates k(T) while building element metadata, before a
    ! current element is set.
    Found = .FALSE.
    DO i = 1, Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      Coeff = GetConstReal(Material, 'Membrane Conductivity Coefficient', Found)
      IF (Found) EXIT
    END DO
    IF (.NOT. Found) CALL Fatal('MembraneConductivity', &
        'No material defines "Membrane Conductivity Coefficient"')
    Expo = GetConstReal(Material, 'Membrane Conductivity Exponent', Found)
    IF (.NOT. Found) CALL Fatal('MembraneConductivity', &
        'Missing "Membrane Conductivity Exponent" in the Membrane material')
    Tfloor = GetConstReal(Material, 'Membrane Conductivity Temperature Floor', Found)
    IF (.NOT. Found) Tfloor = 1.0e-12_dp
    Loaded = .TRUE.
  END IF

  k = Coeff * MAX(T, Tfloor)**Expo
END FUNCTION MembraneConductivity
