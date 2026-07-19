FUNCTION TESHeatSource(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature
  REAL(KIND=dp) :: HeatSource

  REAL(KIND=dp) :: A
  REAL(KIND=dp) :: B
  REAL(KIND=dp) :: Discriminant
  REAL(KIND=dp) :: TesCurrent
  REAL(KIND=dp) :: TesResistance
  REAL(KIND=dp) :: RawResistance
  REAL(KIND=dp) :: TesPower
  REAL(KIND=dp) :: I_BIAS
  REAL(KIND=dp) :: R_SH
  REAL(KIND=dp) :: R0
  REAL(KIND=dp) :: R_MIN
  REAL(KIND=dp) :: ALPHA
  REAL(KIND=dp) :: BETA
  REAL(KIND=dp) :: I0
  REAL(KIND=dp) :: T0
  REAL(KIND=dp) :: TES_VOLUME
  REAL(KIND=dp), SAVE :: AverageTemperature = 0.0_dp
  REAL(KIND=dp), SAVE :: SweepTemperatureSum = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedPower = 0.0_dp
  INTEGER, SAVE :: SweepSampleCount = 0
  INTEGER, SAVE :: LastNodeSeen = -1
  LOGICAL, SAVE :: AverageInitialized = .FALSE.
  LOGICAL :: Found

  ! circuit_layout.R_tes:
  ! R = max(R_min, R0*(1 + alpha*(T_TES-T0)/T0
  !                         + beta*(abs(I_TES)-I0)/I0))
  ! circuit_layout.L1:
  ! L_TES = 1.23e-8 H in series with the TES branch.
  ! For the steady-state parallel circuit, I_TES is the positive root of
  ! I_TES*(R_sh + R(T_TES,I_TES)) = I_bias*R_sh.
  ! All constants are required: silent fallback defaults are a bug source
  ! (redesign plan, Phase 1). Generated case SIFs always provide them.
  I_BIAS = GetConstReal(Model % Constants, 'TES Bias Current', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Bias Current is required')
  R_SH = GetConstReal(Model % Constants, 'TES Shunt Resistance', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Shunt Resistance is required')
  R0 = GetConstReal(Model % Constants, 'TES R0', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES R0 is required')
  R_MIN = GetConstReal(Model % Constants, 'TES Rmin', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Rmin is required')
  ALPHA = GetConstReal(Model % Constants, 'TES Alpha', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Alpha is required')
  BETA = GetConstReal(Model % Constants, 'TES Beta', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Beta is required')
  I0 = GetConstReal(Model % Constants, 'TES I0', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES I0 is required')
  T0 = GetConstReal(Model % Constants, 'TES T0', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES T0 is required')
  TES_VOLUME = GetConstReal(Model % Constants, 'TES Volume', Found)
  IF (.NOT. Found) CALL Fatal('TESHeatSource', 'Constants: TES Volume is required')

  ! Use one TES-average temperature per assembly sweep so the electrical
  ! branch stays lumped, rather than varying point-by-point inside the TES.
  IF (.NOT. AverageInitialized) THEN
    AverageTemperature = T0
    AverageInitialized = .TRUE.
  ELSE IF (LastNodeSeen >= 0 .AND. Node <= LastNodeSeen .AND. SweepSampleCount > 0) THEN
    AverageTemperature = SweepTemperatureSum / REAL(SweepSampleCount, dp)
    SweepTemperatureSum = 0.0_dp
    SweepSampleCount = 0
  END IF

  A = R0 * (1.0_dp + ALPHA * (AverageTemperature - T0) / T0 - BETA)
  B = R0 * BETA / I0
  Discriminant = MAX((R_SH + A)**2 + 4.0_dp * B * I_BIAS * R_SH, 0.0_dp)
  TesCurrent = (SQRT(Discriminant) - (R_SH + A)) / (2.0_dp * B)
  TesCurrent = MAX(MIN(TesCurrent, I_BIAS), 0.0_dp)

  RawResistance = A + B * ABS(TesCurrent)
  IF (RawResistance < R_MIN) THEN
    TesResistance = R_MIN
    TesCurrent = I_BIAS * R_SH / (R_SH + TesResistance)
  ELSE
    TesResistance = RawResistance
  END IF

  TesPower = TesCurrent * TesCurrent * TesResistance
  CachedPower = TesPower

  SweepTemperatureSum = SweepTemperatureSum + Temperature
  SweepSampleCount = SweepSampleCount + 1
  LastNodeSeen = Node

  HeatSource = CachedPower / TES_VOLUME
END FUNCTION TESHeatSource
