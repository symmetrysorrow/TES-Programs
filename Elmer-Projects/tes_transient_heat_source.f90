MODULE TESCircuitModule
  USE DefUtils
  IMPLICIT NONE

  ! Per-circuit state, one element per TES instance (1=unprefixed/legacy,
  ! 2=L, 3=R). Field set and initial values are identical to the SAVE
  ! scalars of the pre-refactor single-instance TESTransientHeatSource.
  TYPE CircuitState
    ! Committed state (end of the last completed timestep)
    INTEGER :: LastTimeStep = -1
    INTEGER :: LastNonlinIter = -1
    LOGICAL :: Initialized = .FALSE.
    LOGICAL :: FileStarted = .FALSE.
    REAL(KIND=dp) :: PreviousCurrent = 0.0_dp
    ! Latest coupling iterates within the current timestep
    REAL(KIND=dp) :: Current = 0.0_dp
    REAL(KIND=dp) :: Resistance = 0.0_dp
    REAL(KIND=dp) :: Power = 0.0_dp
    ! TES average temperature bookkeeping (one assembly sweep = one nonlinear iteration)
    REAL(KIND=dp) :: AverageTemperature = 0.0_dp
    REAL(KIND=dp) :: SweepTemperatureSum = 0.0_dp
    INTEGER :: SweepSampleCount = 0
    ! Aitken/secant relaxation state for the power fixed point within a timestep
    INTEGER :: IterInStep = 0
    REAL(KIND=dp) :: PrevResidual = 0.0_dp
    REAL(KIND=dp) :: Omega = 0.5_dp
    REAL(KIND=dp) :: OmegaCap = 0.5_dp
    REAL(KIND=dp) :: LastDt = -1.0_dp
  END TYPE CircuitState

  TYPE(CircuitState), SAVE :: States(3)

CONTAINS

  ! Common circuit logic shared by all TES instances. InstanceIndex selects
  ! the persistent state slot (States(InstanceIndex)); KeyPrefix selects the
  ! Constants entries (e.g. 'TES ', 'TES L ', 'TES R '); Tag is used as the
  ! function name argument to Info/Fatal; Unit is the Fortran unit number
  ! used for the series CSV (must be distinct across instances that may
  ! write within the same timestep). The numerical formulas and their
  ! execution order are unchanged from the original single-instance
  ! TESTransientHeatSource.
  SUBROUTINE TESCircuitCompute(Model, Node, Temperature, InstanceIndex, KeyPrefix, Tag, Unit, HeatSource)
    TYPE(Model_t) :: Model
    INTEGER :: Node
    REAL(KIND=dp) :: Temperature
    INTEGER :: InstanceIndex
    CHARACTER(LEN=*) :: KeyPrefix
    CHARACTER(LEN=*) :: Tag
    INTEGER :: Unit
    REAL(KIND=dp) :: HeatSource

    INTEGER :: TimeStep
    INTEGER :: NonlinIter
    LOGICAL :: SteadyMode
    CHARACTER(LEN=MAX_NAME_LEN) :: SimulationType
    REAL(KIND=dp) :: Dt
    REAL(KIND=dp) :: A
    REAL(KIND=dp) :: B
    REAL(KIND=dp) :: C
    REAL(KIND=dp) :: Discriminant
    REAL(KIND=dp) :: RawResistance
    REAL(KIND=dp) :: RawCurrent
    REAL(KIND=dp) :: RawPower
    REAL(KIND=dp) :: Residual
    REAL(KIND=dp) :: Denominator
    REAL(KIND=dp) :: CommitTime
    REAL(KIND=dp) :: I_BIAS
    REAL(KIND=dp) :: R_SH
    REAL(KIND=dp) :: L_TES
    REAL(KIND=dp) :: R0
    REAL(KIND=dp) :: R_MIN
    REAL(KIND=dp) :: ALPHA
    REAL(KIND=dp) :: BETA
    REAL(KIND=dp) :: I0
    REAL(KIND=dp) :: T0
    REAL(KIND=dp) :: TES_VOLUME
    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile
    CHARACTER(LEN=256) :: LogMessage
    INTEGER :: IoStatus
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: StateFile
    LOGICAL :: FoundState
    INTEGER :: StateUnit
    LOGICAL :: StateFileOk
    REAL(KIND=dp) :: StateAvgT
    REAL(KIND=dp) :: StateCurrent
    REAL(KIND=dp) :: StateResistance
    REAL(KIND=dp) :: StatePower
    REAL(KIND=dp) :: StatePrevCurrent
    LOGICAL :: IsRoot

    ASSOCIATE (St => States(InstanceIndex))

    ! The circuit is evaluated from this UDF during element assembly.
    ! File output is therefore restricted to rank 0 below.  MPI collectives
    ! must not be called here: ranks enter element callbacks in different
    ! orders, whereas collectives require identical call ordering.
    IsRoot = (ParEnv % MyPE == 0)

    TimeStep = GetTimeStep()
    NonlinIter = GetNonlinIter()
    SteadyMode = .FALSE.
    SimulationType = ListGetString(Model % Simulation, 'Simulation Type', Found)
    IF (Found) SteadyMode = (SimulationType(1:6) == 'steady')

    ! All constants are required: silent fallback defaults are a bug source
    ! (redesign plan, Phase 1). Generated case SIFs always provide them.
    I_BIAS = GetConstReal(Model % Constants, KeyPrefix // 'Bias Current', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Bias Current is required')
    R_SH = GetConstReal(Model % Constants, KeyPrefix // 'Shunt Resistance', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Shunt Resistance is required')
    L_TES = GetConstReal(Model % Constants, KeyPrefix // 'Inductance', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Inductance is required')
    R0 = GetConstReal(Model % Constants, KeyPrefix // 'R0', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'R0 is required')
    R_MIN = GetConstReal(Model % Constants, KeyPrefix // 'Rmin', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Rmin is required')
    ALPHA = GetConstReal(Model % Constants, KeyPrefix // 'Alpha', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Alpha is required')
    BETA = GetConstReal(Model % Constants, KeyPrefix // 'Beta', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Beta is required')
    I0 = GetConstReal(Model % Constants, KeyPrefix // 'I0', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'I0 is required')
    T0 = GetConstReal(Model % Constants, KeyPrefix // 'T0', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'T0 is required')
    TES_VOLUME = GetConstReal(Model % Constants, KeyPrefix // 'Volume', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Volume is required')
    SeriesFile = ListGetString(Model % Constants, KeyPrefix // 'Series File', Found)
    IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Series File is required')

    ! Optional: persisted circuit state (dual-TES cases only; see
    ! docs/dual_tes_plan.md). Absent for every single-pixel case, which keeps
    ! their behavior byte-identical to the pre-existing T0-based
    ! initialization. Steady-state runs WRITE the converged state here (so a
    ! restarted transient/pulse run can seed from it instead of the T0
    ! analytic estimate, which sits ~1.3 mK off the true steady operating
    ! point and otherwise causes a ~20-step Aitken-relaxation transient at
    ! the start of every restarted run). Transient/pulse runs only READ it,
    ! never write, so they cannot clobber the seed a later steady rerun
    ! would need. A distinct unit per instance (Unit+3, i.e. 94/95/96) keeps
    ! it from colliding with the series CSV units (91/92/93).
    StateFile = ListGetString(Model % Constants, KeyPrefix // 'State File', FoundState)
    StateUnit = Unit + 3

    IF (.NOT. St % Initialized) THEN
      StateFileOk = .FALSE.
      IF (.NOT. SteadyMode .AND. FoundState) THEN
        OPEN(UNIT=StateUnit, FILE=TRIM(StateFile), STATUS='OLD', ACTION='READ', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          READ(StateUnit, *, IOSTAT=IoStatus) StateAvgT, StateCurrent, StateResistance, &
            StatePower, StatePrevCurrent
          CLOSE(StateUnit)
          StateFileOk = (IoStatus == 0)
        END IF
      END IF

      IF (StateFileOk) THEN
        ! Seed the circuit from the converged steady-state solution instead
        ! of the T0-based analytic estimate below.
        St % AverageTemperature = StateAvgT
        St % Current = StateCurrent
        St % Resistance = StateResistance
        St % Power = StatePower
        St % PreviousCurrent = StatePrevCurrent
      ELSE
        ! Steady-state circuit solution at T = T0 as the initial electrical state.
        St % AverageTemperature = T0
        A = R0 * (1.0_dp - BETA)
        B = R0 * BETA / I0
        Discriminant = MAX((R_SH + A)**2 + 4.0_dp * B * I_BIAS * R_SH, 0.0_dp)
        St % Current = (SQRT(Discriminant) - (R_SH + A)) / (2.0_dp * B)
        St % Current = MAX(MIN(St % Current, I_BIAS), 0.0_dp)

        RawResistance = A + B * ABS(St % Current)
        IF (RawResistance < R_MIN) THEN
          St % Resistance = R_MIN
          St % Current = I_BIAS * R_SH / (R_SH + St % Resistance)
        ELSE
          St % Resistance = RawResistance
        END IF

        St % PreviousCurrent = St % Current
        St % Power = St % Current * St % Current * St % Resistance
      END IF

      St % LastTimeStep = TimeStep
      St % LastNonlinIter = NonlinIter
      St % IterInStep = 1
      St % PrevResidual = 0.0_dp
      St % Omega = 0.5_dp
      St % SweepTemperatureSum = 0.0_dp
      St % SweepSampleCount = 0
      St % Initialized = .TRUE.

      IF (SteadyMode .AND. FoundState .AND. IsRoot) THEN
        OPEN(UNIT=StateUnit, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          WRITE(StateUnit,'(5ES24.16)') St % AverageTemperature, St % Current, &
            St % Resistance, St % Power, St % PreviousCurrent
          CLOSE(StateUnit)
        END IF
      END IF

    ELSE IF (TimeStep /= St % LastTimeStep .OR. NonlinIter /= St % LastNonlinIter) THEN
      ! A new assembly sweep begins: refresh the TES average temperature from the
      ! sweep that just finished (i.e. the previous nonlinear iterate).
      IF (St % SweepSampleCount > 0) THEN
        St % AverageTemperature = St % SweepTemperatureSum / REAL(St % SweepSampleCount, dp)
      END IF
      St % SweepTemperatureSum = 0.0_dp
      St % SweepSampleCount = 0

      Dt = GetTimeStepSize()
      IF (Dt <= 0.0_dp) Dt = 1.0_dp

      IF (TimeStep /= St % LastTimeStep) THEN
        ! Commit the converged state of the finished timestep and log it.
        CommitTime = GetTime() - Dt
        St % PreviousCurrent = St % Current

        ! Optional transient checkpoint for a restartable single-pixel run.
        ! It is written only after a timestep has converged, so a process loss
        ! can resume from the same thermal field and circuit state.
        IF (.NOT. SteadyMode .AND. FoundState .AND. IsRoot) THEN
          OPEN(UNIT=StateUnit, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
          IF (IoStatus == 0) THEN
            WRITE(StateUnit,'(5ES24.16)') St % AverageTemperature, St % Current, &
              St % Resistance, St % Power, St % PreviousCurrent
            CLOSE(StateUnit)
          END IF
        END IF

        WRITE(LogMessage,'(A,I0,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
          'step=', St % LastTimeStep, ' T=', St % AverageTemperature, ' I_TES=', St % Current, &
          ' R_TES=', St % Resistance, ' P_TES=', St % Power
        IF (IsRoot) CALL Info(Tag, TRIM(LogMessage), Level=4)

        IF (IsRoot) THEN
          IF (.NOT. St % FileStarted) THEN
            OPEN(UNIT=Unit, FILE=TRIM(SeriesFile), STATUS='REPLACE', &
              ACTION='WRITE', IOSTAT=IoStatus)
            IF (IoStatus == 0) THEN
              WRITE(Unit,'(A)') 'time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W'
            END IF
          ELSE
            OPEN(UNIT=Unit, FILE=TRIM(SeriesFile), STATUS='OLD', &
              POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
          END IF
          IF (IoStatus == 0) THEN
            WRITE(Unit,'(ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
              CommitTime, ',', St % AverageTemperature, ',', St % Current, ',', St % Resistance, ',', St % Power
            CLOSE(Unit)
          END IF
        END IF
        St % FileStarted = .TRUE.

        ! Keep Omega: the converged relaxation factor of the previous step is
        ! the best available estimate for the next one.
        St % IterInStep = 0
        St % PrevResidual = 0.0_dp
        St % OmegaCap = 0.5_dp
        IF (St % LastDt > 0.0_dp .AND. ABS(Dt - St % LastDt) > 1.0e-9_dp*St % LastDt) THEN
          ! Timestep size changed: the coupling loop gain scales with
          ! |dP/dT|*dt/(C + G*dt), so the carried relaxation factor is no
          ! longer valid. Restart conservatively to avoid overshoot at
          ! stage boundaries of the timestep ramp.
          St % Omega = 0.1_dp
          St % OmegaCap = 0.25_dp
        END IF
        St % LastDt = Dt
        St % LastTimeStep = TimeStep
      END IF

      ! Re-solve the TES/shunt branch (backward Euler on L*dI/dt) with the latest
      ! temperature iterate. PreviousCurrent stays fixed within the timestep, so
      ! repeating this every nonlinear iteration makes the electrothermal
      ! coupling implicit instead of staggered.
      A = R0 * (1.0_dp + ALPHA * (St % AverageTemperature - T0) / T0 - BETA)
      B = R0 * BETA / I0
      C = R_SH + A + L_TES / Dt

      Discriminant = MAX(C*C + 4.0_dp*B * &
        (I_BIAS*R_SH + L_TES*St % PreviousCurrent/Dt), 0.0_dp)
      RawCurrent = (SQRT(Discriminant) - C) / (2.0_dp*B)
      RawCurrent = MAX(MIN(RawCurrent, I_BIAS), 0.0_dp)

      RawResistance = A + B*ABS(RawCurrent)
      IF (RawResistance < R_MIN) THEN
        RawResistance = R_MIN
        RawCurrent = (I_BIAS*R_SH + L_TES*St % PreviousCurrent/Dt) / &
          (R_SH + RawResistance + L_TES/Dt)
      END IF
      RawPower = RawCurrent*RawCurrent*RawResistance

      ! Aitken (secant) under-relaxation on the Joule power. The plain fixed
      ! point diverges for large dt because |dP/dT|*dt exceeds the local heat
      ! capacity; the secant update finds the stable intersection instead.
      ! Once the power residual is small the electrothermal fixed point is
      ! converged; stop updating so the heat solver can finish its own
      ! conductivity iteration without fresh perturbations. The residual the
      ! secant sees is polluted by that independent Picard drift, so the
      ! adaptive relaxation factor is additionally bounded by a trust-region
      ! cap: any residual growth halves the cap, steady contraction relaxes it.
      Residual = RawPower - St % Power
      IF (ABS(Residual) > 1.0e-6_dp * MAX(ABS(St % Power), 1.0e-30_dp)) THEN
        IF (SteadyMode) THEN
          ! Steady state: no C/dt damping, loop gain |dP/dT|/G ~ 17, and the
          ! conductivity Picard drift corrupts secant slope estimates. A small
          ! fixed factor (stable up to gain ~50) is the robust choice here.
          St % Omega = 0.04_dp
        ELSE IF (St % IterInStep > 0) THEN
          IF (ABS(Residual) > 1.5_dp*ABS(St % PrevResidual)) THEN
            St % OmegaCap = MAX(0.5_dp*St % OmegaCap, 0.02_dp)
          ELSE
            St % OmegaCap = MIN(1.3_dp*St % OmegaCap, 1.0_dp)
          END IF
          Denominator = Residual - St % PrevResidual
          IF (ABS(Denominator) > 1.0e-2_dp*ABS(Residual)) THEN
            St % Omega = -St % Omega * St % PrevResidual / Denominator
          END IF
          St % Omega = MAX(MIN(St % Omega, St % OmegaCap), 0.02_dp)
        ELSE
          ! First update of a timestep: clamp the carried-over factor too.
          St % Omega = MAX(MIN(St % Omega, St % OmegaCap), 0.02_dp)
        END IF
        St % Power = MAX(St % Power + St % Omega*Residual, 0.0_dp)
        St % PrevResidual = Residual
        St % IterInStep = St % IterInStep + 1
      END IF

      St % Current = RawCurrent
      St % Resistance = RawResistance
      St % LastNonlinIter = NonlinIter

      IF (SteadyMode .AND. FoundState .AND. IsRoot) THEN
        OPEN(UNIT=StateUnit, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          WRITE(StateUnit,'(5ES24.16)') St % AverageTemperature, St % Current, &
            St % Resistance, St % Power, St % PreviousCurrent
          CLOSE(StateUnit)
        END IF
      END IF
    END IF

    St % SweepTemperatureSum = St % SweepTemperatureSum + Temperature
    St % SweepSampleCount = St % SweepSampleCount + 1

    HeatSource = St % Power / TES_VOLUME

    END ASSOCIATE
  END SUBROUTINE TESCircuitCompute

END MODULE TESCircuitModule


FUNCTION TESTransientHeatSource(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  USE TESCircuitModule, ONLY: TESCircuitCompute
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature
  REAL(KIND=dp) :: HeatSource

  ! Instance 1: unprefixed, byte-for-byte compatible with the pre-refactor
  ! single-circuit UDF (Constants: 'TES Bias Current', ..., 'TES Series
  ! File'; series CSV on unit 91).
  CALL TESCircuitCompute(Model, Node, Temperature, 1, 'TES ', &
    'TESTransientHeatSource', 91, HeatSource)
END FUNCTION TESTransientHeatSource


FUNCTION TESTransientHeatSourceL(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  USE TESCircuitModule, ONLY: TESCircuitCompute
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature
  REAL(KIND=dp) :: HeatSource

  ! Instance 2: 'TES L ...' Constants, series CSV on unit 92.
  CALL TESCircuitCompute(Model, Node, Temperature, 2, 'TES L ', &
    'TESTransientHeatSourceL', 92, HeatSource)
END FUNCTION TESTransientHeatSourceL


FUNCTION TESTransientHeatSourceR(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  USE TESCircuitModule, ONLY: TESCircuitCompute
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature
  REAL(KIND=dp) :: HeatSource

  ! Instance 3: 'TES R ...' Constants, series CSV on unit 93.
  CALL TESCircuitCompute(Model, Node, Temperature, 3, 'TES R ', &
    'TESTransientHeatSourceR', 93, HeatSource)
END FUNCTION TESTransientHeatSourceR


FUNCTION AbsorberWindowPulseHeatSource(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature
  REAL(KIND=dp) :: HeatSource

  REAL(KIND=dp) :: TimeNow
  REAL(KIND=dp) :: TimePrev
  REAL(KIND=dp) :: Dt
  REAL(KIND=dp) :: Overlap
  REAL(KIND=dp) :: RadiusSquared
  REAL(KIND=dp) :: ENERGY
  REAL(KIND=dp) :: START_TIME
  REAL(KIND=dp) :: DURATION
  REAL(KIND=dp) :: SIGMA
  REAL(KIND=dp) :: X0
  REAL(KIND=dp) :: Y0
  REAL(KIND=dp) :: Z0
  REAL(KIND=dp) :: DISCRETE_NORM
  LOGICAL :: Found

  HeatSource = 0.0_dp

  Dt = GetTimeStepSize()
  IF (Dt <= 0.0_dp) RETURN

  ! All constants are required; the case builder provides them (center and
  ! discrete norm are computed from the mesh at build time).
  ENERGY = GetConstReal(Model % Constants, 'Pulse Energy', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Energy is required')
  START_TIME = GetConstReal(Model % Constants, 'Pulse Start Time', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Start Time is required')
  DURATION = GetConstReal(Model % Constants, 'Pulse Duration', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Duration is required')
  SIGMA = GetConstReal(Model % Constants, 'Pulse Sigma', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Sigma is required')
  X0 = GetConstReal(Model % Constants, 'Pulse Center X', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center X is required')
  Y0 = GetConstReal(Model % Constants, 'Pulse Center Y', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center Y is required')
  Z0 = GetConstReal(Model % Constants, 'Pulse Center Z', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center Z is required')
  ! FE integral of the nodal-sampled Gaussian over the absorber mesh;
  ! normalizing with it deposits exactly ENERGY regardless of how coarsely
  ! the mesh resolves the spatial profile.
  DISCRETE_NORM = GetConstReal(Model % Constants, 'Pulse Discrete Norm', Found)
  IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Discrete Norm is required')

  ! Rectangular pulse in time: distribute ENERGY uniformly over
  ! [START_TIME, START_TIME + DURATION]. Each timestep receives the exact
  ! fraction of the window it overlaps, so the total is preserved for any
  ! timestep staging.
  TimeNow = GetTime()
  TimePrev = TimeNow - Dt
  Overlap = MIN(TimeNow, START_TIME + DURATION) - MAX(TimePrev, START_TIME)
  IF (Overlap <= 0.0_dp) RETURN

  RadiusSquared = (Model % Nodes % x(Node)-X0)**2 + &
                  (Model % Nodes % y(Node)-Y0)**2 + &
                  (Model % Nodes % z(Node)-Z0)**2
  HeatSource = ENERGY * (Overlap/DURATION) * &
               EXP(-RadiusSquared/(2.0_dp*SIGMA**2)) / (DISCRETE_NORM * Dt)
END FUNCTION AbsorberWindowPulseHeatSource
