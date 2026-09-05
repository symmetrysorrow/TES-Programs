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
    LOGICAL :: ParametersInitialized = .FALSE.
    LOGICAL :: HasSeriesFile = .FALSE.
    LOGICAL :: HasStateFile = .FALSE.
    LOGICAL :: SteadyMode = .FALSE.
    REAL(KIND=dp) :: BiasCurrent = 0.0_dp
    REAL(KIND=dp) :: ShuntResistance = 0.0_dp
    REAL(KIND=dp) :: Inductance = 0.0_dp
    REAL(KIND=dp) :: R0 = 0.0_dp
    REAL(KIND=dp) :: RMin = 0.0_dp
    REAL(KIND=dp) :: Alpha = 0.0_dp
    REAL(KIND=dp) :: Beta = 0.0_dp
    REAL(KIND=dp) :: I0 = 0.0_dp
    REAL(KIND=dp) :: T0 = 0.0_dp
    REAL(KIND=dp) :: Volume = 0.0_dp
    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile = ''
    CHARACTER(LEN=MAX_NAME_LEN) :: StateFile = ''
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
    IF (.NOT. St % ParametersInitialized) THEN
      SteadyMode = .FALSE.
      SimulationType = ListGetString(Model % Simulation, 'Simulation Type', Found)
      IF (Found) SteadyMode = (SimulationType(1:6) == 'steady')

      ! All constants are required: silent fallback defaults are a bug source
      ! (redesign plan, Phase 1). Generated case SIFs always provide them.
      St % BiasCurrent = GetConstReal(Model % Constants, KeyPrefix // 'Bias Current', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Bias Current is required')
      St % ShuntResistance = GetConstReal(Model % Constants, KeyPrefix // 'Shunt Resistance', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Shunt Resistance is required')
      St % Inductance = GetConstReal(Model % Constants, KeyPrefix // 'Inductance', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Inductance is required')
      St % R0 = GetConstReal(Model % Constants, KeyPrefix // 'R0', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'R0 is required')
      St % RMin = GetConstReal(Model % Constants, KeyPrefix // 'Rmin', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Rmin is required')
      St % Alpha = GetConstReal(Model % Constants, KeyPrefix // 'Alpha', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Alpha is required')
      St % Beta = GetConstReal(Model % Constants, KeyPrefix // 'Beta', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Beta is required')
      St % I0 = GetConstReal(Model % Constants, KeyPrefix // 'I0', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'I0 is required')
      St % T0 = GetConstReal(Model % Constants, KeyPrefix // 'T0', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'T0 is required')
      St % Volume = GetConstReal(Model % Constants, KeyPrefix // 'Volume', Found)
      IF (.NOT. Found) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Volume is required')
      St % SeriesFile = ListGetString(Model % Constants, KeyPrefix // 'Series File', St % HasSeriesFile)
      IF (.NOT. St % HasSeriesFile) CALL Fatal(Tag, 'Constants: ' // KeyPrefix // 'Series File is required')
      St % StateFile = ListGetString(Model % Constants, KeyPrefix // 'State File', St % HasStateFile)
      St % SteadyMode = SteadyMode
      St % ParametersInitialized = .TRUE.
      CALL Info(Tag, 'cached TES circuit constants and integration metadata', Level=4)
    END IF

    I_BIAS = St % BiasCurrent
    R_SH = St % ShuntResistance
    L_TES = St % Inductance
    R0 = St % R0
    R_MIN = St % RMin
    ALPHA = St % Alpha
    BETA = St % Beta
    I0 = St % I0
    T0 = St % T0
    TES_VOLUME = St % Volume
    SeriesFile = St % SeriesFile
    FoundState = St % HasStateFile
    StateFile = St % StateFile
    SteadyMode = St % SteadyMode

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

  ! TESCircuitCompute defers a timestep's series-CSV row until the *next*
  ! timestep's first assembly call commits it (needed to know the previous
  ! nonlinear sweep has actually finished). The true final timestep of a
  ! transient/pulse run has no such next call, so its row is otherwise never
  ! written even though the solver did compute it. Called once, after the
  ! whole simulation, by TESCircuitFinalizeAll (Exec Solver = After
  ! Simulation) to flush it. A no-op for any instance this case never used
  ! (KeyPrefix's Series File Constant absent) and, via the TransientSimulation
  ! guard in TESCircuitFinalizeAll, for steady-state runs, whose series file
  ! (if any) is intentionally never written here (docs/dual_tes_plan.md).
  SUBROUTINE TESCircuitFlush(Model, InstanceIndex, KeyPrefix, Unit)
    TYPE(Model_t) :: Model
    INTEGER :: InstanceIndex
    CHARACTER(LEN=*) :: KeyPrefix
    INTEGER :: Unit

    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile
    CHARACTER(LEN=256) :: LogMessage
    LOGICAL :: Found
    LOGICAL :: IsRoot
    INTEGER :: IoStatus
    REAL(KIND=dp) :: CommitTime

    ASSOCIATE (St => States(InstanceIndex))

    IF (.NOT. St % Initialized) RETURN

    IsRoot = (ParEnv % MyPE == 0)
    IF (.NOT. IsRoot) RETURN

    SeriesFile = ListGetString(Model % Constants, KeyPrefix // 'Series File', Found)
    IF (.NOT. Found) RETURN

    IF (St % SweepSampleCount > 0) THEN
      St % AverageTemperature = St % SweepTemperatureSum / REAL(St % SweepSampleCount, dp)
    END IF

    CommitTime = GetTime()

    WRITE(LogMessage,'(A,I0,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
      'flush: step=', St % LastTimeStep, ' T=', St % AverageTemperature, ' I_TES=', St % Current, &
      ' R_TES=', St % Resistance, ' P_TES=', St % Power
    CALL Info('TESCircuitFlush', TRIM(LogMessage), Level=4)

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
    St % FileStarted = .TRUE.

    END ASSOCIATE
  END SUBROUTINE TESCircuitFlush

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


! Solver (not a body-force UDF): wired via Exec Solver = After Simulation on a
! standalone Solver section, so it runs exactly once, after the last
! timestep, and flushes whichever of the 3 TESCircuitCompute instances this
! case actually used (see TESCircuitFlush).
SUBROUTINE TESCircuitFinalizeAll(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  USE TESCircuitModule, ONLY: TESCircuitFlush
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  IF (.NOT. TransientSimulation) RETURN

  CALL TESCircuitFlush(Model, 1, 'TES ', 91)
  CALL TESCircuitFlush(Model, 2, 'TES L ', 92)
  CALL TESCircuitFlush(Model, 3, 'TES R ', 93)
END SUBROUTINE TESCircuitFinalizeAll


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
  REAL(KIND=dp) :: TimeIntegral
  REAL(KIND=dp) :: RadiusSquared
  REAL(KIND=dp) :: ENERGY
  REAL(KIND=dp) :: START_TIME
  REAL(KIND=dp) :: DURATION
  REAL(KIND=dp) :: TRANSITION_ZONE
  REAL(KIND=dp) :: SIGMA
  REAL(KIND=dp) :: PULSE_RADIUS
  REAL(KIND=dp) :: X0
  REAL(KIND=dp) :: Y0
  REAL(KIND=dp) :: Z0
  REAL(KIND=dp) :: DISCRETE_NORM
  INTEGER :: TimeStep
  REAL(KIND=dp), SAVE :: CachedEnergy = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedStartTime = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedDuration = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedTransitionZone = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedSigmaSquared = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedPulseRadiusSquared = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedCenterX = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedCenterY = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedCenterZ = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedDiscreteNorm = 0.0_dp
  REAL(KIND=dp), SAVE :: CachedTemporalScale = 0.0_dp
  INTEGER, SAVE :: CachedPulseShape = -1
  INTEGER, SAVE :: CachedPulseTimeStep = -1
  LOGICAL, SAVE :: PulseParametersInitialized = .FALSE.
  INTEGER :: PULSE_SHAPE
  LOGICAL :: Found

  HeatSource = 0.0_dp

  IF (.NOT. PulseParametersInitialized) THEN
    ! Constants, pulse geometry, and the discrete FE normalization are
    ! immutable for a run.  This callback is entered once per absorber
    ! element/node assembly evaluation, so repeatedly looking them up was a
    ! measurable host-side cost rather than useful physics.
    CachedEnergy = GetConstReal(Model % Constants, 'Pulse Energy', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Energy is required')
    CachedStartTime = GetConstReal(Model % Constants, 'Pulse Start Time', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Start Time is required')
    CachedDuration = GetConstReal(Model % Constants, 'Pulse Duration', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Duration is required')
    CachedTransitionZone = GetConstReal(Model % Constants, 'Pulse Transition Zone', Found)
    IF (.NOT. Found) CachedTransitionZone = 0.0_dp
    IF (CachedTransitionZone < 0.0_dp) THEN
      CALL Fatal('AbsorberWindowPulseHeatSource', 'Pulse Transition Zone must be non-negative')
    END IF
    SIGMA = GetConstReal(Model % Constants, 'Pulse Sigma', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Sigma is required')
    CachedSigmaSquared = SIGMA * SIGMA
    CachedPulseShape = NINT(GetConstReal(Model % Constants, 'Pulse Shape', Found))
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Shape is required')
    PULSE_RADIUS = GetConstReal(Model % Constants, 'Pulse Radius', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Radius is required')
    CachedPulseRadiusSquared = PULSE_RADIUS * PULSE_RADIUS
    CachedCenterX = GetConstReal(Model % Constants, 'Pulse Center X', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center X is required')
    CachedCenterY = GetConstReal(Model % Constants, 'Pulse Center Y', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center Y is required')
    CachedCenterZ = GetConstReal(Model % Constants, 'Pulse Center Z', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Center Z is required')
    ! FE integral of the nodal-sampled Gaussian over the absorber mesh.
    CachedDiscreteNorm = GetConstReal(Model % Constants, 'Pulse Discrete Norm', Found)
    IF (.NOT. Found) CALL Fatal('AbsorberWindowPulseHeatSource', 'Constants: Pulse Discrete Norm is required')
    PulseParametersInitialized = .TRUE.
    CALL Info('AbsorberWindowPulseHeatSource', 'cached pulse constants and spatial normalization', Level=4)
  END IF

  Dt = GetTimeStepSize()
  IF (Dt <= 0.0_dp) RETURN
  TimeStep = GetTimeStep()
  IF (TimeStep /= CachedPulseTimeStep) THEN
    ! Temporal rectangular pulse.  The exact interval integral is evaluated
    ! once per timestep; all later callback entries reuse this scale.
    TimeNow = GetTime()
    TimePrev = TimeNow - Dt
    IF (CachedTransitionZone > 0.0_dp) THEN
      TimeIntegral = SmoothStepIntegral(TimeNow - CachedStartTime, 0.5_dp*CachedTransitionZone) - &
                     SmoothStepIntegral(TimePrev - CachedStartTime, 0.5_dp*CachedTransitionZone) - &
                     SmoothStepIntegral(TimeNow - CachedStartTime - CachedDuration, 0.5_dp*CachedTransitionZone) + &
                     SmoothStepIntegral(TimePrev - CachedStartTime - CachedDuration, 0.5_dp*CachedTransitionZone)
    ELSE
      TimeIntegral = MIN(TimeNow, CachedStartTime + CachedDuration) - MAX(TimePrev, CachedStartTime)
    END IF
    IF (TimeIntegral > 0.0_dp) THEN
      CachedTemporalScale = CachedEnergy * (TimeIntegral/CachedDuration) / &
        (CachedDiscreteNorm * Dt)
    ELSE
      CachedTemporalScale = 0.0_dp
    END IF
    CachedPulseTimeStep = TimeStep
  END IF
  IF (CachedTemporalScale <= 0.0_dp) RETURN

  RadiusSquared = (Model % Nodes % x(Node)-CachedCenterX)**2 + &
                  (Model % Nodes % y(Node)-CachedCenterY)**2 + &
                  (Model % Nodes % z(Node)-CachedCenterZ)**2
  SELECT CASE (CachedPulseShape)
  CASE (0)
    HeatSource = EXP(-RadiusSquared/(2.0_dp*CachedSigmaSquared))
  CASE (1)
    IF (RadiusSquared <= CachedPulseRadiusSquared) THEN
      HeatSource = 1.0_dp
    ELSE
      HeatSource = 0.0_dp
    END IF
  CASE DEFAULT
    CALL Fatal('AbsorberWindowPulseHeatSource', 'Pulse Shape must be 0 (Gaussian) or 1 (uniform sphere)')
  END SELECT
  HeatSource = CachedTemporalScale * HeatSource

CONTAINS

  FUNCTION SmoothStepIntegral(X, HalfZone) RESULT(Integral)
    ! Antiderivative of COMSOL flc2hs(X, HalfZone), with Integral=0 for
    ! X <= -HalfZone and Integral=X for X >= HalfZone.  Within the zone,
    ! flc2hs is 1/2 + 15*u/16 - 5*u**3/8 + 3*u**5/16, u=X/HalfZone.
    REAL(KIND=dp), INTENT(IN) :: X
    REAL(KIND=dp), INTENT(IN) :: HalfZone
    REAL(KIND=dp) :: Integral
    REAL(KIND=dp) :: U

    IF (X <= -HalfZone) THEN
      Integral = 0.0_dp
    ELSE IF (X >= HalfZone) THEN
      Integral = X
    ELSE
      U = X / HalfZone
      Integral = HalfZone * (5.0_dp/32.0_dp + U/2.0_dp + 15.0_dp*U**2/32.0_dp - &
        5.0_dp*U**4/32.0_dp + U**6/32.0_dp)
    END IF
  END FUNCTION SmoothStepIntegral
END FUNCTION AbsorberWindowPulseHeatSource
