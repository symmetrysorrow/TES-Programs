! MPI-safe lumped TES circuit coupling.
!
! TESParallelCircuitSolver is a HeatSolver "Nonlinear Pre Solver": Elmer
! calls it before each nonlinear assembly sweep.  Therefore every rank enters
! the reductions below in the same order.  TESParallelHeatSource only reads
! CircuitPower and must remain communication-free because it is called during
! element assembly in rank-dependent orders.
MODULE TESParallelCircuitModule
  USE DefUtils
  IMPLICIT NONE
  LOGICAL, SAVE :: CircuitInitialized = .FALSE.
  LOGICAL, SAVE :: SeriesStarted = .FALSE.
  LOGICAL, SAVE :: IterationSeriesStarted = .FALSE.
  INTEGER, SAVE :: CircuitTimeStep = -1
  INTEGER, SAVE :: CircuitNonlinIter = -1
  INTEGER, SAVE :: CircuitIterInStep = 0
  REAL(KIND=dp), SAVE :: CircuitTemperature = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitPower = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitCurrent = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitResistance = 0.0_dp
  REAL(KIND=dp), SAVE :: PreviousCurrent = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitPrevResidual = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitOmega = 0.5_dp
  REAL(KIND=dp), SAVE :: CircuitOmegaCap = 0.5_dp
  REAL(KIND=dp), SAVE :: CircuitLastDt = 0.0_dp
CONTAINS

  SUBROUTINE TESParallelCircuitSolverCore(Model, Solver, dt, TransientSimulation)
    TYPE(Model_t) :: Model
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: TemperatureVariable
    INTEGER :: t, i, n, body, localCount, p, IoStatus, TimeStep, NonlinIter
    LOGICAL :: Found, FoundSeries, FoundState, FoundIterationSeries, WriteSeries
    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile, StateFile, IterationSeriesFile
    INTEGER :: IterationUnit
    REAL(KIND=dp) :: localSum, globalSum, globalCount, localT
    REAL(KIND=dp) :: I_BIAS, R_SH, R0, R_MIN, ALPHA, BETA, I0, T0, L_TES
    REAL(KIND=dp) :: A, B, C, Discriminant, RawPower, DtLocal, Residual, Denominator
    REAL(KIND=dp) :: StateTemperature, StateCurrent, StateResistance, StatePower, StatePrevious
    REAL(KIND=dp) :: StateLoaded

    CALL Info('TESParallelCircuitSolver', 'entered nonlinear pre-solver', Level=4)
    TemperatureVariable => VariableGet(Model % Mesh % Variables, 'Temperature')
    IF (.NOT. ASSOCIATED(TemperatureVariable)) CALL Fatal('TESParallelCircuitSolver', 'Temperature variable not found')
    CALL Info('TESParallelCircuitSolver', 'temperature variable found', Level=4)

    I_BIAS = GetConstReal(Model % Constants, 'TES Bias Current', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Bias Current is required')
    R_SH = GetConstReal(Model % Constants, 'TES Shunt Resistance', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Shunt Resistance is required')
    R0 = GetConstReal(Model % Constants, 'TES R0', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES R0 is required')
    R_MIN = GetConstReal(Model % Constants, 'TES Rmin', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Rmin is required')
    ALPHA = GetConstReal(Model % Constants, 'TES Alpha', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Alpha is required')
    BETA = GetConstReal(Model % Constants, 'TES Beta', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Beta is required')
    I0 = GetConstReal(Model % Constants, 'TES I0', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES I0 is required')
    T0 = GetConstReal(Model % Constants, 'TES T0', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES T0 is required')
    L_TES = GetConstReal(Model % Constants, 'TES Inductance', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Inductance is required')
    SeriesFile = ListGetString(Model % Constants, 'TES Series File', FoundSeries)
    IterationSeriesFile = ListGetString(Model % Constants, 'TES Iteration Series File', FoundIterationSeries)
    StateFile = ListGetString(Model % Constants, 'TES State File', FoundState)
    body = ListGetInteger(GetSolverParams(), 'TES Body ID', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Body ID is required')
    WriteSeries = GetLogical(GetSolverParams(), 'TES Write Series', Found)
    IF (.NOT. Found) WriteSeries = .TRUE.
    TimeStep = GetTimeStep()
    NonlinIter = GetNonlinIter()
    CALL Info('TESParallelCircuitSolver', 'parameters loaded', Level=4)
    CALL Info('TESParallelCircuitSolver', 'series=' // TRIM(SeriesFile) // &
      ' found=' // MERGE('T', 'F', FoundSeries) // ' write=' // MERGE('T', 'F', WriteSeries), Level=4)
    CALL Info('TESParallelCircuitSolver', 'iteration series=' // TRIM(IterationSeriesFile) // &
      ' found=' // MERGE('T', 'F', FoundIterationSeries), Level=4)

    ! These two collectives are always first, on every rank and every call.
    localSum = 0.0_dp
    localCount = 0
    DO t = 1, Model % Mesh % NumberOfBulkElements
      Element => Model % Mesh % Elements(t)
      IF (Element % BodyId /= body) CYCLE
      n = GetElementNOFNodes(Element)
      IF (n <= 0) CYCLE
      localT = 0.0_dp
      DO i = 1, n
        p = TemperatureVariable % Perm(Element % NodeIndexes(i))
        IF (p > 0) localT = localT + TemperatureVariable % Values(p)
      END DO
      localSum = localSum + localT / REAL(n, dp)
      localCount = localCount + 1
    END DO
    globalSum = ParallelReduction(localSum)
    globalCount = ParallelReduction(REAL(localCount, dp))
    IF (globalCount <= 0.0_dp) CALL Fatal('TESParallelCircuitSolver', 'No TES elements found')
    localT = globalSum / globalCount
    CALL Info('TESParallelCircuitSolver', 'TES temperature reduced', Level=4)

    IF (.NOT. CircuitInitialized) THEN
      ! Read only on root, then broadcast through reductions.  Every rank
      ! nevertheless executes all six reductions in this exact order.
      StateTemperature = 0.0_dp
      StateCurrent = 0.0_dp
      StateResistance = 0.0_dp
      StatePower = 0.0_dp
      StatePrevious = 0.0_dp
      StateLoaded = 0.0_dp
      IF (FoundState .AND. ParEnv % MyPE == 0) THEN
        OPEN(UNIT=98, FILE=TRIM(StateFile), STATUS='OLD', ACTION='READ', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          READ(98, *, IOSTAT=IoStatus) StateTemperature, StateCurrent, StateResistance, &
            StatePower, StatePrevious
          CLOSE(98)
          IF (IoStatus == 0) StateLoaded = 1.0_dp
        END IF
      END IF
      StateLoaded = ParallelReduction(StateLoaded)
      StateTemperature = ParallelReduction(StateTemperature)
      StateCurrent = ParallelReduction(StateCurrent)
      StateResistance = ParallelReduction(StateResistance)
      StatePower = ParallelReduction(StatePower)
      StatePrevious = ParallelReduction(StatePrevious)
      IF (StateLoaded > 0.5_dp) THEN
        CircuitTemperature = StateTemperature
        CircuitCurrent = StateCurrent
        CircuitResistance = StateResistance
        CircuitPower = StatePower
        PreviousCurrent = StatePrevious
      ELSE
        CircuitTemperature = localT
        A = R0 * (1.0_dp + ALPHA * (CircuitTemperature - T0) / T0 - BETA)
        B = R0 * BETA / I0
        Discriminant = MAX((R_SH + A)**2 + 4.0_dp * B * I_BIAS * R_SH, 0.0_dp)
        CircuitCurrent = (SQRT(Discriminant) - (R_SH + A)) / (2.0_dp * B)
        CircuitCurrent = MAX(MIN(CircuitCurrent, I_BIAS), 0.0_dp)
        CircuitResistance = MAX(A + B * ABS(CircuitCurrent), R_MIN)
        IF (CircuitResistance == R_MIN) CircuitCurrent = I_BIAS * R_SH / (R_SH + CircuitResistance)
        CircuitPower = CircuitCurrent * CircuitCurrent * CircuitResistance
        PreviousCurrent = CircuitCurrent
      END IF
      ! A loaded state contains the committed current from the preceding
      ! timestep.  Do not replace it before the first backward-Euler solve.
      CircuitTimeStep = TimeStep
      CircuitNonlinIter = NonlinIter
      CircuitIterInStep = 1
      CircuitPrevResidual = 0.0_dp
      CircuitOmega = 0.5_dp
      CircuitOmegaCap = 0.5_dp
      CircuitInitialized = .TRUE.
      ! A one-step or steady circuit solve can return after initialization;
      ! emit that initial electrical state as well so HYPRE/CPU/GPU cases do
      ! not lose their only machine-readable observable row.
      IF (WriteSeries .AND. FoundSeries .AND. ParEnv % MyPE == 0) THEN
        IF (.NOT. SeriesStarted) THEN
          OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
          IF (IoStatus == 0) WRITE(97,'(A)') &
            'time_s,time_step,nonlinear_iter,tes_temperature_K,tes_current_A,' // &
            'tes_resistance_ohm,tes_power_W,bias_current_A,shunt_resistance_ohm'
          SeriesStarted = .TRUE.
        ELSE
          OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
        END IF
        IF (IoStatus == 0) THEN
          WRITE(97,'(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
            GetTime(), ',', TimeStep, ',', NonlinIter, ',', CircuitTemperature, ',', CircuitCurrent, ',', &
            CircuitResistance, ',', CircuitPower, ',', I_BIAS, ',', R_SH
          CLOSE(97)
        END IF
      END IF
      IF (WriteSeries .AND. FoundIterationSeries .AND. ParEnv % MyPE == 0) THEN
        IterationUnit = 96
        IF (.NOT. IterationSeriesStarted) THEN
          OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
          IF (IoStatus == 0) WRITE(IterationUnit,'(A)') &
            'time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,' // &
            'raw_current_A,tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,' // &
            'relaxed_power_W'
        ELSE
          OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
        END IF
        IF (IoStatus == 0) THEN
          WRITE(IterationUnit,'(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,' // &
            'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
            GetTime(), ',', TimeStep, ',', NonlinIter, ',', CircuitTemperature, ',', PreviousCurrent, ',', &
            CircuitCurrent, ',', CircuitResistance, ',', CircuitPower, ',', 0.0_dp, ',', CircuitOmega, ',', &
            CircuitOmegaCap, ',', CircuitPower
          CLOSE(IterationUnit)
          IterationSeriesStarted = .TRUE.
        END IF
      END IF
      ! Preserve the legacy serial UDF's initial electrical state for the
      ! first assembly.  Its first update happens on the next nonlinear sweep.
      RETURN
    ELSE IF (TimeStep == CircuitTimeStep .AND. NonlinIter == CircuitNonlinIter) THEN
      RETURN
    END IF

    IF (TimeStep /= CircuitTimeStep) THEN
      ! This is the legacy serial UDF's timestep commit: PreviousCurrent is
      ! fixed for every nonlinear sweep of the new BDF step.
      PreviousCurrent = CircuitCurrent
      CircuitIterInStep = 0
      CircuitPrevResidual = 0.0_dp
      DtLocal = GetTimeStepSize()
      CircuitOmegaCap = 0.5_dp
      IF (CircuitLastDt > 0.0_dp .AND. ABS(DtLocal-CircuitLastDt) > 1.0e-9_dp*CircuitLastDt) THEN
        CircuitOmega = 0.1_dp
        CircuitOmegaCap = 0.25_dp
      END IF
      CircuitLastDt = DtLocal
      CircuitTimeStep = TimeStep
    END IF

    ! Same circuit equation and Aitken nonlinear update as the frozen serial
    ! reference.  Only the TES temperature average differs operationally: it
    ! is reduced collectively before this routine is entered.
    CircuitTemperature = localT
    A = R0 * (1.0_dp + ALPHA * (CircuitTemperature - T0) / T0 - BETA)
    B = R0 * BETA / I0
    DtLocal = GetTimeStepSize()
    IF (TransientSimulation .AND. DtLocal > 0.0_dp) THEN
      C = R_SH + A + L_TES / DtLocal
      Discriminant = MAX(C*C + 4.0_dp*B*(I_BIAS*R_SH + L_TES*PreviousCurrent/DtLocal), 0.0_dp)
      CircuitCurrent = (SQRT(Discriminant) - C) / (2.0_dp * B)
    ELSE
      Discriminant = MAX((R_SH + A)**2 + 4.0_dp * B * I_BIAS * R_SH, 0.0_dp)
      CircuitCurrent = (SQRT(Discriminant) - (R_SH + A)) / (2.0_dp * B)
    END IF
    CircuitCurrent = MAX(MIN(CircuitCurrent, I_BIAS), 0.0_dp)
    CircuitResistance = MAX(A + B * ABS(CircuitCurrent), R_MIN)
    IF (CircuitResistance == R_MIN) CircuitCurrent = I_BIAS * R_SH / (R_SH + CircuitResistance)
    RawPower = CircuitCurrent * CircuitCurrent * CircuitResistance
    Residual = RawPower - CircuitPower
    IF (ABS(Residual) > 1.0e-6_dp * MAX(ABS(CircuitPower), 1.0e-30_dp)) THEN
      IF (CircuitIterInStep > 0) THEN
        IF (ABS(Residual) > 1.5_dp*ABS(CircuitPrevResidual)) THEN
          CircuitOmegaCap = MAX(0.5_dp*CircuitOmegaCap, 0.02_dp)
        ELSE
          CircuitOmegaCap = MIN(1.3_dp*CircuitOmegaCap, 1.0_dp)
        END IF
        Denominator = Residual - CircuitPrevResidual
        IF (ABS(Denominator) > 1.0e-2_dp*ABS(Residual)) THEN
          CircuitOmega = -CircuitOmega*CircuitPrevResidual/Denominator
        END IF
        CircuitOmega = MAX(MIN(CircuitOmega, CircuitOmegaCap), 0.02_dp)
      ELSE
        CircuitOmega = MAX(MIN(CircuitOmega, CircuitOmegaCap), 0.02_dp)
      END IF
      CircuitPower = MAX(CircuitPower + CircuitOmega*Residual, 0.0_dp)
      CircuitPrevResidual = Residual
      CircuitIterInStep = CircuitIterInStep + 1
    END IF
    CircuitNonlinIter = NonlinIter

    ! The state format is shared with tes_transient_heat_source.f90.
    ! Steady cases update it on root after every nonlinear circuit update.
    IF (.NOT. TransientSimulation .AND. FoundState .AND. ParEnv % MyPE == 0) THEN
      OPEN(UNIT=98, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
      IF (IoStatus == 0) THEN
        WRITE(98,'(5ES24.16)') CircuitTemperature, CircuitCurrent, CircuitResistance, CircuitPower, PreviousCurrent
        CLOSE(98)
      END IF
    END IF

    ! Emit the canonical electrical row for both steady and transient cases.
    ! The old writer was transient-only, which made a steady electrical
    ! reference impossible to compare and hid missing UDF execution in smoke
    ! cases.  The explicit time-step/nonlinear columns make repeated rows
    ! machine-readable instead of requiring log scraping.
    IF (WriteSeries .AND. FoundSeries .AND. ParEnv % MyPE == 0) THEN
      IF (.NOT. SeriesStarted) THEN
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) WRITE(97,'(A)') &
          'time_s,time_step,nonlinear_iter,tes_temperature_K,tes_current_A,' // &
          'tes_resistance_ohm,tes_power_W,bias_current_A,shunt_resistance_ohm'
        SeriesStarted = .TRUE.
      ELSE
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      END IF
      IF (IoStatus == 0) THEN
        WRITE(97,'(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
          GetTime(), ',', TimeStep, ',', NonlinIter, ',', CircuitTemperature, ',', CircuitCurrent, ',', &
          CircuitResistance, ',', CircuitPower, ',', I_BIAS, ',', R_SH
        CLOSE(97)
      END IF
    END IF
    IF (WriteSeries .AND. FoundIterationSeries .AND. ParEnv % MyPE == 0) THEN
      IterationUnit = 96
      IF (.NOT. IterationSeriesStarted) THEN
        OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) WRITE(IterationUnit,'(A)') &
          'time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,' // &
          'raw_current_A,tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,' // &
          'relaxed_power_W'
      ELSE
        OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      END IF
      IF (IoStatus == 0) THEN
        WRITE(IterationUnit,'(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,' // &
          'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
          GetTime(), ',', TimeStep, ',', NonlinIter, ',', CircuitTemperature, ',', PreviousCurrent, ',', &
          CircuitCurrent, ',', CircuitResistance, ',', RawPower, ',', Residual, ',', CircuitOmega, ',', &
          CircuitOmegaCap, ',', CircuitPower
        CLOSE(IterationUnit)
        IterationSeriesStarted = .TRUE.
      END IF
    END IF
  END SUBROUTINE TESParallelCircuitSolverCore

  FUNCTION TESParallelHeatSourceCore(Model, Node, Temperature) RESULT(HeatSource)
    TYPE(Model_t) :: Model
    INTEGER :: Node
    REAL(KIND=dp) :: Temperature, HeatSource, Volume, Power
    LOGICAL :: Found
    Volume = GetConstReal(Model % Constants, 'TES Volume', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelHeatSource', 'TES Volume is required')
    ! No MPI calls are allowed in this assembly callback.  The intrinsic
    ! HeatSolve circuit hook publishes its collective, rank-consistent power
    ! through Model constants; do not use this module's private state here.
    Power = GetConstReal(Model % Constants, 'TES Parallel Power', Found)
    IF (.NOT. Found) Power = CircuitPower
    HeatSource = Power / Volume
  END FUNCTION TESParallelHeatSourceCore
END MODULE TESParallelCircuitModule

SUBROUTINE TESParallelCircuitSolver(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  USE TESParallelCircuitModule, ONLY: TESParallelCircuitSolverCore
  TYPE(Model_t) :: Model
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  CALL TESParallelCircuitSolverCore(Model, Solver, dt, TransientSimulation)
END SUBROUTINE TESParallelCircuitSolver

FUNCTION TESParallelHeatSource(Model, Node, Temperature) RESULT(HeatSource)
  USE DefUtils
  USE TESParallelCircuitModule, ONLY: TESParallelHeatSourceCore
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature, HeatSource
  HeatSource = TESParallelHeatSourceCore(Model, Node, Temperature)
END FUNCTION TESParallelHeatSource
