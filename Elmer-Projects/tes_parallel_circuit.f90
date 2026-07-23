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
  INTEGER, SAVE :: CircuitTimeStep = -1
  REAL(KIND=dp), SAVE :: CircuitTemperature = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitPower = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitCurrent = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitResistance = 0.0_dp
  REAL(KIND=dp), SAVE :: PreviousCurrent = 0.0_dp
CONTAINS

  SUBROUTINE TESParallelCircuitSolverCore(Model, Solver, dt, TransientSimulation)
    TYPE(Model_t) :: Model
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: TemperatureVariable
    INTEGER :: t, i, n, body, localCount, p, IoStatus, TimeStep
    LOGICAL :: Found, FoundSeries, FoundState, WriteSeries
    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile, StateFile
    REAL(KIND=dp) :: localSum, globalSum, globalCount, localT
    REAL(KIND=dp) :: I_BIAS, R_SH, R0, R_MIN, ALPHA, BETA, I0, T0, L_TES
    REAL(KIND=dp) :: A, B, C, Discriminant, RawPower, DtLocal
    REAL(KIND=dp) :: StateTemperature, StateCurrent, StateResistance, StatePower, StatePrevious
    REAL(KIND=dp) :: StateLoaded

    CALL Info('TESParallelCircuitSolver', 'entered nonlinear pre-solver', Level=4)
    TemperatureVariable => VariableGet(Solver % Mesh % Variables, 'Temperature')
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
    StateFile = ListGetString(Model % Constants, 'TES State File', FoundState)
    body = ListGetInteger(GetSolverParams(), 'TES Body ID', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Body ID is required')
    WriteSeries = GetLogical(GetSolverParams(), 'TES Write Series', Found)
    IF (.NOT. Found) WriteSeries = .TRUE.
    TimeStep = GetTimeStep()
    CALL Info('TESParallelCircuitSolver', 'parameters loaded', Level=4)

    ! These two collectives are always first, on every rank and every call.
    localSum = 0.0_dp
    localCount = 0
    DO t = 1, Solver % NumberOfActiveElements
      Element => GetActiveElement(t, Solver)
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
      CircuitInitialized = .TRUE.
    ELSE IF (TimeStep /= CircuitTimeStep) THEN
      ! Commit precisely once at a timestep boundary.  It must remain fixed
      ! through all nonlinear pre-solver calls within this timestep.
      PreviousCurrent = CircuitCurrent
      CircuitTimeStep = TimeStep
    END IF

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
    CircuitPower = MAX(RawPower, 0.0_dp)

    ! The state format is shared with tes_transient_heat_source.f90.
    ! Steady cases update it on root after every nonlinear circuit update.
    IF (.NOT. TransientSimulation .AND. FoundState .AND. ParEnv % MyPE == 0) THEN
      OPEN(UNIT=98, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
      IF (IoStatus == 0) THEN
        WRITE(98,'(5ES24.16)') CircuitTemperature, CircuitCurrent, CircuitResistance, CircuitPower, PreviousCurrent
        CLOSE(98)
      END IF
    END IF

    IF (TransientSimulation .AND. WriteSeries .AND. FoundSeries .AND. ParEnv % MyPE == 0) THEN
      IF (.NOT. SeriesStarted) THEN
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) WRITE(97,'(A)') 'time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W'
        SeriesStarted = .TRUE.
      ELSE
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      END IF
      IF (IoStatus == 0) THEN
        WRITE(97,'(ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
          GetTime(), ',', CircuitTemperature, ',', CircuitCurrent, ',', CircuitResistance, ',', CircuitPower
        CLOSE(97)
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
