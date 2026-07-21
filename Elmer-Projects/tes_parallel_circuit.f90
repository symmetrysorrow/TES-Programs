! MPI-safe lumped TES circuit prototype.
!
! Unlike a body-force UDF, TESParallelCircuitSolver is invoked once per
! equation sweep.  MPI reductions are therefore collective-safe.  The heat
! source function only reads the resulting common power and never performs
! communication.
MODULE TESParallelCircuitModule
  USE DefUtils
  IMPLICIT NONE
  LOGICAL, SAVE :: CircuitInitialized = .FALSE.
  LOGICAL, SAVE :: SeriesStarted = .FALSE.
  REAL(KIND=dp), SAVE :: CircuitPower = 0.0_dp
  REAL(KIND=dp), SAVE :: CircuitCurrent = 0.0_dp
  REAL(KIND=dp), SAVE :: PreviousCurrent = 0.0_dp
CONTAINS

  SUBROUTINE TESParallelCircuitSolverCore(Model, Solver, dt, TransientSimulation)
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: TemperatureVariable
    INTEGER :: t, i, n, body, localCount, p
    LOGICAL :: Found, FoundSeries, WriteSeries
    CHARACTER(LEN=MAX_NAME_LEN) :: SeriesFile
    INTEGER :: IoStatus
    REAL(KIND=dp) :: localSum, globalSum, globalCount, localT
    REAL(KIND=dp) :: I_BIAS, R_SH, R0, R_MIN, ALPHA, BETA, I0, T0
    REAL(KIND=dp) :: A, B, C, Discriminant, Resistance, RawPower, Omega, DtLocal

    TemperatureVariable => VariableGet(Solver % Mesh % Variables, 'Temperature')
    IF (.NOT. ASSOCIATED(TemperatureVariable)) CALL Fatal('TESParallelCircuitSolver', 'Temperature variable not found')

    I_BIAS = GetConstReal(Model % Constants, 'TES Bias Current', Found)
    R_SH = GetConstReal(Model % Constants, 'TES Shunt Resistance', Found)
    R0 = GetConstReal(Model % Constants, 'TES R0', Found)
    R_MIN = GetConstReal(Model % Constants, 'TES Rmin', Found)
    ALPHA = GetConstReal(Model % Constants, 'TES Alpha', Found)
    BETA = GetConstReal(Model % Constants, 'TES Beta', Found)
    I0 = GetConstReal(Model % Constants, 'TES I0', Found)
    T0 = GetConstReal(Model % Constants, 'TES T0', Found)
    SeriesFile = ListGetString(Model % Constants, 'TES Series File', FoundSeries)
    body = ListGetInteger(GetSolverParams(), 'TES Body ID', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelCircuitSolver', 'TES Body ID is required')
    WriteSeries = GetLogical(GetSolverParams(), 'TES Write Series', Found)
    IF (.NOT. Found) WriteSeries = .TRUE.

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

    A = R0 * (1.0_dp + ALPHA * (localT - T0) / T0 - BETA)
    B = R0 * BETA / I0
    DtLocal = GetTimeStepSize()
    IF (TransientSimulation .AND. DtLocal > 0.0_dp .AND. CircuitInitialized) THEN
      C = R_SH + A + GetConstReal(Model % Constants, 'TES Inductance', Found) / DtLocal
      Discriminant = MAX(C*C + 4.0_dp*B*(I_BIAS*R_SH + &
        GetConstReal(Model % Constants, 'TES Inductance', Found)*PreviousCurrent/DtLocal), 0.0_dp)
      CircuitCurrent = (SQRT(Discriminant) - C) / (2.0_dp * B)
    ELSE
      Discriminant = MAX((R_SH + A)**2 + 4.0_dp * B * I_BIAS * R_SH, 0.0_dp)
      CircuitCurrent = (SQRT(Discriminant) - (R_SH + A)) / (2.0_dp * B)
    END IF
    CircuitCurrent = MAX(MIN(CircuitCurrent, I_BIAS), 0.0_dp)
    Resistance = MAX(A + B * ABS(CircuitCurrent), R_MIN)
    IF (Resistance == R_MIN) CircuitCurrent = I_BIAS * R_SH / (R_SH + Resistance)
    RawPower = CircuitCurrent * CircuitCurrent * Resistance

    IF (.NOT. CircuitInitialized) THEN
      CircuitPower = RawPower
      CircuitInitialized = .TRUE.
    ELSE
      Omega = 0.20_dp
      CircuitPower = MAX(CircuitPower + Omega * (RawPower - CircuitPower), 0.0_dp)
    END IF
    PreviousCurrent = CircuitCurrent
    IF (TransientSimulation .AND. WriteSeries .AND. FoundSeries .AND. ParEnv % MyPE == 0) THEN
      IF (.NOT. SeriesStarted) THEN
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) WRITE(97,'(A)') &
          'time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W'
        SeriesStarted = .TRUE.
      ELSE
        OPEN(UNIT=97, FILE=TRIM(SeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      END IF
      IF (IoStatus == 0) THEN
        WRITE(97,'(ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
          GetTime(), ',', localT, ',', CircuitCurrent, ',', Resistance, ',', CircuitPower
        CLOSE(97)
      END IF
    END IF
    IF (ParEnv % MyPE == 0) WRITE(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
      'TESParallelCircuit: T=', localT, ' I=', CircuitCurrent, ' P=', CircuitPower
  END SUBROUTINE TESParallelCircuitSolverCore

  FUNCTION TESParallelHeatSourceCore(Model, Node, Temperature) RESULT(HeatSource)
    TYPE(Model_t) :: Model
    INTEGER :: Node
    REAL(KIND=dp) :: Temperature, HeatSource, Volume, SharedPower
    LOGICAL :: Found
    Volume = GetConstReal(Model % Constants, 'TES Volume', Found)
    IF (.NOT. Found) CALL Fatal('TESParallelHeatSource', 'TES Volume is required')
    SharedPower = GetConstReal(Model % Constants, 'TES Parallel Power', Found)
    IF (Found) THEN
      HeatSource = SharedPower / Volume
    ELSE
      HeatSource = CircuitPower / Volume
    END IF
  END FUNCTION TESParallelHeatSourceCore
END MODULE TESParallelCircuitModule

SUBROUTINE TESParallelCircuitSolver(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  USE TESParallelCircuitModule, ONLY: TESParallelCircuitSolverCore
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
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
