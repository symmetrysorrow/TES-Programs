!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Module containing a solver for heat equation
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE HeatSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Lists, ONLY : ListAddInteger, ListAddNewInteger
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'HeatSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  
  Params => GetSolverParams()

  ! Default variable name
  CALL ListAddNewString( Params,'Variable','Temperature')

  ! Phase24 is deliberately opt-in.  This selects the existing vectorized
  ! CRS glue path for the legacy HeatSolve route without changing the
  ! physical model or enabling timestep/matrix-value reuse.
  CALL ListAddNewLogical( Params,'Phase24 Vector Assembly',.FALSE.)
  CALL ListAddNewLogical( Params,'Phase24 Wall Profiling',.FALSE.)
  CALL ListAddNewLogical( Params,'Phase24 Static Matrix Reuse',.FALSE.)
  CALL ListAddNewString( Params,'Phase24 Operator Lagging','disabled',.FALSE.)
  ! Disabled-by-default diagnostic for the assembled outer HeatSolve system.
  ! It intentionally excludes auxiliary preconditioner, Schur, and restricted
  ! systems from the capture contract.
  CALL ListAddNewLogical( Params,'Phase24 Outer HeatSolve Capture',.FALSE.)
  CALL ListAddNewString( Params,'Phase24 Outer HeatSolve Capture Prefix', &
       'phase24_outer_capture',.FALSE.)
  CALL ListAddNewInteger( Params,'Phase24 Outer HeatSolve Capture Timestep',1)
  CALL ListAddNewInteger( Params,'Phase24 Outer HeatSolve Capture Max Iterations',3)
  CALL ListAddNewLogical( Params,'Phase24 Full Restriction Capture',.FALSE.)
  CALL ListAddNewString( Params,'Phase24 Full Restriction Capture Prefix', &
       'phase24_full_restriction_capture',.FALSE.)
  CALL ListAddNewInteger( Params,'Phase24 Full Restriction Capture Timestep',1)
  CALL ListAddNewInteger( Params,'Phase24 Full Restriction Capture Max Iterations',3)
  CALL ListAddNewLogical( Params,'Phase24 Restart State Audit',.FALSE.)
  CALL ListAddNewLogical( Params,'Phase24 Restart Audit Fail Fast',.FALSE.)
  CALL ListAddNewString( Params,'Phase24 Restart State Audit Prefix', &
       'phase24_restart_state_audit',.FALSE.)

  ! These use one flag to call library features to compute automatically
  ! a conductivity matrix.
  IF( ListGetLogical(Params,'Calculate Conductance Matrix',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Matrix Symmetric',.TRUE.)
    CALL ListAddNewString( Params,'Constraint Modes Matrix Filename',&
        'ThermalConductanceMatrix.dat',.FALSE.)
    CALL ListRenameAllBC( Model,'Conductivity Body','Constraint Mode Temperature')
  END IF

  ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVITY
  CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
  
END SUBROUTINE HeatSolver_Init


!------------------------------------------------------------------------------
!> Subroutine for solving the energy a.k.a. heat equation in various coordinate systems.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE HeatSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DiffuseConvective
     USE DiffuseConvectiveGeneral
     USE Differentials
     USE Radiation
     USE MaterialModels
     USE Adaptive
     USE DefUtils
     USE SolverUtils, ONLY : SaveLinearSystem
     USE CRSMatrix, ONLY : CRS_ClearPositionCache, CRS_Phase24CacheStats, &
          CRS_PreparePositionCache, CRS_GlueLocalMatrixCached, CRS_AddLocalMatrixCachedTo
     USE ElementDescription, ONLY : ElementInfo

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     INTEGER, PARAMETER :: PHASE_SPATIAL_1 = 1
     INTEGER, PARAMETER :: PHASE_SPATIAL_2 = 2
     INTEGER, PARAMETER :: PHASE_TEMPORAL  = 3
     INTEGER, PARAMETER :: PHASE24_PRISM_MAX_Q = 8
    
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix

     INTEGER :: i,j,k,l,m,n,nd,t,tt,iter,k1,k2,body_id,eq_id,istat,LocalNodes,bf_id
     INTEGER :: Phase24ElementId, Phase24OuterCaptureStep, Phase24OuterCaptureMaxIter

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element,Parent,RadiationElement

     REAL(KIND=dp) :: RelativeChange, &
           Norm,PrevNorm,Text,S,C,C1,Emissivity,StefanBoltzmann, &
           ReferencePressure=0.0d0, SpecificHeatRatio

     CHARACTER(LEN=MAX_NAME_LEN) :: RadiationFlag,ConvectionFlag

     INTEGER :: PhaseChangeModel
     CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel, StabilizeFlag, VarName

     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles,NewtonLinearization = .FALSE., &
         Found, FoundCp, FoundK, GotIt, HeatFluxBC, HeatGapBC, GotMeltPoint, IsRadiation, IsRadiosity, InfBC
! Which compressibility model is used
     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, ConvectionField
     INTEGER :: CompressibilityModel

     LOGICAL :: AllocationsDone = .FALSE.,PhaseSpatial=.FALSE., &
        PhaseChange=.FALSE., CheckLatentHeatRelease=.FALSE., FirstTime, &
        SmartHeaterControl, IntegralHeaterControl, HeaterControlLocal, SmartTolReached=.FALSE., &
        TransientHeaterControl, SmartHeaterAverage, ConstantBulk, SaveBulk, &
	TransientAssembly, Converged, AnyMultiply, NeedFlowSol
     LOGICAL :: Phase24VectorAssembly
     LOGICAL :: Phase24OuterCapture
     LOGICAL :: Phase24Eligible
     LOGICAL :: Phase24WallProfiling
     LOGICAL :: Phase24UseCachedMaterial
     LOGICAL :: Phase24StaticLocal
     LOGICAL :: Phase24StaticAssemblyActive
     LOGICAL :: Phase24PrismAssemblyActive
     LOGICAL :: Phase24StaticMatrixReuse
     LOGICAL :: Phase24StaticCacheHit
     LOGICAL :: Phase24RHSOnly, Phase24MatrixChanged
     LOGICAL :: Phase24BoundaryMatrixInvariantNext
     LOGICAL :: Phase24FastElement
     LOGICAL :: Phase24LaggingEnabled, Phase24LaggingActive, Phase24LagCandidate
     LOGICAL :: Phase24UsedLaggedOperator, Phase24RecoveryThisSolve
     LOGICAL :: Phase24ForceRebuild, Phase24FinalConfirmationDone, Phase24FallbackRebuild
     LOGICAL :: Phase24LHSUnsafe
     LOGICAL, POINTER :: SmartHeaters(:), IntegralHeaters(:)

     TYPE(Variable_t), POINTER :: TempSol,FlowSol,HeatSol,CurrentSol, MeshSol, DensitySol
     TYPE(Variable_t), POINTER :: TimeVar
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER, POINTER :: TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:)

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter,MinNonlinearIter,MDOFs, &
         SmartHeaterBC, SmartHeaterNode, DoneTime=0, bc_elem, nb, NOFactive
     INTEGER :: Phase24LagMode, Phase24LagMaxReuse, Phase24IterationLimit, Phase24LagReuseCount, &
          Phase24LaggedReuseThisTimestep, Phase24KrylovIterations, Phase24KrylovBaseline, &
          Phase24KrylovAbsoluteGrowth
     REAL(KIND=dp) :: TESTempTol, TESMeanNow, FirstSolveTol, SavedLinTol, BDFStepStart
     REAL(KIND=dp), POINTER :: BDFResetTimes(:,:)
     INTEGER :: kReset
     LOGICAL :: FirstSolveResidual, SavedResidualMode, GotSavedTol
     REAL(KIND=dp) :: NonlinearTol,NewtonTol,SmartTol,Relax, &
            SaveRelax,dt,dt0,CumulativeTime, VisibleFraction, PowerScaling=1.0, PrevPowerScaling=1.0, &
            PowerRelax, PowerTimeScale, PowerSensitivity, xave, yave, Normal(3), &
	    dist, mindist, ControlPoint(3), HeatTransferMultiplier
     REAL(KIND=dp) :: Phase24LagThreshold, Phase24LagResidualGrowth, &
          Phase24CoeffRelChange, Phase24PreviousRelativeChange, Phase24KrylovRelativeGrowth
     REAL(KIND=dp) :: Phase24CurrentA0, Phase24CurrentA1, Phase24CurrentA2, &
          Phase24CurrentHOld, Phase24CurrentRatio
     INTEGER :: Phase24CurrentOrder
     LOGICAL :: Phase24AdaptiveDebug
     LOGICAL :: Phase24BDF2Predictor
     LOGICAL :: Phase24BDF2PredictorActive, Phase24PredictorHistoryReady
     REAL(KIND=dp) :: Phase24PredictorDelta, Phase24PredictorMaxDelta
     CHARACTER(LEN=MAX_NAME_LEN) :: Phase24LagModeName, Phase24OuterCapturePrefix

     REAL(KIND=dp), POINTER :: Temperature(:),PrevTemperature(:),FlowSolution(:), &
       ElectricCurrent(:), PhaseChangeIntervals(:,:),ForceVector(:), &
       PrevSolution(:), HC(:), Hwrk(:,:,:),MeshVelocity(:), XX(:), YY(:),ForceHeater(:),&
       RealWork(:,:)

     REAL(KIND=dp), ALLOCATABLE :: vals(:)
     REAL(KIND=dp) :: Jx,Jy,Jz,JAbs, Power, MeltPoint, IntHeatSource

     LOGICAL, ALLOCATABLE, SAVE :: Phase24BodyMaterialKnown(:), Phase24BodyMaterialCached(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24BodyHeatCapacity(:), Phase24BodyHeatConductivity(:)
     INTEGER, SAVE :: Phase24BodyCacheSize = 0

     ! Persistent fixed-mesh metadata for the opt-in scalar P1 tetra path.
     ! The generic solver remains authoritative for every element classified
     ! as dynamic or generic.  These arrays are deliberately contiguous so
     ! the arithmetic-only portion can be consumed by an OpenMP loop later.
     INTEGER, ALLOCATABLE, SAVE :: Phase24ElementClass(:), Phase24ElementBody(:)
     INTEGER, ALLOCATABLE, SAVE :: Phase24ElementPart(:)
     INTEGER, ALLOCATABLE, SAVE :: Phase24ElementNodes(:,:), Phase24ElementEq(:,:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24ElementGrad(:,:,:), Phase24ElementVolume(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24StaticStiff(:,:,:), Phase24SpatialMass(:,:,:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24StaticK(:), Phase24StaticRhoCp(:)
     INTEGER, ALLOCATABLE, SAVE :: Phase24StaticElements(:), Phase24DynamicElements(:), &
          Phase24GenericElements(:), Phase24PrismElements(:), Phase24PrismSlot(:), &
          Phase24PrismNQ(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24PrismGrad(:,:,:,:), Phase24PrismWeight(:,:), &
          Phase24PrismStiff(:,:,:), Phase24PrismMass(:,:,:), Phase24PrismK(:), &
          Phase24PrismRhoCp(:), Phase24PrismVolume(:)
     ! Fast wedges whose scalar conductivity depends on the solution (e.g. a
     ! k(T) UDF): geometry, basis and mass are cached, the stiffness is rebuilt
     ! each assembly from nodal k interpolated at the quadrature points, i.e.
     ! the same form as the generic path.  Any number of such materials: the
     ! nodal k is kept per (material, node) slot, since a node on the interface
     ! of two temperature-dependent materials carries two conductivities.
     ! Slots of material Phase24DynKMats(m) are Phase24DynKStart(m) ..
     ! Phase24DynKStart(m)+Phase24DynKCount(m)-1.
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24PrismBasis(:,:,:), Phase24DynKNodal(:)
     LOGICAL, ALLOCATABLE, SAVE :: Phase24PrismDynamicK(:)
     INTEGER, ALLOCATABLE, SAVE :: Phase24PrismKNodeSlot(:,:), Phase24DynKNodes(:)
     INTEGER, ALLOCATABLE, SAVE :: Phase24DynKMats(:), Phase24DynKStart(:), Phase24DynKCount(:), &
          Phase24ElementDynMat(:)
     INTEGER, SAVE :: Phase24DynKNodeCount = 0, Phase24DynKMatCount = 0, &
          Phase24DynamicPrismCount = 0
     ! Optional preconditioning matrix for thin wedges: the stiffness with the
     ! in-plane/vertical mass-like couplings lumped (an M-matrix-like operator
     ! that classical AMG coarsens well).  Only the preconditioner changes;
     ! the Krylov solve keeps the consistent matrix.  Stored as a delta
     ! parallel to Matrix % Values.
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24PrecDelta(:)
     LOGICAL, SAVE :: Phase24LumpedPrec = .FALSE.
     ! Inner TES circuits.  One TES body per circuit; a single-pixel case uses
     ! the legacy keywords ('TES Body ID', 'TES Bias Current', ...), N circuits
     ! use 'TES Circuit Count = N' (Constants), 'TES Body IDs(N)' (solver) and
     ! numbered constants 'TES <k> Bias Current', ... (see TESKey).  All
     ! per-circuit state lives here; the circuits share the thermal field and
     ! are updated one after another before every assembly.
     TYPE :: TESCircuitState_t
       INTEGER :: Body = 0
       ! nonlinear (Newton / Aitken) state
       REAL(KIND=dp) :: NewtonPrevP = 0.0_dp, NewtonPrevT = 0.0_dp
       REAL(KIND=dp) :: NewtonS = -1.0_dp, NewtonSdt = -1.0_dp
       LOGICAL :: NewtonHasPrev = .FALSE.
       REAL(KIND=dp) :: StepP1 = -1.0_dp, StepP2 = -1.0_dp, StepDt1 = -1.0_dp, StepDt2 = -1.0_dp
       REAL(KIND=dp) :: NewtonPrevStepDt = -1.0_dp
       ! electrical state
       REAL(KIND=dp) :: CircuitPower = 0.0_dp, CircuitCurrent = 0.0_dp
       REAL(KIND=dp) :: CircuitResistance = 0.0_dp, CircuitAverageTemperature = 0.0_dp
       REAL(KIND=dp) :: CircuitSweepTemperature = 0.0_dp
       REAL(KIND=dp) :: PreviousCircuitCurrent = 0.0_dp
       REAL(KIND=dp) :: PreviousCircuitPower = -1.0_dp
       REAL(KIND=dp) :: CircuitPrevResidual = 0.0_dp, CircuitOmega = 0.5_dp, CircuitOmegaCap = 0.5_dp
       REAL(KIND=dp) :: LastCircuitDt = -1.0_dp
       INTEGER :: LastCircuitStep = -HUGE(1), CircuitIterInStep = 0, CircuitCallsInStep = 0
       LOGICAL :: CircuitInitialized = .FALSE., CircuitPowerInitialized = .FALSE.
       LOGICAL :: SeriesStarted = .FALSE., IterationSeriesStarted = .FALSE.
       ! snapshot at an adaptive trial boundary
       REAL(KIND=dp) :: SavedCircuitPower = 0.0_dp, SavedCircuitCurrent = 0.0_dp
       REAL(KIND=dp) :: SavedCircuitResistance = 0.0_dp, SavedCircuitAverageTemperature = 0.0_dp
       REAL(KIND=dp) :: SavedCircuitSweepTemperature = 0.0_dp, SavedPreviousCircuitCurrent = 0.0_dp
       REAL(KIND=dp) :: SavedPreviousCircuitPower = -1.0_dp
       REAL(KIND=dp) :: SavedCircuitPrevResidual = 0.0_dp, SavedCircuitOmega = 0.5_dp
       REAL(KIND=dp) :: SavedCircuitOmegaCap = 0.5_dp, SavedLastCircuitDt = -1.0_dp
       INTEGER :: SavedLastCircuitStep = -HUGE(1), SavedCircuitIterInStep = 0
       INTEGER :: SavedCircuitCallsInStep = 0, SavedTrialId = -1
       LOGICAL :: SavedCircuitInitialized = .FALSE., SavedCircuitPowerInitialized = .FALSE.
       LOGICAL :: SavedAdaptiveState = .FALSE.
       ! circuit law and the Joule power used by the latest thermal solve, for
       ! the coupled residual P_circuit(T_TES) - P_used checked after the solve
       LOGICAL :: ResValid = .FALSE.
       REAL(KIND=dp) :: ResIBias = 0.0_dp, ResRSh = 0.0_dp, ResR0 = 0.0_dp, ResRMin = 0.0_dp
       REAL(KIND=dp) :: ResAlpha = 0.0_dp, ResBeta = 0.0_dp, ResI0 = 1.0_dp, ResT0 = 1.0_dp
       REAL(KIND=dp) :: ResL = 0.0_dp, ResPrevI = 0.0_dp, ResP = 0.0_dp, ResDt = 1.0_dp
       ! TES mean of the previous nonlinear iterate (iterate-change exit tests)
       REAL(KIND=dp) :: ExitPrevT = 0.0_dp
     END TYPE TESCircuitState_t
     TYPE(TESCircuitState_t), ALLOCATABLE, SAVE :: TESC(:)
     INTEGER, SAVE :: TESNCirc = 0
     LOGICAL, SAVE :: TESNumbered = .FALSE.
     INTEGER :: TESk
     LOGICAL :: TESAllValid
     REAL(KIND=dp) :: TESdTOne
     REAL(KIND=dp) :: TESPowerTol, TESPowerRes, TESCoupledTolT, TESdPdT, TESdTCircuit, &
          TESdkRel, TESdkTmax, TESdkTmin, TESdTk
     REAL(KIND=dp), ALLOCATABLE :: TESdkNow(:)
     ! Fast elements whose body force only carries a 'Volumetric Heat Source'
     ! keep the cached operator; the source enters the RHS as
     ! int N_i sum_j N_j q_j (the generic form) from nodal values evaluated
     ! once per assembly for each such body force.
     INTEGER, ALLOCATABLE, SAVE :: Phase24ElementSrc(:), Phase24SrcBF(:), &
          Phase24SrcNodeStart(:), Phase24SrcNodeCount(:), Phase24SrcNodeList(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24SrcLoad(:,:)
     INTEGER, SAVE :: Phase24SrcCount = 0
     INTEGER, SAVE :: Phase24MetadataElements = 0, Phase24StaticCount = 0, &
          Phase24DynamicCount = 0, Phase24GenericCount = 0
     INTEGER, SAVE :: Phase24FastTetCount = 0, Phase24GenericPrismCount = 0, &
          Phase24GenericOtherCount = 0, Phase24PrismCount = 0, Phase24FastPrismCount = 0
     LOGICAL, SAVE :: Phase24MetadataReady = .FALSE.
     LOGICAL, SAVE :: Phase24EpochReady = .FALSE.
     LOGICAL, SAVE :: Phase24BoundaryMatrixInvariant = .FALSE.
     INTEGER, SAVE :: Phase24StructureEpoch = 0, Phase24MatrixEpoch = 0, &
          Phase24RHSEpoch = 0
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24StaticCSRValues(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: Phase24OperatorKReference(:), &
          Phase24OperatorMassReference(:), Phase24OperatorReactionReference(:), &
          Phase24OperatorKCurrent(:), Phase24OperatorMassCurrent(:), &
          Phase24OperatorReactionCurrent(:), Phase24LaggedCSRValues(:)
     LOGICAL, SAVE :: Phase24OperatorReferenceValid = .FALSE., &
          Phase24OperatorSnapshotValid = .FALSE.
     INTEGER, SAVE :: Phase24FullRebuildCount = 0, Phase24CaseACount = 0, &
          Phase24CaseBCount = 0, Phase24LaggedReuseCount = 0, &
          Phase24ForcedRebuildCount = 0, Phase24RecoveryCount = 0
     INTEGER, SAVE :: Phase24RecoveryCooldown = 0
     REAL(KIND=dp), SAVE :: Phase24OperatorDt = -1.0_dp
     INTEGER, SAVE :: Phase24OperatorOrder = -1
     REAL(KIND=dp), SAVE :: Phase24StaticCacheDt = -1.0_dp
     REAL(KIND=dp), SAVE :: Phase24StaticCacheA0 = -1.0_dp
     INTEGER, SAVE :: Phase24StaticCacheOrder = -1, Phase24StaticCacheValueCount = 0
     LOGICAL, SAVE :: Phase24StaticCacheLumped = .FALSE., Phase24StaticCacheValid = .FALSE.

     INTEGER, ALLOCATABLE, SAVE :: Indexes(:), SaveIndexes(:)

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
       FORCE(:), U(:), V(:), W(:), MU(:,:),TimeForce(:), &
       Density(:), LatentHeat(:), HeatTransferCoeff(:), &
       HeatCapacity(:), Enthalpy(:), EnthalpyFraction(:), Viscosity(:), LocalTemperature(:), &
       NodalVal(:), ElectricConductivity(:), Permeability(:), Work(:), C0(:), &
       Pressure(:), dPressuredt(:), GasConstant(:),AText(:), HeaterArea(:), &
       HeaterTarget(:), HeaterScaling(:), HeaterDensity(:), HeaterSource(:), &
       HeatExpansionCoeff(:), ReferenceTemperature(:), PressureCoeff(:), &
       PhaseVelocity(:,:), HeatConductivityIso(:), &
       PerfusionRate(:), PerfusionDensity(:), PerfusionHeatCapacity(:), PerfusionRefTemperature(:)

     REAL(KIND=dp), ALLOCATABLE :: Areas(:), Emiss(:), Reflect(:)
     LOGICAL :: Spectral, Radiosity
     
     SAVE U, V, W, MU, MASS, STIFF, LOAD, PressureCoeff, &
       FORCE, ElementNodes, HeatConductivity, HeatCapacity, HeatTransferCoeff, &
       Enthalpy, EnthalpyFraction, Density, LatentHeat, PhaseVelocity, AllocationsDone, Viscosity, TimeForce, &
       LocalNodes, LocalTemperature, Work, ElectricConductivity, &
       NodalVal, Permeability, C0, dPressuredt, Pressure, &
       GasConstant,AText,Hwrk, XX, YY, ForceHeater, Power, HeaterArea, HeaterTarget, &
       HeaterScaling, HeaterDensity, HeaterSource, SmartHeaters, IntegralHeaters, SmartTolReached,    &
       ReferenceTemperature, HeatExpansionCoeff, PrevPowerScaling, PowerScaling, &
       MeltPoint, DoneTime, SmartHeaterNode, SmartHeaterBC, SmartHeaterAverage, &
       HeatConductivityIso, &
       PerfusionRate, PerfusionDensity, PerfusionHeatCapacity, PerfusionRefTemperature


     INTERFACE
        SUBROUTINE HeatSolver_Boundary_Residual( Model,Edge,Mesh,Quant,Perm,Gnorm,Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END SUBROUTINE HeatSolver_Boundary_Residual

        SUBROUTINE HeatSolver_Edge_Residual( Model,Edge,Mesh,Quant,Perm,Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END SUBROUTINE HeatSolver_Edge_Residual

        SUBROUTINE HeatSolver_Inside_Residual( Model,Element,Mesh,Quant,Perm, Fnorm,Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END SUBROUTINE HeatSolver_Inside_Residual
     END INTERFACE

     REAL(KIND=dp) :: at,at0,totat,st,totst,t1
     REAL(KIND=dp) :: Phase24Tick, Phase24BulkStart, Phase24BoundaryStart, Phase24FinalizeStart, &
          Phase24SolveStart, Phase24IterationStart, Phase24UDF, &
          Phase24Bulk, Phase24Boundary, Phase24Finalize, Phase24Solve, &
          Phase24Material, Phase24Local, Phase24Insert


     CALL Info('HeatSolver','-------------------------------------------',Level=6)
     CALL Info('HeatSolver','Solving the energy equation for temperature',Level=5)

     IF( ListCheckPresentAnyBC( Model,'Heat Gap') ) THEN
       CALL Warn('HeatSolver','The old way of dealing with HeatGap is obsolite!!')
     END IF

     
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     SolverParams => GetSolverParams()
     Phase24OuterCapture = GetLogical(SolverParams, 'Phase24 Outer HeatSolve Capture', Found)
     IF (.NOT. Found) Phase24OuterCapture = .FALSE.
     Phase24OuterCapturePrefix = GetString(SolverParams, &
          'Phase24 Outer HeatSolve Capture Prefix', Found)
     IF (.NOT. Found) Phase24OuterCapturePrefix = 'phase24_outer_capture'
     Phase24OuterCaptureStep = GetInteger(SolverParams, &
          'Phase24 Outer HeatSolve Capture Timestep', Found)
     IF (.NOT. Found .OR. Phase24OuterCaptureStep < 1) Phase24OuterCaptureStep = 1
     Phase24OuterCaptureMaxIter = GetInteger(SolverParams, &
          'Phase24 Outer HeatSolve Capture Max Iterations', Found)
     IF (.NOT. Found .OR. Phase24OuterCaptureMaxIter < 1) Phase24OuterCaptureMaxIter = 3
     Phase24VectorAssembly = GetLogical(SolverParams, 'Phase24 Vector Assembly', Found)
     Phase24LumpedPrec = GetLogical(SolverParams, 'Phase24 Prism Lumped Preconditioner', Found)
     Phase24WallProfiling = GetLogical(SolverParams, 'Phase24 Wall Profiling', Found)
     Phase24StaticMatrixReuse = GetLogical(SolverParams, 'Phase24 Static Matrix Reuse', Found)
     IF (.NOT. Found) Phase24StaticMatrixReuse = .FALSE.
     Phase24LagModeName = GetString(SolverParams, 'Phase24 Operator Lagging', Found)
     IF (.NOT. Found) Phase24LagModeName = 'disabled'
     Phase24LagMode = 0
     SELECT CASE (TRIM(Phase24LagModeName))
     CASE ('conservative')
       Phase24LagMode = 1
     CASE ('adaptive')
       Phase24LagMode = 2
     CASE DEFAULT
       Phase24LagMode = 0
     END SELECT
     Phase24LaggingEnabled = Phase24VectorAssembly .AND. Phase24LagMode > 0
     Phase24LaggingActive = Phase24LaggingEnabled
     IF (Phase24LaggingEnabled .AND. Phase24RecoveryCooldown > 0) &
          Phase24LaggingActive = .FALSE.
     Phase24LagThreshold = GetCReal(SolverParams, 'Phase24 Operator Lag Threshold', Found)
     IF (.NOT. Found .OR. Phase24LagThreshold <= 0.0_dp) THEN
       IF (Phase24LagMode == 1) THEN
         Phase24LagThreshold = 1.0e-4_dp
       ELSE
         Phase24LagThreshold = 5.0e-4_dp
       END IF
     END IF
     Phase24LagMaxReuse = GetInteger(SolverParams, 'Phase24 Operator Lag Max Reuse', Found)
     IF (.NOT. Found .OR. Phase24LagMaxReuse < 1) Phase24LagMaxReuse = 3
     Phase24LagResidualGrowth = GetCReal(SolverParams, &
          'Phase24 Operator Lag Residual Growth', Found)
     IF (.NOT. Found .OR. Phase24LagResidualGrowth <= 1.0_dp) &
          Phase24LagResidualGrowth = 1.25_dp
     Phase24KrylovRelativeGrowth = GetCReal(SolverParams, &
          'Phase24 Operator Lag Krylov Relative Growth', Found)
     IF (.NOT. Found .OR. Phase24KrylovRelativeGrowth <= 0.0_dp) &
          Phase24KrylovRelativeGrowth = 0.50_dp
     Phase24KrylovAbsoluteGrowth = GetInteger(SolverParams, &
          'Phase24 Operator Lag Krylov Absolute Growth', Found)
     IF (.NOT. Found .OR. Phase24KrylovAbsoluteGrowth < 1) &
          Phase24KrylovAbsoluteGrowth = 8
     Phase24LagReuseCount = 0
     Phase24LaggedReuseThisTimestep = 0
     Phase24KrylovBaseline = 0
     Phase24KrylovIterations = 0
     Phase24ForceRebuild = .FALSE.
     Phase24FinalConfirmationDone = .FALSE.
     Phase24FallbackRebuild = .FALSE.
     Phase24LHSUnsafe = .FALSE.
     Phase24PreviousRelativeChange = -1.0_dp
     ! A rejected adaptive trial may have touched numeric HYPRE state during
     ! nonlinear assembly.  Keep the persistent IJ/ParCSR allocation alive,
     ! but invalidate all numeric ownership epochs before the retry.
     IF (ListGetLogical(SolverParams,'Phase24 Adaptive Rejected',Found) .AND. Found) THEN
       Phase24MatrixEpoch = Phase24MatrixEpoch + 1
       Phase24RHSEpoch = 0
       Phase24OperatorReferenceValid = .FALSE.
       Phase24OperatorSnapshotValid = .FALSE.
       Phase24StaticCacheValid = .FALSE.
       Phase24ForceRebuild = .TRUE.
       Phase24RecoveryCooldown = MAX(Phase24RecoveryCooldown,1)
       CALL ListAddLogical(SolverParams,'Phase24 Adaptive Rejected',.FALSE.)
       CALL Info('HeatSolve','Phase24 adaptive retry: numeric lifecycle epochs invalidated',Level=4)
     END IF
     IF (.NOT. Phase24VectorAssembly) THEN
       Phase24OperatorReferenceValid = .FALSE.
       Phase24OperatorSnapshotValid = .FALSE.
     END IF
     Phase24Eligible = Phase24VectorAssembly
     IF (Phase24Eligible) THEN
       IF (CoordinateSystemDimension() /= 3) Phase24Eligible = .FALSE.
       IF (.NOT. ASSOCIATED(Solver % Variable)) Phase24Eligible = .FALSE.
       IF (Phase24Eligible .AND. Solver % Variable % DOFs /= 1) Phase24Eligible = .FALSE.
       IF (.NOT. Phase24Eligible) THEN
         CALL Info('HeatSolve', &
              'Phase24 Vector Assembly disabled: case is outside scalar 3D mixed-element eligibility', &
              Level=4)
       END IF
     END IF
     Phase24VectorAssembly = Phase24Eligible
     IF (.NOT. Phase24VectorAssembly) THEN
       CALL Phase24ClearMetadata()
       IF (ALLOCATED(Phase24BodyMaterialKnown)) DEALLOCATE(Phase24BodyMaterialKnown)
       IF (ALLOCATED(Phase24BodyMaterialCached)) DEALLOCATE(Phase24BodyMaterialCached)
       IF (ALLOCATED(Phase24BodyHeatCapacity)) DEALLOCATE(Phase24BodyHeatCapacity)
       IF (ALLOCATED(Phase24BodyHeatConductivity)) DEALLOCATE(Phase24BodyHeatConductivity)
       Phase24BodyCacheSize = 0
     ELSE IF (Phase24BodyCacheSize /= Model % NumberOfBodies .OR. &
              .NOT. ALLOCATED(Phase24BodyMaterialKnown) .OR. &
              .NOT. AllocationsDone .OR. Solver % MeshChanged) THEN
       IF (ALLOCATED(Phase24BodyMaterialKnown)) DEALLOCATE(Phase24BodyMaterialKnown)
       IF (ALLOCATED(Phase24BodyMaterialCached)) DEALLOCATE(Phase24BodyMaterialCached)
       IF (ALLOCATED(Phase24BodyHeatCapacity)) DEALLOCATE(Phase24BodyHeatCapacity)
       IF (ALLOCATED(Phase24BodyHeatConductivity)) DEALLOCATE(Phase24BodyHeatConductivity)
       IF (Model % NumberOfBodies > 0) THEN
         ALLOCATE(Phase24BodyMaterialKnown(Model % NumberOfBodies), &
              Phase24BodyMaterialCached(Model % NumberOfBodies), &
              Phase24BodyHeatCapacity(Model % NumberOfBodies), &
              Phase24BodyHeatConductivity(Model % NumberOfBodies))
         Phase24BodyMaterialKnown = .FALSE.
         Phase24BodyMaterialCached = .FALSE.
         Phase24BodyHeatCapacity = 0.0_dp
         Phase24BodyHeatConductivity = 0.0_dp
       END IF
       Phase24BodyCacheSize = Model % NumberOfBodies
     END IF
     IF (Phase24VectorAssembly .AND. (.NOT. AllocationsDone .OR. Solver % MeshChanged)) THEN
       CALL CRS_ClearPositionCache()
       CALL Phase24ClearMetadata()
     END IF

     IF (.NOT. Phase24EpochReady .OR. .NOT. AllocationsDone .OR. Solver % MeshChanged) THEN
       Phase24StructureEpoch = Phase24StructureEpoch + 1
       Phase24MatrixEpoch = Phase24MatrixEpoch + 1
       Phase24RHSEpoch = 0
       Phase24EpochReady = .TRUE.
       Phase24BoundaryMatrixInvariant = .FALSE.
     END IF

     Radiosity = GetLogical( SolverParams,'Radiosity Model',Found ) 
     Spectral = GetLogical( SolverParams,'Spectral Model',Found )
     IF( Spectral ) Radiosity = .TRUE.

     IF(.NOT. Radiosity ) THEN
       IF( ListCheckPresent( SolverParams,'Radiator Coordinates' ) .OR. &
           ListCheckPresentAnyBodyForce(Model,'Radiator Coordinates') ) THEN
         CALL Fatal('HeatSolve','For radiative point sources use HeatSolveVec instead!')
       END IF
     END IF

!------------------------------------------------------------------------------
!    The View and Gebhart factors may change. If this is necessary, this is 
!    done within this subroutine. The routine is called in the
!    start as it may affect the matrix topology.
!    Newton lineariarization option is needed only when there is radiation.
!------------------------------------------------------------------------------
     IsRadiation = ListCheckPresentAnyBC( Model,'Radiation')
     
     IF( IsRadiation .AND. .NOT. Radiosity ) THEN
       CALL RadiationFactors( Solver, .FALSE., .FALSE.)       
     END IF

     ! The solver matrix, permutation etc. may change because radiation factors are recomputed
     StiffMatrix => GetMatrix()
     ForceVector => Solver % Matrix % RHS
     TempSol => Solver % Variable
     TempPerm    => TempSol % Perm
     Temperature => TempSol % Values
     VarName = GetVarName( TempSol ) 
  
     LocalNodes = COUNT( TempPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN
     IF(SIZE(Temperature) < LocalNodes) LocalNodes = SIZE(Temperature)

     IF (Phase24VectorAssembly .AND. .NOT. Phase24MetadataReady) THEN
       CALL Phase24BuildMetadata()
     END IF
     
     
     NeedFlowSol = .FALSE.
     DO i=1,Model % NumberOfEquations
       ConvectionFlag = GetString( Model % Equations(i) % Values, 'Convection', Found )
       IF ( ConvectionFlag == 'computed' ) THEN
         NeedFlowSol = .TRUE.
         EXIT
       END IF
     END DO

     FlowSol => NULL()
     IF( NeedFlowSol ) THEN
       ConvectionField = GetString( SolverParams, 'Temperature Convection Field', Found )     
       IF ( Found ) THEN       
         FlowSol => VariableGet( Solver % Mesh % Variables, ConvectionField )
       ELSE
         FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
       END IF
      
       IF ( ASSOCIATED( FlowSol ) ) THEN
         FlowPerm     => FlowSol % Perm
         NSDOFs       =  FlowSol % DOFs
         FlowSolution => FlowSol % Values
       ELSE
         CALL Fatal('HeatSolver','Flow is "computed" but not flow field available!')
       END IF
     END IF
       
     DensitySol => VariableGet( Solver % Mesh % Variables, 'Density' )

     ! Check whether we have some heater controls. This will affect initialization stuff. 
     SmartHeaterControl = ListCheckPresentAnyBodyForce( Model,'Smart Heater Control')
     IntegralHeaterControl = ListCheckPresentAnyBodyForce( Model,'Integral Heat Source')
   
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN
        N = Solver % Mesh % MaxElementDOFs

        IF ( AllocationsDone ) THEN
          DEALLOCATE(  &
                 U, V, W, MU,           &
                 Pressure,              &
                 dPressureDt,           &
                 PressureCoeff,        &
                 ElementNodes % x,      &
                 ElementNodes % y,      &
                 ElementNodes % z,      &
                 Density,Work,          &
                 LatentHeat,            &
                 PhaseVelocity,         &
                 ElectricConductivity,  &
                 Permeability,          &
                 Viscosity,C0,          &
                 HeatTransferCoeff,     &
                 HeatExpansionCoeff,    &
                 ReferenceTemperature,  &
                 MASS,       &
                 LocalTemperature,      &
                 HeatCapacity,Enthalpy, &
                 EnthalpyFraction,      &
                 NodalVal,       &
                 GasConstant, AText,    &
                 HeatConductivity,      &
                 STIFF,LOAD,            &
                 Indexes, SaveIndexes,  &
                 FORCE, TimeForce,      &
                 HeatConductivityIso,   &
                 PerfusionRate,         &
                 PerfusionDensity,      &
                 PerfusionHeatCapacity, &
                 PerfusionRefTemperature )
       END IF

       ALLOCATE( &
                 Indexes(N), SaveIndexes(N),           &
                 U( N ),   V( N ),  W( N ),            &
                 MU( 3,N ),                            &
                 Pressure( N ),                        &
                 dPressureDt( N ),                     &
                 PressureCoeff( N ),                   &
                 ElementNodes % x( N ),                &
                 ElementNodes % y( N ),                &
                 ElementNodes % z( N ),                &
                 Density( N ),Work( N ),               &
                 LatentHeat( N ),                      &
                 PhaseVelocity(3, N),                  &
                 ElectricConductivity(N),              &
                 Permeability(N),                      &
                 Viscosity(N),C0(N),                   &
                 HeatTransferCoeff( N ),               &
                 HeatExpansionCoeff( N ),              &
                 ReferenceTemperature( N ),            &
                 MASS(  2*N,2*N ),                     &
                 LocalTemperature( N ),                &
                 HeatCapacity( N ),Enthalpy( N ),      &
                 EnthalpyFraction( N ),                &
                 NodalVal( N ),                 &
                 GasConstant( N ),AText( N ),          &
                 HeatConductivity( 3,3,N ),            &
                 HeatConductivityIso( N ),             &
                 STIFF( 2*N,2*N ),LOAD( N ), &
                 FORCE( 2*N ), TimeForce(2*N), &
                 PerfusionRate( N ),     &
                 PerfusionDensity( N ),      &
                 PerfusionHeatCapacity( N ), &
                 PerfusionRefTemperature( N ), &
                 STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'HeatSolve', 'Memory allocation error' )
       END IF


       IF ( SmartHeaterControl .OR. IntegralHeaterControl) THEN          
          n = Model % NumberOfBodyForces
          IF ( AllocationsDone ) DEALLOCATE( HeaterArea, HeaterDensity, HeaterSource, &
               HeaterScaling, HeaterTarget, SmartHeaters, IntegralHeaters )
          ALLOCATE( HeaterArea(n), HeaterDensity(n), HeaterSource(n), &
               HeaterScaling(n), HeaterTarget(n), SmartHeaters(n), &
               IntegralHeaters(n) )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          SmartHeaters = .FALSE.
          IntegralHeaters = .FALSE.
       END IF
       
       IF( SmartHeaterControl ) THEN
          IF ( AllocationsDone ) DEALLOCATE( XX, YY, ForceHeater  )
          n = SIZE( Temperature )
          ALLOCATE( XX( n ), YY(n), ForceHeater( n ), STAT=istat )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          XX = 0.0d0 
          YY = 0.0d0
          ForceHeater = 0.0d0
       END IF
       
       NULLIFY( Hwrk )
       AllocationsDone = .TRUE.
    END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     dt = Timestep
     Constants => GetConstants()
     IF( IsRadiation ) THEN
       StefanBoltzmann = ListGetConstReal( Model % Constants, &
                     'Stefan Boltzmann',UnfoundFatal=.TRUE.)
     END IF

!------------------------------------------------------------------------------
     Stabilize = GetLogical( SolverParams,'Stabilize',Found )

     UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
     IF ( .NOT.Found ) UseBubbles = .TRUE.

     StabilizeFlag = GetString( SolverParams, &
          'Stabilization Method',Found )

     SELECT CASE(StabilizeFlag)
     CASE('vms')
       Stabilize = .FALSE.
       UseBubbles= .FALSE.
     CASE('stabilized')
       Stabilize = .TRUE.
       UseBubbles = .FALSE.
     CASE('bubbles')
       Stabilize = .FALSE.
       UseBubbles = .TRUE.
     END SELECT

     NonlinearIter = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 1
     Phase24IterationLimit = NonlinearIter
     IF (Phase24LaggingEnabled) Phase24IterationLimit = NonlinearIter + 2
     MinNonlinearIter = GetInteger( SolverParams, &
                     'Nonlinear System Min Iterations', Found )
     IF ( .NOT.Found ) MinNonlinearIter = 1
     MinNonlinearIter = MAX(1, MIN(MinNonlinearIter, NonlinearIter))

     NonlinearTol  = GetConstReal( SolverParams, &
                     'Nonlinear System Convergence Tolerance',    Found )

     IF( IsRadiation ) THEN
       NewtonTol     = GetConstReal( SolverParams, &
                      'Nonlinear System Newton After Tolerance',  Found )
       NewtonIter    = GetInteger(   SolverParams, &
                      'Nonlinear System Newton After Iterations', Found )
     ELSE
       NewtonTol = 1.0_dp
       NewtonIter =  0
     END IF
     IF ( NewtonIter == 0) NewtonLinearization = .TRUE.

     Relax = GetCReal( SolverParams,'Nonlinear System Relaxation Factor',Found )
     IF ( .NOT.Found ) Relax = 1

     TransientAssembly = TransientSimulation
     dt0 = ListGetConstReal(SolverParams,'Steady State Transition Timestep',Found)
     IF(.NOT. Found) dt0 = ListGetConstReal(SolverParams,'Smart Heater Time Scale',Found)

     IF(Found .AND. dt > dt0) TransientAssembly = .FALSE.

     
     AnyMultiply = ListCheckPresentAnyMaterial( Model, 'Heat Transfer Multiplier' ) 

!------------------------------------------------------------------------------

     TransientHeaterControl = .FALSE.
     IF(SmartHeaterControl) THEN

       ! Mark the smart heaters 
       SmartHeaters = .FALSE.
       bf_id = 0
       DO i = 1,Model % NumberOfBodyForces
         IF( ListGetLogical( Model % BodyForces(i) % Values, &
             'Smart Heater Control', Found ) ) THEN
           SmartHeaters(i) = .TRUE.	     
           bf_id = i
         END IF
       END DO

       ! Find the BC that controls the heater 
       ! If not found assume that smart heater is related to phase change 
       MeltPoint = GetCReal( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Temperature',GotMeltPoint)           
              
       SmartHeaterAverage = .FALSE.
       SmartHeaterNode = ListGetInteger( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Control Node',GotIt) 
       IF(.NOT. GotIt) THEN
         RealWork => ListGetConstRealArray( Model % BodyForces(bf_id) % Values,&
             'Smart Heater Control Point',GotIt) 
         IF( GotIt ) THEN
           ControlPoint(1:3) = RealWork(1:3,1)
           
           mindist = HUGE( mindist )
           DO l=1,Model % NumberOfNodes
             IF( TempPerm(l) == 0 ) CYCLE
             
             jx = Model % Mesh % Nodes % x(l)
             jy = Model % Mesh % Nodes % y(l)
             jz = Model % Mesh % Nodes % z(l)
             
             dist = (ControlPoint(1)-jx)**2 + (ControlPoint(2)-jy)**2 + (ControlPoint(3)-jz)**2
             IF( dist < mindist ) THEN
               mindist = dist
               SmartHeaterNode = l
             END IF
           END DO
         END IF

         WRITE(Message,*) 'Found Control Point at distance:',SQRT(mindist)
         CALL Info('HeatSolve',Message)
         WRITE(Message,*) 'Control Point Index:',SmartHeaterNode
         CALL Info('HeatSolve',Message)        
       END IF
       
       IF( .NOT. GotMeltPoint .OR. SmartHeaterNode == 0) THEN
         GotIt = .FALSE.
         Found = .FALSE.
         SmartHeaterBC = 0
         
         DO i=1,Model % NumberOfBCs
           GotIt = ListGetLogical( Model % BCs(i) % Values,'Smart Heater Boundary', Found ) 
           IF(GotIt) THEN
             SmartHeaterBC = i
             EXIT
           END IF
         END DO
         IF(.NOT. GotIt) THEN
           DO i=1,Model % NumberOfBCs
             GotIt = ListGetLogical( Model % BCs(i) % Values,'Phase Change', Found ) 
             IF(GotIt) THEN
               SmartHeaterBC = i
               EXIT
             END IF
           END DO
         END IF
         IF(SmartHeaterBC == 0) THEN
           CALL Fatal('HeatSolve','Smart Heater Boundary / Phase Change is undefined')
         END IF
         
         MeltPoint = GetCReal( Model % BCs(SmartHeaterBC) % Values,&
             'Smart Heater Temperature',Found)
         IF(.NOT. Found) THEN
           DO k=1, Model % NumberOfMaterials
             MeltPoint = GetCReal( Model % Materials(k) % Values, &
                 'Melting Point', Found )
             IF(Found) EXIT
           END DO
           IF(.NOT. Found) THEN
             CALL Fatal('HeatSolver','Smart Heater Temperature / Melting Point is undefined')
           END IF
         END IF
         
         ! Find the node related to temperature control 
         SmartHeaterAverage = ListGetLogical(Solver % Values,'Smart Heater Average', Found)
         IF(.NOT. SmartHeaterAverage) THEN
           jx = -HUGE(jx)
           DO k = Model % Mesh % NumberOfBulkElements + 1, &
               Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements
             
             Element => Model % Mesh % Elements(k)
             
             IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
               DO l=1,Element % TYPE % NumberOfNodes
                 IF ( Model % Mesh % Nodes % x(Element % NodeIndexes(l)) >= jx ) THEN
                   j = Element % NodeIndexes(l) 
                   jx = Model % Mesh % Nodes % x(Element % NodeIndexes(l))
                 END IF
               END DO
             END IF
           END DO
           SmartHeaterNode = j
         END IF
       END IF

        SmartTol  = GetConstReal( SolverParams, &
             'Smart Heater Control After Tolerance',  Found )
        IF(.NOT. Found) THEN
          SmartTolReached = .TRUE.
          SmartTol = 1.0
        END IF   
     
        PowerTimeScale = ListGetConstReal(Solver % Values, &
             'Smart Heater Time Scale',Found)

        IF(TransientSimulation .AND. dt < PowerTimeScale) THEN
           TransientHeaterControl = .TRUE.
           CALL Info( 'HeatSolve', 'Using Transient Heater Control')
        ELSE
           TransientHeaterControl = .FALSE.
           CALL Info( 'HeatSolve', 'Using Steady-state Heater Control')
        END IF
        
        IF(Solver % DoneTime /= DoneTime) THEN
           PrevPowerScaling = PowerScaling
           DoneTime = Solver % DoneTime
        END IF
     END IF

     IF( IntegralHeaterControl) THEN
        CALL Info( 'HeatSolve', 'Using Integral Heater Control')       
        IntegralHeaters = .FALSE.
        DO i = 1,Model % NumberOfBodyForces
           IntegralHeaters(i) = ListCheckPresent( Model % BodyForces(i) % Values, &
                'Integral Heat Source')
        END DO
     END IF

!------------------------------------------------------------------------------

     ConstantBulk = GetLogical( SolverParams, 'Constant Bulk System', Found )
     SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Save Bulk System', Found )
     SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Calculate Loads', Found )

!------------------------------------------------------------------------------

     SaveRelax = Relax
     CumulativeTime = 0.0d0
     HeaterControlLocal = .FALSE.

     IF(isRadiation) THEN
       nb = Solver % Mesh % NumberOfBoundaryElements
       ALLOCATE(Areas(nb), Emiss(nb), Reflect(nb) )
       Areas=0; Emiss=0; Reflect=0
     END IF
!------------------------------------------------------------------------------
     FirstTime = .TRUE.
     Phase24UDF = 0.0_dp
     Phase24Bulk = 0.0_dp
     Phase24Boundary = 0.0_dp
     Phase24Finalize = 0.0_dp
     Phase24Solve = 0.0_dp
     Phase24Material = 0.0_dp
     Phase24Local = 0.0_dp
     Phase24Insert = 0.0_dp

     ALLOCATE(PrevSolution(LocalNodes))
     
     DO WHILE( CumulativeTime < Timestep-1.0d-12 .OR. .NOT. TransientSimulation )
!------------------------------------------------------------------------------
!    The first time around this has been done by the caller...
!------------------------------------------------------------------------------
     IF ( TransientSimulation .AND. .NOT.FirstTime ) &
       CALL InitializeTimestep(Solver)
     FirstTime = .FALSE.
!------------------------------------------------------------------------------
!    Save current solution
!------------------------------------------------------------------------------
     PrevSolution = Temperature(1:LocalNodes)
     IF ( TransientSimulation ) THEN
       PrevTemperature => Solver % Variable % PrevValues(:,1)
     END IF
     ! Apply event-boundary order resets before choosing the nonlinear initial
     ! field and circuit state for this timestep.
     IF (TransientSimulation .AND. Solver % Order > 1) THEN
       BDFResetTimes => ListGetConstRealArray(SolverParams, &
            'Phase24 BDF1 Step Start Times', Found)
       IF (Found .AND. ASSOCIATED(BDFResetTimes)) THEN
         BDFStepStart = GetTime() - Solver % dt
         DO kReset = 1, SIZE(BDFResetTimes,1)
           IF (ABS(BDFStepStart - BDFResetTimes(kReset,1)) <= 1.0e-3_dp*Solver % dt) THEN
             IF (Solver % DoneTime /= 1) CALL Info('HeatSolve', &
                  'Phase24 BDF order reset to 1 for this step',Level=4)
             Solver % DoneTime = 1
           END IF
         END DO
       END IF
     END IF
     Phase24BDF2Predictor = ListGetLogical(Model % Simulation, &
          'Phase24 BDF2 Predictor',Found)
     Phase24BDF2PredictorActive = .FALSE.
     IF (Phase24BDF2Predictor .AND. TransientSimulation .AND. &
          MIN(Solver % DoneTime,Solver % Order) == 2 .AND. Solver % dt > 0.0_dp) THEN
       Phase24PredictorMaxDelta = ListGetConstReal(Model % Simulation, &
            'Phase24 BDF2 Predictor Max Delta',Found)
       IF (.NOT. Found) Phase24PredictorMaxDelta = 1.0e-5_dp
       Phase24PredictorHistoryReady = .FALSE.
       Phase24PredictorDelta = HUGE(1.0_dp)
       IF (ASSOCIATED(Solver % Variable % PrevValues)) THEN
         IF (SIZE(Solver % Variable % PrevValues,2) >= 2) THEN
           CALL Phase24GetBDFCoefficients(2,Solver % dt,Phase24CurrentHOld, &
                Phase24CurrentA0,Phase24CurrentA1,Phase24CurrentA2)
           Phase24CurrentRatio = Solver % dt / Phase24CurrentHOld
           Phase24PredictorDelta = Phase24CurrentRatio * MAXVAL(ABS( &
                Solver % Variable % PrevValues(1:LocalNodes,1) - &
                Solver % Variable % PrevValues(1:LocalNodes,2)))
           Phase24PredictorHistoryReady = .TRUE.
         END IF
       END IF
       ! A large extrapolated jump can move the electrothermal Picard solve
       ! to a different branch before its usual two-iteration exit.  Require
       ! every rank to stay within the trust region, particularly at pulses
       ! and when the variable timestep grows sharply.
       Phase24BDF2PredictorActive = ParallelReduction(MERGE(0.0_dp,1.0_dp, &
            Phase24PredictorHistoryReady .AND. &
            Phase24PredictorDelta <= Phase24PredictorMaxDelta)) < 0.5_dp
       IF (Phase24BDF2PredictorActive) THEN
           Temperature(1:LocalNodes) = Solver % Variable % PrevValues(1:LocalNodes,1) + &
                Phase24CurrentRatio * (Solver % Variable % PrevValues(1:LocalNodes,1) - &
                Solver % Variable % PrevValues(1:LocalNodes,2))
           PrevSolution = Temperature(1:LocalNodes)
       END IF
     END IF
     CALL ListAddLogical(SolverParams,'Phase24 BDF2 Predictor Active', &
          Phase24BDF2PredictorActive)
!------------------------------------------------------------------------------
     
     totat = 0.0d0
     totst = 0.0d0

     Norm = Solver % Variable % Norm

     CALL DefaultStart()
     
     
     DO iter=1,Phase24IterationLimit
       at  = CPUTime()
       at0 = RealTime()
       Phase24IterationStart = at0
       Phase24UDF = 0.0_dp
       Phase24Bulk = 0.0_dp
       Phase24Boundary = 0.0_dp
       Phase24Finalize = 0.0_dp
       Phase24Solve = 0.0_dp
       Phase24Material = 0.0_dp
       Phase24Local = 0.0_dp
       Phase24Insert = 0.0_dp
       Phase24RHSOnly = .FALSE.
       Phase24UsedLaggedOperator = .FALSE.
       Phase24RecoveryThisSolve = .FALSE.
       Phase24MatrixChanged = .TRUE.
       IF (Phase24LaggingEnabled) CALL ListAddLogical(SolverParams, &
            'Phase24 HYPRE Recovery Active',.FALSE.)
       IF (Phase24VectorAssembly .AND. Phase24MetadataReady) THEN
         IF (ALLOCATED(Phase24OperatorKCurrent)) THEN
           IF (Phase24OperatorReferenceValid) THEN
             Phase24OperatorKCurrent = Phase24OperatorKReference
             Phase24OperatorMassCurrent = Phase24OperatorMassReference
             Phase24OperatorReactionCurrent = Phase24OperatorReactionReference
           ELSE
             Phase24OperatorKCurrent = 0.0_dp
             Phase24OperatorMassCurrent = 0.0_dp
             Phase24OperatorReactionCurrent = 0.0_dp
           END IF
         END IF
         IF (Phase24LaggingEnabled .AND. Phase24OperatorReferenceValid) THEN
           IF (ABS(Phase24OperatorDt - Solver % dt) > &
                1.0e-14_dp * MAX(1.0_dp,ABS(Solver % dt)) .OR. &
               Phase24OperatorOrder /= MIN(Solver % DoneTime, Solver % Order)) THEN
             Phase24ForceRebuild = .TRUE.
           END IF
         END IF
       END IF

       ! Optional MPI-safe lumped TES circuit update. This point is reached
       ! collectively by all ranks once per HeatSolve nonlinear iteration.
       IF (Phase24WallProfiling) CALL ResetTimer('Phase24 UDF Circuit')
       Phase24Tick = RealTime()
       CALL TESInnerCircuitUpdate()
       IF (Phase24WallProfiling) Phase24UDF = Phase24UDF + RealTime() - Phase24Tick
       IF (Phase24WallProfiling) CALL CheckTimer('Phase24 UDF Circuit',Level=4,Delete=.TRUE.)

       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       WRITE( Message,* ) 'TEMPERATURE ITERATION', iter
       CALL Info( 'HeatSolve', Message, Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', 'Starting Assembly...', Level=4 )

500    IF ( ConstantBulk .AND. ASSOCIATED(Solver % Matrix % BulkValues) ) THEN
         Solver % Matrix % Values = Solver % Matrix % BulkValues
         Solver % Matrix % RHS = Solver % Matrix % BulkRHS
         GOTO 1000
       END IF
            
       IF(Radiosity) THEN
         CALL RadiationFactors( Solver, .FALSE., NewtonLinearization)         
         CALL TabulateBoundaryAverages(Solver % Mesh, Emiss, Reflect)         
       ELSE
         CALL TabulateBoundaryAverages(Solver % Mesh, Emiss) 
       END IF
         
!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
 
       IF ( SmartHeaterControl .OR. IntegralHeaterControl ) THEN
          IF( SmartHeaterControl) ForceHeater = 0.0d0
          HeaterArea = 0.0d0
          HeaterSource = 0.0d0
          HeaterScaling = 1.0d0
          HeaterDensity = 0.0d0
          HeaterTarget = 0.0d0
          HeaterControlLocal = .FALSE.

          DO t=1,Solver % NumberOfActiveElements             
             Element => GetActiveElement(t)             
             BodyForce => GetBodyForce()
             
             IF ( .NOT. ASSOCIATED( BodyForce ) ) CYCLE
             bf_id = GetBodyForceId()
             
             IF( .NOT. (SmartHeaters(bf_id) .OR. IntegralHeaters(bf_id) ) ) CYCLE

             n = GetElementNOFNodes()

             Material => GetMaterial()

             Density(1:n) = GetReal( Material, 'Density' )
             Load(1:n) = GetReal( BodyForce, 'Heat Source' )

             s = ElementArea( Solver % Mesh, Element, n )

             IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
                  CurrentCoordinateSystem() == CylindricSymmetric ) s = 2 * PI * s

             HeaterSource(bf_id) = HeaterSource(bf_id) + s * SUM(Density(1:n) * Load(1:n)) / n
             HeaterArea(bf_id) = HeaterArea(bf_id) + s
             HeaterDensity(bf_id) = HeaterDensity(bf_id) + s * SUM( Density(1:n) ) / n
          END DO

          DO i = 1,Model % NumberOfBodyForces
             IF( IntegralHeaters(i) .OR. SmartHeaters(i) ) THEN
                HeaterDensity(i) = HeaterDensity(i) / HeaterArea(i)
             END IF
             IF(IntegralHeaters(i)) THEN
                HeaterTarget(i) = GetCReal(  Model % BodyForces(i) % Values, &
                     'Integral Heat Source', Found )
                HeaterScaling(i) = HeaterTarget(i) / HeaterSource(i)
             END IF
          END DO
       END IF

!------------------------------------------------------------------------------
       body_id = -1
       NULLIFY(Material)
!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       CALL StartAdvanceOutput( 'HeatSolve', 'Assembly:' )
       NofActive = GetNOFActive()

     IF (Phase24WallProfiling) CALL ResetTimer('Phase24 Bulk Assembly')
     Phase24BulkStart = RealTime()
     Phase24StaticAssemblyActive = .FALSE.
     Phase24PrismAssemblyActive = .FALSE.
     Phase24StaticCacheHit = .FALSE.
     Phase24CurrentOrder = MIN(Solver % DoneTime, Solver % Order)
     CALL Phase24GetBDFCoefficients(Phase24CurrentOrder,Solver % dt, &
          Phase24CurrentHOld,Phase24CurrentA0,Phase24CurrentA1,Phase24CurrentA2)
     Phase24CurrentRatio = Solver % dt / MAX(Phase24CurrentHOld,1.0e-30_dp)
     Phase24AdaptiveDebug = ListGetLogical(Model % Simulation,'Adaptive Debug',Found)
     IF (Phase24AdaptiveDebug .AND. TransientAssembly .AND. &
          Phase24CurrentOrder >= 1) THEN
       WRITE(Message,'(A,A,A,ES16.8,A,ES16.8,A,ES16.8,A,I0,A,3ES16.8)') &
            'Phase24 adaptive coefficients: kind=', &
            TRIM(ListGetString(SolverParams,'Phase24 Adaptive Trial Kind',GotIt)), &
            ', h_new=',Solver % dt,', h_old=',Phase24CurrentHOld,', r=', &
            Phase24CurrentRatio,', order=',Phase24CurrentOrder,', a0/a1/a2=', &
            Phase24CurrentA0,Phase24CurrentA1,Phase24CurrentA2
       CALL Info('HeatSolve',TRIM(Message),Level=4)
     END IF
     Phase24BoundaryMatrixInvariantNext = .TRUE.
     IF (Phase24VectorAssembly .AND. Phase24MetadataReady) THEN
       ! The static matrix is invariant over nonlinear iterations for a fixed
       ! mesh, timestep/order, and mass formulation.  Keep the transient RHS
       ! dynamic: on a cache hit the two routines below visit the compact
       ! static lists only for RHS history terms.
       IF (Phase24StaticMatrixReuse .AND. Phase24StaticCacheValid .AND. &
            ALLOCATED(Phase24StaticCSRValues)) THEN
         ! Keep these checks nested: Fortran does not guarantee short-circuit
         ! evaluation of .AND., so SIZE on an unallocated cache is unsafe.
         IF (Phase24StaticCacheValueCount == SIZE(Solver % Matrix % Values)) THEN
           IF (Phase24StaticCacheOrder == Phase24CurrentOrder) THEN
             IF (Phase24StaticCacheLumped .EQV. Solver % Matrix % Lumped) THEN
               IF (ABS(Phase24StaticCacheDt - Solver % dt) <= &
                    1.0e-14_dp * MAX(1.0_dp,ABS(Solver % dt))) THEN
                 IF (ABS(Phase24StaticCacheA0 - Phase24CurrentA0) <= &
                      1.0e-14_dp * MAX(1.0_dp,ABS(Phase24CurrentA0))) THEN
                   Solver % Matrix % Values = Phase24StaticCSRValues
                   Phase24StaticCacheHit = .TRUE.
                 END IF
               END IF
             END IF
           END IF
       END IF
      END IF
      IF (Phase24AdaptiveDebug .AND. TransientAssembly .AND. &
           Phase24CurrentOrder >= 1 .AND. Phase24StaticMatrixReuse) THEN
        WRITE(Message,'(A,L1,A,L1,A,ES16.8,A,ES16.8,A,I0)') &
             'Phase24 static cache decision: hit=',Phase24StaticCacheHit, &
             ', valid=',Phase24StaticCacheValid,', cached_a0=', &
             Phase24StaticCacheA0,', current_a0=',Phase24CurrentA0, &
             ', epoch=',Phase24MatrixEpoch
        CALL Info('HeatSolve',TRIM(Message),Level=4)
      END IF
      IF (Phase24SrcCount > 0) CALL Phase24EvaluateSources()
      IF (Phase24StaticCacheHit) THEN
         CALL Phase24AssembleStaticElements(Phase24StaticAssemblyActive,RHSOnly=.TRUE.)
         CALL Phase24AssemblePrismElements(Phase24PrismAssemblyActive,RHSOnly=.TRUE.)
         CALL Info('HeatSolve','Phase24 static matrix cache hit; RHS history refreshed',Level=4)
       ELSE
         CALL Phase24AssembleStaticElements(Phase24StaticAssemblyActive)
         CALL Phase24AssemblePrismElements(Phase24PrismAssemblyActive)
         IF (Phase24StaticMatrixReuse .AND. &
              (Phase24StaticAssemblyActive .OR. Phase24PrismAssemblyActive)) THEN
           IF (.NOT. ALLOCATED(Phase24StaticCSRValues) .OR. &
                SIZE(Phase24StaticCSRValues) /= SIZE(Solver % Matrix % Values)) THEN
             IF (ALLOCATED(Phase24StaticCSRValues)) DEALLOCATE(Phase24StaticCSRValues)
             ALLOCATE(Phase24StaticCSRValues(SIZE(Solver % Matrix % Values)))
           END IF
           Phase24StaticCSRValues = Solver % Matrix % Values
           Phase24StaticCacheDt = Solver % dt
           Phase24StaticCacheA0 = Phase24CurrentA0
           Phase24StaticCacheOrder = Phase24CurrentOrder
           Phase24StaticCacheValueCount = SIZE(Solver % Matrix % Values)
           Phase24StaticCacheLumped = Solver % Matrix % Lumped
           Phase24StaticCacheValid = .TRUE.
           CALL Info('HeatSolve','Phase24 static matrix cache built',Level=4)
         END IF
       END IF
       ! Exact RHS-only assembly is enabled only when every active bulk
       ! element is represented by the fixed geometry/material cache.  Any
       ! dynamic or generic element is a possible operator dependency and
       ! therefore keeps the conservative matrix-dirty path.
       IF (Phase24StaticCacheHit .AND. Phase24DynamicCount == 0 .AND. &
            Phase24DynamicPrismCount == 0 .AND. &
            Phase24GenericCount == 0 .AND. Phase24BoundaryMatrixInvariant) THEN
         Phase24RHSOnly = .TRUE.
         Phase24MatrixChanged = .FALSE.
       END IF
     END IF
     ! Dynamic eligible tetrahedra use the cached CRS destinations below.  If
     ! static assembly did not prepare the cache (for example when there are
     ! no static elements), prepare it once before the mixed traversal.
     IF (Phase24VectorAssembly .AND. Phase24MetadataReady .AND. &
          .NOT. Phase24StaticAssemblyActive .AND. Phase24DynamicCount > 0 .AND. &
          Solver % Matrix % Format == MATRIX_CRS .AND. &
          Solver % LocalSystemMode == 0 .AND. .NOT. Solver % PeriodicFlipActive) THEN
       CALL CRS_PreparePositionCache(Solver % Matrix, Phase24MetadataElements)
     END IF
     DO t=1,NofActive
         
         CALL AdvanceOutput(t,NofActive)
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where temperature 
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)

         ! Skip elements already assembled by the Phase24 fast paths before
         ! any per-element material/equation lookup (the original skip sat
         ! after that work, costing seconds per sweep on large meshes).
         IF (Phase24StaticAssemblyActive .OR. Phase24PrismAssemblyActive) THEN
           IF (Element % ElementIndex > 0 .AND. &
                Element % ElementIndex <= Phase24MetadataElements) THEN
             IF (Phase24StaticAssemblyActive .AND. &
                  Phase24ElementClass(Element % ElementIndex) == 1) CYCLE
             IF (Phase24PrismAssemblyActive .AND. &
                  Phase24ElementClass(Element % ElementIndex) == 4) CYCLE
           END IF
         END IF

!------------------------------------------------------------------------------
         IF ( Element % BodyId /= body_id ) THEN
!------------------------------------------------------------------------------
           Equation => GetEquation()
           ConvectionFlag = GetString( Equation, 'Convection', Found )

           Material => GetMaterial()
!------------------------------------------------------------------------------
           CompressibilityFlag = GetString( Material, &
                 'Compressibility Model', Found)
           IF ( .NOT.Found ) CompressibilityModel = Incompressible

           SELECT CASE( CompressibilityFlag )

             CASE( 'incompressible' )
               CompressibilityModel = Incompressible

             CASE( 'user defined' )
               CompressibilityModel = UserDefined1

             CASE( 'perfect gas', 'perfect gas equation 1' )
               CompressibilityModel = PerfectGas1

             CASE( 'thermal' )
               CompressibilityModel = Thermal

             CASE DEFAULT
               CompressibilityModel = Incompressible
           END SELECT
!------------------------------------------------------------------------------

           PhaseModel = GetString( Equation, 'Phase Change Model',Found )
           IF(.NOT. Found) PhaseModel = GetString( Material, 'Phase Change Model',Found )

           PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')
           IF ( PhaseChange ) THEN
              CheckLatentHeatRelease = GetLogical( Equation, &
                   'Check Latent Heat Release',Found )
           END IF
         END IF
!------------------------------------------------------------------------------

         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )
         Phase24ElementId = Element % ElementIndex
         Phase24FastElement = .FALSE.

         ! Static, source-free P1 tetrahedra were assembled from the
         ! persistent geometry/material cache above.  Leave all other
         ! elements in this original loop so dynamic properties, sources,
         ! boundaries, and unsupported element behavior remain unchanged.
         ! Do not combine these tests into one .AND. expression.  Fortran
         ! does not require short-circuit evaluation, so the fast-path array
         ! must only be indexed after its activation and bounds are known.
         IF (Phase24StaticAssemblyActive) THEN
           IF (Phase24ElementId > 0 .AND. Phase24ElementId <= Phase24MetadataElements) THEN
             IF (Phase24ElementClass(Phase24ElementId) == 1) CYCLE
           END IF
         END IF
         IF (Phase24PrismAssemblyActive) THEN
           IF (Phase24ElementId > 0 .AND. Phase24ElementId <= Phase24MetadataElements) THEN
             IF (Phase24ElementClass(Phase24ElementId) == 4) CYCLE
           END IF
         END IF

         ! Class 2 is the dynamic-but-eligible P1 tetra bucket.  Class 3
         ! includes prisms and every other unsupported element; those must use
         ! the original generic insertion path, even when the opt-in is on.
         IF (Phase24VectorAssembly .AND. Phase24MetadataReady .AND. &
              Phase24ElementId > 0 .AND. Phase24ElementId <= Phase24MetadataElements) THEN
           IF (Phase24ElementClass(Phase24ElementId) == 2 .AND. &
                Solver % Matrix % Format == MATRIX_CRS .AND. &
                Solver % LocalSystemMode == 0 .AND. .NOT. Solver % PeriodicFlipActive) THEN
             Phase24FastElement = .TRUE.
           END IF
         END IF

         IF (Phase24WallProfiling) Phase24Tick = RealTime()
         CALL GetScalarLocalSolution( LocalTemperature )
!------------------------------------------------------------------------------
!        Get element material parameters
!------------------------------------------------------------------------------
         Phase24UseCachedMaterial = .FALSE.
         IF (Phase24VectorAssembly) THEN
           IF (body_id > 0 .AND. body_id <= Phase24BodyCacheSize) THEN
             IF (.NOT. Phase24BodyMaterialKnown(body_id)) THEN
               Phase24BodyHeatCapacity(body_id) = GetCReal(Material, 'Heat Capacity', FoundCp)
               Phase24BodyHeatConductivity(body_id) = GetCReal(Material, 'Heat Conductivity', FoundK)
               Phase24BodyMaterialCached(body_id) = FoundCp .AND. FoundK
               Phase24BodyMaterialKnown(body_id) = .TRUE.
             END IF
             Phase24UseCachedMaterial = Phase24BodyMaterialCached(body_id)
           END IF
         END IF

         IF (Phase24UseCachedMaterial) THEN
           ! Only constant scalar properties enter this cache.  MATC/UDF and
           ! anisotropic properties stay on the generic per-element path.
           HeatCapacity(1:n) = Phase24BodyHeatCapacity(body_id)
           HeatConductivity = 0.0d0
           DO i=1,3
             HeatConductivity(i,i,1:n) = Phase24BodyHeatConductivity(body_id)
           END DO
         ELSE
           HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity', Found )

           CALL ListGetRealArray( Material,'Heat Conductivity',Hwrk,n, &
                        Element % NodeIndexes )
           HeatConductivity = 0.0d0
           IF ( SIZE(Hwrk,1) == 1 ) THEN
             DO i=1,3
               HeatConductivity( i,i,1:n ) = Hwrk( 1,1,1:n )
             END DO
           ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
             DO i=1,MIN(3,SIZE(Hwrk,1))
               HeatConductivity(i,i,1:n) = Hwrk(i,1,1:n)
             END DO
           ELSE
             DO i=1,MIN(3,SIZE(Hwrk,1))
               DO j=1,MIN(3,SIZE(Hwrk,2))
                 HeatConductivity( i,j,1:n ) = Hwrk(i,j,1:n)
               END DO
             END DO
           END IF
         END IF
!------------------------------------------------------------------------------

         IF ( CompressibilityModel == PerfectGas1 ) THEN

           ! Read Specific Heat Ratio:
           !--------------------------
           SpecificHeatRatio = GetConstReal( Material, &
               'Specific Heat Ratio', Found )
           IF ( .NOT.Found ) SpecificHeatRatio = 5.d0/3.d0

           ! For an ideal gas, \gamma, c_p and R are really a constant
           ! GasConstant is an array only since HeatCapacity formally is:
           !-------------------------------------------------------------
           GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) * &
               HeatCapacity(1:n) / SpecificHeatRatio

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 1.0_dp
         ELSE IF ( CompressibilityModel == Thermal ) THEN
           ReferenceTemperature(1:n) = GetReal( Material, 'Reference Temperature' )
           HeatExpansionCoeff(1:n) = GetReal( Material, 'Heat Expansion Coefficient' )

           Density(1:n) = GetReal( Material,'Density' )
           Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                ( LocalTemperature(1:n) - ReferenceTemperature(1:n) ) )

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) &
             PressureCoeff(1:n) = LocalTemperature(1:n) * HeatExpansionCoeff(1:n) / &
                   ( 1-HeatExpansionCoeff(1:n)*( &
                               LocalTemperature(1:n)-ReferenceTemperature(1:n)) )
         ELSE IF ( CompressibilityModel == UserDefined1 ) THEN
           IF ( ASSOCIATED( DensitySol ) ) THEN
             CALL GetScalarLocalSolution( Density, 'Density' ) 
           ELSE
             Density(1:n) = GetReal( Material,'Density' )
           END IF
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
         ELSE
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
           Density(1:n) = GetReal( Material, 'Density' )
         END IF

!------------------------------------------------------------------------------
! Take pressure deviation p_d as the dependent variable, p = p_0 + p_d
! for PerfectGas, read p_0
!------------------------------------------------------------------------------
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', Found)
           IF ( .NOT.Found ) ReferencePressure = 0.0d0
         END IF
!------------------------------------------------------------------------------

         HeaterControlLocal = .FALSE.
         Load = 0.0D0
         Pressure = 0.0d0
         dPressuredt = 0.0d0
!------------------------------------------------------------------------------
!        Check for convection model
!------------------------------------------------------------------------------
         C1 = 1.0D0
         U = 0._dp
         V = 0._dp
         W = 0._dp

         MU = 0.0d0
         CALL GetVectorLocalSolution( MU, 'Mesh Velocity' )

         IF ( ConvectionFlag == 'constant' ) THEN

           U(1:n) = GetReal( Material, 'Convection Velocity 1', Found )
           IF ( .NOT. Found ) &
              U(1:n) = GetReal( Equation, 'Convection Velocity 1', Found )
           V(1:n) = GetReal( Material, 'Convection Velocity 2', Found )
           IF ( .NOT. Found ) &
             V(1:n) = GetReal( Equation, 'Convection Velocity 2', Found )
           W(1:n) = GetReal( Material, 'Convection Velocity 3', Found )
           IF ( .NOT. Found ) &
             W(1:n) = GetReal( Equation, 'Convection Velocity 3', Found )

         ELSE IF ( ConvectionFlag == 'computed' ) THEN
           DO i=1,n
             k = FlowPerm(Element % NodeIndexes(i))
             IF ( k > 0 ) THEN
!------------------------------------------------------------------------------
               Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
               SELECT CASE( CompressibilityModel )
                 CASE( PerfectGas1 )
                   Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalTemperature(i) )
               END SELECT
               IF ( TransientSimulation ) THEN
                 dPressureDt(i) = ( FlowSolution(NSDOFs*k) - &
                     FlowSol % PrevValues(NSDOFs*k,1) ) / dt
               END IF
!------------------------------------------------------------------------------

               SELECT CASE( NSDOFs )
               CASE(3)
                 U(i) = FlowSolution( NSDOFs*k-2 )
                 V(i) = FlowSolution( NSDOFs*k-1 )
                 W(i) = 0.0D0

               CASE(4)
                 U(i) = FlowSolution( NSDOFs*k-3 )
                 V(i) = FlowSolution( NSDOFs*k-2 )
                 W(i) = FlowSolution( NSDOFs*k-1 )
               END SELECT
             ELSE
               U(i) = 0.0d0
               V(i) = 0.0d0
               W(i) = 0.0d0
             END IF
           END DO
         ELSE
           IF ( ALL(MU==0) ) C1 = 0.0D0 
         END IF

         IF (Phase24WallProfiling) Phase24Material = Phase24Material + RealTime() - Phase24Tick
         HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
 
!------------------------------------------------------------------------------
!        Check if modelling Phase Change with Eulerian approach 
!------------------------------------------------------------------------------
         PhaseSpatial = .FALSE.
         IF (  PhaseChange ) THEN
           CALL EffectiveHeatCapacity()
         END IF

         Viscosity = 0.0d0
!------------------------------------------------------------------------------
!        Add body forces, if any
!------------------------------------------------------------------------------
         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
!------------------------------------------------------------------------------
!          Frictional viscous heating
!------------------------------------------------------------------------------
           IF ( GetLogical( BodyForce, 'Friction Heat',Found) ) THEN
              Viscosity(1:n) = GetReal( Material,'Viscosity' )
           END IF
!------------------------------------------------------------------------------
!          Given heat source
!------------------------------------------------------------------------------
           Load(1:n) = GetReal( BodyForce, 'Volumetric Heat Source', Found )
           IF(.NOT. Found ) THEN
             Load(1:n) = Density(1:n) *  GetReal( BodyForce, 'Heat Source', Found )
           END IF
             
           IF ( SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
              IF(  SmartHeaters(bf_id) ) THEN
               HeaterControlLocal = .TRUE.
               IF( TransientHeaterControl ) THEN
                 Load(1:n) = PrevPowerScaling * Load(1:n)
                 HeaterScaling(bf_id) = PrevPowerScaling
               END IF
             END IF
           END IF

           IF ( IntegralHeaterControl ) THEN
              IF( IntegralHeaters(bf_id) ) THEN
                 Load(1:n) = Load(1:n) * HeaterScaling(bf_id) 
              END IF
           END IF
         END IF
	
         C0 = 0.0_dp
!------------------------------------------------------------------------------
! Note at this point HeatCapacity = \rho * c_p OR \rho * (c_p - R)
! and C1 = 0 (diffusion) or 1 (convection)
!------------------------------------------------------------------------------
	 
!------------------------------------------------------------------------------
!          Perfusion (added as suggested by Matthias Zenker)
!------------------------------------------------------------------------------
         IF( ASSOCIATED(BodyForce) ) THEN
           PerfusionRate(1:n) = GetReal( BodyForce, 'Perfusion Rate', Found )
         
           IF ( Found ) THEN
             PerfusionRefTemperature(1:n) = GetReal( BodyForce, 'Perfusion Reference Temperature' )
             PerfusionDensity(1:n) = GetReal( BodyForce, 'Perfusion Density' )
             PerfusionHeatCapacity(1:n) = GetReal( BodyForce, 'Perfusion Heat Capacity' )
             C0(1:n) = PerfusionHeatCapacity(1:n) * PerfusionRate(1:n) * PerfusionDensity(1:n) 
             Load(1:n) = Load(1:n) + C0(1:n) * PerfusionRefTemperature(1:n)           
           END IF
         END IF

         ! Phase24 tracks only coefficients that can enter the operator.  Heat
         ! sources and circuit terms remain RHS-only; convection, perfusion,
         ! and geometry multipliers are conservative operator dependencies.
         IF (Phase24VectorAssembly .AND. Phase24MetadataReady) THEN
           IF (ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' .OR. &
                AnyMultiply) Phase24LHSUnsafe = .TRUE.
           IF (ASSOCIATED(BodyForce)) THEN
             IF (MAXVAL(ABS(C0(1:n))) > 0.0_dp) Phase24LHSUnsafe = .TRUE.
           END IF
           IF (Phase24ElementId > 0 .AND. Phase24ElementId <= Phase24MetadataElements) THEN
             IF (ALLOCATED(Phase24OperatorKCurrent)) THEN
               Phase24OperatorKCurrent(Phase24ElementId) = &
                    MAXVAL(ABS(HeatConductivity(:,:,1:n)))
               Phase24OperatorMassCurrent(Phase24ElementId) = &
                    MAXVAL(ABS(HeatCapacity(1:n)))
               Phase24OperatorReactionCurrent(Phase24ElementId) = &
                    MAXVAL(ABS(C0(1:n)))
             END IF
           END IF
         END IF

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         Phase24StaticLocal = .FALSE.
         IF (Phase24VectorAssembly .AND. Phase24MetadataReady) THEN
           IF (Phase24ElementId > 0 .AND. Phase24ElementId <= Phase24MetadataElements) THEN
             IF (Phase24ElementClass(Phase24ElementId) == 1 .AND. &
                  CurrentCoordinateSystem() == Cartesian .AND. &
                  CompressibilityModel == Incompressible .AND. .NOT. PhaseChange .AND. &
                  .NOT. AnyMultiply .AND. &
                  .NOT. ASSOCIATED(BodyForce) .AND. ConvectionFlag /= 'computed' .AND. &
                  ConvectionFlag /= 'constant') THEN
           ! Strict eligibility is intentional: no source, convection,
           ! stabilization, phase change, anisotropy, or dynamic density may
           ! enter this shortcut.  Everything else remains generic.
               MASS = 0.0_dp
               STIFF = 0.0_dp
               FORCE = 0.0_dp
               MASS(1:4,1:4) = Phase24SpatialMass(:,:,Phase24ElementId)
               STIFF(1:4,1:4) = Phase24StaticStiff(:,:,Phase24ElementId)
               Phase24StaticLocal = .TRUE.
             END IF
           END IF
         END IF
         IF (Phase24WallProfiling) Phase24Tick = RealTime()
         IF ( .NOT. Phase24StaticLocal .AND. CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, UseBubbles, Element, n, ElementNodes )

!------------------------------------------------------------------------------
         ELSE IF ( .NOT. Phase24StaticLocal ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, Element, n, ElementNodes )
!------------------------------------------------------------------------------
         END IF
         IF (Phase24WallProfiling) Phase24Local = Phase24Local + RealTime() - Phase24Tick
!------------------------------------------------------------------------------

         ! The heat equation may have lower dimensional elements active also.
         ! For example, heat transfer through a pipe could be expressed by 1d elements.
         ! Then the multiplier should be the area of the pipe when included in 3D mesh.
         IF( AnyMultiply ) THEN
           HeatTransferMultiplier = GetCReal( Material, 'Heat Transfer Multiplier', Found )
           IF( Found ) THEN
             MASS = HeatTransferMultiplier * MASS
             STIFF = HeatTransferMultiplier * STIFF
             FORCE = HeatTransferMultiplier * FORCE
           END IF
         END IF


         IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN

           IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
             CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF

           IF (Phase24WallProfiling) Phase24Tick = RealTime()
           CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
               ForceHeater, FORCE, n, 1, TempPerm(Element % NodeIndexes) )
           IF (Phase24WallProfiling) Phase24Insert = Phase24Insert + RealTime() - Phase24Tick
         ELSE
            Bubbles = UseBubbles .AND. .NOT.Stabilize .AND. &
            ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )
            
!------------------------------------------------------------------------------
!           If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
            TimeForce  = 0.0_dp
            IF ( TransientAssembly ) THEN
               IF ( ConstantBulk ) THEN
                 CALL DefaultUpdateMass( MASS )
               ELSE
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
               END IF
            ELSE IF ( Solver % NOFEigenValues>0 ) THEN
              CALL DefaultUpdateDamp(MASS)
            END IF
!------------------------------------------------------------------------------
!           Update global matrices from local matrices
!------------------------------------------------------------------------------
            IF (  Bubbles ) THEN
               CALL Condensate( N, STIFF, FORCE, TimeForce )
            END IF

            IF (Phase24WallProfiling) Phase24Tick = RealTime()
            IF (Phase24FastElement) THEN
              CALL DefaultUpdateEquations( STIFF, FORCE, VecAssembly=.TRUE., &
                   Phase24PositionCache=.TRUE., RHSOnly=Phase24RHSOnly )
            ELSE
              CALL DefaultUpdateEquations( STIFF, FORCE, RHSOnly=Phase24RHSOnly )
            END IF
            IF (Phase24WallProfiling) Phase24Insert = Phase24Insert + RealTime() - Phase24Tick
         END IF
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
!------------------------------------------------------------------------------
      
      CALL DefaultFinishBulkAssembly()
      IF (Phase24WallProfiling) Phase24Bulk = RealTime() - Phase24BulkStart
      IF (Phase24WallProfiling) CALL CheckTimer('Phase24 Bulk Assembly',Level=4,Delete=.TRUE.)


1000  CONTINUE

     

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      IF (Phase24WallProfiling) CALL ResetTimer('Phase24 Boundary Assembly')
      Phase24BoundaryStart = RealTime()
      DO bc_elem = 1, Solver % Mesh % NumberOfBoundaryElements
        
        Element => GetBoundaryElement(bc_elem)
        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE
        
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()


        ! This checks whether there are any Dirichlet conditions on the 
        ! smart heater boundary. If there are the r.h.s. must be zero as 
        ! there can possibly not be any effect on temperature.
        !-----------------------------------------------------------------
	IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN
          IF( ListCheckPresent(BC, Varname) ) THEN
             nd = GetElementDOFs(Indexes)
             ForceHeater(TempPerm(Indexes(1:nd))) = 0.0_dp
          END IF
        END IF

        HeatFluxBC = GetLogical( BC, 'Heat Flux BC', Found )
        ! A boundary carrying only a fixed Dirichlet value does not change the
        ! finalized operator between nonlinear iterations.  Only active
        ! Neumann/Newton data, a gap, or radiation can introduce a changing
        ! boundary operator here; fixed Dirichlet elimination is handled by
        ! DefaultDirichletBCs and remains part of the stored CSR snapshot.
        IF ( Found .AND. .NOT. HeatFluxBC ) CYCLE

        HeatGapBC = ListGetLogical( BC, 'Heat Gap', Found )
        IF ( HeatGapBC .OR. ListCheckPresent(BC,'Radiation') ) &
             Phase24BoundaryMatrixInvariantNext = .FALSE.
        CALL AddHeatFluxBC()

        IF ( HeatGapBC ) THEN
          CALL FindGapIndexes( Element, Indexes, n )
          SaveIndexes(1:n) = Element % NodeIndexes
          Element % NodeIndexes = Indexes(1:n)
          CALL AddHeatFluxBC()
          Element % NodeIndexes = SaveIndexes(1:n)
        END IF

      END DO   ! Neumann & Newton BCs
      Phase24BoundaryMatrixInvariant = Phase24BoundaryMatrixInvariantNext
      IF (Phase24WallProfiling) Phase24Boundary = RealTime() - Phase24BoundaryStart
      IF (Phase24WallProfiling) CALL CheckTimer('Phase24 Boundary Assembly',Level=4,Delete=.TRUE.)
!------------------------------------------------------------------------------


      IF ( TransientSimulation .AND. ConstantBulk ) CALL AddGlobalTime()

      IF (Phase24WallProfiling) CALL ResetTimer('Phase24 Matrix Finalization')
      Phase24FinalizeStart = RealTime()
      CALL DefaultFinishBoundaryAssembly()
      CALL DefaultFinishAssembly()
      IF (Phase24WallProfiling) Phase24Finalize = RealTime() - Phase24FinalizeStart
      IF (Phase24WallProfiling) CALL CheckTimer('Phase24 Matrix Finalization',Level=4,Delete=.TRUE.)
      CALL Info( 'HeatSolve', 'Assembly done', Level=4 )

      CALL DefaultDirichletBCs()
      IF (Phase24LumpedPrec) CALL Phase24BuildLumpedPrec()

      IF (Phase24VectorAssembly .AND. Phase24MetadataReady) &
           CALL Phase24FinalizeOperatorLag()

      IF (Phase24VectorAssembly .AND. Phase24EpochReady) THEN
        IF (Phase24MatrixChanged) Phase24MatrixEpoch = Phase24MatrixEpoch + 1
        Phase24RHSEpoch = Phase24RHSEpoch + 1
        CALL ListAddInteger(SolverParams,'Phase24 Structure Epoch',Phase24StructureEpoch)
        CALL ListAddInteger(SolverParams,'Phase24 Matrix Epoch',Phase24MatrixEpoch)
        CALL ListAddInteger(SolverParams,'Phase24 RHS Epoch',Phase24RHSEpoch)
        CALL ListAddInteger(SolverParams,'Phase24 Full Rebuilds',Phase24FullRebuildCount)
        CALL ListAddInteger(SolverParams,'Phase24 Case A Reuses',Phase24CaseACount)
        CALL ListAddInteger(SolverParams,'Phase24 Case B Refreshes',Phase24CaseBCount)
        CALL ListAddInteger(SolverParams,'Phase24 Lagged Reuses',Phase24LaggedReuseCount)
        CALL ListAddInteger(SolverParams,'Phase24 Forced Rebuilds',Phase24ForcedRebuildCount)
        CALL ListAddInteger(SolverParams,'Phase24 Recovery Events',Phase24RecoveryCount)
        IF (Phase24RHSOnly) CALL Info('HeatSolve', &
             'Phase24 exact RHS-only assembly; matrix epoch unchanged',Level=4)
      END IF

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      at = CPUTime() - at
      st = CPUTime()

      PrevNorm = Norm

      IF(SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
      
        IF(.NOT. TransientHeaterControl) THEN
          ! These are control loops. Do use them to check convergence or advance the
          ! nonlinear iteration flag.
          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
          CALL ListAddLogical(SolverParams,'Skip Advance Nonlinear Iter',.TRUE.)

          Relax = GetCReal( SolverParams,'Nonlinear System Relaxation Factor', Found )
          
          IF ( Found .AND. Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,&
                'Nonlinear System Relaxation Factor', 1.0d0 )
          ELSE
            Relax = 1.0d0
          END IF          

          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              ForceHeater, XX, Norm, 1, Solver )         
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, YY, Norm, 1, Solver )

          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.FALSE.)
          CALL ListAddLogical(SolverParams,'Skip Advance Nonlinear Iter',.FALSE.)
        ELSE                    
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, Temperature, Norm, 1, Solver )
          YY = Temperature
        END IF

        IF(.NOT. SmartHeaterAverage) THEN
          xave = XX(TempPerm(SmartHeaterNode))
          yave = YY(TempPerm(SmartHeaterNode))
        ELSE          
          xave = 0.0d0
          yave = 0.0d0
          j = 0
          
          DO k = Model % Mesh % NumberOfBulkElements + 1, &
              Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements            

            Element => Model % Mesh % Elements(k)            
            IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
              l = Element % TYPE % NumberOfNodes
              j = j + l
              xave = xave + SUM( XX(TempPerm(Element % NodeIndexes)) )
              yave = yave + SUM( YY(TempPerm(Element % NodeIndexes)) )
            END IF
          END DO
          xave = xave / j
          yave = yave / j 
          CALL ListAddConstReal(Model % Simulation,'res: Smart Heater Temperature',yave)
        END IF

        IF(.NOT. TransientHeaterControl) THEN
          IF ( ASSOCIATED(Solver % Variable % NonlinValues) ) THEN
            Solver % Variable % NonlinValues = Temperature
          END IF

          PowerScaling = (MeltPoint - yave) / xave 
          Temperature = YY + PowerScaling * XX

          ! The change is computed separately for the controlled temperature field
          !-----------------------------------------------------------------------
          CALL ComputeChange(Solver,.FALSE.,LocalNodes,Temperature)
          Norm = Solver % Variable % Norm

        END IF

        IF(dt > PowerTimeScale) THEN
          IF ( Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', Relax )
          END IF
        END IF
      ELSE
!------------------------------------------------------------------------------
!     Check stepsize for nonlinear iteration
!------------------------------------------------------------------------------
        IF( DefaultLinesearch( Converged ) ) GOTO 500
        IF( Converged ) THEN
          IF (Phase24LaggingEnabled .AND. Phase24LaggedReuseThisTimestep > 0 .AND. &
               .NOT. Phase24FinalConfirmationDone) THEN
            Phase24FinalConfirmationDone = .TRUE.
            Phase24ForceRebuild = .TRUE.
            Converged = .FALSE.
            CALL Info('HeatSolve', &
                 'Phase24 lagging: scheduling fresh-matrix convergence confirmation',Level=4)
          ELSE
            EXIT
          END IF
        END IF

        IF (Phase24WallProfiling) CALL ResetTimer('Phase24 Linear Solve')
        IF (ListGetLogical(SolverParams,'Phase24 HYPRE Diagnostic',Found)) THEN
          TimeVar => VariableGet(Model % Mesh % Variables,'Time')
          IF (ASSOCIATED(TimeVar)) CALL ListAddConstReal(SolverParams, &
               'Phase24 Diagnostic Physical Time',TimeVar % Values(1))
          CALL ListAddConstReal(SolverParams,'Phase24 Diagnostic Dt',Solver % dt)
          CALL ListAddInteger(SolverParams,'Phase24 Diagnostic BDF Order', &
               MIN(Solver % DoneTime,Solver % Order))
          CALL ListAddInteger(SolverParams,'Phase24 Diagnostic Nonlinear Iteration',iter)
        END IF
        Phase24SolveStart = RealTime()
        IF (Phase24OuterCapture .AND. GetTimeStep() == Phase24OuterCaptureStep .AND. &
             iter <= Phase24OuterCaptureMaxIter) CALL Phase24CaptureOuterSystem('before')
        ! Opt-in: the first solve of a timestep starts from the previous
        ! step's field, where a tolerance relative to ||b|| (dominated by the
        ! bath/mass rows) stops before the TES/pulse update is resolved.  Solve
        ! it in residual form (tolerance relative to ||b-Ax0||) instead.
        FirstSolveResidual = .FALSE.
        IF (TransientSimulation .AND. iter == 1) THEN
          FirstSolveTol = GetConstReal(SolverParams, &
               'Phase24 First Solve Residual Tolerance', Found)
          IF (Found .AND. FirstSolveTol > 0.0_dp) THEN
            FirstSolveResidual = .TRUE.
            SavedResidualMode = GetLogical(SolverParams,'Linear System Residual Mode',Found)
            SavedLinTol = GetConstReal(SolverParams,'Linear System Convergence Tolerance',GotSavedTol)
            CALL ListAddLogical(SolverParams,'Linear System Residual Mode',.TRUE.)
            CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',FirstSolveTol)
          END IF
        END IF
        Norm = DefaultSolve()
        IF (FirstSolveResidual) THEN
          CALL ListAddLogical(SolverParams,'Linear System Residual Mode',SavedResidualMode)
          IF (GotSavedTol) CALL ListAddConstReal(SolverParams, &
               'Linear System Convergence Tolerance',SavedLinTol)
        END IF
        IF (Phase24OuterCapture .AND. GetTimeStep() == Phase24OuterCaptureStep .AND. &
             iter <= Phase24OuterCaptureMaxIter) CALL Phase24CaptureOuterSystem('after')
        IF (Phase24WallProfiling) Phase24Solve = RealTime() - Phase24SolveStart
        IF (Phase24WallProfiling) CALL CheckTimer('Phase24 Linear Solve',Level=4,Delete=.TRUE.)
      END IF

      IF (Phase24LaggingEnabled .AND. Phase24LaggingActive) THEN
        Phase24RecoveryThisSolve = ListGetLogical(SolverParams, &
             'Phase24 HYPRE Recovery Active',Found)
        IF (Phase24RecoveryThisSolve) THEN
          Phase24LaggingActive = .FALSE.
          Phase24FallbackRebuild = .TRUE.
          Phase24ForceRebuild = .TRUE.
          Phase24ForcedRebuildCount = Phase24ForcedRebuildCount + 1
          Phase24RecoveryCount = Phase24RecoveryCount + 1
          Phase24RecoveryCooldown = 1
          CALL Info('HeatSolve', &
               'Phase24 lagging fallback: HYPRE recovery requested rebuild-every-iteration',Level=4)
        ELSE IF (Phase24UsedLaggedOperator) THEN
          Phase24KrylovIterations = ListGetInteger(SolverParams, &
               'Phase24 HYPRE Native Last Krylov Iterations',Found)
          IF (Found .AND. Phase24KrylovBaseline > 0 .AND. &
               REAL(Phase24KrylovIterations,dp) > MAX( &
                 REAL(Phase24KrylovBaseline + Phase24KrylovAbsoluteGrowth,dp), &
                 (1.0_dp + Phase24KrylovRelativeGrowth) * &
                    REAL(Phase24KrylovBaseline,dp))) THEN
            Phase24LaggingActive = .FALSE.
            Phase24FallbackRebuild = .TRUE.
            Phase24ForceRebuild = .TRUE.
            Phase24ForcedRebuildCount = Phase24ForcedRebuildCount + 1
            CALL Info('HeatSolve', &
                 'Phase24 lagging fallback: Krylov iteration growth requested fresh operator',Level=4)
          END IF
        END IF
      END IF

      IF (Phase24LaggingEnabled) THEN
        Phase24KrylovIterations = ListGetInteger(SolverParams, &
             'Phase24 HYPRE Native Last Krylov Iterations',Found)
        IF (Found .AND. Phase24KrylovIterations > 0 .AND. &
             .NOT. Phase24UsedLaggedOperator .AND. .NOT. Phase24RecoveryThisSolve) &
          Phase24KrylovBaseline = Phase24KrylovIterations
        IF (.NOT. Phase24UsedLaggedOperator .AND. .NOT. Phase24RecoveryThisSolve .AND. &
             Phase24RecoveryCooldown > 0) THEN
          Phase24RecoveryCooldown = 0
          Phase24LaggingActive = .TRUE.
          CALL Info('HeatSolve', &
               'Phase24 lagging cooldown complete after fresh state',Level=4)
        END IF
      END IF

      IF (Phase24WallProfiling) CALL Phase24Report()


      IF( SmartHeaterControl .OR. IntegralHeaterControl) THEN
         CALL Info( 'HeatSolve', 'Heater Control Information', Level=4 )
         DO i=1,Model % NumberOfBodyForces
            IF( .NOT. (SmartHeaters(i) .OR. IntegralHeaters(i))) CYCLE
            IF( SmartHeaters(i) )  HeaterScaling(i) = PowerScaling

            WRITE( Message, '(A,T35,I15)' ) 'Heater for body: ', i
            CALL Info( 'HeatSolve', Message, Level=4 )
            IF(SmartHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Smart heater'
            IF(IntegralHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Integral heater'
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Volume (m^3): ', HeaterArea(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            s = HeaterSource(i) * HeaterScaling(i)
            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Power (W): ', s
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater scaling: ', HeaterScaling(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            WRITE( Message, '(A,T35,ES15.4)' ) 'Heater Power Density (W/kg): ', s/(HeaterDensity(i) * HeaterArea(i))
            CALL Info( 'HeatSolve', Message, Level=4 )

            CALL ListAddConstReal(Model % Simulation,'res: Heater Power Scaling '//I2S(i),HeaterScaling(i))
            CALL ListAddConstReal(Model % Simulation,'res: Heater Power Density '//I2S(i),&
                 s/(HeaterDensity(i) * HeaterArea(i)))
         END DO
      END IF


      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'HeatSolve', Message, Level=4 )
!------------------------------------------------------------------------------
!     If modelling phase change (and if requested by the user), check if any
!     node has jumped over the phase change interval, and if so, reduce
!     timestep and or relaxation and recompute.
!------------------------------------------------------------------------------
      IF (PhaseChange .AND. CheckLatentHeatRelease ) THEN
!------------------------------------------------------------------------------
        IF ( CheckLatentHeat() ) THEN
          Temperature(1:LocalNodes) = PrevSolution
          Norm = PrevNorm

          IF ( TransientSimulation ) THEN
            dt = dt / 2
            Solver % dt = dt
            WRITE( Message, * ) &
                  'Latent heat release check: reducing timestep to: ',dt
            CALL Info( 'HeatSolve', Message, Level=4 )
          ELSE
            Relax = Relax / 2
            CALL  ListAddConstReal( Solver % Values,  &
                 'Nonlinear System Relaxation Factor', Relax )
            WRITE( Message, * ) &
                 'Latent heat release check: reducing relaxation to: ',Relax
            CALL Info( 'HeatSolve', Message, Level=4 )
          END IF

          CYCLE
        END IF
        IF ( .NOT.TransientSimulation ) PrevSolution=Temperature(1:LocalNodes)
      END IF
!------------------------------------------------------------------------------
     
      RelativeChange = Solver % Variable % NonlinChange
      IF (Phase24VectorAssembly) CALL ListAddConstReal(SolverParams, &
           'Phase24 Nonlinear Relative Change',RelativeChange)

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'HeatSolve', Message, Level=4 )

      IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
               NewtonLinearization = .TRUE.
      Converged =  ( Solver % Variable % NonlinConverged > 0 ) .AND. &
          ( .NOT. SmartHeaterControl .OR. SmartTolReached )
      ! Opt-in: also require the TES mean temperature to have settled between
      ! the last two iterates of this step (field-norm changes are dominated
      ! by the bath/absorber and can stop while the TES circuit still moves).
      ! Opt-in: converge on the coupled residual of the final state instead of
      ! on iterate changes.  The thermal solve is exact for the power it used,
      ! so the remaining electrothermal inconsistency is P_circuit(T_TES) -
      ! P_used.  Unlike increment tests this cannot be fooled by a good
      ! (extrapolated) initial guess, and it allows a one-solve step when the
      ! predictor is already consistent.
      TESPowerTol = -1.0_dp
      TESPowerTol = GetConstReal(SolverParams, &
           'TES Inner Circuit Power Residual Tolerance', Found)
      IF (.NOT. Found) TESPowerTol = -1.0_dp
      TESCoupledTolT = -1.0_dp
      TESCoupledTolT = GetConstReal(SolverParams, &
           'Phase24 Coupled Temperature Tolerance', Found)
      TESAllValid = ALLOCATED(TESC)
      IF (TESAllValid) TESAllValid = ALL(TESC(:) % ResValid)
      IF (TESCoupledTolT > 0.0_dp .AND. TESAllValid) THEN
        ! Both remaining nonlinearities, expressed as a bound on the TES
        ! temperature error of accepting this iterate:
        !  circuit: dT = F/(1/S + |dP/dT|) <= |F|/|dP/dT| for the negative
        !           electrothermal feedback (dP/dT < 0), without knowing the
        !           step's thermal response S;
        !  k(T):    the membrane conductivity was evaluated at the iterate the
        !           solve started from; dT <= max|dk/k| * membrane T drop.
        ! This replaces the field-increment test, which measures predictor
        ! error rather than inconsistency.
        ! Worst circuit (each TES body has its own circuit).
        TESdTCircuit = 0.0_dp
        DO TESk = 1, TESNCirc
          TESMeanNow = TESMeanTemperatureNow(TESk)
          TESPowerRes = TESCircuitPowerAt(TESk, TESMeanNow) - TESC(TESk) % ResP
          TESdPdT = (TESCircuitPowerAt(TESk, TESMeanNow+1.0e-7_dp) - &
               TESCircuitPowerAt(TESk, TESMeanNow-1.0e-7_dp)) / 2.0e-7_dp
          IF (TESdPdT < 0.0_dp) THEN
            TESdTOne = ABS(TESPowerRes) / ABS(TESdPdT)
          ELSE IF (iter > 1) THEN
            ! Steps much shorter than L/R hold the current (dP/dT > 0): no
            ! S-free bound exists, so use the TES change between iterates.
            TESdTOne = ABS(TESMeanNow - TESC(TESk) % ExitPrevT)
          ELSE
            TESdTOne = HUGE(1.0_dp)
          END IF
          TESC(TESk) % ExitPrevT = TESMeanNow
          TESdTCircuit = MAX(TESdTCircuit, TESdTOne)
        END DO
        TESdkRel = 0.0_dp
        TESdkTmax = -HUGE(1.0_dp)
        TESdkTmin = HUGE(1.0_dp)
        IF (Phase24DynKNodeCount > 0) THEN
          ALLOCATE(TESdkNow(Phase24DynKNodeCount))
          CALL Phase24EvaluateDynK(TESdkNow)
          DO i = 1, Phase24DynKNodeCount
            IF (Phase24DynKNodal(i) > 0.0_dp) TESdkRel = MAX(TESdkRel, &
                 ABS(TESdkNow(i) - Phase24DynKNodal(i)) / Phase24DynKNodal(i))
            j = TempPerm(Phase24DynKNodes(i))
            IF (j > 0) THEN
              TESdkTmax = MAX(TESdkTmax, Temperature(j))
              TESdkTmin = MIN(TESdkTmin, Temperature(j))
            END IF
          END DO
          DEALLOCATE(TESdkNow)
        END IF
        TESdkRel = ParallelReduction(TESdkRel, 2)
        TESdkTmax = ParallelReduction(TESdkTmax, 2)
        TESdkTmin = ParallelReduction(TESdkTmin, 1)
        IF (TESdkTmax >= TESdkTmin) THEN
          TESdTk = TESdkRel * (TESdkTmax - TESdkTmin)
          Converged = TESdTCircuit <= TESCoupledTolT .AND. TESdTk <= TESCoupledTolT .AND. &
               ( .NOT. SmartHeaterControl .OR. SmartTolReached )
        ELSE
          ! No temperature-dependent conductivity: only the circuit couples.
          TESdTk = 0.0_dp
          IF (TESdTCircuit > TESCoupledTolT) Converged = .FALSE.
        END IF
        WRITE(Message,'(A,I0,A,ES11.3,A,ES11.3,A,ES11.3)') 'TES coupled error: iter=',iter, &
             ', dT_circuit=',TESdTCircuit,', dT_k=',TESdTk,', field change=',RelativeChange
        CALL Info('HeatSolve',TRIM(Message),Level=4)
      ELSE IF (TESPowerTol > 0.0_dp .AND. TESAllValid) THEN
        TESdTOne = 0.0_dp
        DO TESk = 1, TESNCirc
          TESMeanNow = TESMeanTemperatureNow(TESk)
          TESPowerRes = TESCircuitPowerAt(TESk, TESMeanNow) - TESC(TESk) % ResP
          TESdTOne = MAX(TESdTOne, ABS(TESPowerRes)/MAX(ABS(TESC(TESk) % ResP),1.0e-30_dp))
        END DO
        IF (TESdTOne > TESPowerTol) Converged = .FALSE.
        WRITE(Message,'(A,I0,A,ES12.4,A,ES12.4)') 'TES coupled residual: iter=',iter, &
             ', |dP|/P=',TESdTOne,', field change=',RelativeChange
        CALL Info('HeatSolve',TRIM(Message),Level=4)
      ELSE IF (TransientSimulation) THEN
        TESTempTol = GetConstReal(SolverParams, &
             'TES Inner Circuit Temperature Tolerance', Found)
        IF (Found .AND. TESTempTol > 0.0_dp .AND. ALLOCATED(TESC)) THEN
          DO TESk = 1, TESNCirc
            TESMeanNow = TESMeanTemperatureNow(TESk)
            IF (iter <= 1) THEN
              Converged = .FALSE.
            ELSE IF (ABS(TESMeanNow - TESC(TESk) % ExitPrevT) > TESTempTol) THEN
              Converged = .FALSE.
            END IF
            TESC(TESk) % ExitPrevT = TESMeanNow
          END DO
        END IF
      END IF
      IF (Phase24LaggingEnabled .AND. Phase24LaggingActive .AND. &
           Phase24LagReuseCount > 0 .AND. Phase24PreviousRelativeChange >= 0.0_dp .AND. &
           RelativeChange > MAX(Phase24PreviousRelativeChange,1.0e-30_dp) * &
                Phase24LagResidualGrowth) THEN
        Phase24LaggingActive = .FALSE.
        Phase24FallbackRebuild = .TRUE.
        Phase24ForceRebuild = .TRUE.
        Phase24ForcedRebuildCount = Phase24ForcedRebuildCount + 1
        CALL Info('HeatSolve', &
             'Phase24 lagging fallback: nonlinear change increased; rebuilding every iteration',Level=4)
      END IF
      Phase24PreviousRelativeChange = RelativeChange
      IF( Converged .AND. iter >= MinNonlinearIter ) THEN
        IF (Phase24LaggingEnabled .AND. Phase24LaggedReuseThisTimestep > 0 .AND. &
             .NOT. Phase24FinalConfirmationDone) THEN
          Phase24FinalConfirmationDone = .TRUE.
          Phase24ForceRebuild = .TRUE.
          Converged = .FALSE.
          CALL Info('HeatSolve', &
               'Phase24 lagging: scheduling fresh-matrix convergence confirmation',Level=4)
        ELSE
          EXIT
        END IF
      END IF

      IF(SmartHeaterControl) THEN
        IF ( RelativeChange < SmartTol ) THEN
          SmartTolReached = .TRUE.
          YY = Temperature
        END IF
      END IF
      
!------------------------------------------------------------------------------
    END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------
    IF(TransientHeaterControl) THEN
      PowerRelax = GetCReal(Solver % Values,'Smart Heater Relaxation Factor', GotIt)
      IF(.NOT. GotIt) PowerRelax = 1.0_dp
      PowerSensitivity = ListGetConstReal(Solver % Values,'Smart Heater Power Sensivity',GotIt)
      IF(.NOT. GotIt) PowerSensitivity = 4.0_dp
      PowerScaling = PowerScaling * (1 + PowerSensitivity * PowerRelax * (MeltPoint/yave - 1.0d0) ) 

      IF( ListGetLogical( Solver % Values,'Smart Heater Transient Speedup',GotIt ) ) THEN
        Temperature = Temperature * (1 + PowerRelax * (MeltPoint/yave - 1.0d0)   )     
      END IF
      YY = Temperature
    END IF

!------------------------------------------------------------------------------
!   Compute cumulative time done by now and time remaining
!------------------------------------------------------------------------------
    IF ( .NOT. TransientSimulation ) EXIT
    CumulativeTime = CumulativeTime + dt
    dt = Timestep - CumulativeTime

   END DO ! time interval
   Solver % dt = Timestep

   IF (TransientSimulation .AND. ListGetLogical(Model % Simulation, &
        'Phase24 Step Controller',Found)) CALL Phase24StepController()

   CALL DefaultFinish()
   
!------------------------------------------------------------------------------
   CALL  ListAddConstReal( Solver % Values,  &
        'Nonlinear System Relaxation Factor', SaveRelax )
!------------------------------------------------------------------------------

   DEALLOCATE( PrevSolution )

   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) THEN
     IF(.NOT. ListGetLogical( Solver % Values,'Library Adaptivity',Found )) THEN
       CALL RefineMesh( Model,Solver,Temperature,TempPerm, &
           HeatSolver_Inside_Residual, HeatSolver_Edge_Residual, &
           HeatSolver_Boundary_Residual )
     END IF
   END IF
     
CONTAINS

  ! Save exactly the outer HeatSolve system around DefaultSolve.  The regular
  ! save hook in DefaultSolve records the pre-solve state; the explicit second
  ! call below records the post-solve state under a distinct prefix.  This
  ! routine is intentionally adjacent to DefaultSolve, so captures cannot be
  ! mistaken for an inner preconditioner, Schur, or restricted solve.
  SUBROUTINE Phase24CaptureOuterSystem(Stage)
    CHARACTER(LEN=*), INTENT(IN) :: Stage
    CHARACTER(LEN=MAX_NAME_LEN) :: CaptureBase, SavedPrefix
    LOGICAL :: PrefixFound

    WRITE(CaptureBase,'(A,"/ts",I4.4,"_nl",I4.4,"/outer_")') &
         TRIM(Phase24OuterCapturePrefix), GetTimeStep(), iter
    SavedPrefix = GetString(SolverParams,'Linear System Save Prefix',PrefixFound)
    IF (TRIM(Stage) == 'before') THEN
      CALL ListAddString(SolverParams,'Linear System Save Prefix', &
           TRIM(CaptureBase)//'before')
      CALL SaveLinearSystem(Solver)
    ELSE IF (TRIM(Stage) == 'after') THEN
      CALL ListAddString(SolverParams,'Linear System Save Prefix', &
           TRIM(CaptureBase)//'after')
      CALL SaveLinearSystem(Solver)
    ELSE
      CALL Fatal('Phase24CaptureOuterSystem','Unknown capture stage')
    END IF
    IF (PrefixFound) THEN
      CALL ListAddString(SolverParams,'Linear System Save Prefix',TRIM(SavedPrefix))
    END IF

    WRITE(Message,'(A,A,A,I0,A,I0,A,ES16.8,A,ES16.8,A,I0,A,I0,A,I0)') &
         'PHASE24_OUTER_CAPTURE stage=',TRIM(Stage),', timestep=',GetTimeStep(), &
         ', nonlinear_iteration=',iter,', time=',GetTime(),', dt=',Solver % dt, &
         ', bdf_order=',MIN(Solver % DoneTime,Solver % Order),', matrix_epoch=', &
         Phase24MatrixEpoch,', rhs_epoch=',Phase24RHSEpoch
    CALL Info('HeatSolve',TRIM(Message),Level=4)
  END SUBROUTINE Phase24CaptureOuterSystem

  SUBROUTINE Phase24GetBDFCoefficients(Order,HNew,HOld,A0,A1,A2)
    INTEGER, INTENT(IN) :: Order
    REAL(KIND=dp), INTENT(IN) :: HNew
    REAL(KIND=dp), INTENT(OUT) :: HOld,A0,A1,A2
    TYPE(Variable_t), POINTER :: DtVar

    HOld = HNew
    A0 = 0.0_dp
    A1 = 0.0_dp
    A2 = 0.0_dp
    IF (HNew <= 0.0_dp) RETURN

    SELECT CASE(MIN(MAX(Order,0),5))
    CASE(1)
      A0 = 1.0_dp / HNew
      A1 = -A0
    CASE(2)
      DtVar => VariableGet(Solver % Mesh % Variables,'Timestep size')
      IF (ASSOCIATED(DtVar)) THEN
        IF (ASSOCIATED(DtVar % PrevValues)) THEN
          HOld = DtVar % PrevValues(1,1)
        END IF
      END IF
      IF (HOld <= 0.0_dp) HOld = HNew
      A0 = 1.0_dp / HNew + 1.0_dp / (HNew + HOld)
      A1 = -(1.0_dp + HNew / HOld) / HNew
      A2 = (HNew / HOld) / (HNew + HOld)
    CASE(3)
      A0 = (11.0_dp/6.0_dp) / HNew
      A1 = -3.0_dp / HNew
      A2 = (3.0_dp/2.0_dp) / HNew
    CASE(4)
      A0 = (25.0_dp/12.0_dp) / HNew
      A1 = -4.0_dp / HNew
      A2 = 3.0_dp / HNew
    CASE(5)
      A0 = (137.0_dp/60.0_dp) / HNew
      A1 = -5.0_dp / HNew
      A2 = 5.0_dp / HNew
    END SELECT
  END SUBROUTINE Phase24GetBDFCoefficients

  SUBROUTINE Phase24FinalizeOperatorLag()
    INTEGER :: q
    REAL(KIND=dp) :: denom,delta

    IF (.NOT. Phase24MetadataReady) RETURN
    IF (.NOT. ALLOCATED(Phase24OperatorKCurrent)) RETURN

    Phase24CoeffRelChange = 0.0_dp
    IF (Phase24OperatorReferenceValid) THEN
      DO q=1,Phase24MetadataElements
        IF (Phase24ElementClass(q) /= 2 .AND. Phase24ElementClass(q) /= 3) CYCLE
        denom = MAX(MAX(ABS(Phase24OperatorKReference(q)), &
             ABS(Phase24OperatorKCurrent(q))),1.0e-30_dp)
        delta = ABS(Phase24OperatorKCurrent(q)-Phase24OperatorKReference(q)) / denom
        Phase24CoeffRelChange = MAX(Phase24CoeffRelChange,delta)
        denom = MAX(MAX(ABS(Phase24OperatorMassReference(q)), &
             ABS(Phase24OperatorMassCurrent(q))),1.0e-30_dp)
        delta = ABS(Phase24OperatorMassCurrent(q)-Phase24OperatorMassReference(q)) / denom
        Phase24CoeffRelChange = MAX(Phase24CoeffRelChange,delta)
        denom = MAX(MAX(ABS(Phase24OperatorReactionReference(q)), &
             ABS(Phase24OperatorReactionCurrent(q))),1.0e-30_dp)
        delta = ABS(Phase24OperatorReactionCurrent(q)- &
             Phase24OperatorReactionReference(q)) / denom
        Phase24CoeffRelChange = MAX(Phase24CoeffRelChange,delta)
      END DO
    END IF

    Phase24LagCandidate = Phase24LaggingActive .AND. &
         Phase24OperatorReferenceValid .AND. Phase24OperatorSnapshotValid .AND. &
         Phase24BoundaryMatrixInvariant .AND. .NOT. Phase24LHSUnsafe .AND. &
         .NOT. Phase24ForceRebuild .AND. &
         Phase24CoeffRelChange <= Phase24LagThreshold .AND. &
         Phase24LagReuseCount < Phase24LagMaxReuse
    IF (Phase24LagCandidate .AND. Phase24LagMode == 1 .AND. &
         Phase24PreviousRelativeChange >= 0.0_dp) THEN
      ! Do not lag while the nonlinear state is still moving materially.  The
      ! adaptive mode uses the same safety gate; its additional freedom is the
      ! coefficient threshold and bounded reuse, not an undocumented solver
      ! tolerance change.
      IF (Phase24PreviousRelativeChange > MAX(1.0e-3_dp, &
           100.0_dp*MAX(NonlinearTol,1.0e-12_dp))) Phase24LagCandidate = .FALSE.
    END IF
    WRITE(Message,'(A,ES12.4,A,L1,A,L1,A,L1,A,L1)') &
         'Phase24 operator decision: coeff_change=',Phase24CoeffRelChange, &
         ', boundary_invariant=',Phase24BoundaryMatrixInvariant, &
         ', lhs_safe=',.NOT. Phase24LHSUnsafe, ', force=',Phase24ForceRebuild, &
         ', candidate=',Phase24LagCandidate
    CALL Info('HeatSolve',Message,Level=5)

    IF (Phase24LagCandidate) THEN
      Solver % Matrix % Values = Phase24LaggedCSRValues
      Phase24RHSOnly = .TRUE.
      Phase24UsedLaggedOperator = .TRUE.
      Phase24MatrixChanged = .FALSE.
      Phase24LagReuseCount = Phase24LagReuseCount + 1
      Phase24LaggedReuseThisTimestep = Phase24LaggedReuseThisTimestep + 1
      Phase24LaggedReuseCount = Phase24LaggedReuseCount + 1
      WRITE(Message,'(A,ES12.4,A,I0)') 'Phase24 adaptive operator lag reuse; coeff_change=', &
           Phase24CoeffRelChange, ', reuse=',Phase24LagReuseCount
      CALL Info('HeatSolve',Message,Level=4)
    ELSE
      IF (Phase24ForceRebuild) Phase24ForcedRebuildCount = Phase24ForcedRebuildCount + 1
      IF (Phase24OperatorReferenceValid .AND. Phase24LaggingActive .AND. &
           Phase24CoeffRelChange > Phase24LagThreshold) &
           Phase24LagReuseCount = 0
      IF (.NOT. Phase24RHSOnly) THEN
        IF (.NOT. ALLOCATED(Phase24LaggedCSRValues) .OR. &
             SIZE(Phase24LaggedCSRValues) /= SIZE(Solver % Matrix % Values)) THEN
          IF (ALLOCATED(Phase24LaggedCSRValues)) DEALLOCATE(Phase24LaggedCSRValues)
          ALLOCATE(Phase24LaggedCSRValues(SIZE(Solver % Matrix % Values)))
        END IF
        Phase24LaggedCSRValues = Solver % Matrix % Values
        Phase24OperatorKReference = Phase24OperatorKCurrent
        Phase24OperatorMassReference = Phase24OperatorMassCurrent
        Phase24OperatorReactionReference = Phase24OperatorReactionCurrent
        Phase24OperatorReferenceValid = .TRUE.
        Phase24OperatorSnapshotValid = .TRUE.
        Phase24OperatorDt = Solver % dt
        Phase24OperatorOrder = MIN(Solver % DoneTime, Solver % Order)
        Phase24FullRebuildCount = Phase24FullRebuildCount + 1
        CALL Info('HeatSolve','Phase24 fresh operator stored for adaptive lagging',Level=4)
      ELSE IF (.NOT. Phase24OperatorReferenceValid) THEN
        ! Static-only RHS reuse has no dynamic coefficient samples, but its
        ! matrix is still a valid reference for the lifecycle.
        Phase24OperatorKReference = Phase24OperatorKCurrent
        Phase24OperatorMassReference = Phase24OperatorMassCurrent
        Phase24OperatorReactionReference = Phase24OperatorReactionCurrent
        Phase24OperatorReferenceValid = .TRUE.
        Phase24OperatorDt = Solver % dt
        Phase24OperatorOrder = MIN(Solver % DoneTime, Solver % Order)
      END IF
      Phase24LagReuseCount = 0
    END IF
    Phase24ForceRebuild = .FALSE.
    CALL ListAddConstReal(SolverParams,'Phase24 Coefficient Relative Change',Phase24CoeffRelChange)
  END SUBROUTINE Phase24FinalizeOperatorLag

  SUBROUTINE Phase24ClearMetadata()
    IF (ALLOCATED(Phase24ElementClass)) DEALLOCATE(Phase24ElementClass)
    IF (ALLOCATED(Phase24ElementBody)) DEALLOCATE(Phase24ElementBody)
    IF (ALLOCATED(Phase24ElementPart)) DEALLOCATE(Phase24ElementPart)
    IF (ALLOCATED(Phase24ElementNodes)) DEALLOCATE(Phase24ElementNodes)
    IF (ALLOCATED(Phase24ElementEq)) DEALLOCATE(Phase24ElementEq)
    IF (ALLOCATED(Phase24ElementGrad)) DEALLOCATE(Phase24ElementGrad)
    IF (ALLOCATED(Phase24ElementVolume)) DEALLOCATE(Phase24ElementVolume)
    IF (ALLOCATED(Phase24StaticStiff)) DEALLOCATE(Phase24StaticStiff)
    IF (ALLOCATED(Phase24SpatialMass)) DEALLOCATE(Phase24SpatialMass)
    IF (ALLOCATED(Phase24StaticK)) DEALLOCATE(Phase24StaticK)
    IF (ALLOCATED(Phase24StaticRhoCp)) DEALLOCATE(Phase24StaticRhoCp)
    IF (ALLOCATED(Phase24StaticElements)) DEALLOCATE(Phase24StaticElements)
    IF (ALLOCATED(Phase24DynamicElements)) DEALLOCATE(Phase24DynamicElements)
    IF (ALLOCATED(Phase24GenericElements)) DEALLOCATE(Phase24GenericElements)
    IF (ALLOCATED(Phase24PrismElements)) DEALLOCATE(Phase24PrismElements)
    IF (ALLOCATED(Phase24PrismSlot)) DEALLOCATE(Phase24PrismSlot)
    IF (ALLOCATED(Phase24PrismNQ)) DEALLOCATE(Phase24PrismNQ)
    IF (ALLOCATED(Phase24PrismGrad)) DEALLOCATE(Phase24PrismGrad)
    IF (ALLOCATED(Phase24PrismWeight)) DEALLOCATE(Phase24PrismWeight)
    IF (ALLOCATED(Phase24PrismStiff)) DEALLOCATE(Phase24PrismStiff)
    IF (ALLOCATED(Phase24PrismMass)) DEALLOCATE(Phase24PrismMass)
    IF (ALLOCATED(Phase24PrismK)) DEALLOCATE(Phase24PrismK)
    IF (ALLOCATED(Phase24PrismRhoCp)) DEALLOCATE(Phase24PrismRhoCp)
    IF (ALLOCATED(Phase24PrismVolume)) DEALLOCATE(Phase24PrismVolume)
    IF (ALLOCATED(Phase24PrismBasis)) DEALLOCATE(Phase24PrismBasis)
    IF (ALLOCATED(Phase24PrismDynamicK)) DEALLOCATE(Phase24PrismDynamicK)
    IF (ALLOCATED(Phase24PrismKNodeSlot)) DEALLOCATE(Phase24PrismKNodeSlot)
    IF (ALLOCATED(Phase24DynKNodes)) DEALLOCATE(Phase24DynKNodes)
    IF (ALLOCATED(Phase24DynKNodal)) DEALLOCATE(Phase24DynKNodal)
    IF (ALLOCATED(Phase24DynKMats)) DEALLOCATE(Phase24DynKMats)
    IF (ALLOCATED(Phase24DynKStart)) DEALLOCATE(Phase24DynKStart)
    IF (ALLOCATED(Phase24DynKCount)) DEALLOCATE(Phase24DynKCount)
    IF (ALLOCATED(Phase24ElementDynMat)) DEALLOCATE(Phase24ElementDynMat)
    Phase24DynKNodeCount = 0
    Phase24DynKMatCount = 0
    Phase24DynamicPrismCount = 0
    IF (ALLOCATED(Phase24ElementSrc)) DEALLOCATE(Phase24ElementSrc)
    IF (ALLOCATED(Phase24SrcBF)) DEALLOCATE(Phase24SrcBF)
    IF (ALLOCATED(Phase24SrcNodeStart)) DEALLOCATE(Phase24SrcNodeStart)
    IF (ALLOCATED(Phase24SrcNodeCount)) DEALLOCATE(Phase24SrcNodeCount)
    IF (ALLOCATED(Phase24SrcNodeList)) DEALLOCATE(Phase24SrcNodeList)
    IF (ALLOCATED(Phase24SrcLoad)) DEALLOCATE(Phase24SrcLoad)
    Phase24SrcCount = 0
    IF (ALLOCATED(Phase24OperatorKReference)) DEALLOCATE(Phase24OperatorKReference)
    IF (ALLOCATED(Phase24OperatorMassReference)) DEALLOCATE(Phase24OperatorMassReference)
    IF (ALLOCATED(Phase24OperatorReactionReference)) DEALLOCATE(Phase24OperatorReactionReference)
    IF (ALLOCATED(Phase24OperatorKCurrent)) DEALLOCATE(Phase24OperatorKCurrent)
    IF (ALLOCATED(Phase24OperatorMassCurrent)) DEALLOCATE(Phase24OperatorMassCurrent)
    IF (ALLOCATED(Phase24OperatorReactionCurrent)) DEALLOCATE(Phase24OperatorReactionCurrent)
    IF (ALLOCATED(Phase24LaggedCSRValues)) DEALLOCATE(Phase24LaggedCSRValues)
    Phase24MetadataElements = 0
    Phase24StaticCount = 0
    Phase24DynamicCount = 0
    Phase24GenericCount = 0
    Phase24FastTetCount = 0
    Phase24GenericPrismCount = 0
    Phase24GenericOtherCount = 0
    Phase24PrismCount = 0
    Phase24FastPrismCount = 0
    Phase24MetadataReady = .FALSE.
    Phase24OperatorReferenceValid = .FALSE.
    Phase24OperatorSnapshotValid = .FALSE.
    Phase24OperatorDt = -1.0_dp
    Phase24OperatorOrder = -1
    IF (ALLOCATED(Phase24StaticCSRValues)) DEALLOCATE(Phase24StaticCSRValues)
    Phase24StaticCacheDt = -1.0_dp
    Phase24StaticCacheA0 = -1.0_dp
    Phase24StaticCacheOrder = -1
    Phase24StaticCacheValueCount = 0
    Phase24StaticCacheLumped = .FALSE.
    Phase24StaticCacheValid = .FALSE.
  END SUBROUTINE Phase24ClearMetadata

  SUBROUTINE Phase24BuildMetadata()
    INTEGER :: i,j,k,nBulk,idx,matId,q,nPrism
    INTEGER :: nStatic,nDynamic,nGeneric
    INTEGER :: nodeIds(4),eqIds(4),prismNodeIds(6),prismEqIds(6)
    REAL(KIND=dp) :: detJ,volume,kappa,cp,rho
    REAL(KIND=dp) :: Basis4(4),grad4(4,3),Basis6(6),grad6(6,3),wq
    REAL(KIND=dp) :: StaticK,StaticRhoCp
    REAL(KIND=dp) :: fastPercent
    LOGICAL :: stat,found,foundK,foundCp,foundKappa,foundRho,eligible
    LOGICAL :: foundEquation,foundBodyForce,staticEligible,prismEligible
    CHARACTER(LEN=MAX_NAME_LEN) :: metadataConvection,metadataPhase
    TYPE(Element_t), POINTER :: MetadataElement
    TYPE(ValueList_t), POINTER :: MetadataMaterial
    TYPE(ValueList_t), POINTER :: MetadataEquation,MetadataBodyForce
    TYPE(Nodes_t) :: MetadataNodes
    TYPE(GaussIntegrationPoints_t) :: PrismIntegStuff
    LOGICAL :: kConstant,dynK
    LOGICAL, ALLOCATABLE :: DynFlag(:)
    INTEGER, ALLOCATABLE :: NodeMark(:)
    INTEGER :: bodyMat,nDynNodes,bfId,srcBF(64),nSrc,srcTotal,m
    TYPE(ValueListEntry_t), POINTER :: kEntry

    CALL Phase24ClearMetadata()
    nBulk = Solver % Mesh % NumberOfBulkElements
    IF (nBulk <= 0) THEN
      Phase24MetadataReady = .TRUE.
      RETURN
    END IF

    ALLOCATE(Phase24ElementClass(nBulk), Phase24ElementBody(nBulk), &
         Phase24ElementPart(nBulk), &
         Phase24ElementNodes(6,nBulk), Phase24ElementEq(6,nBulk), &
         Phase24ElementGrad(4,3,nBulk), Phase24ElementVolume(nBulk), &
         Phase24StaticStiff(4,4,nBulk), Phase24SpatialMass(4,4,nBulk), &
         Phase24StaticK(nBulk), Phase24StaticRhoCp(nBulk), &
         Phase24StaticElements(nBulk), Phase24DynamicElements(nBulk), &
         Phase24GenericElements(nBulk), Phase24PrismElements(nBulk), &
         Phase24OperatorKReference(nBulk), Phase24OperatorMassReference(nBulk), &
         Phase24OperatorReactionReference(nBulk), Phase24OperatorKCurrent(nBulk), &
         Phase24OperatorMassCurrent(nBulk), Phase24OperatorReactionCurrent(nBulk))
    ALLOCATE(DynFlag(nBulk), Phase24ElementSrc(nBulk), Phase24ElementDynMat(nBulk))
    DynFlag = .FALSE.
    Phase24ElementDynMat = 0
    Phase24ElementSrc = 0
    Phase24ElementClass = 3
    Phase24ElementBody = 0
    Phase24ElementPart = -1
    Phase24ElementNodes = 0
    Phase24ElementEq = 0
    Phase24ElementGrad = 0.0_dp
    Phase24ElementVolume = 0.0_dp
    Phase24StaticStiff = 0.0_dp
    Phase24SpatialMass = 0.0_dp
    Phase24StaticK = 0.0_dp
    Phase24StaticRhoCp = 0.0_dp
    Phase24OperatorKReference = 0.0_dp
    Phase24OperatorMassReference = 0.0_dp
    Phase24OperatorReactionReference = 0.0_dp
    Phase24OperatorKCurrent = 0.0_dp
    Phase24OperatorMassCurrent = 0.0_dp
    Phase24OperatorReactionCurrent = 0.0_dp
    nStatic = 0
    nDynamic = 0
    nGeneric = 0
    nPrism = 0

    DO i=1,nBulk
      MetadataElement => Solver % Mesh % Elements(i)
      Phase24ElementBody(i) = MetadataElement % BodyId
      Phase24ElementPart(i) = MetadataElement % PartIndex

      ! Production hybrid meshes use code 706 for linear six-node wedges.
      ! Classify a wedge as fast only when its material and equation reduce to
      ! the same constant scalar diffusion/mass form as the tetra shortcut.
      prismEligible = ASSOCIATED(MetadataElement % TYPE) .AND. &
           MetadataElement % TYPE % ElementCode == 706 .AND. &
           SIZE(MetadataElement % NodeIndexes) >= 6
      IF (prismEligible) THEN
        prismNodeIds = MetadataElement % NodeIndexes(1:6)
        prismEqIds = 0
        IF (ALL(prismNodeIds > 0) .AND. MAXVAL(prismNodeIds) <= SIZE(TempPerm)) &
             prismEqIds = TempPerm(prismNodeIds)
        PrismIntegStuff = GaussPoints(MetadataElement)
        prismEligible = ALL(prismEqIds > 0) .AND. &
             PrismIntegStuff % n <= PHASE24_PRISM_MAX_Q
        MetadataMaterial => GetMaterial(MetadataElement, found)
        IF (.NOT. found .OR. .NOT. ASSOCIATED(MetadataMaterial)) prismEligible = .FALSE.
        IF (prismEligible) THEN
          kappa = GetCReal(MetadataMaterial, 'Heat Conductivity', foundK)
          cp = GetCReal(MetadataMaterial, 'Heat Capacity', foundCp)
          rho = GetCReal(MetadataMaterial, 'Density', foundRho)
          prismEligible = foundK .AND. foundCp .AND. foundRho
          ! Same rule as the tetra path: a Variable/Procedure/MATC property
          ! must not be cached (it would be frozen at its metadata-time value).
          ! Density and heat capacity must be constant (cached mass); a
          ! solution-dependent *scalar* conductivity of a single material is
          ! handled as a dynamic-k wedge whose stiffness is rebuilt each time.
          dynK = .FALSE.
          IF (prismEligible) prismEligible = &
               ListCheckIsConstant(MetadataMaterial, 'Heat Capacity') .AND. &
               ListCheckIsConstant(MetadataMaterial, 'Density')
          IF (prismEligible) THEN
            kConstant = ListCheckIsConstant(MetadataMaterial, 'Heat Conductivity')
            IF (.NOT. kConstant) THEN
              prismEligible = .FALSE.
              kEntry => ListFind(MetadataMaterial, 'Heat Conductivity', found)
              bodyMat = ListGetInteger(CurrentModel % Bodies(MetadataElement % BodyId) % Values, &
                   'Material', found)
              IF (ASSOCIATED(kEntry) .AND. found .AND. bodyMat > 0) THEN
                IF (kEntry % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
                     kEntry % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR .OR. &
                     kEntry % TYPE == LIST_TYPE_CONSTANT_SCALAR_PROC) THEN
                  prismEligible = .TRUE.
                  dynK = .TRUE.
                  Phase24ElementDynMat(i) = bodyMat
                END IF
              END IF
            END IF
          END IF
          MetadataEquation => GetEquation(MetadataElement, foundEquation)
          MetadataBodyForce => GetBodyForce(MetadataElement, foundBodyForce)
          IF (.NOT. foundEquation) prismEligible = .FALSE.
          bfId = 0
          IF (foundBodyForce) THEN
            IF (Phase24SourceOnlyBodyForce(MetadataBodyForce)) THEN
              bfId = ListGetInteger(CurrentModel % Bodies(MetadataElement % BodyId) % Values, &
                   'Body Force', found)
              IF (.NOT. found .OR. bfId <= 0) prismEligible = .FALSE.
            ELSE
              prismEligible = .FALSE.
            END IF
          END IF
          IF (prismEligible) THEN
            metadataConvection = GetString(MetadataEquation, 'Convection', found)
            IF (found .AND. (metadataConvection == 'computed' .OR. &
                 metadataConvection == 'constant')) prismEligible = .FALSE.
            metadataPhase = GetString(MetadataEquation, 'Phase Change Model', found)
            IF (found .AND. metadataPhase(1:4) /= 'none') prismEligible = .FALSE.
          END IF
        END IF
        IF (prismEligible) THEN
          nPrism = nPrism + 1
          Phase24PrismElements(nPrism) = i
          Phase24ElementClass(i) = 4
          Phase24ElementNodes(:,i) = 0
          Phase24ElementEq(:,i) = 0
          Phase24ElementNodes(1:6,i) = prismNodeIds
          Phase24ElementEq(1:6,i) = prismEqIds
          Phase24StaticK(i) = kappa
          Phase24StaticRhoCp(i) = rho * cp
          DynFlag(i) = dynK
          Phase24ElementSrc(i) = bfId
          CYCLE
        END IF
      END IF

      eligible = ASSOCIATED(MetadataElement % TYPE) .AND. &
           MetadataElement % TYPE % ElementCode == 504 .AND. &
           SIZE(MetadataElement % NodeIndexes) >= 4
      IF (.NOT. eligible) THEN
        nGeneric = nGeneric + 1
        Phase24GenericElements(nGeneric) = i
        CYCLE
      END IF

      nodeIds = MetadataElement % NodeIndexes(1:4)
      eqIds = 0
      IF (ALL(nodeIds > 0) .AND. MAXVAL(nodeIds) <= SIZE(TempPerm)) THEN
        eqIds = TempPerm(nodeIds)
      END IF
      Phase24ElementNodes(:,i) = nodeIds
      Phase24ElementEq(:,i) = eqIds

      CALL GetElementNodes(MetadataNodes, UElement=MetadataElement)
      stat = ElementInfo(MetadataElement, MetadataNodes, 0.25_dp, 0.25_dp, &
           0.25_dp, detJ, Basis4, grad4)
      IF (.NOT. stat .OR. detJ <= 0.0_dp .OR. .NOT. ALL(eqIds > 0)) THEN
        nGeneric = nGeneric + 1
        Phase24GenericElements(nGeneric) = i
        CYCLE
      END IF

      Phase24ElementGrad(:,:,i) = grad4
      Phase24ElementVolume(i) = detJ / 6.0_dp
      MetadataMaterial => GetMaterial(MetadataElement, found)
      IF (.NOT. found .OR. .NOT. ASSOCIATED(MetadataMaterial)) THEN
        nDynamic = nDynamic + 1
        Phase24DynamicElements(nDynamic) = i
        CYCLE
      END IF

      kappa = GetCReal(MetadataMaterial, 'Heat Conductivity', foundK)
      cp = GetCReal(MetadataMaterial, 'Heat Capacity', foundCp)
      rho = GetCReal(MetadataMaterial, 'Density', foundRho)
      staticEligible = foundK .AND. foundCp .AND. foundRho .AND. &
           .NOT. Solver % PeriodicFlipActive
      ! GetCReal also "finds" a Variable/Procedure/MATC property (evaluated at
      ! the current state); caching it would freeze e.g. a k(T) material at its
      ! initial temperature.  Only truly constant properties may be static.
      IF (staticEligible) staticEligible = &
           ListCheckIsConstant(MetadataMaterial, 'Heat Conductivity') .AND. &
           ListCheckIsConstant(MetadataMaterial, 'Heat Capacity') .AND. &
           ListCheckIsConstant(MetadataMaterial, 'Density')
      MetadataEquation => GetEquation(MetadataElement, foundEquation)
      MetadataBodyForce => GetBodyForce(MetadataElement, foundBodyForce)
      IF (.NOT. foundEquation) staticEligible = .FALSE.
      bfId = 0
      IF (foundBodyForce) THEN
        IF (Phase24SourceOnlyBodyForce(MetadataBodyForce)) THEN
          bfId = ListGetInteger(CurrentModel % Bodies(MetadataElement % BodyId) % Values, &
               'Body Force', found)
          IF (.NOT. found .OR. bfId <= 0) staticEligible = .FALSE.
        ELSE
          staticEligible = .FALSE.
        END IF
      END IF
      IF (staticEligible) THEN
        metadataConvection = GetString(MetadataEquation, 'Convection', found)
        IF (found .AND. (metadataConvection == 'computed' .OR. &
             metadataConvection == 'constant')) staticEligible = .FALSE.
        metadataPhase = GetString(MetadataEquation, 'Phase Change Model', found)
        IF (found .AND. metadataPhase(1:4) /= 'none') staticEligible = .FALSE.
      END IF
      IF (staticEligible) THEN
        nStatic = nStatic + 1
        Phase24StaticElements(nStatic) = i
        Phase24ElementClass(i) = 1
        Phase24ElementSrc(i) = bfId
        StaticK = kappa
        StaticRhoCp = rho * cp
        Phase24StaticK(i) = StaticK
        Phase24StaticRhoCp(i) = StaticRhoCp
      ELSE
        nDynamic = nDynamic + 1
        Phase24DynamicElements(nDynamic) = i
        Phase24ElementClass(i) = 2
      END IF
    END DO

    Phase24StaticCount = nStatic
    Phase24DynamicCount = nDynamic
    Phase24GenericCount = nGeneric
    Phase24FastTetCount = nStatic + nDynamic
    Phase24PrismCount = nPrism
    Phase24FastPrismCount = nPrism
    IF (nPrism > 0) THEN
      ALLOCATE(Phase24PrismSlot(nBulk), Phase24PrismNQ(nPrism), &
           Phase24PrismGrad(6,3,PHASE24_PRISM_MAX_Q,nPrism), &
           Phase24PrismWeight(PHASE24_PRISM_MAX_Q,nPrism), &
           Phase24PrismStiff(6,6,nPrism), Phase24PrismMass(6,6,nPrism), &
           Phase24PrismK(nPrism), Phase24PrismRhoCp(nPrism), &
           Phase24PrismVolume(nPrism), &
           Phase24PrismBasis(6,PHASE24_PRISM_MAX_Q,nPrism), &
           Phase24PrismDynamicK(nPrism), Phase24PrismKNodeSlot(6,nPrism))
      Phase24PrismBasis = 0.0_dp
      Phase24PrismDynamicK = .FALSE.
      Phase24PrismKNodeSlot = 0
      Phase24PrismSlot = 0
      Phase24PrismNQ = 0
      Phase24PrismGrad = 0.0_dp
      Phase24PrismWeight = 0.0_dp
      Phase24PrismStiff = 0.0_dp
      Phase24PrismMass = 0.0_dp
      Phase24PrismK = 0.0_dp
      Phase24PrismRhoCp = 0.0_dp
      Phase24PrismVolume = 0.0_dp
      DO i=1,nPrism
        idx = Phase24PrismElements(i)
        Phase24PrismSlot(idx) = i
        Phase24PrismDynamicK(i) = DynFlag(idx)
        Phase24PrismK(i) = Phase24StaticK(idx)
        ! A dynamic-k wedge keeps a zero cached stiffness; it is rebuilt at
        ! every assembly from the current nodal conductivity.
        IF (Phase24PrismDynamicK(i)) Phase24PrismK(i) = 0.0_dp
        Phase24PrismRhoCp(i) = Phase24StaticRhoCp(idx)
        MetadataElement => Solver % Mesh % Elements(idx)
        CALL GetElementNodes(MetadataNodes, UElement=MetadataElement)
        PrismIntegStuff = GaussPoints(MetadataElement)
        Phase24PrismNQ(i) = PrismIntegStuff % n
        DO q=1,PrismIntegStuff % n
          stat = ElementInfo(MetadataElement, MetadataNodes, &
               PrismIntegStuff % u(q), PrismIntegStuff % v(q), &
               PrismIntegStuff % w(q), detJ, Basis6, grad6)
          IF (.NOT. stat .OR. detJ <= 0.0_dp) THEN
            CALL Fatal('HeatSolve','Phase24 prism geometry cache failed')
          END IF
          wq = PrismIntegStuff % s(q) * detJ
          Phase24PrismGrad(:,:,q,i) = grad6
          Phase24PrismBasis(:,q,i) = Basis6
          Phase24PrismWeight(q,i) = wq
          Phase24PrismVolume(i) = Phase24PrismVolume(i) + wq
          DO j=1,6
            DO k=1,6
              Phase24PrismStiff(j,k,i) = Phase24PrismStiff(j,k,i) + &
                   wq * Phase24PrismK(i) * SUM(grad6(j,:) * grad6(k,:))
              Phase24PrismMass(j,k,i) = Phase24PrismMass(j,k,i) + &
                   wq * Phase24PrismRhoCp(i) * Basis6(j) * Basis6(k)
            END DO
          END DO
        END DO
      END DO
      ! Unique node list of the dynamic-k wedges: their conductivity is
      ! evaluated once per assembly with a single ListGetReal call.
      Phase24DynamicPrismCount = COUNT(Phase24PrismDynamicK)
      IF (Phase24DynamicPrismCount > 0) THEN
        ! Distinct temperature-dependent materials, in order of appearance.
        ALLOCATE(Phase24DynKMats(CurrentModel % NumberOfMaterials), &
             Phase24DynKStart(CurrentModel % NumberOfMaterials), &
             Phase24DynKCount(CurrentModel % NumberOfMaterials))
        Phase24DynKMatCount = 0
        DO i=1,nPrism
          IF (.NOT. Phase24PrismDynamicK(i)) CYCLE
          k = Phase24ElementDynMat(Phase24PrismElements(i))
          IF (.NOT. ANY(Phase24DynKMats(1:Phase24DynKMatCount) == k)) THEN
            Phase24DynKMatCount = Phase24DynKMatCount + 1
            Phase24DynKMats(Phase24DynKMatCount) = k
          END IF
        END DO
        ! One slot per (material, node); slots of a material are contiguous.
        ALLOCATE(NodeMark(Solver % Mesh % NumberOfNodes))
        nDynNodes = 0
        DO m=1,Phase24DynKMatCount
          NodeMark = 0
          Phase24DynKStart(m) = nDynNodes + 1
          DO i=1,nPrism
            IF (.NOT. Phase24PrismDynamicK(i)) CYCLE
            idx = Phase24PrismElements(i)
            IF (Phase24ElementDynMat(idx) /= Phase24DynKMats(m)) CYCLE
            DO j=1,6
              k = Solver % Mesh % Elements(idx) % NodeIndexes(j)
              IF (NodeMark(k) == 0) THEN
                nDynNodes = nDynNodes + 1
                NodeMark(k) = nDynNodes
              END IF
              Phase24PrismKNodeSlot(j,i) = NodeMark(k)
            END DO
          END DO
          Phase24DynKCount(m) = nDynNodes - Phase24DynKStart(m) + 1
        END DO
        ALLOCATE(Phase24DynKNodes(nDynNodes), Phase24DynKNodal(nDynNodes))
        DO i=1,nPrism
          IF (.NOT. Phase24PrismDynamicK(i)) CYCLE
          idx = Phase24PrismElements(i)
          DO j=1,6
            Phase24DynKNodes(Phase24PrismKNodeSlot(j,i)) = Solver % Mesh % Elements(idx) % NodeIndexes(j)
          END DO
        END DO
        Phase24DynKNodal = 0.0_dp
        Phase24DynKNodeCount = nDynNodes
        DEALLOCATE(NodeMark)
        WRITE(Message,'(A,I0,A,I0,A,I0)') 'Phase24 dynamic-k prisms: elements=', &
             Phase24DynamicPrismCount, ', node slots=', nDynNodes, ', materials=', Phase24DynKMatCount
        CALL Info('HeatSolve',TRIM(Message),Level=4)
      END IF
    END IF
    ! Compact list of source-only body forces used by fast elements, and the
    ! unique node set of each (for one ListGetReal per body force).
    nSrc = 0
    DO i=1,nBulk
      bfId = Phase24ElementSrc(i)
      IF (bfId <= 0) CYCLE
      IF (Phase24ElementClass(i) /= 1 .AND. Phase24ElementClass(i) /= 4) THEN
        Phase24ElementSrc(i) = 0
        CYCLE
      END IF
      k = 0
      DO j=1,nSrc
        IF (srcBF(j) == bfId) k = j
      END DO
      IF (k == 0) THEN
        nSrc = nSrc + 1
        IF (nSrc > SIZE(srcBF)) CALL Fatal('HeatSolve','Phase24: too many source body forces')
        srcBF(nSrc) = bfId
        k = nSrc
      END IF
      Phase24ElementSrc(i) = k
    END DO
    Phase24SrcCount = nSrc
    IF (nSrc > 0) THEN
      ALLOCATE(Phase24SrcBF(nSrc), Phase24SrcNodeStart(nSrc), Phase24SrcNodeCount(nSrc), &
           NodeMark(Solver % Mesh % NumberOfNodes))
      Phase24SrcBF = srcBF(1:nSrc)
      srcTotal = 0
      DO k=1,nSrc
        NodeMark = 0
        Phase24SrcNodeStart(k) = srcTotal + 1
        Phase24SrcNodeCount(k) = 0
        DO i=1,nBulk
          IF (Phase24ElementSrc(i) /= k) CYCLE
          DO j=1,Solver % Mesh % Elements(i) % TYPE % NumberOfNodes
            q = Solver % Mesh % Elements(i) % NodeIndexes(j)
            IF (NodeMark(q) == 0) THEN
              NodeMark(q) = 1
              Phase24SrcNodeCount(k) = Phase24SrcNodeCount(k) + 1
            END IF
          END DO
        END DO
        srcTotal = srcTotal + Phase24SrcNodeCount(k)
      END DO
      ALLOCATE(Phase24SrcNodeList(srcTotal), &
           Phase24SrcLoad(Solver % Mesh % NumberOfNodes, nSrc))
      Phase24SrcLoad = 0.0_dp
      DO k=1,nSrc
        NodeMark = 0
        j = Phase24SrcNodeStart(k) - 1
        DO i=1,nBulk
          IF (Phase24ElementSrc(i) /= k) CYCLE
          DO idx=1,Solver % Mesh % Elements(i) % TYPE % NumberOfNodes
            q = Solver % Mesh % Elements(i) % NodeIndexes(idx)
            IF (NodeMark(q) == 0) THEN
              NodeMark(q) = 1
              j = j + 1
              Phase24SrcNodeList(j) = q
            END IF
          END DO
        END DO
        WRITE(Message,'(A,I0,A,I0,A,I0)') 'Phase24 fast source: body force=', &
             Phase24SrcBF(k), ', elements=', COUNT(Phase24ElementSrc == k), &
             ', nodes=', Phase24SrcNodeCount(k)
        CALL Info('HeatSolve',TRIM(Message),Level=4)
      END DO
      DEALLOCATE(NodeMark)
    END IF
    Phase24GenericPrismCount = 0
    DO i=1,nGeneric
      idx = Phase24GenericElements(i)
      IF (ASSOCIATED(Solver % Mesh % Elements(idx) % TYPE)) THEN
        IF (Solver % Mesh % Elements(idx) % TYPE % ElementCode == 706) &
             Phase24GenericPrismCount = Phase24GenericPrismCount + 1
      END IF
    END DO
    Phase24GenericOtherCount = nGeneric - Phase24GenericPrismCount

    ! Geometry is now fixed.  Generate reusable element matrices in parallel;
    ! this is intentionally arithmetic-only and has no sparse or ValueList
    ! side effects, so it is safe for OpenMP and suitable for later hot-loop
    ! parallel assembly.
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(Phase24StaticElements,Phase24StaticCount,Phase24ElementGrad, &
    !$OMP& Phase24ElementVolume,Phase24StaticK,Phase24StaticRhoCp, &
    !$OMP& Phase24StaticStiff,Phase24SpatialMass) &
    !$OMP PRIVATE(i,j,k,idx)
    DO i=1,Phase24StaticCount
      idx = Phase24StaticElements(i)
      DO j=1,4
        DO k=1,4
          Phase24StaticStiff(j,k,idx) = Phase24StaticK(idx) * Phase24ElementVolume(idx) * &
               SUM(Phase24ElementGrad(j,:,idx) * Phase24ElementGrad(k,:,idx))
          IF (j == k) THEN
            Phase24SpatialMass(j,k,idx) = Phase24StaticRhoCp(idx) * &
                 Phase24ElementVolume(idx) / 10.0_dp
          ELSE
            Phase24SpatialMass(j,k,idx) = Phase24StaticRhoCp(idx) * &
                 Phase24ElementVolume(idx) / 20.0_dp
          END IF
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

    Phase24MetadataElements = nBulk
    Phase24MetadataReady = .TRUE.
    WRITE(Message,'(A,3(I0,A))') 'Phase24 metadata: static=',nStatic, &
         ', dynamic=',nDynamic, ', generic=',nGeneric
    CALL Info('HeatSolve',TRIM(Message),Level=4)
    fastPercent = 100.0_dp * REAL(Phase24FastTetCount+Phase24FastPrismCount,dp) / &
         MAX(1.0_dp,REAL(nBulk,dp))
    WRITE(Message,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F8.3,A)') &
         'Phase24 element accounting: total=',nBulk, &
         ', fast_tet=',Phase24FastTetCount, &
         ', fast_prism=',Phase24FastPrismCount, &
         ', generic_prism=',Phase24GenericPrismCount, &
         ', generic_other=',Phase24GenericOtherCount, &
         ', accounted=',Phase24FastTetCount+Phase24FastPrismCount+Phase24GenericCount, &
         ', fast_percent=',fastPercent,'%' 
    CALL Info('HeatSolve',TRIM(Message),Level=4)
  END SUBROUTINE Phase24BuildMetadata

  SUBROUTINE Phase24AssembleStaticElements(Active,RHSOnly)
    LOGICAL, INTENT(OUT) :: Active
    LOGICAL, INTENT(IN), OPTIONAL :: RHSOnly
    INTEGER :: i,j,k,idx,eq,order,srcK
    REAL(KIND=dp) :: LocalMass(4,4),LocalStiff(4,4),LocalForce(4)
    REAL(KIND=dp) :: load4(4),wsrc
    REAL(KIND=dp) :: timeScale,sourceTerm,historyScale1,historyScale2,previousDt
    TYPE(Variable_t), POINTER :: DtVar
    LOGICAL :: useBdf,OnlyRHS

    OnlyRHS = .FALSE.
    IF (PRESENT(RHSOnly)) OnlyRHS = RHSOnly
    Active = .FALSE.
    IF (Phase24StaticCount <= 0) RETURN
    IF (Solver % Matrix % Format /= MATRIX_CRS) RETURN
    IF (Solver % PeriodicFlipActive) RETURN
    IF (Solver % LocalSystemMode > 0) RETURN
    IF (ConstantBulk .OR. Solver % NOFEigenValues > 0) RETURN

    useBdf = .FALSE.
    IF (TransientAssembly) THEN
      IF (GetString(Solver % Values, 'Timestepping Method', Found) == 'bdf') useBdf = .TRUE.
      IF (.NOT. useBdf) RETURN
    END IF

    IF (.NOT. OnlyRHS) CALL CRS_PreparePositionCache(Solver % Matrix, Phase24MetadataElements)
    Active = .TRUE.
    WRITE(Message,'(A,I0,A,I0)') 'Phase24 fast static assembly entered: elements=', &
         Phase24StaticCount, ', order=', MIN(Solver % DoneTime, Solver % Order)
    CALL Info('HeatSolve',TRIM(Message),Level=4)
    order = MIN(Solver % DoneTime, Solver % Order)
    CALL Phase24GetBDFCoefficients(order,Solver % dt,previousDt,timeScale, &
         historyScale1,historyScale2)
    historyScale1 = -historyScale1
    historyScale2 = -historyScale2
    IF (TransientAssembly .AND. order > 5) RETURN

    ! The matrix cache is fixed before the loop.  Each element writes a
    ! private local contribution; cached CSR values and RHS entries use
    ! atomics, so no global critical region surrounds the loop.
    !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(static) &
    !$OMP SHARED(Phase24StaticElements,Phase24StaticCount,Phase24ElementEq, &
    !$OMP& Phase24StaticStiff,Phase24SpatialMass,Phase24StaticRhoCp, &
    !$OMP& Phase24ElementVolume,Phase24ElementPart,ParEnv, &
    !$OMP& Phase24ElementSrc,Phase24SrcLoad,Phase24ElementNodes, &
         !$OMP& Phase24MetadataElements,TransientAssembly,order,timeScale,historyScale1,historyScale2,Solver,OnlyRHS) &
    !$OMP PRIVATE(i,j,k,idx,eq,LocalMass,LocalStiff,LocalForce,sourceTerm,srcK,load4,wsrc)
    DO i=1,Phase24StaticCount
      idx = Phase24StaticElements(i)
      IF (ParEnv % PEs > 1 .AND. Phase24ElementPart(idx) /= ParEnv % MyPE) CYCLE
      IF (Solver % Matrix % Lumped) THEN
        LocalMass = 0.0_dp
        DO j=1,4
          LocalMass(j,j) = Phase24StaticRhoCp(idx) * &
               Phase24ElementVolume(idx) / 4.0_dp
        END DO
      ELSE
        LocalMass = Phase24SpatialMass(:,:,idx)
      END IF
      LocalStiff = Phase24StaticStiff(:,:,idx)
      LocalForce = 0.0_dp
      IF (TransientAssembly .AND. order > 0) THEN
        DO j=1,4
          sourceTerm = 0.0_dp
          DO k=1,4
            SELECT CASE(order)
            CASE(1)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   historyScale1 * Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1)
            CASE(2)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (historyScale1*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) + &
                   historyScale2*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2))
            CASE(3)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (3.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   1.5_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (1.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3))
            CASE(4)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (4.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   3.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (4.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3) - &
                   0.25_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),4))
            CASE(5)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (5.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   5.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (10.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3) - &
                   1.25_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),4) + &
                   0.2_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),5))
            END SELECT
          END DO
          IF (order <= 2) THEN
            LocalForce(j) = sourceTerm
          ELSE
            LocalForce(j) = sourceTerm * timeScale
          END IF
        END DO
        LocalStiff = LocalStiff + timeScale * LocalMass
      END IF
      ! Volumetric source: consistent P1 tetra integral, V/20*(sum(q)+q_j).
      srcK = Phase24ElementSrc(idx)
      IF (srcK > 0) THEN
        load4 = Phase24SrcLoad(Phase24ElementNodes(1:4,idx),srcK)
        wsrc = Phase24ElementVolume(idx) / 20.0_dp
        DO j=1,4
          LocalForce(j) = LocalForce(j) + wsrc * (SUM(load4) + load4(j))
        END DO
      END IF
      IF (.NOT. OnlyRHS) THEN
        CALL CRS_GlueLocalMatrixCached(Solver % Matrix,4,1, &
             Phase24ElementEq(:,idx),LocalStiff,idx,.FALSE.,.FALSE.)
      END IF
      DO j=1,4
        eq = Phase24ElementEq(j,idx)
        IF (eq > 0) THEN
          !$OMP ATOMIC
          Solver % Matrix % RHS(eq) = Solver % Matrix % RHS(eq) + LocalForce(j)
        END IF
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE Phase24AssembleStaticElements

  FUNCTION TESMeanTemperatureNow(kc) RESULT(Tmean)
    INTEGER, INTENT(IN) :: kc
    ! Element-averaged TES temperature, the same definition as the circuit.
    REAL(KIND=dp) :: Tmean
    INTEGER :: t, i, n, p, body, cnt
    LOGICAL :: gotBody
    REAL(KIND=dp) :: acc
    TYPE(Element_t), POINTER :: El

    body = TESC(kc) % Body
    acc = 0.0_dp
    cnt = 0
    DO t = 1, Solver % NumberOfActiveElements
      El => GetActiveElement(t)
      IF (ParEnv % PEs > 1 .AND. El % PartIndex /= ParEnv % MyPE) CYCLE
      IF (El % BodyId /= body) CYCLE
      n = GetElementNOFNodes(El)
      IF (n <= 0) CYCLE
      DO i = 1, n
        p = TempPerm(El % NodeIndexes(i))
        IF (p > 0) acc = acc + Temperature(p) / REAL(n, dp)
      END DO
      cnt = cnt + 1
    END DO
    acc = ParallelReduction(acc)
    Tmean = acc / MAX(ParallelReduction(REAL(cnt, dp)), 1.0_dp)
  END FUNCTION TESMeanTemperatureNow

  ! Error-controlled step size for variable-step BDF2 (no rejections).  The
  ! local error of the converged TES temperature is estimated from the
  ! quadratic predictor through the three previous steps:
  !   p(t_n+1) - x(t_n+1) = -x'''/6 h(h+h1)(h+h1+h2)
  !   LTE_BDF2 ~ (2/9) h^3 x'''   =>  LTE ~ 4/3 h^2 |x - p| / ((h+h1)(h+h1+h2))
  ! The next step is dt*0.9*(tol/LTE)^(1/3), limited to [0.2, 2] x dt (BDF2
  ! zero stability needs ratios below 1+sqrt(2)).  ElmerSolver uses the value
  ! in the last timestep interval only ('Phase24 Step Controller'), so the
  ! pulse and its nanosecond ramp keep their fixed steps.
  SUBROUTINE Phase24StepController()
    REAL(KIND=dp), ALLOCATABLE, SAVE :: HistT(:,:)
    REAL(KIND=dp), SAVE :: HistTime(4) = 0.0_dp
    INTEGER, SAVE :: HistN = 0
    INTEGER :: kc
    REAL(KIND=dp) :: tnow, tol, errk, tEnd, maxdt, h, h1, h2, p, err, fac, dtnext, remaining
    REAL(KIND=dp) :: t0, t1, t2, L0, L1, L2
    LOGICAL :: gotEnd, gotMax

    tnow = GetTime()
    IF (.NOT. ALLOCATED(TESC)) RETURN
    IF (.NOT. ALLOCATED(HistT)) THEN
      ALLOCATE(HistT(4,TESNCirc)); HistT = 0.0_dp
    END IF
    ! One TES temperature history per circuit; the worst circuit sets dt.
    HistT(1:3,:) = HistT(2:4,:); HistTime(1:3) = HistTime(2:4)
    DO kc = 1, TESNCirc
      HistT(4,kc) = TESMeanTemperatureNow(kc)
    END DO
    HistTime(4) = tnow
    HistN = MIN(HistN + 1, 4)

    tol = ListGetCReal(Model % Simulation,'Phase24 Step Controller Tolerance',Found)
    IF (.NOT. Found .OR. tol <= 0.0_dp) tol = 2.0e-8_dp
    tEnd = ListGetCReal(Model % Simulation,'Phase24 Step Controller End Time',gotEnd)
    maxdt = ListGetCReal(Model % Simulation,'Phase24 Step Controller Max Timestep',gotMax)

    dtnext = Solver % dt
    err = -1.0_dp
    IF (HistN == 4) THEN
      h = HistTime(4) - HistTime(3)
      h1 = HistTime(3) - HistTime(2)
      h2 = HistTime(2) - HistTime(1)
      IF (h > 0.0_dp .AND. h1 > 0.0_dp .AND. h2 > 0.0_dp) THEN
        t0 = HistTime(1); t1 = HistTime(2); t2 = HistTime(3)
        L0 = (tnow-t1)*(tnow-t2)/((t0-t1)*(t0-t2))
        L1 = (tnow-t0)*(tnow-t2)/((t1-t0)*(t1-t2))
        L2 = (tnow-t0)*(tnow-t1)/((t2-t0)*(t2-t1))
        err = 0.0_dp
        DO kc = 1, TESNCirc
          p = L0*HistT(1,kc) + L1*HistT(2,kc) + L2*HistT(3,kc)
          errk = (4.0_dp/3.0_dp) * h*h * ABS(HistT(4,kc) - p) / ((h+h1)*(h+h1+h2))
          err = MAX(err, errk)
        END DO
        fac = 0.9_dp * (tol / MAX(err, 1.0e-300_dp))**(1.0_dp/3.0_dp)
        ! Hysteresis: keep the step unless it must shrink or can grow by a
        ! useful amount, so equal-dt runs keep the circuit Newton slope and the
        ! step history smooth; growth is limited to 1.5x per step.
        IF (fac >= 1.0_dp .AND. fac < 1.2_dp) fac = 1.0_dp
        fac = MIN(MAX(fac, 0.2_dp), 1.5_dp)
        dtnext = Solver % dt * fac
      END IF
    END IF
    IF (gotMax .AND. maxdt > 0.0_dp) dtnext = MIN(dtnext, maxdt)
    IF (gotEnd) THEN
      remaining = tEnd - tnow
      IF (remaining <= 1.0e-9_dp * Solver % dt) THEN
        CALL ListAddConstReal(Model % Simulation,'Exit Condition',1.0_dp)
      ELSE IF (remaining <= dtnext) THEN
        dtnext = remaining
      ELSE IF (remaining < 1.5_dp * dtnext) THEN
        ! Avoid a tiny final step: split the rest in two.
        dtnext = 0.5_dp * remaining
      END IF
    END IF
    CALL ListAddConstReal(Model % Simulation,'Phase24 Controller Next Timestep',dtnext)
    WRITE(Message,'(A,ES12.5,A,ES10.3,A,ES10.3,A,ES10.3)') 'Phase24 step controller: t=',tnow, &
         ', dt=',Solver % dt,', lte=',err,', next dt=',dtnext
    CALL Info('HeatSolve',TRIM(Message),Level=4)
  END SUBROUTINE Phase24StepController

  ! Element-averaged TES temperature of the converged previous step (history
  ! slot 1).  A predictor may already have overwritten the current field.
  FUNCTION TESMeanPreviousStep(kc) RESULT(Tmean)
    INTEGER, INTENT(IN) :: kc
    REAL(KIND=dp) :: Tmean
    INTEGER :: t, i, n, p, body, cnt
    LOGICAL :: gotBody
    REAL(KIND=dp) :: acc
    TYPE(Element_t), POINTER :: El

    body = TESC(kc) % Body
    acc = 0.0_dp
    cnt = 0
    DO t = 1, Solver % NumberOfActiveElements
      El => GetActiveElement(t)
      IF (ParEnv % PEs > 1 .AND. El % PartIndex /= ParEnv % MyPE) CYCLE
      IF (El % BodyId /= body) CYCLE
      n = GetElementNOFNodes(El)
      IF (n <= 0) CYCLE
      DO i = 1, n
        p = TempPerm(El % NodeIndexes(i))
        IF (p > 0) acc = acc + Solver % Variable % PrevValues(p,1) / REAL(n, dp)
      END DO
      cnt = cnt + 1
    END DO
    acc = ParallelReduction(acc)
    Tmean = acc / MAX(ParallelReduction(REAL(cnt, dp)), 1.0_dp)
  END FUNCTION TESMeanPreviousStep

  ! Joule power of the implicit-Euler circuit step at TES temperature T, with
  ! the coefficients of the latest TESInnerCircuitUpdate call (same law).
  FUNCTION TESCircuitPowerAt(kc, T) RESULT(P)
    INTEGER, INTENT(IN) :: kc
    REAL(KIND=dp), INTENT(IN) :: T
    REAL(KIND=dp) :: P, Ac, Bc, Cc, Dc, Ic, Rc, Ld
    Ld = 0.0_dp
    IF (TESC(kc) % ResDt > 0.0_dp) Ld = TESC(kc) % ResL / TESC(kc) % ResDt
    Ac = TESC(kc) % ResR0 * (1.0_dp + TESC(kc) % ResAlpha * (T - TESC(kc) % ResT0) / TESC(kc) % ResT0 - TESC(kc) % ResBeta)
    Bc = TESC(kc) % ResR0 * TESC(kc) % ResBeta / TESC(kc) % ResI0
    Cc = TESC(kc) % ResRSh + Ac + Ld
    Dc = MAX(Cc*Cc + 4.0_dp*Bc*(TESC(kc) % ResIBias*TESC(kc) % ResRSh + Ld*TESC(kc) % ResPrevI), 0.0_dp)
    Ic = MAX(MIN((SQRT(Dc) - Cc) / (2.0_dp*Bc), TESC(kc) % ResIBias), 0.0_dp)
    Rc = Ac + Bc*ABS(Ic)
    IF (Rc < TESC(kc) % ResRMin) THEN
      Rc = TESC(kc) % ResRMin
      Ic = (TESC(kc) % ResIBias*TESC(kc) % ResRSh + Ld*TESC(kc) % ResPrevI) / &
           (TESC(kc) % ResRSh + Rc + Ld)
    END IF
    P = Ic*Ic*Rc
  END FUNCTION TESCircuitPowerAt

  FUNCTION Phase24SourceOnlyBodyForce(List) RESULT(SourceOnly)
    ! True when the body force carries nothing but a volumetric heat source
    ! (and its name), i.e. it touches the RHS but never the operator.
    TYPE(ValueList_t), POINTER :: List
    LOGICAL :: SourceOnly
    TYPE(ValueListEntry_t), POINTER :: entry

    SourceOnly = .FALSE.
    IF (.NOT. ASSOCIATED(List)) RETURN
    entry => List % Head
    DO WHILE (ASSOCIATED(entry))
      SELECT CASE (entry % Name)
      CASE ('name')
      CASE ('volumetric heat source')
        SourceOnly = .TRUE.
      CASE DEFAULT
        SourceOnly = .FALSE.
        RETURN
      END SELECT
      entry => entry % Next
    END DO
  END FUNCTION Phase24SourceOnlyBodyForce

  SUBROUTINE Phase24EvaluateSources()
    ! Nodal 'Volumetric Heat Source' of every fast source body force, one
    ! serial ListGetReal call each (UDFs/ValueLists are not thread-safe).
    INTEGER :: k,a,b
    LOGICAL :: gotSrc
    REAL(KIND=dp), ALLOCATABLE :: vals(:)

    DO k=1,Phase24SrcCount
      a = Phase24SrcNodeStart(k)
      b = a + Phase24SrcNodeCount(k) - 1
      IF (b < a) CYCLE
      ALLOCATE(vals(b-a+1))
      vals = ListGetReal(CurrentModel % BodyForces(Phase24SrcBF(k)) % Values, &
           'Volumetric Heat Source', b-a+1, Phase24SrcNodeList(a:b), gotSrc)
      IF (.NOT. gotSrc) vals = 0.0_dp
      Phase24SrcLoad(Phase24SrcNodeList(a:b),k) = vals
      DEALLOCATE(vals)
    END DO
  END SUBROUTINE Phase24EvaluateSources

  ! Nodal conductivity of every temperature-dependent material at its slots
  ! (serial: ValueList/UDF calls are not thread-safe), one ListGetReal per
  ! material with the current Temperature.
  SUBROUTINE Phase24EvaluateDynK(kNodal)
    REAL(KIND=dp), INTENT(OUT) :: kNodal(:)
    INTEGER :: m, a, c
    LOGICAL :: gotK
    DO m = 1, Phase24DynKMatCount
      a = Phase24DynKStart(m); c = Phase24DynKCount(m)
      IF (c <= 0) CYCLE
      kNodal(a:a+c-1) = ListGetReal(CurrentModel % Materials(Phase24DynKMats(m)) % Values, &
           'Heat Conductivity', c, Phase24DynKNodes(a:a+c-1), gotK)
      IF (.NOT. gotK) CALL Fatal('HeatSolve','Phase24 dynamic-k conductivity not found')
    END DO
  END SUBROUTINE Phase24EvaluateDynK

  SUBROUTINE Phase24AssemblePrismElements(Active,RHSOnly)
    LOGICAL, INTENT(OUT) :: Active
    LOGICAL, INTENT(IN), OPTIONAL :: RHSOnly
    INTEGER :: i,j,k,idx,slot,eq,order,q
    REAL(KIND=dp) :: LocalMass(6,6),LocalStiff(6,6),LocalForce(6)
    REAL(KIND=dp) :: nodalK(6),kq,load6(6),lq,LumpDelta(6,6)
    INTEGER :: srcK
    LOGICAL :: foundDynK,doLump
    REAL(KIND=dp) :: timeScale,sourceTerm,totalMass,historyScale1,historyScale2,previousDt
    TYPE(Variable_t), POINTER :: DtVar
    LOGICAL :: useBdf,OnlyRHS

    OnlyRHS = .FALSE.
    IF (PRESENT(RHSOnly)) OnlyRHS = RHSOnly
    Active = .FALSE.
    IF (Phase24PrismCount <= 0) RETURN
    IF (Solver % Matrix % Format /= MATRIX_CRS) RETURN
    IF (Solver % PeriodicFlipActive) RETURN
    IF (Solver % LocalSystemMode > 0) RETURN
    IF (ConstantBulk .OR. Solver % NOFEigenValues > 0) RETURN

    useBdf = .FALSE.
    IF (TransientAssembly) THEN
      IF (GetString(Solver % Values, 'Timestepping Method', Found) == 'bdf') useBdf = .TRUE.
      IF (.NOT. useBdf) RETURN
    END IF

    order = MIN(Solver % DoneTime, Solver % Order)
    IF (TransientAssembly .AND. order > 5) RETURN
    IF (.NOT. OnlyRHS) CALL CRS_PreparePositionCache(Solver % Matrix, Phase24MetadataElements)
    Active = .TRUE.
    WRITE(Message,'(A,I0,A,I0)') 'Phase24 fast prism assembly entered: elements=', &
         Phase24PrismCount, ', order=', order
    CALL Info('HeatSolve',TRIM(Message),Level=4)

    CALL Phase24GetBDFCoefficients(order,Solver % dt,previousDt,timeScale, &
         historyScale1,historyScale2)
    historyScale1 = -historyScale1
    historyScale2 = -historyScale2

    ! On an RHS-only pass the matrix is unchanged, so the delta is kept.
    doLump = Phase24LumpedPrec .AND. .NOT. OnlyRHS
    IF (doLump) THEN
      IF (ALLOCATED(Phase24PrecDelta)) THEN
        IF (SIZE(Phase24PrecDelta) /= SIZE(Solver % Matrix % Values)) DEALLOCATE(Phase24PrecDelta)
      END IF
      IF (.NOT. ALLOCATED(Phase24PrecDelta)) ALLOCATE(Phase24PrecDelta(SIZE(Solver % Matrix % Values)))
      Phase24PrecDelta = 0.0_dp
    END IF

    ! Current nodal conductivity of the dynamic-k wedges, evaluated serially
    ! (ValueList/UDF calls are not thread-safe) in one call for all nodes.
    IF (Phase24DynamicPrismCount > 0 .AND. .NOT. OnlyRHS) THEN
      CALL Phase24EvaluateDynK(Phase24DynKNodal)
    END IF

    !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(static) &
    !$OMP SHARED(Phase24PrismElements,Phase24PrismCount,Phase24PrismSlot, &
    !$OMP& Phase24ElementEq,Phase24PrismMass,Phase24PrismStiff, &
    !$OMP& Phase24PrismRhoCp,Phase24PrismVolume,Phase24ElementPart,ParEnv, &
    !$OMP& Phase24PrismDynamicK,Phase24PrismKNodeSlot,Phase24DynKNodal, &
    !$OMP& Phase24PrismNQ,Phase24PrismBasis,Phase24PrismGrad,Phase24PrismWeight, &
    !$OMP& Phase24ElementSrc,Phase24SrcLoad,Phase24ElementNodes, &
         !$OMP& Phase24MetadataElements,TransientAssembly,order,timeScale,historyScale1,historyScale2,Solver,OnlyRHS, &
    !$OMP& doLump,Phase24PrecDelta) &
    !$OMP PRIVATE(i,j,k,q,idx,slot,eq,LocalMass,LocalStiff,LocalForce,sourceTerm,totalMass,nodalK,kq, &
    !$OMP& srcK,load6,lq,LumpDelta)
    DO i=1,Phase24PrismCount
      idx = Phase24PrismElements(i)
      IF (ParEnv % PEs > 1 .AND. Phase24ElementPart(idx) /= ParEnv % MyPE) CYCLE
      slot = Phase24PrismSlot(idx)
      IF (Solver % Matrix % Lumped) THEN
        LocalMass = 0.0_dp
        totalMass = Phase24PrismRhoCp(slot) * Phase24PrismVolume(slot)
        DO j=1,6
          LocalMass(j,j) = totalMass / 6.0_dp
        END DO
      ELSE
        LocalMass = Phase24PrismMass(:,:,slot)
      END IF
      LocalStiff = Phase24PrismStiff(:,:,slot)
      IF (Phase24PrismDynamicK(slot) .AND. .NOT. OnlyRHS) THEN
        ! Generic-path form: nodal k interpolated with the basis at each
        ! quadrature point.
        nodalK = Phase24DynKNodal(Phase24PrismKNodeSlot(:,slot))
        LocalStiff = 0.0_dp
        DO q=1,Phase24PrismNQ(slot)
          kq = SUM(Phase24PrismBasis(:,q,slot) * nodalK)
          DO j=1,6
            DO k=1,6
              LocalStiff(j,k) = LocalStiff(j,k) + Phase24PrismWeight(q,slot) * kq * &
                   SUM(Phase24PrismGrad(j,:,q,slot) * Phase24PrismGrad(k,:,q,slot))
            END DO
          END DO
        END DO
      END IF
      IF (doLump) CALL Phase24WedgeLumpDelta(LocalStiff, LumpDelta)
      LocalForce = 0.0_dp
      IF (TransientAssembly .AND. order > 0) THEN
        DO j=1,6
          sourceTerm = 0.0_dp
          DO k=1,6
            SELECT CASE(order)
            CASE(1)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   historyScale1 * Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1)
            CASE(2)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (historyScale1*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) + &
                   historyScale2*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2))
            CASE(3)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (3.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   1.5_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (1.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3))
            CASE(4)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (4.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   3.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (4.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3) - &
                   0.25_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),4))
            CASE(5)
              sourceTerm = sourceTerm + LocalMass(j,k) * &
                   (5.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),1) - &
                   5.0_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),2) + &
                   (10.0_dp/3.0_dp)*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),3) - &
                   1.25_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),4) + &
                   0.2_dp*Solver % Variable % PrevValues(Phase24ElementEq(k,idx),5))
            END SELECT
          END DO
          IF (order <= 2) THEN
            LocalForce(j) = sourceTerm
          ELSE
            LocalForce(j) = sourceTerm * timeScale
          END IF
        END DO
        LocalStiff = LocalStiff + timeScale * LocalMass
      END IF
      ! Volumetric source with the cached wedge quadrature (generic form).
      srcK = Phase24ElementSrc(idx)
      IF (srcK > 0) THEN
        load6 = Phase24SrcLoad(Phase24ElementNodes(1:6,idx),srcK)
        DO q=1,Phase24PrismNQ(slot)
          lq = SUM(Phase24PrismBasis(:,q,slot) * load6) * Phase24PrismWeight(q,slot)
          DO j=1,6
            LocalForce(j) = LocalForce(j) + lq * Phase24PrismBasis(j,q,slot)
          END DO
        END DO
      END IF
      IF (.NOT. OnlyRHS) THEN
        CALL CRS_GlueLocalMatrixCached(Solver % Matrix,6,1, &
             Phase24ElementEq(:,idx),LocalStiff,idx,.FALSE.,.FALSE.)
        IF (doLump) CALL CRS_AddLocalMatrixCachedTo(6,idx,Phase24ElementEq(:,idx), &
             LumpDelta,Phase24PrecDelta)
      END IF
      DO j=1,6
        eq = Phase24ElementEq(j,idx)
        IF (eq > 0) THEN
          !$OMP ATOMIC
          Solver % Matrix % RHS(eq) = Solver % Matrix % RHS(eq) + LocalForce(j)
        END IF
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE Phase24AssemblePrismElements

  ! Lumped-minus-consistent stiffness of a wedge whose nodes j and j+3 are
  ! the vertical pairs.  For a right wedge K = Kt(x)Mz + Mt(x)Kz; the lumped
  ! form replaces Mz by diag(h/2) and Mt by diag(A/3), which keeps row sums
  ! and removes the positive mass-like couplings of thin layers.
  SUBROUTINE Phase24WedgeLumpDelta(K, Delta)
    REAL(KIND=dp), INTENT(IN) :: K(6,6)
    REAL(KIND=dp), INTENT(OUT) :: Delta(6,6)
    REAL(KIND=dp) :: Sb(3,3), St(3,3), m(3)
    INTEGER :: a
    Sb = K(1:3,1:3) + K(1:3,4:6)
    St = K(4:6,4:6) + K(4:6,1:3)
    Sb = 0.5_dp * (Sb + TRANSPOSE(Sb))
    St = 0.5_dp * (St + TRANSPOSE(St))
    m = (SUM(K(1:3,1:3) - 2.0_dp*K(1:3,4:6), DIM=2) + &
         SUM(K(4:6,4:6) - 2.0_dp*K(4:6,1:3), DIM=2)) / 6.0_dp
    Delta = 0.0_dp
    Delta(1:3,1:3) = Sb
    Delta(4:6,4:6) = St
    DO a=1,3
      Delta(a,a) = Delta(a,a) + m(a)
      Delta(a+3,a+3) = Delta(a+3,a+3) + m(a)
      Delta(a,a+3) = -m(a)
      Delta(a+3,a) = -m(a)
    END DO
    Delta = Delta - K
  END SUBROUTINE Phase24WedgeLumpDelta

  ! PrecValues = Values + lumping delta, with Dirichlet rows/columns taken
  ! from the eliminated matrix unchanged.
  SUBROUTINE Phase24BuildLumpedPrec()
    TYPE(Matrix_t), POINTER :: Am
    LOGICAL, ALLOCATABLE :: isD(:)
    INTEGER :: i, j, nr
    Am => Solver % Matrix
    IF (.NOT. ALLOCATED(Phase24PrecDelta)) RETURN
    IF (SIZE(Phase24PrecDelta) /= SIZE(Am % Values)) RETURN
    nr = Am % NumberOfRows
    ALLOCATE(isD(nr))
    IF (ALLOCATED(Am % ConstrainedDOF)) THEN
      isD = Am % ConstrainedDOF(1:nr)
    ELSE
      DO i=1,nr
        isD(i) = .TRUE.
        DO j=Am % Rows(i),Am % Rows(i+1)-1
          IF (Am % Cols(j) /= i .AND. Am % Values(j) /= 0.0_dp) THEN
            isD(i) = .FALSE.
            EXIT
          END IF
        END DO
      END DO
    END IF
    IF (.NOT. ASSOCIATED(Am % PrecValues)) THEN
      ALLOCATE(Am % PrecValues(SIZE(Am % Values)))
    ELSE IF (SIZE(Am % PrecValues) /= SIZE(Am % Values)) THEN
      DEALLOCATE(Am % PrecValues)
      ALLOCATE(Am % PrecValues(SIZE(Am % Values)))
    END IF
    DO i=1,nr
      DO j=Am % Rows(i),Am % Rows(i+1)-1
        IF (isD(i) .OR. isD(Am % Cols(j))) THEN
          Am % PrecValues(j) = Am % Values(j)
        ELSE
          Am % PrecValues(j) = Am % Values(j) + Phase24PrecDelta(j)
        END IF
      END DO
    END DO
    DEALLOCATE(isD)
  END SUBROUTINE Phase24BuildLumpedPrec

  SUBROUTINE Phase24Report()
    REAL(KIND=dp) :: Traversal
    INTEGER :: CacheBuilds, CacheHits, CacheInvalidations

    IF (.NOT. Phase24WallProfiling) RETURN

    CALL CRS_Phase24CacheStats(CacheBuilds, CacheHits, CacheInvalidations)
    WRITE(Message,'(A,3(I0,A))') 'Phase24 CSR cache: builds=',CacheBuilds, &
         ', hits=',CacheHits, ', invalidations=',CacheInvalidations
    CALL Info('HeatSolve',TRIM(Message),Level=4)

    ! Bulk is timed as a whole.  The remainder is deliberately reported as
    ! traversal/other bulk work so the fine-grained intervals remain additive.
    Traversal = MAX(Phase24Bulk - Phase24Material - Phase24Local - Phase24Insert, 0.0_dp)
    WRITE(Message,'(A,ES24.16)') 'component=element_traversal wall_seconds=', Traversal
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=material_property_evaluation wall_seconds=', Phase24Material
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=local_element_calculation wall_seconds=', Phase24Local
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=global_sparse_insertion wall_seconds=', Phase24Insert
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=boundary_element_assembly wall_seconds=', Phase24Boundary
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=matrix_format_conversion wall_seconds=', Phase24Finalize
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=hypre_solve wall_seconds=', Phase24Solve
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,ES24.16)') 'component=udf_circuit wall_seconds=', Phase24UDF
    CALL Info('PHASE24_EVENT',TRIM(Message),Level=4)
    WRITE(Message,'(A,4(I0,1X))') &
         'Phase24 operator counters: full_matrix lagged_reuse forced_rebuild recovery=', &
         Phase24FullRebuildCount,Phase24LaggedReuseCount,Phase24ForcedRebuildCount,Phase24RecoveryCount
    CALL Info('HeatSolve',TRIM(Message),Level=4)
  END SUBROUTINE Phase24Report

  ! Update every inner TES circuit (one TES body each) before an assembly.
  SUBROUTINE TESInnerCircuitUpdate()
    INTEGER :: kc
    LOGICAL :: gotU
    IF (.NOT. ListGetLogical(Solver % Values, 'TES Inner Circuit Update', gotU)) RETURN
    CALL TESCircuitSetup()
    DO kc = 1, TESNCirc
      CALL TESInnerCircuitUpdateOne(kc)
    END DO
  END SUBROUTINE TESInnerCircuitUpdate

  ! Circuit count, TES bodies and keyword style, once.
  SUBROUTINE TESCircuitSetup()
    INTEGER :: n
    LOGICAL :: gotN, gotB
    INTEGER, POINTER :: bodies(:)
    IF (ALLOCATED(TESC)) RETURN
    n = ListGetInteger(Model % Constants, 'TES Circuit Count', gotN)
    IF (gotN) THEN
      IF (n < 1) CALL Fatal('HeatSolve', 'TES Circuit Count must be positive')
      TESNumbered = .TRUE.
      bodies => ListGetIntegerArray(Solver % Values, 'TES Body IDs', gotB)
      IF (.NOT. gotB) CALL Fatal('HeatSolve', 'TES Circuit Count requires TES Body IDs(n) in the solver')
      IF (SIZE(bodies,1) /= n) CALL Fatal('HeatSolve', 'TES Body IDs must list TES Circuit Count bodies')
      ALLOCATE(TESC(n))
      TESC(:) % Body = bodies(1:n)
    ELSE
      n = 1
      TESNumbered = .FALSE.
      ALLOCATE(TESC(1))
      TESC(1) % Body = ListGetInteger(Solver % Values, 'TES Body ID', gotB)
      IF (.NOT. gotB) CALL Fatal('HeatSolver', 'TES Body ID is required for inner circuit update')
    END IF
    TESNCirc = n
    WRITE(Message,'(A,I0,A,L1)') 'TES inner circuits: ', n, ', numbered keywords=', TESNumbered
    CALL Info('HeatSolve', TRIM(Message), Level=3)
  END SUBROUTINE TESCircuitSetup

  ! Constants key of circuit kc: 'TES <suffix>' (single legacy circuit) or
  ! 'TES <kc> <suffix>'.
  FUNCTION TESKey(kc, suffix) RESULT(key)
    INTEGER, INTENT(IN) :: kc
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    CHARACTER(:), ALLOCATABLE :: key
    IF (TESNumbered) THEN
      key = 'TES '//I2S(kc)//' '//suffix
    ELSE
      key = 'TES '//suffix
    END IF
  END FUNCTION TESKey

  FUNCTION TESTag(kc) RESULT(tag)
    INTEGER, INTENT(IN) :: kc
    CHARACTER(:), ALLOCATABLE :: tag
    IF (TESNumbered) THEN
      tag = 'TESInnerCircuit'//I2S(kc)
    ELSE
      tag = 'TESInnerCircuit'
    END IF
  END FUNCTION TESTag

  SUBROUTINE TESInnerCircuitUpdateOne(kc)
    INTEGER, INTENT(IN) :: kc
    LOGICAL :: Enabled, GotIt, FoundState, FoundSeries, FoundIterationSeries, FreezePower, StepCommit
    LOGICAL :: RestartAudit, RestartAuditFatal, StateExists, StateOpenSuccess, StateParseSuccess
    LOGICAL :: AdaptiveRejected, AdaptiveRestoreTrialState, BDF2Predictor
    INTEGER :: CircuitBody, t, i, n, p, LocalCount, TimeStepNow, IoStatus, SeriesUnit, IterationUnit
    INTEGER :: TrialId
    TYPE(Element_t), POINTER :: CircuitElement
    REAL(KIND=dp) :: LocalSum, GlobalSum, GlobalCount, AverageT, UpdateTemperature
    REAL(KIND=dp) :: AbsLocalSum, AbsGlobalSum, AbsGlobalCount, AbsorberAverageT
    REAL(KIND=dp) :: StyLocalSum, StyGlobalSum, StyGlobalCount, StycastAverageT
    REAL(KIND=dp) :: IBias, RSh, R0c, RMin, AlphaC, BetaC, I0c, T0c, Lc
    REAL(KIND=dp) :: A, B, Cc, Disc, ResistanceC, RawCurrent, RawPower, Residual, Denominator
    REAL(KIND=dp) :: StateTemperature, StateCurrent, StateResistance, StatePower, StatePrevious, StateLoaded
    REAL(KIND=dp) :: MidTemperature, ElectricalTau, CommitCurrent, CommitResistance, CommitPower, StepDt
    REAL(KIND=dp) :: PredictionRatio, PredictedPower, PreviousStepCurrent, SaveUpdateT
    LOGICAL :: UseNewton
    INTEGER :: sgn
    REAL(KIND=dp) :: Pk, Tk, Fk, dPdT, Snew, Pnew, hT, Tpm, Apm, Cpm, Dpm, Ipm, Rpm, Ppm(2)
    REAL(KIND=dp) :: NewtonSecantMinRel, Ldt
    ! Converged powers of the last two steps (for a first-solve predictor).
    CHARACTER(LEN=MAX_NAME_LEN) :: StateFile, SeriesFile, IterationSeriesFile, RestartAuditPrefix

    ASSOCIATE ( &
      NewtonPrevP => TESC(kc) % NewtonPrevP,&
      NewtonPrevT => TESC(kc) % NewtonPrevT,&
      NewtonS => TESC(kc) % NewtonS,&
      NewtonSdt => TESC(kc) % NewtonSdt,&
      NewtonHasPrev => TESC(kc) % NewtonHasPrev,&
      StepP1 => TESC(kc) % StepP1,&
      StepP2 => TESC(kc) % StepP2,&
      StepDt1 => TESC(kc) % StepDt1,&
      StepDt2 => TESC(kc) % StepDt2,&
      NewtonPrevStepDt => TESC(kc) % NewtonPrevStepDt,&
      CircuitPower => TESC(kc) % CircuitPower,&
      CircuitCurrent => TESC(kc) % CircuitCurrent,&
      CircuitResistance => TESC(kc) % CircuitResistance,&
      CircuitAverageTemperature => TESC(kc) % CircuitAverageTemperature,&
      CircuitSweepTemperature => TESC(kc) % CircuitSweepTemperature,&
      PreviousCircuitCurrent => TESC(kc) % PreviousCircuitCurrent,&
      PreviousCircuitPower => TESC(kc) % PreviousCircuitPower,&
      CircuitPrevResidual => TESC(kc) % CircuitPrevResidual,&
      CircuitOmega => TESC(kc) % CircuitOmega,&
      CircuitOmegaCap => TESC(kc) % CircuitOmegaCap,&
      LastCircuitDt => TESC(kc) % LastCircuitDt,&
      LastCircuitStep => TESC(kc) % LastCircuitStep,&
      CircuitIterInStep => TESC(kc) % CircuitIterInStep,&
      CircuitCallsInStep => TESC(kc) % CircuitCallsInStep,&
      CircuitInitialized => TESC(kc) % CircuitInitialized,&
      CircuitPowerInitialized => TESC(kc) % CircuitPowerInitialized,&
      SeriesStarted => TESC(kc) % SeriesStarted,&
      IterationSeriesStarted => TESC(kc) % IterationSeriesStarted,&
      SavedCircuitPower => TESC(kc) % SavedCircuitPower,&
      SavedCircuitCurrent => TESC(kc) % SavedCircuitCurrent,&
      SavedCircuitResistance => TESC(kc) % SavedCircuitResistance,&
      SavedCircuitAverageTemperature => TESC(kc) % SavedCircuitAverageTemperature,&
      SavedCircuitSweepTemperature => TESC(kc) % SavedCircuitSweepTemperature,&
      SavedPreviousCircuitCurrent => TESC(kc) % SavedPreviousCircuitCurrent,&
      SavedPreviousCircuitPower => TESC(kc) % SavedPreviousCircuitPower,&
      SavedCircuitPrevResidual => TESC(kc) % SavedCircuitPrevResidual,&
      SavedCircuitOmega => TESC(kc) % SavedCircuitOmega,&
      SavedCircuitOmegaCap => TESC(kc) % SavedCircuitOmegaCap,&
      SavedLastCircuitDt => TESC(kc) % SavedLastCircuitDt,&
      SavedLastCircuitStep => TESC(kc) % SavedLastCircuitStep,&
      SavedCircuitIterInStep => TESC(kc) % SavedCircuitIterInStep,&
      SavedCircuitCallsInStep => TESC(kc) % SavedCircuitCallsInStep,&
      SavedTrialId => TESC(kc) % SavedTrialId,&
      SavedCircuitInitialized => TESC(kc) % SavedCircuitInitialized,&
      SavedCircuitPowerInitialized => TESC(kc) % SavedCircuitPowerInitialized,&
      SavedAdaptiveState => TESC(kc) % SavedAdaptiveState )
    Enabled = ListGetLogical(Solver % Values, 'TES Inner Circuit Update', GotIt)
    IF (.NOT. GotIt .OR. .NOT. Enabled) RETURN
    BDF2Predictor = ListGetLogical(Solver % Values, 'Phase24 BDF2 Predictor Active', GotIt)
    IF (.NOT. GotIt) BDF2Predictor = .FALSE.
    RestartAudit = ListGetLogical(Solver % Values, 'Phase24 Restart State Audit', GotIt)
    IF (.NOT. GotIt) RestartAudit = .FALSE.
    RestartAuditFatal = ListGetLogical(Solver % Values, 'Phase24 Restart Audit Fail Fast', GotIt)
    IF (.NOT. GotIt) RestartAuditFatal = .FALSE.
    RestartAuditPrefix = ListGetString(Solver % Values, 'Phase24 Restart State Audit Prefix', GotIt)
    IF (.NOT. GotIt) RestartAuditPrefix = 'phase24_restart_state_audit'
    StepCommit = ListGetLogical(Solver % Values, 'TES Inner Circuit Step Commit', GotIt)
    IF (.NOT. GotIt) StepCommit = .FALSE.
    AdaptiveRejected = ListGetLogical(Solver % Values, 'Adaptive Step Rejected', GotIt)
    IF (.NOT. GotIt) AdaptiveRejected = .FALSE.
    AdaptiveRestoreTrialState = ListGetLogical(Solver % Values, 'Adaptive Restore Trial State', GotIt)
    IF (.NOT. GotIt) AdaptiveRestoreTrialState = .FALSE.
    TrialId = ListGetInteger(Solver % Values, 'Adaptive Trial ID', GotIt)
    IF (.NOT. GotIt) TrialId = -1
    FreezePower = ListGetLogical(Solver % Values, 'TES Inner Circuit Freeze Power', GotIt)
    StateFile = ListGetString(Model % Constants, TESKey(kc,'State File'), FoundState)
    INQUIRE(FILE=TRIM(StateFile), EXIST=StateExists)
    StateOpenSuccess = .FALSE.
    StateParseSuccess = .FALSE.
    SeriesFile = ListGetString(Model % Constants, TESKey(kc,'Series File'), FoundSeries)
    IterationSeriesFile = ListGetString(Model % Constants, TESKey(kc,'Iteration Series File'), FoundIterationSeries)
    CircuitBody = TESC(kc) % Body
    LocalSum = 0.0_dp
    LocalCount = 0
    DO t = 1, Solver % NumberOfActiveElements
      CircuitElement => GetActiveElement(t)
      ! Parallel meshes include halo elements on neighbouring ranks.  Count
      ! each physical element only on its owning rank before the all-reduce.
      IF (ParEnv % PEs > 1 .AND. CircuitElement % PartIndex /= ParEnv % MyPE) CYCLE
      IF (CircuitElement % BodyId /= CircuitBody) CYCLE
      n = GetElementNOFNodes(CircuitElement)
      IF (n <= 0) CYCLE
      DO i = 1, n
        p = TempPerm(CircuitElement % NodeIndexes(i))
        IF (p > 0) LocalSum = LocalSum + Temperature(p) / REAL(n, dp)
      END DO
      LocalCount = LocalCount + 1
    END DO
    GlobalSum = ParallelReduction(LocalSum)
    GlobalCount = ParallelReduction(REAL(LocalCount, dp))
    IF (GlobalCount <= 0.0_dp) CALL Fatal('HeatSolver', 'No TES elements found for inner circuit update')
    IF (.NOT. CircuitInitialized) THEN
      WRITE(Message,'(A,I0,A,I0,A,ES14.6,A,I0)') 'TES circuit elements: local=', LocalCount, &
           ', global=', NINT(GlobalCount), ', mean T=', GlobalSum/GlobalCount, &
           ', active elements=', Solver % NumberOfActiveElements
      CALL Info('HeatSolve', TRIM(Message), Level=3)
    END IF

    ! The driver can retry a trial after this routine has already advanced
    ! the lumped circuit.  Restore the snapshot made at the trial boundary
    ! before evaluating the retry; output interpolation never reaches here.
    IF ((AdaptiveRejected .OR. AdaptiveRestoreTrialState) .AND. SavedAdaptiveState) THEN
      CircuitPower = SavedCircuitPower
      CircuitCurrent = SavedCircuitCurrent
      CircuitResistance = SavedCircuitResistance
      CircuitAverageTemperature = SavedCircuitAverageTemperature
      CircuitSweepTemperature = SavedCircuitSweepTemperature
      PreviousCircuitCurrent = SavedPreviousCircuitCurrent
      PreviousCircuitPower = SavedPreviousCircuitPower
      CircuitPrevResidual = SavedCircuitPrevResidual
      CircuitOmega = SavedCircuitOmega
      CircuitOmegaCap = SavedCircuitOmegaCap
      LastCircuitDt = SavedLastCircuitDt
      LastCircuitStep = SavedLastCircuitStep
      CircuitIterInStep = SavedCircuitIterInStep
      CircuitCallsInStep = SavedCircuitCallsInStep
      CircuitInitialized = SavedCircuitInitialized
      CircuitPowerInitialized = SavedCircuitPowerInitialized
      CALL ListAddLogical(Solver % Values,'Adaptive Step Rejected',.FALSE.)
      CALL ListAddLogical(Solver % Values,'Adaptive Restore Trial State',.FALSE.)
    END IF
    IF (TrialId /= SavedTrialId) THEN
      SavedCircuitPower = CircuitPower
      SavedCircuitCurrent = CircuitCurrent
      SavedCircuitResistance = CircuitResistance
      SavedCircuitAverageTemperature = CircuitAverageTemperature
      SavedCircuitSweepTemperature = CircuitSweepTemperature
      SavedPreviousCircuitCurrent = PreviousCircuitCurrent
      SavedPreviousCircuitPower = PreviousCircuitPower
      SavedCircuitPrevResidual = CircuitPrevResidual
      SavedCircuitOmega = CircuitOmega
      SavedCircuitOmegaCap = CircuitOmegaCap
      SavedLastCircuitDt = LastCircuitDt
      SavedLastCircuitStep = LastCircuitStep
      SavedCircuitIterInStep = CircuitIterInStep
      SavedCircuitCallsInStep = CircuitCallsInStep
      SavedCircuitInitialized = CircuitInitialized
      SavedCircuitPowerInitialized = CircuitPowerInitialized
      SavedTrialId = TrialId
      SavedAdaptiveState = .TRUE.
    END IF
    AverageT = GlobalSum / GlobalCount
    ! Collect the absorber average separately from active elements.  VTU
    ! partitions contain halo values, so a raw pointwise min/max is not a
    ! reliable pulse-energy diagnostic in MPI runs.
    AbsLocalSum = 0.0_dp
    LocalCount = 0
    DO t = 1, Solver % NumberOfActiveElements
      CircuitElement => GetActiveElement(t)
      IF (ParEnv % PEs > 1 .AND. CircuitElement % PartIndex /= ParEnv % MyPE) CYCLE
      IF (CircuitElement % BodyId /= 1) CYCLE
      n = GetElementNOFNodes(CircuitElement)
      IF (n <= 0) CYCLE
      DO i = 1, n
        p = TempPerm(CircuitElement % NodeIndexes(i))
        IF (p > 0) AbsLocalSum = AbsLocalSum + Temperature(p) / REAL(n, dp)
      END DO
      LocalCount = LocalCount + 1
    END DO
    AbsGlobalSum = ParallelReduction(AbsLocalSum)
    AbsGlobalCount = ParallelReduction(REAL(LocalCount, dp))
    IF (AbsGlobalCount > 0.0_dp) THEN
      AbsorberAverageT = AbsGlobalSum / AbsGlobalCount
    ELSE
      AbsorberAverageT = 0.0_dp
    END IF
    ! Record the intermediate Stycast temperature as a cheap thermal-path
    ! diagnostic.  The existing TESInnerCircuit log already reports the TES
    ! average and the absorber average; adding Stycast lets short-window
    ! runs distinguish interface/mesh transport from circuit-temperature
    ! extraction without writing large VTU fields.
    StyLocalSum = 0.0_dp
    LocalCount = 0
    DO t = 1, Solver % NumberOfActiveElements
      CircuitElement => GetActiveElement(t)
      IF (ParEnv % PEs > 1 .AND. CircuitElement % PartIndex /= ParEnv % MyPE) CYCLE
      IF (CircuitElement % BodyId /= 3) CYCLE
      n = GetElementNOFNodes(CircuitElement)
      IF (n <= 0) CYCLE
      DO i = 1, n
        p = TempPerm(CircuitElement % NodeIndexes(i))
        IF (p > 0) StyLocalSum = StyLocalSum + Temperature(p) / REAL(n, dp)
      END DO
      LocalCount = LocalCount + 1
    END DO
    StyGlobalSum = ParallelReduction(StyLocalSum)
    StyGlobalCount = ParallelReduction(REAL(LocalCount, dp))
    IF (StyGlobalCount > 0.0_dp) THEN
      StycastAverageT = StyGlobalSum / StyGlobalCount
    ELSE
      StycastAverageT = 0.0_dp
    END IF
    IBias = GetConstReal(Model % Constants, TESKey(kc,'Bias Current'), GotIt)
    RSh = GetConstReal(Model % Constants, TESKey(kc,'Shunt Resistance'), GotIt)
    R0c = GetConstReal(Model % Constants, TESKey(kc,'R0'), GotIt)
    RMin = GetConstReal(Model % Constants, TESKey(kc,'Rmin'), GotIt)
    AlphaC = GetConstReal(Model % Constants, TESKey(kc,'Alpha'), GotIt)
    BetaC = GetConstReal(Model % Constants, TESKey(kc,'Beta'), GotIt)
    I0c = GetConstReal(Model % Constants, TESKey(kc,'I0'), GotIt)
    T0c = GetConstReal(Model % Constants, TESKey(kc,'T0'), GotIt)
    Lc = GetConstReal(Model % Constants, TESKey(kc,'Inductance'), GotIt)
    IF (.NOT. CircuitInitialized) THEN
      ! Restore on root and distribute through collective reductions.  Every
      ! rank executes the same reductions, including when no state file is set.
      StateTemperature = 0.0_dp
      StateCurrent = 0.0_dp
      StateResistance = 0.0_dp
      StatePower = 0.0_dp
      StatePrevious = 0.0_dp
      StateLoaded = 0.0_dp
      ! A steady solve can also be a restart of a previously converged Gate3
      ! field.  In that case the thermal .result alone is not the augmented
      ! state: restore the circuit checkpoint before the first assembly too.
      ! Fresh steady cases still fall back exactly as before when the file is
      ! absent or the Constants entry is not present.
      IF (FoundState .AND. ParEnv % MyPE == 0) THEN
        OPEN(UNIT=98, FILE=TRIM(StateFile), STATUS='OLD', ACTION='READ', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          StateOpenSuccess = .TRUE.
          READ(98, *, IOSTAT=IoStatus) StateTemperature, StateCurrent, StateResistance, StatePower, StatePrevious
          CLOSE(98)
          IF (IoStatus == 0) THEN
            StateParseSuccess = .TRUE.
            StateLoaded = 1.0_dp
          END IF
        END IF
      END IF
      StateLoaded = ParallelReduction(StateLoaded)
      StateTemperature = ParallelReduction(StateTemperature)
      StateCurrent = ParallelReduction(StateCurrent)
      StateResistance = ParallelReduction(StateResistance)
      StatePower = ParallelReduction(StatePower)
      StatePrevious = ParallelReduction(StatePrevious)
      IF (RestartAudit .AND. ParEnv % MyPE == 0) THEN
        CALL Phase24WriteRestartAudit(TRIM(RestartAuditPrefix), 'load_after', StateFile, &
             StateExists, StateOpenSuccess, StateParseSuccess, StateLoaded > 0.5_dp, &
             StateTemperature, StateCurrent, StateResistance, StatePower, StatePrevious, &
             CircuitAverageTemperature, CircuitSweepTemperature, CircuitOmega, CircuitOmegaCap, &
             CircuitPrevResidual, LastCircuitDt, LastCircuitStep, CircuitIterInStep, CircuitCallsInStep, &
             CircuitInitialized, CircuitPowerInitialized, GetTime(), GetTimeStep(), dt, &
             MIN(Solver % DoneTime, Solver % Order), 'checkpoint_or_default')
      END IF
      IF (RestartAuditFatal .AND. TransientSimulation .AND. StateLoaded <= 0.5_dp) THEN
        CALL Fatal('HeatSolve', 'Phase24 restart audit: checkpoint load failed')
      END IF
      IF (StateLoaded > 0.5_dp) THEN
        CircuitPower = StatePower
        CircuitCurrent = StateCurrent
        CircuitResistance = StateResistance
        CircuitAverageTemperature = StateTemperature
        PreviousCircuitCurrent = StatePrevious
        CircuitPowerInitialized = .TRUE.
      ELSE
        ! Match the frozen serial UDF: initialize the electrical branch at
        ! the analytic T0 operating point before the first transient step.
        A = R0c * (1.0_dp - BetaC)
        B = R0c * BetaC / I0c
        Disc = MAX((RSh + A)**2 + 4.0_dp*B*IBias*RSh, 0.0_dp)
        CircuitCurrent = MAX(MIN((SQRT(Disc) - (RSh + A))/(2.0_dp*B), IBias), 0.0_dp)
        ResistanceC = MAX(A + B*ABS(CircuitCurrent), RMin)
        IF (ResistanceC == RMin) CircuitCurrent = IBias*RSh/(RSh + ResistanceC)
        CircuitPower = CircuitCurrent*CircuitCurrent*ResistanceC
        CircuitResistance = ResistanceC
        CircuitAverageTemperature = T0c
        PreviousCircuitCurrent = CircuitCurrent
        CircuitPowerInitialized = .TRUE.
      END IF
      CircuitInitialized = .TRUE.
      PreviousCircuitPower = CircuitPower
      LastCircuitStep = GetTimeStep()
      IF (StepCommit .AND. TransientSimulation .AND. dt > 0.0_dp) LastCircuitDt = dt
      ! The serial body-force UDF samples Temperature during assembly and
      ! consumes that average only when the following sweep begins.  Preserve
      ! the same one-sweep delay instead of using the newly solved field one
      ! nonlinear iteration too early.
      CircuitSweepTemperature = AverageT
      ! The serial UDF initializes IterInStep to one.  Its next nonlinear
      ! sweep therefore follows the carried-state secant path, whereas true
      ! timestep boundaries below explicitly restart at zero.
      CircuitIterInStep = 1
      CircuitCallsInStep = 1
      ! The serial UDF only initializes on its first assembly call.  In
      ! particular, do not turn that initialization into a first Aitken
      ! update, because its temperature is not a completed prior sweep.
      CALL ListAddConstReal(Model % Constants, TESKey(kc,'Parallel Power'), CircuitPower)
      IF (RestartAudit .AND. ParEnv % MyPE == 0) THEN
        CALL Phase24WriteRestartAudit(TRIM(RestartAuditPrefix), 'circuit_init_after', StateFile, &
             StateExists, StateOpenSuccess, StateParseSuccess, StateLoaded > 0.5_dp, &
             CircuitAverageTemperature, CircuitCurrent, CircuitResistance, CircuitPower, PreviousCircuitCurrent, &
             CircuitAverageTemperature, CircuitSweepTemperature, CircuitOmega, CircuitOmegaCap, &
             CircuitPrevResidual, LastCircuitDt, LastCircuitStep, CircuitIterInStep, CircuitCallsInStep, &
             CircuitInitialized, CircuitPowerInitialized, GetTime(), GetTimeStep(), dt, &
             MIN(Solver % DoneTime, Solver % Order), 'loaded_or_default_then_init')
        CALL Phase24WriteRestartAudit(TRIM(RestartAuditPrefix), 'pre_first_assembly', StateFile, &
             StateExists, StateOpenSuccess, StateParseSuccess, StateLoaded > 0.5_dp, &
             CircuitAverageTemperature, CircuitCurrent, CircuitResistance, CircuitPower, PreviousCircuitCurrent, &
             CircuitAverageTemperature, CircuitSweepTemperature, CircuitOmega, CircuitOmegaCap, &
             CircuitPrevResidual, LastCircuitDt, LastCircuitStep, CircuitIterInStep, CircuitCallsInStep, &
             CircuitInitialized, CircuitPowerInitialized, GetTime(), GetTimeStep(), dt, &
             MIN(Solver % DoneTime, Solver % Order), 'preserved_before_first_assembly')
      END IF
      IF (FoundIterationSeries .AND. ParEnv % MyPE == 0) THEN
        IterationUnit = 96
        OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          WRITE(IterationUnit,'(A)') &
            'time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,' // &
            'raw_current_A,tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,relaxed_power_W'
          WRITE(IterationUnit, &
            '(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,' // &
            'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
            GetTime(), ',', GetTimeStep(), ',', iter, ',', CircuitAverageTemperature, ',', &
            PreviousCircuitCurrent, ',', CircuitCurrent, ',', CircuitResistance, ',', CircuitPower, ',', &
            0.0_dp, ',', CircuitOmega, ',', CircuitOmegaCap, ',', CircuitPower
          CLOSE(IterationUnit)
          IterationSeriesStarted = .TRUE.
        END IF
      END IF
      RETURN
    END IF
    UpdateTemperature = CircuitSweepTemperature
    ! Default: the circuit sees the TES temperature of the previous call, a
    ! one-iteration lag on top of the Picard lag.  Opt-in: use the current
    ! field (standard Picard); the converged fixed point is unchanged, only
    ! the number of coupled iterations is.
    IF (BDF2Predictor .OR. &
         ListGetLogical(Solver % Values, 'TES Inner Circuit Current Temperature', GotIt) .OR. &
         ListGetLogical(Solver % Values, 'TES Inner Circuit Newton', GotIt)) &
         UpdateTemperature = AverageT
    A = R0c * (1.0_dp + AlphaC * (UpdateTemperature - T0c) / T0c - BetaC)
    B = R0c * BetaC / I0c
    TimeStepNow = GetTimeStep()
    IF (TransientSimulation .AND. TimeStepNow /= LastCircuitStep) THEN
      ! This call starts a new timestep, so the values retained here belong
      ! to the converged preceding timestep.  Commit and emit that preceding
      ! state before computing any new-step circuit iterate.
      ! The commit must see the converged field of that step: a timestep
      ! predictor has already replaced the current field by its extrapolation,
      ! which would otherwise leak into the emitted current and into the
      ! inductor history (PreviousCircuitCurrent).
      SaveUpdateT = UpdateTemperature
      IF (ASSOCIATED(Solver % Variable % PrevValues)) UpdateTemperature = TESMeanPreviousStep(kc)
      CommitCurrent = CircuitCurrent
      CommitResistance = CircuitResistance
      CommitPower = CircuitPower
      PreviousStepCurrent = PreviousCircuitCurrent
      IF (StepCommit) THEN
        StepDt = LastCircuitDt
        IF (StepDt <= 0.0_dp) StepDt = dt
        IF (StepDt <= 0.0_dp) CALL Fatal('HeatSolve', 'TES step commit requires positive timestep')
        A = R0c * (1.0_dp + AlphaC * (UpdateTemperature - T0c) / T0c - BetaC)
        B = R0c * BetaC / I0c
        Cc = RSh + A + Lc / StepDt
        Disc = MAX(Cc*Cc + 4.0_dp*B*(IBias*RSh + Lc*PreviousCircuitCurrent/StepDt), 0.0_dp)
        CommitCurrent = MAX(MIN((SQRT(Disc)-Cc)/(2.0_dp*B), IBias), 0.0_dp)
        CommitResistance = A + B*ABS(CommitCurrent)
        IF (CommitResistance < RMin) THEN
          CommitResistance = RMin
          CommitCurrent = (IBias*RSh + Lc*PreviousCircuitCurrent/StepDt) / &
            (RSh + CommitResistance + Lc/StepDt)
        END IF
        CommitPower = CommitCurrent*CommitCurrent*CommitResistance
      ELSE
        ElectricalTau = Lc / MAX(RSh + CircuitResistance, 1.0e-30_dp)
        IF (CircuitCallsInStep <= 2 .AND. LastCircuitDt >= 100.0_dp*ElectricalTau .AND. &
          CircuitCurrent < I0c .AND. &
          (UpdateTemperature-CircuitAverageTemperature) * &
          (CircuitCurrent-PreviousCircuitCurrent) > 0.0_dp) THEN
        ! A two-sweep thermal convergence can leave the circuit evaluated
        ! with the preceding sweep even though the just-converged temperature
        ! moves in the opposite electrothermal direction.  On timesteps much
        ! longer than the electrical time constant, reconcile that stale
        ! endpoint at the midpoint of the two thermal sweeps.  This removes a
        ! partition-dependent one-point spike without changing the Joule-power
        ! fixed point used by the completed thermal solve.
        MidTemperature = 0.5_dp * (CircuitAverageTemperature + UpdateTemperature)
        A = R0c * (1.0_dp + AlphaC * (MidTemperature - T0c) / T0c - BetaC)
        B = R0c * BetaC / I0c
        Cc = RSh + A + Lc / dt
        Disc = MAX(Cc*Cc + 4.0_dp*B*(IBias*RSh + Lc*PreviousCircuitCurrent/dt), 0.0_dp)
        CommitCurrent = MAX(MIN((SQRT(Disc) - Cc) / (2.0_dp*B), IBias), 0.0_dp)
        CommitResistance = A + B*ABS(CommitCurrent)
        IF (CommitResistance < RMin) THEN
          CommitResistance = RMin
          CommitCurrent = (IBias*RSh + Lc*PreviousCircuitCurrent/dt) / &
            (RSh + CommitResistance + Lc/dt)
        END IF
        END IF
      END IF
      IF (FoundSeries .AND. ParEnv % MyPE == 0) THEN
        SeriesUnit = 97
        IF (.NOT. SeriesStarted) THEN
          OPEN(UNIT=SeriesUnit, FILE=TRIM(SeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
          IF (IoStatus == 0) THEN
            WRITE(SeriesUnit,'(A)') 'time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W'
            SeriesStarted = .TRUE.
          END IF
        ELSE
          OPEN(UNIT=SeriesUnit, FILE=TRIM(SeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
        END IF
        IF (IoStatus == 0) THEN
          WRITE(SeriesUnit,'(ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
            GetTime()-dt, ',', UpdateTemperature, ',', CommitCurrent, ',', CommitResistance, ',', CommitPower
          CLOSE(SeriesUnit)
        END IF
      END IF
      IF (StepCommit) THEN
        PreviousCircuitCurrent = CommitCurrent
        CircuitCurrent = CommitCurrent
        CircuitResistance = CommitResistance
        CircuitPower = CommitPower
        CircuitAverageTemperature = UpdateTemperature
      ELSE
        PreviousCircuitCurrent = CircuitCurrent
      END IF
      IF (BDF2Predictor .AND. .NOT. FreezePower .AND. &
           PreviousCircuitPower >= 0.0_dp .AND. dt > 0.0_dp) THEN
        CALL Phase24GetBDFCoefficients(2,dt,StepDt, &
             Phase24CurrentA0,Phase24CurrentA1,Phase24CurrentA2)
        PredictionRatio = dt / StepDt
        PredictedPower = CommitPower + PredictionRatio*(CommitPower-PreviousCircuitPower)
        CircuitPower = MIN(MAX(PredictedPower,0.5_dp*CommitPower),1.5_dp*CommitPower)
        CircuitCurrent = MIN(MAX(CommitCurrent + PredictionRatio* &
             (CommitCurrent-PreviousStepCurrent),0.0_dp),IBias)
        IF (ListGetLogical(Solver % Values,'TES Inner Circuit Predict From Temperature',GotIt)) THEN
          ! Start from the circuit fixed point at the predicted TES
          ! temperature (the field predictor has already been applied), so the
          ! first solve begins from a consistent (T, P) pair.
          Apm = R0c * (1.0_dp + AlphaC * (SaveUpdateT - T0c) / T0c - BetaC)
          Cpm = RSh + Apm + Lc / dt
          Dpm = MAX(Cpm*Cpm + 4.0_dp*B*(IBias*RSh + Lc*PreviousCircuitCurrent/dt), 0.0_dp)
          Ipm = MAX(MIN((SQRT(Dpm) - Cpm) / (2.0_dp*B), IBias), 0.0_dp)
          Rpm = Apm + B*ABS(Ipm)
          IF (Rpm < RMin) THEN
            Rpm = RMin
            Ipm = (IBias*RSh + Lc*PreviousCircuitCurrent/dt) / (RSh + Rpm + Lc/dt)
          END IF
          CircuitPower = MIN(MAX(Ipm*Ipm*Rpm,0.5_dp*CommitPower),1.5_dp*CommitPower)
          CircuitCurrent = Ipm
        END IF
      END IF
      PreviousCircuitPower = CommitPower
      CircuitPrevResidual = 0.0_dp
      CircuitOmegaCap = 0.5_dp
      CircuitIterInStep = 0
      CircuitCallsInStep = 0
      IF (LastCircuitDt > 0.0_dp .AND. ABS(dt-LastCircuitDt) > 1.0e-9_dp*LastCircuitDt) THEN
        CircuitOmega = 0.1_dp
        CircuitOmegaCap = 0.25_dp
      END IF
      NewtonPrevStepDt = LastCircuitDt
      LastCircuitDt = dt
      LastCircuitStep = TimeStepNow
      UpdateTemperature = SaveUpdateT
    END IF
    ! The optional commit-only reconciliation above uses midpoint
    ! coefficients.  Restore the actual sweep coefficients before advancing
    ! the circuit state for the new timestep.
    A = R0c * (1.0_dp + AlphaC * (UpdateTemperature - T0c) / T0c - BetaC)
    B = R0c * BetaC / I0c
    IF (TransientSimulation .AND. CircuitInitialized .AND. dt > 0.0_dp) THEN
      Cc = RSh + A + Lc / dt
      Disc = MAX(Cc*Cc + 4.0_dp*B*(IBias*RSh + Lc*PreviousCircuitCurrent/dt), 0.0_dp)
      RawCurrent = (SQRT(Disc) - Cc) / (2.0_dp*B)
    ELSE
      Disc = MAX((RSh + A)**2 + 4.0_dp*B*IBias*RSh, 0.0_dp)
      RawCurrent = (SQRT(Disc) - (RSh + A)) / (2.0_dp*B)
    END IF
    RawCurrent = MAX(MIN(RawCurrent, IBias), 0.0_dp)
    ResistanceC = A + B*ABS(RawCurrent)
    IF (ResistanceC < RMin) THEN
      ResistanceC = RMin
      IF (TransientSimulation .AND. dt > 0.0_dp) THEN
        RawCurrent = (IBias*RSh + Lc*PreviousCircuitCurrent/dt) / &
          (RSh + ResistanceC + Lc/dt)
      ELSE
        RawCurrent = IBias*RSh/(RSh + ResistanceC)
      END IF
    END IF
    RawPower = RawCurrent*RawCurrent*ResistanceC
    Residual = 0.0_dp
    ! Numerical Aitken update of the nonlinear electrothermal fixed point.
    ! This is the same nonlinear solve procedure as the validated serial UDF;
    ! it is not a physical time relaxation or an added circuit parameter.
    UseNewton = ListGetLogical(Solver % Values, 'TES Inner Circuit Newton', GotIt)
    ! Steady solves use the same Newton update without the inductor (L/dt).
    Ldt = 0.0_dp
    IF (TransientSimulation .AND. dt > 0.0_dp) Ldt = Lc / dt
    IF (UseNewton .AND. CircuitPowerInitialized .AND. .NOT. FreezePower .AND. &
         (.NOT. TransientSimulation .OR. dt > 0.0_dp)) THEN
      ! Opt-in scalar Newton on F(P) = Pcircuit(T_TES(P)) - P.  Within a
      ! timestep the heat problem is linear in the TES power, so dT/dP is the
      ! secant of the last two (P,T) pairs of this step (kept across steps of
      ! equal dt); dPcircuit/dT is a central difference of the implicit
      ! circuit law.  The converged fixed point is the same as the Aitken path.
      Pk = CircuitPower
      Tk = UpdateTemperature
      Fk = RawPower - Pk
      NewtonSecantMinRel = 1.0e-10_dp
      IF (ListGetLogical(Solver % Values,'TES Inner Circuit Newton Slope Reuse',GotIt)) &
           NewtonSecantMinRel = 1.0e-7_dp
      IF (CircuitCallsInStep == 0) THEN
        ! The field still belongs to the previous step.  Predict P for the
        ! first solve by linear extrapolation of the last converged step
        ! powers when the step size is unchanged; otherwise keep P.
        NewtonHasPrev = .FALSE.
        StepP2 = StepP1
        StepDt2 = StepDt1
        StepP1 = Pk
        StepDt1 = NewtonPrevStepDt
        IF (.NOT. BDF2Predictor .AND. StepP2 > 0.0_dp .AND. ABS(StepDt1 - dt) <= 1.0e-9_dp*dt .AND. &
             ABS(StepDt2 - dt) <= 1.0e-9_dp*dt) THEN
          CircuitPower = MIN(MAX(2.0_dp*StepP1 - StepP2, 0.5_dp*Pk), 1.5_dp*Pk)
        END IF
      ELSE
        IF (NewtonHasPrev .AND. ABS(Pk - NewtonPrevP) > NewtonSecantMinRel*MAX(Pk,1.0e-30_dp)) THEN
          Snew = (Tk - NewtonPrevT) / (Pk - NewtonPrevP)
          ! A secant over a power change near the linear-solve noise can be
          ! wildly large and then freezes the Newton update; dT/dP of a step
          ! can never exceed the static T/P.
          IF (Snew > 0.0_dp .AND. Snew <= Tk / MAX(Pk,1.0e-30_dp)) THEN
            NewtonS = Snew
            NewtonSdt = dt
          END IF
        END IF
        hT = 1.0e-7_dp
        DO sgn = 1, 2
          Tpm = Tk + MERGE(hT, -hT, sgn == 1)
          Apm = R0c * (1.0_dp + AlphaC * (Tpm - T0c) / T0c - BetaC)
          Cpm = RSh + Apm + Ldt
          Dpm = MAX(Cpm*Cpm + 4.0_dp*B*(IBias*RSh + Ldt*PreviousCircuitCurrent), 0.0_dp)
          Ipm = MAX(MIN((SQRT(Dpm) - Cpm) / (2.0_dp*B), IBias), 0.0_dp)
          Rpm = Apm + B*ABS(Ipm)
          IF (Rpm < RMin) THEN
            Rpm = RMin
            Ipm = (IBias*RSh + Ldt*PreviousCircuitCurrent) / (RSh + Rpm + Ldt)
          END IF
          Ppm(sgn) = Ipm*Ipm*Rpm
        END DO
        dPdT = (Ppm(1) - Ppm(2)) / (2.0_dp*hT)
        IF (NewtonS > 0.0_dp .AND. ABS(NewtonSdt - dt) <= 1.0e-9_dp*dt) THEN
          Pnew = Pk + Fk / (1.0_dp - dPdT*NewtonS)
        ELSE IF (NewtonS > 0.0_dp .AND. NewtonSdt > 0.0_dp .AND. &
             dt >= 0.5_dp*NewtonSdt .AND. dt <= 2.0_dp*NewtonSdt .AND. &
             ListGetLogical(Solver % Values,'TES Inner Circuit Newton Slope Reuse',GotIt)) THEN
          ! dT/dP of a step grows with dt (between dt/C and 1/G); scale the
          ! slope of the previous step size geometrically.  Only the rate of
          ! the Newton iteration depends on this estimate, not its fixed point.
          Pnew = Pk + Fk / (1.0_dp - dPdT*NewtonS*SQRT(dt/NewtonSdt))
        ELSE
          Pnew = Pk + 0.5_dp*Fk
        END IF
        Pnew = MIN(MAX(Pnew, 0.5_dp*Pk), 1.5_dp*Pk)
        NewtonPrevP = Pk
        NewtonPrevT = Tk
        NewtonHasPrev = .TRUE.
        CircuitPower = Pnew
      END IF
      Residual = Fk
      CircuitPrevResidual = Fk
      CircuitIterInStep = CircuitIterInStep + 1
    ELSE IF (.NOT. FreezePower .OR. .NOT. CircuitPowerInitialized) THEN
      Residual = RawPower - CircuitPower
      IF (ABS(Residual) > 1.0e-6_dp * MAX(ABS(CircuitPower), 1.0e-30_dp)) THEN
        IF (.NOT. TransientSimulation) THEN
          ! Match the frozen serial UDF's steady solve.  Its large undamped
          ! electrothermal loop gain makes the adaptive secant unreliable.
          CircuitOmega = 0.04_dp
        ELSE IF (CircuitIterInStep > 0) THEN
          IF (ABS(Residual) > 1.5_dp*ABS(CircuitPrevResidual)) THEN
            CircuitOmegaCap = MAX(0.5_dp*CircuitOmegaCap, 0.02_dp)
          ELSE
            CircuitOmegaCap = MIN(1.3_dp*CircuitOmegaCap, 1.0_dp)
          END IF
          Denominator = Residual - CircuitPrevResidual
          IF (ABS(Denominator) > 1.0e-2_dp*ABS(Residual)) CircuitOmega = -CircuitOmega*CircuitPrevResidual/Denominator
        END IF
        CircuitOmega = MAX(MIN(CircuitOmega, CircuitOmegaCap), 0.02_dp)
        CircuitPower = MAX(CircuitPower + CircuitOmega*Residual, 0.0_dp)
        CircuitPrevResidual = Residual
        CircuitIterInStep = CircuitIterInStep + 1
      END IF
    END IF
    CircuitPowerInitialized = .TRUE.
    TESC(kc) % ResValid = .NOT. FreezePower .AND. (.NOT. TransientSimulation .OR. dt > 0.0_dp)
    TESC(kc) % ResIBias = IBias; TESC(kc) % ResRSh = RSh; TESC(kc) % ResR0 = R0c; TESC(kc) % ResRMin = RMin
    TESC(kc) % ResAlpha = AlphaC; TESC(kc) % ResBeta = BetaC; TESC(kc) % ResI0 = I0c; TESC(kc) % ResT0 = T0c
    TESC(kc) % ResL = Lc; TESC(kc) % ResPrevI = PreviousCircuitCurrent; TESC(kc) % ResP = CircuitPower
    TESC(kc) % ResDt = 0.0_dp
    IF (TransientSimulation) TESC(kc) % ResDt = dt
    CircuitCurrent = RawCurrent
    CircuitResistance = ResistanceC
    CircuitAverageTemperature = UpdateTemperature
    CircuitSweepTemperature = AverageT
    CircuitCallsInStep = CircuitCallsInStep + 1
    IF (.NOT. TransientSimulation .AND. FoundState .AND. ParEnv % MyPE == 0) THEN
      OPEN(UNIT=98, FILE=TRIM(StateFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
      IF (IoStatus == 0) THEN
        ! A steady solution has no earlier transient step.  Seed the
        ! subsequent backward-Euler solve with its converged current rather
        ! than the module's startup value (which is zero on a fresh run).
        WRITE(98,'(5ES24.16)') UpdateTemperature, CircuitCurrent, ResistanceC, CircuitPower, CircuitCurrent
        CLOSE(98)
      END IF
    END IF
    CALL ListAddConstReal(Model % Constants, TESKey(kc,'Parallel Power'), CircuitPower)
    IF (FoundIterationSeries .AND. ParEnv % MyPE == 0) THEN
      IterationUnit = 96
      IF (.NOT. IterationSeriesStarted) THEN
        OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
        IF (IoStatus == 0) THEN
          WRITE(IterationUnit,'(A)') &
            'time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,' // &
            'raw_current_A,tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,relaxed_power_W'
          IterationSeriesStarted = .TRUE.
        END IF
      ELSE
        OPEN(UNIT=IterationUnit, FILE=TRIM(IterationSeriesFile), STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      END IF
      IF (IoStatus == 0) THEN
        WRITE(IterationUnit, &
          '(ES24.16,A,I0,A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,' // &
          'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16)') &
          GetTime(), ',', TimeStepNow, ',', iter, ',', UpdateTemperature, ',', PreviousCircuitCurrent, ',', &
          RawCurrent, ',', ResistanceC, ',', RawPower, ',', Residual, ',', CircuitOmega, ',', &
          CircuitOmegaCap, ',', CircuitPower
        CLOSE(IterationUnit)
      END IF
    END IF
    IF (ParEnv % MyPE == 0) WRITE(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
      TESTag(kc)//': T=', UpdateTemperature, ' I=', CircuitCurrent, ' P=', CircuitPower
    IF (ParEnv % MyPE == 0 .AND. AbsGlobalCount > 0.0_dp) WRITE(*,'(A,ES14.6)') &
      TESTag(kc)//': AbsorberT=', AbsorberAverageT
    IF (ParEnv % MyPE == 0 .AND. StyGlobalCount > 0.0_dp) WRITE(*,'(A,ES14.6)') &
      TESTag(kc)//': StycastT=', StycastAverageT
    END ASSOCIATE
  END SUBROUTINE TESInnerCircuitUpdateOne

  SUBROUTINE Phase24WriteRestartAudit(Prefix, Stage, StateFile, StateExists, StateOpenSuccess, &
       StateParseSuccess, StateLoaded, Temperature, Current, Resistance, Power, PreviousCurrent, &
       AverageTemperature, SweepTemperature, Omega, OmegaCap, PrevResidual, LastDt, LastStep, &
       IterInStep, CallsInStep, Initialized, PowerInitialized, PhysicalTime, TimeStep, Dt, BdfOrder, Source)
    CHARACTER(LEN=*), INTENT(IN) :: Prefix, Stage, StateFile, Source
    LOGICAL, INTENT(IN) :: StateExists, StateOpenSuccess, StateParseSuccess, StateLoaded
    LOGICAL, INTENT(IN) :: Initialized, PowerInitialized
    REAL(KIND=dp), INTENT(IN) :: Temperature, Current, Resistance, Power, PreviousCurrent
    REAL(KIND=dp), INTENT(IN) :: AverageTemperature, SweepTemperature, Omega, OmegaCap, PrevResidual, LastDt
    REAL(KIND=dp), INTENT(IN) :: PhysicalTime, Dt
    INTEGER, INTENT(IN) :: LastStep, IterInStep, CallsInStep, TimeStep, BdfOrder
    INTEGER, SAVE :: Unit = 0
    INTEGER :: IoStatus
    LOGICAL, SAVE :: Started = .FALSE.

    IF (ParEnv % MyPE /= 0) RETURN
    IF (.NOT. Started) THEN
      Unit = 95
      OPEN(UNIT=Unit, FILE=TRIM(Prefix)//'.csv', STATUS='REPLACE', ACTION='WRITE', IOSTAT=IoStatus)
      IF (IoStatus /= 0) CALL Fatal('HeatSolve', 'Phase24 restart audit: cannot open audit output')
      WRITE(Unit,'(A)') 'stage,physical_time,timestep,dt,bdf_order,state_file,state_exists,state_open_success,state_parse_success,state_loaded,current,previous_current,resistance,power,average_temperature,sweep_temperature,omega,omega_cap,prev_residual,last_circuit_dt,last_circuit_step,circuit_iter_in_step,circuit_calls_in_step,circuit_initialized,circuit_power_initialized,source'
      Started = .TRUE.
    ELSE
      OPEN(UNIT=Unit, FILE=TRIM(Prefix)//'.csv', STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=IoStatus)
      IF (IoStatus /= 0) CALL Fatal('HeatSolve', 'Phase24 restart audit: cannot append audit output')
    END IF
    WRITE(Unit,'(A,A,ES24.16,A,I0,A,ES24.16,A,I0,A,A,A,L1,A,L1,A,L1,A,L1,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,I0,A,I0,A,I0,A,L1,A,L1,A,A)') &
      TRIM(Stage),',',PhysicalTime,',',TimeStep,',',Dt,',',BdfOrder,',',TRIM(StateFile),',', &
      StateExists,',',StateOpenSuccess,',',StateParseSuccess,',',StateLoaded,',',Current,',',PreviousCurrent,',', &
      Resistance,',',Power,',',AverageTemperature,',',SweepTemperature,',',Omega,',',OmegaCap,',',PrevResidual,',', &
      LastDt,',',LastStep,',',IterInStep,',',CallsInStep,',',Initialized,',',PowerInitialized,',',TRIM(Source)
    CLOSE(Unit)
  END SUBROUTINE Phase24WriteRestartAudit


!------------------------------------------------------------------------------
! To save some time tabulate data needed for the diffuse gray radiation. 
!------------------------------------------------------------------------------
  SUBROUTINE TabulateBoundaryAverages( Mesh, Emiss, Reflect )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), ALLOCATABLE :: Emiss(:)
    REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: Reflect(:)
 !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: bindex, nb, n, j, noactive
    INTEGER :: ElemInds(12)
    REAL(KIND=dp) :: NodalVal(12), Absorp
    
     nb = Mesh % NumberOfBoundaryElements
     NoActive = 0
     
     DO j=1,nb
!      bindex = j + Mesh % NumberOfBulkElements
!      Element => Mesh % Elements(bindex)
       Element => GetBoundaryElement(j)

       BC => GetBC(Element)
       IF(.NOT. ASSOCIATED( BC ) ) CYCLE

       IF( .NOT. ListCheckPresent( BC,'Radiation' ) ) CYCLE

       IF(.NOT. ALLOCATED( Emiss ) ) THEN
         ALLOCATE( Emiss(nb) )
         Emiss = 0.0_dp
         IF( PRESENT( Reflect ) ) THEN
           ALLOCATE( Reflect(nb) )
           Reflect = 0.0_dp
         END IF
       END IF

       n = GetElementNOFNodes(Element)

       NodalVal(1:n) = GetReal(BC,'Emissivity',Found)
       IF (Found) THEN
         Emiss(j) = SUM(NodalVal(1:n)) / n
         IF( PRESENT(Reflect)) THEN
           NodalVal(1:n) = GetReal(BC,'Absorptivity',Found)
           IF(Found) THEN
             Absorp = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp = Emiss(j)
           END IF
           NodalVal(1:n) = GetReal(BC,'Reflectivity',Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp
           END IF
         END IF
       ELSE
         NodalVal(1:n) = GetParentMatProp('Emissivity',Element)
         Emiss(j) = SUM(NodalVal(1:n)) / n
         IF( PRESENT( Reflect ) ) THEN
           NodalVal(1:n) = GetParentMatProp('Absorptivity',Element, Found)
           IF(Found) THEN
             Absorp = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp = Emiss(j)
           END IF
           NodalVal(1:n) = GetParentMatProp('Reflectivity',Element, Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp
           END IF
         END IF
       END IF
     END DO

   END SUBROUTINE TabulateBoundaryAverages
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE AddHeatFluxBC()
!------------------------------------------------------------------------------
      LOGICAL :: ElementForThisPartition
      CALL GetElementNodes( ElementNodes )

      HeatTransferCoeff = 0.0_dp
      LOAD  = 0.0_dp
      AText = 0.0_dp
!------------------------------------------------------------------------------
!     BC: -k@T/@n = \epsilon\sigma(T^4 - Text^4)
!------------------------------------------------------------------------------
      RadiationFlag = GetString( BC, 'Radiation', Found )

      IF ( Found .AND. RadiationFlag /= 'none' ) THEN
        Emissivity = Emiss(bc_elem)

!------------------------------------------------------------------------------
        IsRadiosity = .FALSE.

        ElementForThisPartition = Element % PartIndex == ParEnv % myPE

        IF (  RadiationFlag == 'idealized' ) THEN
          AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
          IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
        ELSE IF (ElementForThisPartition) THEN
          IF( Radiosity ) THEN
            CALL RadiosityRadiation( Model, Solver, Element, & 
                n, Temperature, TempPerm, ForceVector )
            IsRadiosity = .TRUE.
          ELSE          
            CALL DiffuseGrayRadiation( Model, Solver, Element, & 
                Temperature, TempPerm, ForceVector, VisibleFraction, Text)
            
            IF( GetLogical( BC, 'Radiation Boundary Open', Found) ) THEN
              AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
              IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
              IF( VisibleFraction >= 1.0_dp ) THEN
                Atext(1:n) = Text
              ELSE
                Atext(1:n) = ( (1 - VisibleFraction) * Atext(1:n)**4 + &
                    VisibleFraction * Text**4 ) ** 0.25_dp
              END IF
            ELSE
              AText(1:n) = Text
            END IF
          END IF
        END IF
!------------------------------------------------------------------------------
!       Add our own contribution to surface temperature (and external
!       if using linear type iteration or idealized radiation)
!------------------------------------------------------------------------------
        IF(.NOT.IsRadiosity.AND.(RadiationFlag=='idealized'.OR.ElementForThisPartition)) THEN
          DO j=1,n
            k = TempPerm(Element % NodeIndexes(j))
            Text = AText(j)
            IF ( .NOT. HeatGapBC .AND. NewtonLinearization ) THEN
              HeatTransferCoeff(j) = Emissivity * 4*Temperature(k)**3 * &
                  StefanBoltzmann
              LOAD(j) = Emissivity*(3*Temperature(k)**4+Text**4) * &
                  StefanBoltzmann
            ELSE
              HeatTransferCoeff(j) = Emissivity * (Temperature(k)**3 +   &
                  Temperature(k)**2*Text+Temperature(k)*Text**2 + Text**3) * &
                  StefanBoltzmann 
              LOAD(j) = HeatTransferCoeff(j) * Text
            END IF
          END DO
        END IF
      END IF  ! of radition
!------------------------------------------------------------------------------

      Work(1:n)  = GetReal( BC, 'Heat Transfer Coefficient',Found )
      IF ( Found ) THEN
       AText(1:n) = GetReal( BC, 'External Temperature',Found )
       DO j=1,n
!------------------------------------------------------------------------------
!         BC: -k@T/@n = \alpha(T - Text)
!------------------------------------------------------------------------------
          k = TempPerm(Element % NodeIndexes(j))
          LOAD(j) = LOAD(j) + Work(j) * AText(j)
          HeatTransferCoeff(j) = HeatTransferCoeff(j) + Work(j)
        END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = (rho*L)*v.n 
!     Heating related to pulling is possible only in ss cases where pull velocity
!     is desrcibed.
!------------------------------------------------------------------------------

      IF( GetLogical( BC, 'Phase Change',Found ) ) THEN
         PhaseVelocity(1,1:n) = GetReal( BC,'Phase Velocity 1', Found  )
         PhaseVelocity(2,1:n) = GetReal( BC,'Phase Velocity 2', Found  )
         PhaseVelocity(3,1:n) = GetReal( BC,'Phase Velocity 3', Found  )
  
         ! Ensure that the latent heat and density come from the same side
         LatentHeat(1:n) = GetParentMatProp( 'Latent Heat', &
              UElement = Element, UParent = Parent )
         IF(.NOT. ASSOCIATED(Parent) ) THEN
           CALL Warn('HeatSolve','Parent not associated')
         ELSE
           k = GetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
           Density(1:n) = GetReal( Model % Materials(k) % Values, 'Density' )
         END IF

         ! This could be rather put as a new type of BC into the assembly routine and 
         ! then the Normal could be taken at the proper Gaussian integration points. 
         Normal = NormalVector( Element, ElementNodes, 0.0_dp, 0.0_dp, .TRUE. )

         DO i=1,n
            LOAD(i) = LOAD(i) + &
                 LatentHeat(i) * Density(i) * SUM( Normal(1:3) * PhaseVelocity(1:3,i))
         END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = g
!------------------------------------------------------------------------------
      LOAD(1:n) = LOAD(1:n) +  GetReal( BC, 'Heat Flux', Found )

      InfBC = ListGetLogical( BC,'Infinity BC '//TRIM(VarName),GotIt)
      IF( InfBC ) THEN
        AText(1:n) = GetReal( BC,'Infinity BC '//TRIM(VarName)//' Offset',GotIt)
        ! currently only isotropic heat conductivity supported
        HeatConductivityIso(1:n) = GetParentMatProp('Heat Conductivity',Element,GotIt)
        IF(.NOT. GotIt) THEN
          CALL Fatal( 'HeatSolver','Could not find > Heat Conductivity < for parent!' )           
        END IF
      END IF


!------------------------------------------------------------------------------
!     Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() == Cartesian ) THEN
        CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
            LOAD,HeatTransferCoeff,InfBC,HeatConductivityIso,AText(1:n),&
            Element,n,ElementNodes )
      ELSE
        IF( InfBC ) THEN
          CALL Fatal('HeatSolver','Infinity BC not implemented only for cartersian case!')
        END IF
        CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
            LOAD,HeatTransferCoeff,Element,n,ElementNodes ) 
      END IF

!------------------------------------------------------------------------------
!     Update global matrices from local matrices
!------------------------------------------------------------------------------
      IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
        MASS = 0.d0
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF

      IF ( HeatGapBC ) &
        CALL AddHeatGap( Solver, Element, STIFF, TempPerm)

      CALL DefaultUpdateEquations( STIFF, FORCE, VecAssembly=Phase24VectorAssembly, &
           Phase24PositionCache=Phase24VectorAssembly, RHSOnly=Phase24RHSOnly )
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatFluxBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddGlobalTime()
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,n
      REAL(KIND=dp) :: FORCE(1)
      REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
      SAVE STIFF, MASS, X
      REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),MASS(:,:), X(:,:)

      IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
         IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,X )
         n = 0
         DO i=1,Solver % Matrix % NumberOfRows
           n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
         END DO
         k = SIZE(Solver % Variable % PrevValues,2)
         ALLOCATE( STIFF(1,n),MASS(1,n),X(n,k) )
 
         SaveValues => Solver % Variable % Values
      END IF

      DO i=1,Solver % Matrix % NumberOFRows
        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
          n=n+1
          STIFF(1,n) = Solver % Matrix % Values(j)
          MASS(1,n)  = Solver % Matrix % MassValues(j)
          X(n,:) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),:)
        END DO
        FORCE(1) = Solver % Matrix % RHS(i)
        Solver % Matrix % Force(i,1) = FORCE(1)
        k = MIN( Solver % DoneTime, Solver % Order )
        CALL BDFLocal( n, dt, MASS, STIFF, FORCE, X, k )

        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
           n=n+1
          Solver % Matrix % Values(j) = STIFF(1,n)
        END DO
        Solver % Matrix % RHS(i) = FORCE(1)
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddGlobalTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE DiffuseGrayRadiation( Model, Solver, Element,  &
      Temperature, TempPerm, ForceVector,AngleFraction, Text)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: TempPerm(:)
      REAL(KIND=dp) :: Temperature(:), ForceVector(:)
      REAL(KIND=dp) :: AngleFraction, Text
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Area, Asum, gEmissivity, Base(12), Load, TransCoeff
      REAL(KIND=dp), POINTER :: Fact(:)
      INTEGER :: i,j,k,l,m,ImplicitFactors, nf,nr, bindex, nb
      INTEGER, POINTER :: ElementList(:)
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------

      Asum = 0.0_dp
      IF ( .NOT. NewtonLinearization ) THEN
        Text = ComputeRadiationLoad( Model, Solver % Mesh, Element, &
           Temperature, TempPerm, Emissivity, AngleFraction, Areas, Emiss )

      ELSE   !  Full Newton-Raphson solver
!------------------------------------------------------------------------------
!       Go through surfaces (j) this surface (i) is getting
!       radiated from.
!------------------------------------------------------------------------------

        ElementList => Element % BoundaryInfo % RadiationFactors % Elements

        nf = Element % BoundaryInfo % RadiationFactors % NumberOfFactors

        bindex = Element % ElementIndex - Solver % Mesh % NumberOfBulkElements
        Area  = Areas(bindex)
        CALL GetBase( Base, Element, n, ElementNodes )

        Fact => Element % BoundaryInfo % RadiationFactors % Factors

        DO j=1,nf
          RadiationElement => Solver % Mesh % Elements(ElementList(j))
          Text = Fact(j)
          Asum = Asum + Text

!------------------------------------------------------------------------------
!         Gebhart factors are given elementwise at the center
!         of the element, so take average of nodal temperatures
!         (or integrate over surface j)
!------------------------------------------------------------------------------

          k = RadiationElement % TYPE % NumberOfNodes
          ImplicitFactors = Element % BoundaryInfo % RadiationFactors % NumberOfImplicitFactors
          IF(ImplicitFactors == 0) &
              ImplicitFactors = Element % BoundaryInfo % RadiationFactors % NumberOfFactors

          IF(j <= ImplicitFactors) THEN
            
            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )**(1._dp/4._dp)
!------------------------------------------------------------------------------
!          Linearization of the G_jiT^4_j term
!------------------------------------------------------------------------------
           LOAD = -3 * Text * S**4 * StefanBoltzmann
           TransCoeff = -4 * Text * S**3 * StefanBoltzmann
!------------------------------------------------------------------------------
!          Integrate the contribution of surface j over surface i
!          and add to global matrix
!------------------------------------------------------------------------------
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              DO l=1,k
                k2 = TempPerm( RadiationElement % NodeIndexes(l) )
                CALL AddToMatrixElement( StiffMatrix,k1,k2,TransCoeff*Base(m)/k )
              END DO
              ForceVector(k1) = ForceVector(k1) + Load*Base(m)
            END DO

          ELSE
            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )
            
            LOAD = Text * S * StefanBoltzmann
            
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              ForceVector(k1) = ForceVector(k1) + LOAD*Base(m)
            END DO            
          END IF 

        END DO

!------------------------------------------------------------------------------
!       We have already added all external temperature contributions
!       to the matrix for the Newton type iteration
!------------------------------------------------------------------------------
        AngleFraction = Asum / Emissivity
        Text = 0.0

      END IF  !  of newton-raphson

    END SUBROUTINE DiffuseGrayRadiation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE RadiosityRadiation( Model, Solver, Element,  &
      n, Temperature, TempPerm, ForceVector)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: n
      INTEGER :: TempPerm(:)
      REAL(KIND=dp) :: Temperature(:), ForceVector(:)      
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Emis1, Refl1, RadCoeffAtIp, RadLoadAtIp, TempAtIp, s, x, y, z, DetJ
      REAL(KIND=dp) :: Basis(27)
      REAL(KIND=dp), POINTER :: Fact(:)
      INTEGER :: t,p,q,k1,k2
      INTEGER, POINTER :: pIndexes(:)
      TYPE(GaussIntegrationPoints_t), TARGET :: IP
      LOGICAL :: stat
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------
      
      pIndexes => Element % NodeIndexes     
      
      Emis1 = Emissivity            
      IF(.NOT.Spectral ) THEN
        Refl1 = Reflect(bc_elem)
        Emis1 = Emis1 / Refl1
      END IF
      Fact => Element % BoundaryInfo % RadiationFactors % Factors
      
      TempAtIp = SUM(Temperature(TempPerm(pIndexes(1:n))))/n

      IF(NewtonLinearization) THEN
        RadLoadAtIp =  3 * Emis1 * TempAtIp**4 * StefanBoltzmann + &
             (Fact(1) - Fact(2)*TempAtIp)
        RadCoeffAtIp = 4 * Emis1 * TempAtIp**3 * StefanBoltzmann - Fact(2)
      ELSE
        RadCoeffAtIp = Emis1 * StefanBoltzmann * TempAtIp**3
        RadLoadAtIp = Fact(1)
      END IF

      IP = GaussPoints( Element )
      
      DO t=1,IP % n
        stat = ElementInfo( Element,ElementNodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
        s = detJ * IP % s(t)
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( ElementNodes % x(1:n)*Basis )
          y = SUM( ElementNodes % y(1:n)*Basis )
          z = SUM( ElementNodes % z(1:n)*Basis )
          s = s * CoordinateSqrtMetric( x,y,z )
        END IF
        DO p=1,n
          k1 = TempPerm( pIndexes(p) )
          DO q=1,n
            k2 = TempPerm( pIndexes(q) )
            CALL AddToMatrixElement( Solver % Matrix,k1,k2,s*Basis(p)*Basis(q)*RadCoeffAtIP)
          END DO
          ForceVector(k1) = ForceVector(k1) + s *Basis(p) * RadLoadAtIp
        END DO
      END DO
      
    END SUBROUTINE RadiosityRadiation
!------------------------------------------------------------------------------


    
!------------------------------------------------------------------------------
    SUBROUTINE EffectiveHeatCapacity()
      LOGICAL :: Found, Specific, GotFraction
      REAL(KIND=dp), ALLOCATABLE :: dT(:)
      REAL(KIND=dp) :: dT0

!------------------------------------------------------------------------------
!     See if temperature gradient indside the element is large enough 
!     to use  the c_p = SQRT( (dH/dx)^2 / (dT/dx)^2 ), otherwise
!     use c_p = dH/dT, or if in time dependent simulation, use
!     c_p = (dH/dt) / (dT/dt), if requested. 
!------------------------------------------------------------------------------

      SELECT CASE(PhaseModel)
!------------------------------------------------------------------------------
        CASE( 'spatial 1' )
          PhaseChangeModel = PHASE_SPATIAL_1
!------------------------------------------------------------------------------

        CASE( 'spatial 2' )
!------------------------------------------------------------------------------
! Check if local variation of temperature is large enough to actually use the
! Spatial 2 model. Should perhaps be scaled to element size (or actually
! compute the gradient, but this will do for now...).
!------------------------------------------------------------------------------
          s = MAXVAL(LocalTemperature(1:n))-MINVAL(LocalTemperature(1:n))
          IF ( s < AEPS ) THEN
            PhaseChangeModel = PHASE_SPATIAL_1
          ELSE
            PhaseChangeModel = PHASE_SPATIAL_2
          END IF

!------------------------------------------------------------------------------
! Note that here HeatCapacity is miused for saving dT.
!------------------------------------------------------------------------------
        CASE('temporal')
          IF ( TransientSimulation )  THEN
            ALLOCATE( dT(n) )
            dT(1:n) = Temperature(TempPerm(Element % NodeIndexes)) - &
                     PrevTemperature(TempPerm(Element % NodeIndexes))

            IF ( ANY(ABS(dT(1:n)) < AEPS) ) THEN
              PhaseChangeModel = PHASE_SPATIAL_1
            ELSE
              PhaseChangeModel = PHASE_TEMPORAL
            END IF
          ELSE
             PhaseChangeModel = PHASE_SPATIAL_1
          END IF

!------------------------------------------------------------------------------
        CASE DEFAULT
          PhaseChangeModel = PHASE_SPATIAL_1

      END SELECT
!------------------------------------------------------------------------------

      PhaseSpatial = ( PhaseChangeModel == PHASE_SPATIAL_2 )
      Specific = ListCheckPresent( Material,'Specific Enthalpy')

      EnthalpyFraction(1:n) = ListGetReal(Material,'Enthalpy Fraction',&
          n,Element % NodeIndexes,GotFraction)          
 

!-----------------------------------------------------------------------------
      SELECT CASE( PhaseChangeModel )

!------------------------------------------------------------------------------
! This phase change model is available only for some type of real entries 
! that have an implemented analytical derivation rule.
!-----------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_1 )

        Work(1:n) = ListGetReal( Material, &
            'Effective Heat Capacity', n,Element % NodeIndexes, Found )
        IF ( .NOT. Found ) THEN
          dT0 = ListGetCReal( Material,'Enthalpy Temperature Differential',Found )
          IF(.NOT. Found) dT0 = 1.0d-3
          IF( Specific ) THEN
            Work(1:n) = ListGetDerivValue( Material, &
                'Specific Enthalpy', n,Element % NodeIndexes, dT0 )
            Work(1:n) = Density(1:n) * Work(1:n)
          ELSE
            Work(1:n) = ListGetDerivValue( Material, &
                'Enthalpy', n,Element % NodeIndexes, dT0 )
          END IF
        END IF

        IF( GotFraction ) THEN
          HeatCapacity(1:n) = HeatCapacity(1:n) + EnthalpyFraction(1:n) * Work(1:n) 
        ELSE
          HeatCapacity(1:n) = HeatCapacity(1:n) + Work(1:n) 
        END IF
          
!---------------------------------------------------------------------------------------
! Note that for the 'spatial 2' model the evaluation of c_p is done in each integration
! point and thus Enthalphy and PhaseSpatial flag are used instead of HeatCapacity directly.
!-----------------------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_2 )
        IF( Specific ) THEN
          Enthalpy(1:n) = ListGetReal(Material,'Specific Enthalpy',n,Element % NodeIndexes)
          Enthalpy(1:n) = Density(1:n) * Enthalpy(1:n)
        ELSE
          Enthalpy(1:n) = ListGetReal(Material,'Enthalpy',n,Element % NodeIndexes)          
        END IF
        
        IF( GotFraction ) THEN
          CALL Warn('EffectiveHeatCapacity',&
              '> Enthalpy Fraction < not treated yet by spatial 2 phase change')
        END IF

          
!------------------------------------------------------------------------------
      CASE( PHASE_TEMPORAL )
        ! When retrieving the value of enthalphy on the previous timestep 
        ! the relevant entries of the Temperature solution in the global vector
        ! are tampered in order to make the ListGetReal command work as wanted. 
        ! 1) Values at current temperature     
        !------------------------------------------------------------------------
        IF( Specific ) THEN
          Work(1:n) = ListGetReal( Material,'Specific Enthalpy',n,Element % NodeIndexes )
        ELSE
          Work(1:n) = ListGetReal( Material,'Enthalpy',n,Element % NodeIndexes )
        END IF

        ! 2) Values at previous temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) 

        IF( Specific ) THEN
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Specific Enthalpy', &
              n,Element % NodeIndexes )          
          Work(1:n) = Density(1:n) * Work(1:n) / dT(1:n)
        ELSE
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Enthalpy', &
              n,Element % NodeIndexes )
          Work(1:n) = Work(1:n) / dT(1:n)
        END IF

        IF( GotFraction ) THEN
          HeatCapacity(1:n) = HeatCapacity(1:n) + EnthalpyFraction(1:n) * Work(1:n) 
        ELSE
          HeatCapacity(1:n) = HeatCapacity(1:n) + Work(1:n) 
        END IF


        ! Revert to current temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) + dT(1:n)        

!------------------------------------------------------------------------------
      END SELECT

!------------------------------------------------------------------------------
    END SUBROUTINE EffectiveHeatCapacity
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION CheckLatentHeat() RESULT(Failure)
!------------------------------------------------------------------------------
      LOGICAL :: Failure, PhaseChange, CheckLatentHeatRelease
      INTEGER :: t, eq_id, body_id
      CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel
!------------------------------------------------------------------------------

      Failure = .FALSE.
!------------------------------------------------------------------------------
      DO t=1,Solver % Mesh % NumberOfBulkElements
!------------------------------------------------------------------------------
!       Check if this element belongs to a body where temperature 
!       has been calculated
!------------------------------------------------------------------------------
        Element => Solver % Mesh % Elements(t)

        NodeIndexes => Element % NodeIndexes
        IF ( ANY( TempPerm( NodeIndexes ) <= 0 ) ) CYCLE

        body_id = Element % Bodyid
        eq_id = ListGetInteger( Model % Bodies(body_id) % Values, &
            'Equation', minv=1, maxv=Model % NumberOfEquations )

        PhaseModel = ListGetString( Model % Equations(eq_id) % Values, &
                          'Phase Change Model',Found )

        PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')

        IF ( PhaseChange ) THEN
          CheckLatentHeatRelease = ListGetLogical(Model % Equations(eq_id) % &
                    Values, 'Check Latent Heat Release',Found )
        END IF
        IF ( .NOT. ( PhaseChange .AND. CheckLatentHeatRelease ) ) CYCLE

        n = Element % TYPE % NumberOfNodes
!------------------------------------------------------------------------------
!       Set the current element pointer in the model structure to
!       reflect the element being processed
!------------------------------------------------------------------------------
        Model % CurrentElement => Element
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!       Get element material parameters
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies(body_id) % Values,'Material', &
                minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(k) % Values

        PhaseChangeIntervals => ListGetConstRealArray( Material, &
                        'Phase Change Intervals' )

        DO k=1,n
          i = TempPerm( NodeIndexes(k) )
          DO j=1,SIZE(PhaseChangeIntervals,2)
            IF ( ( Temperature(i)  < PhaseChangeIntervals(1,j) .AND. &
                   PrevSolution(i) > PhaseChangeIntervals(2,j) ).OR. &
                 ( Temperature(i)  > PhaseChangeIntervals(2,j) .AND. &
                   PrevSolution(i) < PhaseChangeIntervals(1,j) )  ) THEN
              Failure = .TRUE.
              EXIT
            END IF
          END DO
          IF ( Failure ) EXIT
        END DO
        IF ( Failure ) EXIT
      END DO
!------------------------------------------------------------------------------
    END FUNCTION CheckLatentHeat
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
   SUBROUTINE GetBase( Base, Element, n, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Base(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n,  m

     REAL(KIND=dp) :: Basis(n), DetJ

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     Base = 0._dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
!------------------------------------------------------------------------------
       DO p=1,N
         Base(p) = Base(p) + s * Basis(p)
       END DO
     END DO
   END SUBROUTINE GetBase
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE IntegOverA( BoundaryMatrix, BoundaryVector, &
     LOAD, NodalAlpha, Element, n, m, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), LOAD,NodalAlpha

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n,  m

     REAL(KIND=dp) :: Basis(n), DetJ

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0_dp
     BoundaryMatrix = 0.0_dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     Force = LOAD
     Alpha = NodalAlpha / m

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
!------------------------------------------------------------------------------
!      Force = SUM( LOAD(1:n) * Basis )
!      Alpha = SUM( NodalAlpha(1:n) * Basis )

       DO p=1,N
         DO q=1,M
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + s * Alpha * Basis(p)
         END DO
       END DO

       DO p=1,N
         BoundaryVector(p) = BoundaryVector(p) + s * Force * Basis(p)
       END DO
     END DO
   END SUBROUTINE IntegOverA
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindGapIndexes( Element, Indexes, n )
!------------------------------------------------------------------------------
      TYPE(Element_t) :: Element
      INTEGER :: n,Indexes(:)
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      Left  => Element % BoundaryInfo % Left
      Right => Element % BoundaryInfo % Right

      IF ( .NOT.ASSOCIATED(Left) .OR. .NOT.ASSOCIATED(Right) ) RETURN

      l = 0
      DO i=1,n
        Parent => Left
        k = Element % NodeIndexes(i)

        IF ( ANY( Parent % NodeIndexes == k ) ) &
          Parent => Right

        x0 = ElementNodes % x(i)
        y0 = ElementNodes % y(i)
        z0 = ElementNodes % z(i)
        DO j=1,Parent % TYPE % NumberOfNodes
          k = Parent % NodeIndexes(j)
          x = Solver % Mesh % Nodes % x(k) - x0
          y = Solver % Mesh % Nodes % y(k) - y0
          z = Solver % Mesh % Nodes % z(k) - z0
          IF ( x**2 + y**2 + z**2 < AEPS ) EXIT
        END DO
        Indexes(i) = k
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE FindGapIndexes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddHeatGap( Solver, Element, STIFF, TempPerm )
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: TempPerm(:)
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l, Ind(n)
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      CALL FindGapIndexes( Element, Ind, n )
      DO i=1,n
        DO j=1,n
          k = TempPerm( Element % NodeIndexes(i) )
          l = TempPerm( Ind(j) )
          IF ( k > 0 .AND. l > 0 ) THEN
            CALL AddToMatrixElement( Solver % Matrix,k,l,-STIFF(i,j) )
          END IF
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatGap
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE HeatSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE HeatSolver_Boundary_Residual( Model, Edge, Mesh, Quant, Perm,Gnorm, Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Radiation

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,n,l,t,dim,Pn,En,nd
     LOGICAL :: stat, Found
     INTEGER, ALLOCATABLE :: Indexes(:)

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), ExtTemperature(:), &
       TransferCoeff(:), EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
       dBasisdx(:,:), Temperature(:), Flux(:), NodalEmissivity(:)

     REAL(KIND=dp) :: Conductivity, Emissivity, StefanBoltzmann

     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, gx, gy, gz

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE., Dirichlet
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Indicator = 0.0d0
     Gnorm     = 0.0d0

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT
!    
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     en = Edge % TYPE % NumberOfNodes
     pn = Element % TYPE % NumberOfNodes

     ALLOCATE(EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     nd = GetElementNOFDOFs(Element)
     ALLOCATE( Temperature(nd), Basis(nd), ExtTemperature(nd), &
        TransferCoeff(en), x(en), y(en), z(en), EdgeBasis(nd), &
        dBasisdx(nd,3), NodalConductivity(nd), Flux(nd), &
        NodalEmissivity(nd), Indexes(nd) ) 

     nd = GetElementDOFs(Indexes,Element)

     ALLOCATE(Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )
     Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))


     DO l = 1,en
       DO k = 1,pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO
!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------

     Indicator    = 0.0d0
     EdgeLength   = 0.0d0
     ResidualNorm = 0.0d0

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                 'Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        Dirichlet = ListCheckPresent( Model % BCs(j) % Values,'Temperature')       
        
!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:en) = ListGetReal( Model % BCs(j) % Values, &
          'Heat Flux', en, Edge % NodeIndexes, Found )

!       ...convective heat transfer:
!       ----------------------------
        TransferCoeff(1:en) =  ListGetReal( Model % BCs(j) % Values, &
          'Heat Transfer Coefficient', en, Edge % NodeIndexes, Found )

        ExtTemperature(1:en) = ListGetReal( Model % BCs(j) % Values, &
          'External Temperature', en, Edge % NodeIndexes, Found )

!       ...black body radiation:
!       ------------------------
        Emissivity      = 0.0d0
        StefanBoltzmann = 0.0d0

        SELECT CASE(ListGetString(Model % BCs(j) % Values,'Radiation',Found))
           !------------------
           CASE( 'idealized' )
           !------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                  ListGetConstReal( Model % Constants,'Stefan Boltzmann',UnfoundFatal=.TRUE. )

           !---------------------
           CASE( 'diffuse gray' )
           !---------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

              ExtTemperature(1:en) =  ComputeRadiationLoad( Model, &
                      Mesh, Edge, Quant, Perm, Emissivity )
        END SELECT

!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Heat Conductivity', Hwrk, en, Edge % NodeIndexes )
        NodalConductivity(1:en) = Hwrk(1,1,1:en)

!       elementwise nodal solution:
!       ---------------------------
        nd = GetElementDOFs(Indexes,Element)
        Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0

        IntegStuff = GaussPoints(Edge)

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
                    EdgeBasis, dBasisdx )
           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ
           ELSE
              gx = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
              gy = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
              gz = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, gx, gy, gz )

              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:en) * x(1:en) )
           v = SUM( EdgeBasis(1:en) * y(1:en) )
           w = SUM( EdgeBasis(1:en) * z(1:en) )
           stat = ElementInfo(Element,Nodes, u, v, w, detJ,Basis,dBasisdx )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Conductivity = SUM( NodalConductivity(1:en) * EdgeBasis(1:en) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:en) * EdgeBasis(1:en) )

!          convective ...:
!          ----------------
           Residual = Residual + SUM(TransferCoeff(1:en) * EdgeBasis(1:en)) * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) )

!          black body radiation...:
!          -------------------------
           Residual = Residual + &
                Emissivity * StefanBoltzmann * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) ** 4 - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) ** 4 )

!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,dim
                 Residual = Residual + Conductivity  * &
                    SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(k)

                 Gnorm = Gnorm + s * (Conductivity * &
                       SUM(dBasisdx(1:nd,k) * Temperature(1:nd)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,dim
                 DO l=1,dim
                    Residual = Residual + Metric(k,l) * Conductivity  * &
                       SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
                      SUM(dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l))**2
                 END DO
              END DO
           END IF

           EdgeLength   = EdgeLength + s
           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * Residual ** 2
           END IF
        END DO
        EXIT
     END DO

     IF ( CoordinateSystemDimension() == 3 ) EdgeLength = SQRT(EdgeLength)

!    Gnorm = EdgeLength * Gnorm
     Indicator = EdgeLength * ResidualNorm
!------------------------------------------------------------------------------
   END SUBROUTINE HeatSolver_Boundary_Residual
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE HeatSolver_Edge_Residual(Model,Edge,Mesh,Quant,Perm,Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,dim,En,Pn,nd
     INTEGER, ALLOCATABLE :: Indexes(:)
     LOGICAL :: stat, Found
     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), &
            EdgeBasis(:), Basis(:), dBasisdx(:,:), Temperature(:)

     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump, Conductivity

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Residual, ResidualNorm, Area

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

     !    Initialize:
     !    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
     END DO

     Grad = 0.0d0
!
!    ---------------------------------------------

     n = Mesh % MaxElementDOFs
     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     en = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( NodalConductivity(en), EdgeBasis(en), Basis(n), &
        dBasisdx(n,3), x(en), y(en), z(en), Temperature(n), Indexes(n) )

!    Integrate square of jump over edge:
!    -----------------------------------
     ResidualNorm = 0.0d0
     EdgeLength   = 0.0d0
     Indicator    = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
           v = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
           w = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! 
        ! Compute flux over the edge as seen by elements
        ! on both sides of the edge:
        ! ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           pn = Element % TYPE % NumberOfNodes

           DO j = 1,en
              DO k = 1,pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM(EdgeBasis(1:en) * x(1:en))
           v = SUM(EdgeBasis(1:en) * y(1:en))
           w = SUM(EdgeBasis(1:en) * z(1:en))
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           nd = GetElementDOFs(Indexes,Element)
           Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
           Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
           Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))

           stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx)
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                   'Heat Conductivity', Hwrk,en, Edge % NodeIndexes )

           NodalConductivity(1:en) = Hwrk( 1,1,1:en )
           Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
!
!          Temperature at element nodal points:
!          ------------------------------------
           Temperature(1:nd) = Quant( Perm(Indexes(1:nd)) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,dim
              Grad(j,i) = Conductivity * SUM( dBasisdx(1:nd,j) * Temperature(1:nd) )
           END DO
        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0d0
        DO k=1,dim
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,dim
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
     END DO

     IF (dim==3) EdgeLength = SQRT(EdgeLength)
     Indicator = EdgeLength * ResidualNorm

!------------------------------------------------------------------------------
   END SUBROUTINE HeatSolver_Edge_Residual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE HeatSolver_Inside_Residual( Model, Element, Mesh, &
        Quant, Perm, Fnorm, Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,n,t,dim,nd
     INTEGER, ALLOCATABLE :: Indexes(:)

     LOGICAL :: stat, Found, Compressible, VolSource
     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalCapacity(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
     REAL(KIND=dp), ALLOCATABLE :: Velo(:,:), Pressure(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Temperature(:), PrevTemp(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:)

     REAL(KIND=dp) :: u, v, w, s, detJ, Density, Capacity

     REAL(KIND=dp) :: SpecificHeatRatio, ReferencePressure, dt
     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area, Conductivity

     TYPE( ValueList_t ), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Fnorm     = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
     IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN

     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

!    Alllocate local arrays
!    ----------------------
     nd = GetElementNOFDOFs(Element)
     n = GetElementNOFNodes(Element)
     ALLOCATE( NodalDensity(nd), NodalCapacity(nd), NodalConductivity(nd),      &
         Velo(3,nd), Pressure(nd), NodalSource(nd), Temperature(nd), PrevTemp(nd), &
         Basis(nd), dBasisdx(nd,3), ddBasisddx(nd,3,3), Indexes(nd) )
!
!    Element nodal points:
!    ---------------------
     ALLOCATE( Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )

     nd = GetElementDOFs(Indexes,Element)
     Nodes % x = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z = Mesh % Nodes % z(Indexes(1:nd))
!
!    Elementwise nodal solution:
!    ---------------------------
     Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))
!
!    Check for time dep.
!    -------------------
     PrevTemp(1:nd) = Temperature(1:nd)
     dt = Model % Solver % dt
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Temperature', .TRUE. )
        PrevTemp(1:nd) = Var % PrevValues(Var % Perm(Indexes(1:nd)),1)
     END IF
!
!    Material parameters: conductivity, heat capacity and density
!    -------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL ListGetRealArray(Material,'Heat Conductivity',Hwrk,n,Element % NodeIndexes)
     NodalConductivity(1:n) = Hwrk( 1,1,1:n)
     NodalDensity(1:n) = GetReal(Material, 'Density', Found )
     NodalCapacity(1:n) = GetReal(Material, 'Heat Capacity', Found )
!
!    Check for compressible flow equations:
!    --------------------------------------
     Compressible = .FALSE.

     IF (  ListGetString( Material, 'Compressibility Model', Found ) == &
                 'perfect gas equation 1' ) THEN

        Compressible = .TRUE.

        Pressure = 0.0d0
        Var => VariableGet( Mesh % Variables, 'Pressure', .TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Pressure(1:n) = Var % Values(Var % Perm(Indexes(1:n)))
        END IF

        ReferencePressure = GetConstReal( Material, 'Reference Pressure' )
        SpecificHeatRatio = GetConstReal( Material, 'Specific Heat Ratio' )
        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * NodalCapacity(1:n) * Temperature(1:n) )
     END IF
!
!    Get (possible) convection velocity at the nodes of the element:
!    ----------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                minv=1, maxv=Model % NumberOFEquations )

     Velo = 0.0d0
     SELECT CASE( ListGetString( Model % Equations(k) % Values, &
                         'Convection', Found ) )

        !-----------------
        CASE( 'constant' )
        !-----------------

           Velo(1,1:n) = GetReal( Material, 'Convection Velocity 1',Found )
           Velo(2,1:n) = GetReal( Material, 'Convection Velocity 2',Found )
           Velo(3,1:n) = GetReal( Material, 'Convection Velocity 3',Found )

        !-----------------
        CASE( 'computed' )
        !-----------------

           Var => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              IF ( ALL(Var % Perm(Indexes(1:n)) > 0 ) ) THEN
                 Velo(1,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
                 Var => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(2,1:n) = Var % Values(Var % Perm(Indexes(1:n)) )
                 Var => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(3,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
              END IF
           END IF

     END SELECT

!
!    Heat source:
!    ------------
!
     k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',Found, &
                 1, Model % NumberOFBodyForces)

     NodalSource = 0.0d0
     IF( k > 0 ) THEN
       NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
           'Volumetric Heat Source',VolSource ) 
       IF( .NOT. VolSource ) THEN
         NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
             'Heat Source',  Found )
       END IF
     END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0d0
     Area = 0.0d0

     IntegStuff = GaussPoints( Element )
     ddBasisddx = 0

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM(Basis(1:nd) * Nodes % x(1:nd))
           v = SUM(Basis(1:nd) * Nodes % y(1:nd))
           w = SUM(Basis(1:nd) * Nodes % z(1:nd))

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Capacity     = SUM(NodalCapacity(1:n) * Basis(1:n))
        Density      = SUM(NodalDensity(1:n) * Basis(1:n))
        Conductivity = SUM(NodalConductivity(1:n) * Basis(1:n))
!
!       Residual of the convection-diffusion (heat) equation:
!        R = \rho * c_p * (@T/@t + u.grad(T)) - &
!            div(C grad(T)) + p div(u) - h,
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = \rho * c_p * (@T/@t + u^j T_{,j}) - &
!          g^{jk} (C T_{,j}}_{,k} + p div(u) - h
!       ---------------------------------------------------
!

        IF( VolSource ) THEN
          Residual = -SUM( NodalSource(1:n) * Basis(1:n) )
        ELSE
          Residual = -Density * SUM( NodalSource(1:n) * Basis(1:n) )
        END IF
          
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,dim
!
!             - grad(C).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                 SUM( NodalConductivity(1:n) * dBasisdx(1:n,j) )

!
!             - C div(grad(T)):
!             -------------------
!
              Residual = Residual - Conductivity * &
                 SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,j) )
           END DO
        ELSE
           DO j=1,dim
              DO k=1,dim
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Conductivity * &
                    SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,dim
                    Residual = Residual + Metric(j,k) * Conductivity * &
                      Symb(j,k,l) * SUM( Temperature(1:nd) * dBasisdx(1:nd,l) )
                 END DO
              END DO
           END DO
        END IF

!       + \rho * c_p * (@T/@t + u.grad(T)):
!       -----------------------------------
        Residual = Residual + Density * Capacity *  &
           SUM((Temperature(1:nd)-PrevTemp(1:nd))*Basis(1:nd)) / dt

        DO j=1,dim
           Residual = Residual + &
              Density * Capacity * SUM(Velo(j,1:n) * Basis(1:n)) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) )
        END DO


        IF ( Compressible ) THEN
!
!          + p div(u) or p u^j_{,j}:
!          -------------------------
!
           DO j=1,dim
              Residual = Residual + &
                 SUM( Pressure(1:n) * Basis(1:n) ) * &
                      SUM( Velo(j,1:n) * dBasisdx(1:n,j) )

              IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO k=1,dim
                    Residual = Residual + &
                       SUM( Pressure(1:n) * Basis(1:n) ) * &
                           Symb(j,k,j) * SUM( Velo(k,1:n) * Basis(1:n) )
                 END DO
              END IF
           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
        DO i=1,dim
           Fnorm = Fnorm + s * (Density *SUM(NodalSource(1:n)*Basis(1:n)))**2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
     END DO

!    Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm
!------------------------------------------------------------------------------
   END SUBROUTINE HeatSolver_Inside_Residual
!------------------------------------------------------------------------------
