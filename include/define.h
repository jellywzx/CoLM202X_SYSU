! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land subgrid type classification:
!    Select one of the following options.
#undef LULC_USGS
#undef LULC_IGBP
#undef LULC_IGBP_PFT
#define LULC_IGBP_PC

! 2.1 3D Urban model (put it temporarily here):
#undef URBAN_MODEL
!    Dependence:  only LULC_IGBP subgrid type for
!    single point URBAN_MODEL right now.
#if (defined URBAN_MODEL && defined SinglePoint)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 3. If defined, debug information is output.
#define CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#define RangeCheck
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. If defined, MPI parallelization is enabled.
#define USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef   Campbell_SOIL_MODEL
#define  vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define CatchLateralFlow
!    Conflicts :
#ifndef CATCHMENT
#undef CatchLateralFlow
#endif

! 6. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood
#if (defined SinglePoint)
#undef CaMa_Flood
#endif
#ifndef USEMPI
#undef CaMa_Flood
#endif

#define GridRiverLakeFlow
!    Conflicts :
#if (defined CATCHMENT || defined SinglePoint || defined CaMa_Flood)
#undef GridRiverLakeFlow
#endif

! 7. If defined, BGC model is used.
#undef BGC

!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#ifndef LULC_IGBP_PC
#undef BGC
#endif
#endif
! 7.1 If defined, CROP model is used
#undef CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif

! 8. If defined, open Land use and land cover change mode.
#undef LULCC

! 9. If defined, data assimilation is used.
#undef DataAssimilation
#if (defined DataAssimilation)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 10. Interface to AI model.
#undef USESplitAI

! 11. External lake models.
#undef EXTERNAL_LAKE

! 12. Hyperspectral scheme.
#undef HYPERSPECTRAL

! 12b. If defined, extended canopy interception schemes are enabled.
#define extend_interception

! 13. If defined, the tracer subsystem is enabled (isotope, solute,
!     particle, and gas families).
!     This repository template currently enables TRACER; switch to #undef
!     TRACER only for builds that intentionally exclude all tracer species.
#define TRACER

! The extended canopy implementation requires tracer phase-change outputs.
! Without TRACER, retain the standard interception and thermal interfaces.
#ifndef TRACER
#undef extend_interception
#endif
! Match the extended THERMAL dummy arguments for coupled flood evaporation.
#if defined(extend_interception) && (defined(CaMa_Flood) || defined(GridRiverLakeFlow))
#define CANOPY_FLOOD_FEEDBACK
#else
#undef CANOPY_FLOOD_FEEDBACK
#endif
!    Conflicts: TRACER requires VariablySaturatedFlow soil hydrology
!    (vanGenuchten_Mualem_SOIL_MODEL). Campbell_SOIL_MODEL cannot silently
!    disable TRACER because that changes the requested physics at compile time.
#if (defined TRACER) && (defined Campbell_SOIL_MODEL)
#error "TRACER requires vanGenuchten_Mualem_SOIL_MODEL"
#endif
#if (defined TRACER) && (defined BGC)
#if (!defined LULC_IGBP_PFT && !defined LULC_IGBP_PC)
#error "Methane requires BGC with LULC_IGBP_PFT or LULC_IGBP_PC"
#endif
#endif
