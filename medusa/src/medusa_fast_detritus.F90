#include "fabm_driver.h"

!
!**************************************************************
!            FABM-MEDUSA fast-sinking detritus mineralisation
!**************************************************************
! Note that fast-sinking detritus creation and nutrient regeneration are part of medusa_pelagic module
!
module medusa_fast_detritus

   use fabm_types
   use fabm_particle

   implicit none

   private

  type,extends(type_particle_model),public :: type_medusa_fast_detritus
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZDIC,id_ZDIN,id_ZSIL,id_ZOXY,id_ZFER,id_ZDET,id_ZDTC,id_ZALK
      type (type_bottom_state_variable_id) :: id_ZSEDSI,id_ZSEDC,id_ZSEDN,id_ZSEDCA,id_ZSEDFE, id_ZSEDP
      type (type_dependency_id)            :: id_dz
      type (type_dependency_id)            :: id_ftempc,id_ftempn,id_ftempsi,id_ftempfe,id_ftempca
      type (type_dependency_id)            :: id_freminc1,id_freminn1,id_freminsi1,id_freminfe1,id_freminca1
      type (type_dependency_id)            :: id_om_cal
      type (type_diagnostic_variable_id)   :: id_freminc,id_freminn,id_freminsi,id_freminfe,id_freminca
      type (type_diagnostic_variable_id)   :: id_ffastc_loc,id_ffastn_loc,id_ffastca_loc,id_ffastsi_loc
      type (type_horizontal_diagnostic_variable_id) :: id_ffastc,id_ffastn,id_ffastsi,id_ffastfe,id_ffastca
      type (type_horizontal_diagnostic_variable_id) :: id_SEAFLRN,id_SEAFLRSI,id_SEAFLRFE,id_SEAFLRC,id_SEAFLRCA
      type (type_horizontal_dependency_id) :: id_ffastc1,id_ffastn1,id_ffastfe1,id_ffastsi1,id_ffastca1,id_CAL_CCD
      type (type_diagnostic_variable_id)   :: id_tempc,id_tempn,id_tempsi,id_tempfe,id_tempca,id_freminn2d
      type (type_horizontal_diagnostic_variable_id) :: id_OCAL_LVL
      ! Parameters
      real(rk) :: xthetanit,xthetarem,xo2min,xrfn
      integer :: seafloor

   contains

      procedure :: initialize
      procedure :: get_light => do_fast_detritus
      procedure :: do
      procedure :: do_bottom

  end type type_medusa_fast_detritus

contains

   subroutine initialize(self,configunit)

   class(type_medusa_fast_detritus),intent(inout),target :: self
   integer,                         intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   call self%register_implemented_routines((/source_do, source_do_bottom, source_do_column/))
   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)

   ! Create diagnostics for fast-sinking detritus that acts as state variables, so they can receive sources.
   call self%register_diagnostic_variable(self%id_tempc,'tempc','mmol C m-3','fast-sinking detritus (C)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempn,'tempn','mmol N m-3','fast-sinking detritus (N)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempfe,'tempfe','-','fast-sinking detritus (Fe)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempsi,'tempsi','mmol Si m-3','fast-sinking detritus (Si)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempca,'tempca','mmol CaCO3 m-3','fast-sinking detritus (CaCO3)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)

   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_tempc)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_tempn)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_tempsi)

   ! Retrieve the sources of fast-sinking detritus.
   call self%register_dependency(self%id_ftempc, 'ftempc', 'mmol C m-3 s-1', 'production of fast-sinking detritus (C)')
   call self%register_dependency(self%id_ftempn, 'ftempn', 'mmol N m-3 s-1', 'production of fast-sinking detritus (N)')
   call self%register_dependency(self%id_ftempfe,'ftempfe', 'mmol Fe m-3 s-1', 'production of fast-sinking detritus (Fe)')
   call self%register_dependency(self%id_ftempsi,'ftempsi', 'mmol Si m-3 s-1', 'production of fast-sinking detritus (Si)')
   call self%register_dependency(self%id_ftempca,'ftempca', 'mmol CaCO3 m-3 s-1', 'production of fast-sinking detritus (CaCO3)')

   call self%register_dependency(self%id_CAL_CCD,'CAL_CCD','m','calcite CCD depth')

   call self%request_coupling(self%id_ftempca,'tempca_sms_tot')
   call self%request_coupling(self%id_ftempc, 'tempc_sms_tot')
   call self%request_coupling(self%id_ftempn, 'tempn_sms_tot')
   call self%request_coupling(self%id_ftempfe, 'tempfe_sms_tot')
   call self%request_coupling(self%id_ftempsi, 'tempsi_sms_tot')

   call self%register_state_dependency(self%id_ZOXY,'OXY','mmol O_2 m-3', 'dissolved oxygen')
   call self%register_state_dependency(self%id_ZDIN,'DIN','mmol N m-3', 'nitrogen nutrient')
   call self%register_state_dependency(self%id_ZSIL,'SIL','mmol Si m-3', 'silicic acid')
   call self%register_state_dependency(self%id_ZFER,'FER','mmol Fe m-3', 'iron nutrient')
   call self%register_state_dependency(self%id_ZDIC,'DiC','mmol C m-3', 'dissolved inorganic carbon')
   call self%register_state_dependency(self%id_ZDET,'DET','mmol N m-3', 'detritus nitrogen')
   call self%register_state_dependency(self%id_ZDTC,'DTC','mmol C m-3', 'detritus carbon')
   call self%register_state_dependency(self%id_ZALK,'ALK','meq/m**3', 'total alkalinity')

   call self%register_diagnostic_variable(self%id_freminc,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_freminn,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_freminsi,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_freminfe,'freminfe','mmol Fe m-3 s-1','remineralisation of detritus (Fe)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_freminca,'freminca','mmol CaCO3 m-3 s-1','remineralisation of calcite (CaCO3)',missing_value=0.0_rk,source=source_do_column)

   call self%register_dependency(self%id_freminc1,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)')
   call self%register_dependency(self%id_freminn1,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)')
   call self%register_dependency(self%id_freminsi1,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)')
   call self%register_dependency(self%id_freminfe1,'freminfe','mmol Fe m-3 s-1','remineralisation of detritus (Fe)')
   call self%register_dependency(self%id_freminca1,'freminca','mmol CaCO3 m-3 s-1','remineralisation of calcite (CaCO3)')

   call self%register_diagnostic_variable(self%id_ffastc_loc,'ffastc_loc','mmol C m-2 s-1','local remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_ffastn_loc,'ffastn_loc','mmol C m-2 s-1','local remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column)
call self%register_diagnostic_variable(self%id_ffastca_loc,'ffastca_loc','mmol Ca m-2 s-1','local remineralisation of detritus (Ca)',missing_value=0.0_rk,source=source_do_column)
call self%register_diagnostic_variable(self%id_ffastsi_loc,'ffastsi_loc','mmol Si m-2 s-1','local remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column)

   call self%register_diagnostic_variable(self%id_ffastc,'ffastc','mmol C m-2 s-1','remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_ffastn,'ffastn','mmol N m-2 s-1','remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_ffastfe,'ffastfe','mmol Fe m-2 s-1','remineralisation of detritus (Fe)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_ffastsi,'ffastsi','mmol Si m-2 s-1','remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column)
   call self%register_diagnostic_variable(self%id_ffastca,'ffastca','mmol CaCO3 m-2 s-1','remineralisation of calcite (CaCO3)',missing_value=0.0_rk,source=source_do_column,output=output_none)

   call self%register_diagnostic_variable(self%id_SEAFLRN,'SEAFLRN','mmolN/m2/s','Seafloor flux of N',source=source_do_column)
   call self%register_diagnostic_variable(self%id_SEAFLRSI,'SEAFLRSI','mmolSi/m2/s','Seafloor flux of Si',source=source_do_column)
   call self%register_diagnostic_variable(self%id_SEAFLRFE,'SEAFLRFE','mmolFe/m2/s','Seafloor flux of Fe',source=source_do_column)
   call self%register_diagnostic_variable(self%id_SEAFLRC,'SEAFLRC','mmolC/m2/s','Seafloor flux of C',source=source_do_column)
   call self%register_diagnostic_variable(self%id_SEAFLRCA,'SEAFLRCA','mmolCaCO3/m2/s','Seafloor flux of CaCO3',source=source_do_column)

   call self%register_dependency(self%id_om_cal,'OM_CAL3','-','calcite saturation')
   call self%register_diagnostic_variable(self%id_OCAL_LVL,'OCAL_LVL','m','Calcite CCD level',source=source_do_column)
   call self%register_diagnostic_variable(self%id_freminn2d,'freminn2d','mmol N m-2 s-1','remineralisation of detritus (N) * dz',source=source_do_column)
   call self%register_horizontal_dependency(self%id_ffastc1,'ffastc','mmol C m-2 s-1','remineralisation of detritus (C)')
   call self%register_horizontal_dependency(self%id_ffastn1,'ffastn','mmol N m-2 s-1','remineralisation of detritus (N)')
   call self%register_horizontal_dependency(self%id_ffastfe1,'ffastfe','mmol Fe m-2 s-1','remineralisation of detritus (Fe)')
   call self%register_horizontal_dependency(self%id_ffastsi1,'ffastsi','mmol Si m-2 s-1','remineralisation of detritus (Si)')
   call self%register_horizontal_dependency(self%id_ffastca1,'ffastca','mmol CaCO3 m-2 s-1','remineralisation of calcite (CaCO3)')

   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)

   call self%get_parameter(self%seafloor,'seafloor','-','seafloor parameterisation: 1-inorganic returns, 2-organic returns, 3-coupled benthic model', default = 3)
   call self%get_parameter(self%xrfn,'xrfn','umol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=0.03_rk)

   if (self%seafloor .eq. 3) then
         call self%register_state_dependency(self%id_ZSEDSI,'BEN_SI','mmol Si m-2', 'sediment (Si)')
         call self%register_state_dependency(self%id_ZSEDCA,'BEN_CA','mmol Ca m-2', 'sediment (Ca)')
         call self%register_state_dependency(self%id_ZSEDC,'BEN_C','mmol C m-2', 'sediment (C)')
         call self%register_state_dependency(self%id_ZSEDN,'BEN_N','mmol N m-2', 'sediment (N)')
         call self%register_state_dependency(self%id_ZSEDP,'BEN_P','mmol P m-2', 'sediment (P)')
         call self%register_state_dependency(self%id_ZSEDFE,'BEN_FE','mmol Fe m-2', 'sediment (Fe)')
         call self%request_coupling_to_model(self%id_ZSEDC, 'BEN', standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_ZSEDN, 'BEN', standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_ZSEDP, 'BEN', standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_ZSEDSI, 'BEN', standard_variables%total_silicate)
   end if

   end subroutine initialize

   subroutine do_fast_detritus(self,_ARGUMENTS_VERTICAL_)
   class(type_medusa_fast_detritus),intent(in) :: self
      !----------------------------------------------------------------------------------------------------------
      ! This version of MEDUSA includes Yool et al (2011) ballast model (iball=1 in the original implementation):
      !  1.  Fast detritus is stored as a 2D array  [ ffastX  ]
      !  2.  Fast detritus is added level-by-level  [ ftempX  ]
      !  3.  Fast detritus is not remineralised in the top box [ freminX ]
      !  4.  Remaining fast detritus is remineralised in the bottom  [ fsedX ] box
      !
      ! The method couples C, N and Fe remineralisation to the remineralisation of particulate Si and CaCO3
      !
   _DECLARE_ARGUMENTS_VERTICAL_

   real(rk) :: dz,fq0,fq1,fq2,fq3,fq4,fq5,fq6,fq7,fq8,fprotf
   real(rk) :: xmassc = 12.011_rk
   real(rk) :: xmassca = 100.086_rk
   real(rk) :: xmasssi = 60.084_rk
   real(rk) :: xprotca = 0.07_rk
   real(rk) :: xprotsi = 0.026_rk
   real(rk) :: xfastc = 188._rk, xfastsi = 2000._rk, xfastca = 3500._rk

   real(rk) :: ffastc,ffastn,ffastca,ffastsi,ffastfe
   real(rk) :: ftempc,ftempn,ftempfe,ftempsi,ftempca
   real(rk) :: freminc,freminn,freminfe,freminsi,freminca
   real(rk) :: om_cal,fdep1,cal_ccd
   integer ::  ccd_count
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      !==========================================================================================
      ! BALLAST SUBMODEL
      !==========================================================================================
      ! 
      !------------------------------------------------------------------------------------------
      ! Fast-sinking detritus fluxes, pt. 1: REMINERALISATION
      ! aside from explicitly modelled, slow-sinking detritus, the model includes an implicit representation of detrital
      ! particles that sink too quickly to be modelled with explicit state variables; this sinking flux is instead
      ! instantaneously remineralised down the water column using the version of Armstrong et al. (2002)'s ballast model
      ! used by Dunne et al. (2007); the version of this model here considers silicon and calcium carbonate ballast
      ! minerals; this section of the code redistributes the fast sinking material generated locally down the water column;
      ! this differs from Dunne et al. (2007) in that fast sinking material is distributed at *every* level below that it is
      ! generated, rather than at every level below some fixed depth; this scheme is also different in that sinking 
      ! material generated in one level is aggregated with that generated by shallower levels; this should make the 
      ! ballast model more self-consistent (famous last words)
      !-----------------------------------------------------------------------------------------
      !
   fdep1 = 0._rk
   ccd_count = 0
      ! this is the SURFACE OCEAN BOX (no remineralisation on first iteration of the vertical loop)
   ffastc=0._rk
   ffastn=0._rk
   ffastca=0._rk
   ffastsi=0._rk
   ffastfe=0._rk

   _VERTICAL_LOOP_BEGIN_

    _GET_(self%id_dz,dz)
   ! oranic carbon
   fq0      = ffastc                              ! how much organic C enters this box        (mol)
   fq1      = (fq0 * xmassc)                      ! how much it weighs                        (mass)
   fq2      = (ffastca * xmassca)                 ! how much CaCO3 enters this box            (mass)
   fq3      = (ffastsi * xmasssi)                 ! how much opal enters this box             (mass)
   fq4      = (fq2 * xprotca) + (fq3 * xprotsi)   ! total protected organic C                 (mass)
   !
   !   This next term is calculated for C but used for N and Fe as well
   !   it needs to be protected in case ALL C is protected
   !
   if (fq4.lt.fq1) then
      fprotf   = (fq4 / (fq1 + tiny(fq1)))        ! protected fraction of total organic C     (non-dim)
   else
      fprotf   = 1._rk                            ! all organic C is protected                (non-dim)
   endif
   fq5      = (1._rk - fprotf)                    ! unprotected fraction of total organic C   (non-dim)
   fq6      = (fq0 * fq5)                         ! how much organic C is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))         ! how much unprotected C leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))              ! how much total C leaves this box          (mol)
   freminc  = (fq0 - fq8) / dz                    ! C remineralisation in this box            (mol)
   ffastc = fq8

  _SET_DIAGNOSTIC_(self%id_freminc,freminc)

   ! organic nitrogen
   fq0      = ffastn                              ! how much organic N enters this box        (mol)
   fq5      = (1._rk - fprotf)                    ! unprotected fraction of total organic N   (non-dim)
   fq6      = (fq0 * fq5)                         ! how much organic N is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))         ! how much unprotected N leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))              ! how much total N leaves this box          (mol)
   freminn  = (fq0 - fq8) / dz                    ! N remineralisation in this box            (mol)
   ffastn = fq8
  
  _SET_DIAGNOSTIC_(self%id_freminn2d,freminn * dz / d_per_s)
  _SET_DIAGNOSTIC_(self%id_freminn,freminn / d_per_s)

   ! organic iron
   fq0      = ffastfe                             ! how much organic Fe enters this box       (mol)
   fq5      = (1._rk - fprotf)                    ! unprotected fraction of total organic Fe  (non-dim)
   fq6      = (fq0 * fq5)                         ! how much organic Fe is unprotected        (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))         ! how much unprotected Fe leaves this box   (mol)
   fq8      = (fq7 + (fq0 * fprotf))              ! how much total Fe leaves this box         (mol)
   freminfe = (fq0 - fq8) / dz                    ! Fe remineralisation in this box           (mol)
   ffastfe = fq8

  _SET_DIAGNOSTIC_(self%id_freminfe,freminfe / d_per_s)

   ! biogenic silicon
   fq0      = ffastsi                             ! how much  opal centers this box           (mol) 
   fq1      = fq0 * exp(-(dz / xfastsi))          ! how much  opal leaves this box            (mol)
   freminsi = (fq0 - fq1) / dz                    ! Si remineralisation in this box           (mol)
   _SET_DIAGNOSTIC_(self%id_freminsi,freminsi / d_per_s)
   ffastsi = fq1

   ! biogenic calcium carbonate

    fdep1 = fdep1 + dz
                       
   _GET_HORIZONTAL_(self%id_CAL_CCD,cal_ccd)      ! use calculated CCD field
   _GET_(self%id_om_cal,om_cal)
 
   fq0      = ffastca                             ! how much CaCO3 enters this box            (mol)
   if (fdep1-dz.le.cal_ccd) then                  ! whole grid cell above CCD
    fq1      = fq0                                ! above lysocline, no Ca dissolves          (mol)
    freminca = 0._rk                              ! above lysocline, no Ca dissolves          (mol)
    ccd_count = ccd_count + 1
   elseif (fdep1-dz.ge.cal_ccd) then              ! whole grid cell below CCD
    fq1      = fq0 * exp(-(dz / xfastca))         ! how much CaCO3 leaves this box            (mol)
    freminca = (fq0 - fq1) / dz                   ! Ca remineralisation in this box           (mol)
   else                                           ! partial grid cell below CCD
    fq2      = fdep1 - cal_ccd                    ! amount of grid cell below CCD             (m)
    fq1      = fq0 * exp (-(fq2 / xfastca))       ! how much CaCO3 leaves this box            (mol)
    freminca = (fq0 - fq1) / dz                   ! Ca remineralisation in this box           (mol)
   endif

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_OCAL_LVL,ccd_count)
   _SET_DIAGNOSTIC_(self%id_freminca,freminca)
   ffastca = fq1

   !-------------------------------------------------------------------------------------------------
   ! Fast-sinking detritus fluxes, pt. 2: UPDATE FAST FLUXES
   ! here locally calculated additions to the fast-sinking flux are added to the total fast-sinking flux;
   !this is done here such that material produced in a particular layer is only remineralised below this layer
   !-------------------------------------------------------------------------------------------------
   !add sinking material generated in this layer to running totals
    _GET_(self%id_ftempc,ftempc)                   ! === organic carbon ===
    _GET_(self%id_ftempn,ftempn)                   ! === organic nitrogen ===
    _GET_(self%id_ftempfe,ftempfe)                 ! === organic iron ===
    _GET_(self%id_ftempsi,ftempsi)                 ! === biogenic silicon ===
    _GET_(self%id_ftempca,ftempca)                 ! === biogenic calcium carbonate ===
    ffastc  = ffastc + ftempc * dz                 ! diatom and mesozooplankton mortality
    ffastn  = ffastn  + ftempn * dz                ! diatom and mesozooplankton mortality
    ffastfe = ffastfe + ftempfe * dz               ! diatom and mesozooplankton mortality
    ffastsi = ffastsi + ftempsi * dz               ! diatom mortality and grazed diatoms
    ffastca = ffastca + ftempca * dz               ! latitudinally-based fraction of total primary production

   _SET_DIAGNOSTIC_(self%id_ffastc_loc,ffastc)
   _SET_DIAGNOSTIC_(self%id_ffastn_loc,ffastn)
   _SET_DIAGNOSTIC_(self%id_ffastsi_loc,ffastsi)
   _SET_DIAGNOSTIC_(self%id_ffastca_loc,ffastca)

   _VERTICAL_LOOP_END_

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastc,ffastc)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastn,ffastn)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastfe,ffastfe)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastsi,ffastsi)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastca,ffastca)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_SEAFLRN,ffastn / d_per_s)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_SEAFLRSI,ffastsi / d_per_s)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_SEAFLRFE,ffastfe / d_per_s)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_SEAFLRC,ffastc / d_per_s)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_SEAFLRCA,ffastca / d_per_s)

   end subroutine do_fast_detritus

   subroutine do(self,_ARGUMENTS_DO_)

     class(type_medusa_fast_detritus), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_

     real(rk) :: ZOXY
     real(rk) :: freminc,freminn,freminsi,freminfe,freminca
     real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _LOOP_BEGIN_
       _GET_(self%id_ZOXY,ZOXY)
       _GET_(self%id_freminc1,freminc)
       _GET_(self%id_freminn1,freminn)
       _GET_(self%id_freminsi1,freminsi)
       _GET_(self%id_freminfe1,freminfe)
       _GET_(self%id_freminca1,freminca)
       _SET_ODE_(self%id_ZDIC, + freminc + freminca)
       _SET_ODE_(self%id_ZDIN, + freminn * d_per_s)
       _SET_ODE_(self%id_ZSIL, + freminsi * d_per_s)
       _SET_ODE_(self%id_ZFER, + freminfe * d_per_s)

       if (ZOXY .ge. self%xo2min) _SET_ODE_(self%id_ZOXY, - self%xthetarem * freminc - self%xthetanit * freminn * d_per_s)

      _SET_ODE_(self%id_ZALK, 2._rk * freminca)

   _LOOP_END_

   end subroutine do

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

     class(type_medusa_fast_detritus), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: ffastc,ffastn,ffastsi,ffastca,ffastfe
     real(rk) :: ZOXY
   !----------------------------------------------------------
   ! Fast-sinking detritus fluxes, pt. 3: SEAFLOOR
   ! remineralise all remaining fast-sinking detritus to dissolved nutrients;
   ! Three options for seafloor handling recreate original MEDUSA options 
   ! Option 1: (jfdfate.eq.0 .and. jorgben.eq.0) - immediate remineralisation
   ! Option 2: (jfdfate.eq.1 .and. jorgben.eq.0) - conversion of fast C into slow C
   ! Option 3: (jfdfate.eq.0 .and. jorgben.eq.1) - conversion of fast C into benthic C
   ! Note that Option 3 requires coupling of benthic module to allow mineralisation (and transformation) of benthic material
   !
    _HORIZONTAL_LOOP_BEGIN_
    _GET_(self%id_ZOXY,ZOXY)
    _GET_HORIZONTAL_(self%id_ffastc1,ffastc)
    _GET_HORIZONTAL_(self%id_ffastn1,ffastn)
    _GET_HORIZONTAL_(self%id_ffastsi1,ffastsi)
    _GET_HORIZONTAL_(self%id_ffastca1,ffastca)
    _GET_HORIZONTAL_(self%id_ffastfe1,ffastfe)

     if (self%seafloor .eq. 1) then ! C remineralisation

    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, + ffastc)
     if (ZOXY .ge. self%xo2min) _SET_BOTTOM_EXCHANGE_(self%id_ZOXY, - self%xthetarem * ffastc - self%xthetanit * ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, + ffastca)
    _SET_BOTTOM_EXCHANGE_(self%id_ZALK, 2._rk * ffastca)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIN, + ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + ffastsi)
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, + ffastn * self%xrfn)

     elseif (self%seafloor .eq. 2) then ! fast C -> slow C
     ! not tested!
     if (ZOXY .ge. self%xo2min) _SET_BOTTOM_EXCHANGE_(self%id_ZOXY, - self%xthetarem * ffastc - self%xthetanit * ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZALK, 2._rk * ffastca)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDTC, + ffastc + ffastca)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDET, + ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + ffastsi)
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, + ffastn * self%xrfn)

     elseif (self%seafloor .eq. 3) then ! fast C -> benthic C. Did we couple benthic module?

    _SET_BOTTOM_ODE_(self%id_ZSEDC, + ffastc)
    _SET_BOTTOM_ODE_(self%id_ZSEDN, + ffastn)
    _SET_BOTTOM_ODE_(self%id_ZSEDSI, + ffastsi)
    _SET_BOTTOM_ODE_(self%id_ZSEDCA, + ffastca)
    _SET_BOTTOM_ODE_(self%id_ZSEDFE, + ffastfe)

     _SET_BOTTOM_ODE_(self%id_ZSEDP, + ffastc/106._rk)

     end if
 
    _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module medusa_fast_detritus
