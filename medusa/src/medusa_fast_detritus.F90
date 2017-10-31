#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA fast detritus mineralisation
!                    and implicit sinking
!*********************************************************

module medusa_fast_detritus

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_fast_detritus
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZDIC,id_ZDIN,id_ZSIL,id_ZOXY
      type (type_dependency_id)            :: id_dz
      type (type_dependency_id)            :: id_ftempc,id_ftempn,id_ftempsi,id_ftempfe,id_ftempca,id_freminc1,id_freminn1,id_freminsi1 
      type (type_diagnostic_variable_id)   :: id_freminc,id_freminn,id_freminsi  
      type (type_diagnostic_variable_id)   :: id_tempc,id_tempn,id_tempsi,id_tempfe,id_tempca
      ! Parameters
      real(rk) :: xthetanit,xthetarem,xo2min

   contains

      procedure :: initialize
      procedure :: get_light => do_fast_detritus
      procedure :: do

  end type type_medusa_fast_detritus

contains

   subroutine initialize(self,configunit)

   class(type_medusa_fast_detritus),intent(inout),target :: self
   integer,                         intent(in)           :: configunit

   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)

   ! Create diagnostics for fast-sinking detritus that acts as state variables, so they can receive sources.
   call self%register_diagnostic_variable(self%id_tempc,'tempc','mmol C m-3','fast-sinking detritus (C)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempn,'tempn','mmol N m-3','fast-sinking detritus (N)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempfe,'tempfe','-','fast-sinking detritus (Fe)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%register_diagnostic_variable(self%id_tempsi,'tempsi','mmol Si m-3','fast-sinking detritus (Si)', act_as_state_variable=.true., missing_value=0.0_rk,source=source_none, output=output_none)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_tempc)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_tempn)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_tempsi)

   ! Retrieve the sources of fast-sinking detritus.
   call self%register_dependency(self%id_ftempc, 'ftempc', 'mmol C m-3 s-1', 'production of fast-sinking detritus (C)')
   call self%register_dependency(self%id_ftempn, 'ftempn', 'mmol N m-3 s-1', 'production of fast-sinking detritus (N)')
   call self%register_dependency(self%id_ftempfe,'ftempfe', 'mmol Fe m-3 s-1', 'production of fast-sinking detritus (Fe)')
   call self%register_dependency(self%id_ftempsi,'ftempsi', 'mmol Si m-3 s-1', 'production of fast-sinking detritus (Si)')
  ! call self%register_dependency(self%id_ftempca,'ftempca', '', 'production of fast-sinking detritus (CaCO3)')
   call self%request_coupling(self%id_ftempc, 'tempc_sms_tot')
   call self%request_coupling(self%id_ftempn, 'tempn_sms_tot')
   call self%request_coupling(self%id_ftempfe, 'tempfe_sms_tot')
   call self%request_coupling(self%id_ftempsi, 'tempsi_sms_tot')

   call self%register_state_dependency(self%id_ZOXY,'ZOXY','mmol O_2 m-3', 'dissolved oxygen')
   call self%register_state_dependency(self%id_ZDIN,'ZDIN','mmol N m-3', 'nitrogen nutrient')
   call self%register_state_dependency(self%id_ZSIL,'ZSIL','mmol Si m-3', 'silicic acid')
   call self%register_state_dependency(self%id_ZDIC,'ZDIC','mmol C m-3', 'dissolved inorganic carbon')

   call self%register_diagnostic_variable(self%id_freminc,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminn,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminsi,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column,output=output_none)

   call self%register_dependency(self%id_freminc1,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)')
   call self%register_dependency(self%id_freminn1,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)')
   call self%register_dependency(self%id_freminsi1,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)')

   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)

   end subroutine initialize

   subroutine do_fast_detritus(self,_ARGUMENTS_VERTICAL_)
   class(type_medusa_fast_detritus),intent(in) :: self
   
   _DECLARE_ARGUMENTS_VERTICAL_

   real(rk) :: dz,fq0,fq1,fq2,fq3,fq4,fq5,fq6,fq7,fq8,fprotf
   real(rk) :: xmassc = 12.011_rk
   real(rk) :: xmassca = 100.086_rk
   real(rk) :: xmasssi = 60.084_rk
   real(rk) :: xprotca = 0.07_rk
   real(rk) :: xprotsi = 0.026_rk
   real(rk) :: xfastc = 188_rk, xfastsi = 2000_rk

   real(rk) :: ffastc,ffastn,ffastca,ffastsi,ffastfe
   real(rk) :: ftempc,ftempn,ftempfe,ftempsi,ftempca
   real(rk) :: freminc,freminn,freminfe,freminsi,freminca
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   ffastc=0._rk
   ffastn=0._rk
   ffastca=0._rk
   ffastsi=0._rk
   ffastfe=0._rk

   _VERTICAL_LOOP_BEGIN_

    _GET_(self%id_dz,dz)

!   !Carbon
   fq0      = ffastc                            !! how much organic C enters this box        (mol)
   fq1      = (fq0 * xmassc)                    !! how much it weighs                        (mass)
   fq2      = 0._rk !(ffastca * xmassca)        !! how much CaCO3 enters this box            (mass)
   fq3      = ffastsi * xmasssi                 !! how much opal enters this box            (mass)
   fq4      = (fq2 * xprotca) + (fq3 * xprotsi) !! total protected organic C                 (mass)
!
!   !! this next term is calculated for C but used for N and Fe as well
!   !! it needs to be protected in case ALL C is protected
!
   if (fq4.lt.fq1) then
     fprotf   = (fq4 / (fq1 + tiny(fq1)))     !! protected fraction of total organic C     (non-dim)
   else
     fprotf   = 1._rk                         !! all organic C is protected                (non-dim)
   endif
   fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic C   (non-dim)
   fq6      = (fq0 * fq5)                     !! how much organic C is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected C leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))          !! how much total C leaves this box          (mol)
   freminc  = (fq0 - fq8) / dz                !! C remineralisation in this box            (mol)
   _SET_DIAGNOSTIC_(self%id_freminc,freminc)
   ffastc = fq8

   !Nitrogen
   fq0      = ffastn                          !! how much organic N enters this box        (mol)
   fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic N   (non-dim)
   fq6      = (fq0 * fq5)                     !! how much organic N is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected N leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))          !! how much total N leaves this box          (mol)
   freminn  = (fq0 - fq8) / dz                !! N remineralisation in this box            (mol)
   _SET_DIAGNOSTIC_(self%id_freminn,freminn)
   ffastn = fq8

   !Iron
   fq0      = ffastfe                         !! how much organic Fe enters this box       (mol)
   fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic Fe  (non-dim)
   fq6      = (fq0 * fq5)                     !! how much organic Fe is unprotected        (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected Fe leaves this box   (mol)
   fq8      = (fq7 + (fq0 * fprotf))          !! how much total Fe leaves this box         (mol)            
   freminfe = (fq0 - fq8) / dz                !! Fe remineralisation in this box           (mol)
   ffastfe = fq8

   !biogenic silicon
   fq0      = ffastsi                         !! how much  opal centers this box           (mol) 
   fq1      = fq0 * exp(-(dz / xfastsi))      !! how much  opal leaves this box            (mol)
   freminsi = (fq0 - fq1) / dz                !! Si remineralisation in this box           (mol)
   _SET_DIAGNOSTIC_(self%id_freminsi,freminsi)
   ffastsi = fq1

   !biogenic calcium carbonate
  ! fq0      = ffastca                           !! how much CaCO3 enters this box            (mol)
  ! if (fdep.le.fccd_dep) then
  ! !! whole grid cell above CCD
  ! fq1      = fq0                               !! above lysocline, no Ca dissolves          (mol)
  ! freminca = 0.0                               !! above lysocline, no Ca dissolves          (mol)
  ! fccd = real(jk)                              !! which is the last level above the CCD?    (#)
  ! elseif (fdep.ge.fccd_dep) then
  ! !! whole grid cell below CCD
  ! fq1      = fq0 * exp(-(dz / xfastca))        !! how much CaCO3 leaves this box            (mol)
  ! freminca = (fq0 - fq1) / dz                  !! Ca remineralisation in this box           (mol)
  ! else
  ! !! partial grid cell below CCD
  ! fq2      = fdep1 - fccd_dep                  !! amount of grid cell below CCD             (m)
  ! fq1      = fq0 * exp(-(fq2 / xfastca))       !! how much CaCO3 leaves this box            (mol)
  ! freminca = (fq0 - fq1) / dz                  !! Ca remineralisation in this box           (mol)
  ! endif
  ! ffastca = fq1

    _GET_(self%id_ftempc,ftempc)
    _GET_(self%id_ftempn,ftempn)
    _GET_(self%id_ftempfe,ftempfe)
    _GET_(self%id_ftempsi,ftempsi)
   !_GET_(self%id_ftempca,ftempca)

    ffastc  = ffastc + ftempc * dz
    ffastn  = ffastn  + ftempn * dz
    ffastfe = ffastfe + ftempfe * dz
    ffastsi = ffastsi + ftempsi * dz
  !  ffastca = ffastca + ftempca

   _VERTICAL_LOOP_END_

   end subroutine do_fast_detritus

   subroutine do(self,_ARGUMENTS_DO_)

     class(type_medusa_fast_detritus), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_

     real(rk) :: ZOXY
     real(rk) :: freminc,freminn,freminsi

    _LOOP_BEGIN_

       _GET_(self%id_ZOXY,ZOXY)
       _GET_(self%id_freminc1,freminc)
       _GET_(self%id_freminn1,freminn)
       _GET_(self%id_freminsi1,freminsi)

       _SET_ODE_(self%id_ZDIC, + freminc)
       _SET_ODE_(self%id_ZDIN, + freminn)
       _SET_ODE_(self%id_ZSIL, + freminsi)
       if (ZOXY .ge. self%xo2min) _SET_ODE_(self%id_ZOXY, - self%xthetarem * freminc - self%xthetanit * freminn)

   _LOOP_END_

   end subroutine do

end module medusa_fast_detritus
