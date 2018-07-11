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
      type (type_state_variable_id)        :: id_ZDIC,id_ZDIN,id_ZSIL,id_ZOXY,id_ZFER,id_ZDET,id_ZDTC,id_ZALK
      type (type_bottom_state_variable_id) :: id_ZSEDSI,id_ZSEDC,id_ZSEDN,id_ZSEDCA
      type (type_dependency_id)            :: id_dz,id_depth
      type (type_dependency_id)            :: id_ftempc,id_ftempn,id_ftempsi,id_ftempfe,id_ftempca
      type (type_dependency_id)            :: id_freminc1,id_freminn1,id_freminsi1,id_freminfe1,id_freminca1
      type (type_dependency_id)            :: id_om_cal
      type (type_diagnostic_variable_id)   :: id_freminc,id_freminn,id_freminsi,id_freminfe,id_freminca
      type (type_diagnostic_variable_id)   :: id_ffastc_loc,id_ffastca_loc,id_ffastsi_loc
      type (type_horizontal_diagnostic_variable_id) :: id_ffastc,id_ffastn,id_ffastsi,id_ffastfe,id_ffastca
      type (type_horizontal_dependency_id) :: id_ffastc1,id_ffastn1,id_ffastfe1,id_ffastsi1,id_ffastca1
      type (type_diagnostic_variable_id)   :: id_tempc,id_tempn,id_tempsi,id_tempfe,id_tempca
      ! Parameters
      real(rk) :: xthetanit,xthetarem,xo2min,xrfn,xfe_sed
      integer :: seafloor
      integer :: iball

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

   call self%get_parameter(self%xfe_sed,'xfe_sed','mmol Fe m-2 d-1','sedimentary flux of iron',default=0.000228_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)
   call self%get_parameter(self%iball,'iball','ballast model formulation (1- ballast model (Yool et al., 2011), 2- ballast-sans-ballast model)', default=1)

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

   call self%request_coupling(self%id_ftempca,'tempca_sms_tot')
   call self%request_coupling(self%id_ftempc, 'tempc_sms_tot')
   call self%request_coupling(self%id_ftempn, 'tempn_sms_tot')
   call self%request_coupling(self%id_ftempfe, 'tempfe_sms_tot')
   call self%request_coupling(self%id_ftempsi, 'tempsi_sms_tot')

   call self%register_state_dependency(self%id_ZOXY,'ZOXY','mmol O_2 m-3', 'dissolved oxygen')
   call self%register_state_dependency(self%id_ZDIN,'ZDIN','mmol N m-3', 'nitrogen nutrient')
   call self%register_state_dependency(self%id_ZSIL,'ZSIL','mmol Si m-3', 'silicic acid')
   call self%register_state_dependency(self%id_ZFER,'ZFER','mmol Fe m-3', 'iron nutrient')
   call self%register_state_dependency(self%id_ZDIC,'ZDIC','mmol C m-3', 'dissolved inorganic carbon')
   call self%register_state_dependency(self%id_ZDET,'ZDET','mmol N m-3', 'detritus nitrogen')
   call self%register_state_dependency(self%id_ZDTC,'ZDTC','mmol C m-3', 'detritus carbon')
   call self%register_state_dependency(self%id_ZALK,'ZALK','meq/m**3', 'total alkalinity')

   call self%register_diagnostic_variable(self%id_freminc,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminn,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminsi,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminfe,'freminfe','mmol Fe m-3 s-1','remineralisation of detritus (Fe)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_freminca,'freminca','mmol CaCO3 m-3 s-1','remineralisation of calcite (CaCO3)',missing_value=0.0_rk,source=source_do_column,output=output_none)

   call self%register_dependency(self%id_freminc1,'freminc','mmol C m-3 s-1','remineralisation of detritus (C)')
   call self%register_dependency(self%id_freminn1,'freminn','mmol N m-3 s-1','remineralisation of detritus (N)')
   call self%register_dependency(self%id_freminsi1,'freminsi','mmol Si m-3 s-1','remineralisation of detritus (Si)')
   call self%register_dependency(self%id_freminfe1,'freminfe','mmol Fe m-3 s-1','remineralisation of detritus (Fe)')
   call self%register_dependency(self%id_freminca1,'freminca','mmol CaCO3 m-3 s-1','remineralisation of calcite (CaCO3)')

   call self%register_diagnostic_variable(self%id_ffastc_loc,'ffastc_loc','mmol C m-2 s-1','local remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column)
call self%register_diagnostic_variable(self%id_ffastca_loc,'ffastca_loc','mmol Ca m-2 s-1','local remineralisation of detritus (Ca)',missing_value=0.0_rk,source=source_do_column)
call self%register_diagnostic_variable(self%id_ffastsi_loc,'ffastsi_loc','mmol Si m-2 s-1','local remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column)

   call self%register_diagnostic_variable(self%id_ffastc,'ffastc','mmol C m-2 s-1','remineralisation of detritus (C)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_ffastn,'ffastn','mmol N m-2 s-1','remineralisation of detritus (N)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_ffastfe,'ffastfe','mmol Fe m-2 s-1','remineralisation of detritus (Fe)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_ffastsi,'ffastsi','mmol Si m-2 s-1','remineralisation of detritus (Si)',missing_value=0.0_rk,source=source_do_column,output=output_none)
   call self%register_diagnostic_variable(self%id_ffastca,'ffastca','mmol CaCO3 m-2 s-1','remineralisation of calcite (CaCO3)',missing_value=0.0_rk,source=source_do_column,output=output_none)

   call self%register_dependency(self%id_om_cal,'om_cal','-','calcite saturation')

   call self%register_horizontal_dependency(self%id_ffastc1,'ffastc','mmol C m-2 s-1','remineralisation of detritus (C)')
   call self%register_horizontal_dependency(self%id_ffastn1,'ffastn','mmol N m-2 s-1','remineralisation of detritus (N)')
   call self%register_horizontal_dependency(self%id_ffastfe1,'ffastfe','mmol Fe m-2 s-1','remineralisation of detritus (Si)')
   call self%register_horizontal_dependency(self%id_ffastsi1,'ffastsi','mmol Si m-2 s-1','remineralisation of detritus (Si)')
   call self%register_horizontal_dependency(self%id_ffastca1,'ffastca','mmol CaCO3 m-2 s-1','remineralisation of calcite (CaCO3)')

   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   call self%register_dependency(self%id_depth, standard_variables%depth)

   call self%get_parameter(self%seafloor,'seafloor','-','seafloor handling: 1-inorganic returns, 2-organic returns, 3-coupled benthic model', default = 1)
   call self%get_parameter(self%xrfn,'xrfn','umol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=0.03_rk) !worth to double-check
   if (self%seafloor .eq. 3)  call self%register_state_dependency(self%id_ZSEDSI,'ZSEDSI','mmol Si m-2', 'sediment (Si)')
   if (self%seafloor .eq. 3)  call self%register_state_dependency(self%id_ZSEDC,'ZSEDC','mmol C m-2', 'sediment (C)')
   if (self%seafloor .eq. 3)  call self%register_state_dependency(self%id_ZSEDN,'ZSEDN','mmol N m-2', 'sediment (N)')
   if (self%seafloor .eq. 3)  call self%register_state_dependency(self%id_ZSEDCA,'ZSEDCA','mmol N m-2', 'sediment (CaCO3)')

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
   real(rk) :: xfastc = 188._rk, xfastsi = 2000._rk, xfastca = 3500._rk

   real(rk) :: ffastc,ffastn,ffastca,ffastsi,ffastfe
   real(rk) :: ftempc,ftempn,ftempfe,ftempsi,ftempca
   real(rk) :: freminc,freminn,freminfe,freminsi,freminca
   real(rk) :: om_cal

   ffastc=0._rk
   ffastn=0._rk
   ffastca=0._rk
   ffastsi=0._rk
   ffastfe=0._rk

   _VERTICAL_LOOP_BEGIN_

    _GET_(self%id_dz,dz)

!   !Carbon
   fq0      = ffastc                            !! how much organic C enters this box        (mol)
   if (self%iball .eq. 1) then
     fq1      = (fq0 * xmassc)                    !! how much it weighs                        (mass)
     fq2      = (ffastca * xmassca)               !! how much CaCO3 enters this box            (mass)
     fq3      = (ffastsi * xmasssi)               !! how much opal enters this box             (mass)
     fq4      = (fq2 * xprotca) + (fq3 * xprotsi) !! total protected organic C                 (mass)
!
!   !! this next term is calculated for C but used for N and Fe as well
!   !! it needs to be protected in case ALL C is protected
!
     if (fq4.lt.fq1) then
        fprotf   = (fq4 / (fq1 + tiny(fq1)))      !! protected fraction of total organic C     (non-dim)
     else
        fprotf   = 1._rk                         !! all organic C is protected                (non-dim)
     endif
     fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic C   (non-dim)
     fq6      = (fq0 * fq5)                     !! how much organic C is unprotected         (mol)
     fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected C leaves this box    (mol)
     fq8      = (fq7 + (fq0 * fprotf))          !! how much total C leaves this box          (mol)
     freminc  = (fq0 - fq8) / dz                !! C remineralisation in this box            (mol)
     ffastc = fq8
   else
     fq1=fq0 * exp(-(dz/xfastc))
     freminc = (fq0 - fq1) / dz
     ffastc = fq1
   end if
  _SET_DIAGNOSTIC_(self%id_freminc,freminc)
  _SET_DIAGNOSTIC_(self%id_ffastc_loc,ffastc)

   !Nitrogen
   fq0      = ffastn                            !! how much organic N enters this box        (mol)
   if (self%iball .eq. 1) then
     fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic N   (non-dim)
     fq6      = (fq0 * fq5)                     !! how much organic N is unprotected         (mol)
     fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected N leaves this box    (mol)
     fq8      = (fq7 + (fq0 * fprotf))          !! how much total N leaves this box          (mol)
     freminn  = (fq0 - fq8) / dz                !! N remineralisation in this box            (mol)
     ffastn = fq8
   else
     fq1 = fq0 * exp(-(dz / xfastc))
     freminn = (fq0 - fq1) / dz
     ffastn = fq1
   end if
   _SET_DIAGNOSTIC_(self%id_freminn,freminn)


   !Iron
   fq0      = ffastfe                           !! how much organic Fe enters this box       (mol)
   if (self%iball .eq. 1) then
     fq5      = (1._rk - fprotf)                !! unprotected fraction of total organic Fe  (non-dim)
     fq6      = (fq0 * fq5)                     !! how much organic Fe is unprotected        (mol)
     fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected Fe leaves this box   (mol)
     fq8      = (fq7 + (fq0 * fprotf))          !! how much total Fe leaves this box         (mol)            
     freminfe = (fq0 - fq8) / dz                !! Fe remineralisation in this box           (mol)
     ffastfe = fq8
   else
     fq1 = fq0 * exp(-(dz / xfastc))
     freminfe = (fq0 - fq1) / dz
     ffastfe = fq1
   end if
  _SET_DIAGNOSTIC_(self%id_freminfe,freminfe)

   !biogenic silicon
   fq0      = ffastsi                         !! how much  opal centers this box           (mol) 
   fq1      = fq0 * exp(-(dz / xfastsi))      !! how much  opal leaves this box            (mol)
   freminsi = (fq0 - fq1) / dz                !! Si remineralisation in this box           (mol)
   _SET_DIAGNOSTIC_(self%id_freminsi,freminsi)
   ffastsi = fq1
   _SET_DIAGNOSTIC_(self%id_ffastsi_loc,ffastsi)

   !biogenic calcium carbonate
  
   _GET_(self%id_om_cal,om_cal)
   fq0      = ffastca                           !! how much CaCO3 enters this box            (mol)
   if (om_cal .ge. 1._rk) then
   fq1      = fq0                               !! above lysocline, no Ca dissolves          (mol)
   freminca = 0._rk                               !! above lysocline, no Ca dissolves          (mol)
   elseif (om_cal.lt. 1._rk) then
   fq1      = fq0 * exp(-(dz / xfastca))        !! how much CaCO3 leaves this box            (mol)
   freminca = (fq0 - fq1) / dz                  !! Ca remineralisation in this box           (mol)
   endif
   _SET_DIAGNOSTIC_(self%id_freminca,freminca)
   ffastca = fq1
   _SET_DIAGNOSTIC_(self%id_ffastca_loc,ffastca)

    _GET_(self%id_ftempc,ftempc)
    _GET_(self%id_ftempn,ftempn)
    _GET_(self%id_ftempfe,ftempfe)
    _GET_(self%id_ftempsi,ftempsi)
    _GET_(self%id_ftempca,ftempca)

    ffastc  = ffastc + ftempc * dz
    ffastn  = ffastn  + ftempn * dz
   ! ffastfe = ffastfe + ftempfe * dz
    ffastsi = ffastsi + ftempsi * dz
    ffastca = ffastca + ftempca * dz
   _VERTICAL_LOOP_END_

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastc,ffastc)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastn,ffastn)
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastfe,ffastfe)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastsi,ffastsi)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ffastca,ffastca)
   end subroutine do_fast_detritus

   subroutine do(self,_ARGUMENTS_DO_)

     class(type_medusa_fast_detritus), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_

     real(rk) :: ZOXY
     real(rk) :: freminc,freminn,freminsi,freminfe,freminca

    _LOOP_BEGIN_
       _GET_(self%id_ZOXY,ZOXY)
       _GET_(self%id_freminc1,freminc)
       _GET_(self%id_freminn1,freminn)
       _GET_(self%id_freminsi1,freminsi)
       _GET_(self%id_freminfe1,freminfe)
       _GET_(self%id_freminca1,freminca)
       _SET_ODE_(self%id_ZDIC, + freminc + freminca)
       _SET_ODE_(self%id_ZDIN, + freminn)
       _SET_ODE_(self%id_ZSIL, + freminsi)
       _SET_ODE_(self%id_ZFER, + freminfe)

       if (ZOXY .ge. self%xo2min) _SET_ODE_(self%id_ZOXY, - self%xthetarem * freminc - self%xthetanit * freminn)

      _SET_ODE_(self%id_ZALK, 2._rk * freminca)

   _LOOP_END_

   end subroutine do

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

     class(type_medusa_fast_detritus), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: ffastc,ffastn,ffastsi,ffastca !ffastfe
     real(rk) :: depth

     !TO-DO: oxygen consumption

    _HORIZONTAL_LOOP_BEGIN_

    _GET_HORIZONTAL_(self%id_ffastc1,ffastc)
    _GET_HORIZONTAL_(self%id_ffastn1,ffastn)
    _GET_HORIZONTAL_(self%id_ffastsi1,ffastsi)
    _GET_HORIZONTAL_(self%id_ffastca1,ffastca)

 !   _GET_HORIZONTAL_(self%id_ffastfe1,ffastfe)

    _GET_(self%id_depth,depth)
    if (depth .lt. 500._rk) then 
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, self%xfe_sed) ! sit here for now...
    end if

     if (self%seafloor .eq. 1) then

    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, + ffastc)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIN, + ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + ffastsi)
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, + ffastn * self%xrfn)

     elseif (self%seafloor .eq. 2) then

    _SET_BOTTOM_EXCHANGE_(self%id_ZDTC, + ffastc + ffastca)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDET, + ffastn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + ffastsi)
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, + ffastn * self%xrfn)

     elseif (self%seafloor .eq. 3) then !medusa_benthic is coupled

    _SET_BOTTOM_ODE_(self%id_ZSEDC, + ffastc)
    _SET_BOTTOM_ODE_(self%id_ZSEDN, + ffastn)
    _SET_BOTTOM_ODE_(self%id_ZSEDSI, + ffastsi)
    _SET_BOTTOM_ODE_(self%id_ZSEDCA, + ffastca)

     end if
 
    _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module medusa_fast_detritus
