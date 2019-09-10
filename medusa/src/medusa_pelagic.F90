#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA pelagic biogeochemistry
!*********************************************************

module medusa_pelagic

   use fabm_types
   use fabm_particle

   implicit none

   private

  type,extends(type_particle_model),public :: type_medusa_pelagic
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZCHN,id_ZCHD,id_ZPHN,id_ZPHD,id_ZPDS,id_ZDIN,id_ZFER,id_ZSIL,id_ZDET,id_ZDTC,id_ZZMI,id_ZZME,id_ZDIC,id_ZALK,id_ZOXY
      type (type_dependency_id)            :: id_temp,id_salt, id_om_cal, id_xpar
      type (type_state_variable_id)   :: id_tempc,id_tempn,id_tempsi,id_tempfe,id_tempca
      type (type_diagnostic_variable_id) :: id_prn,id_prd,id_mpn,id_mpd,id_OPAL,id_OPALDISS,id_detn,id_detc,id_MDET,id_MDETC
      type (type_diagnostic_variable_id) :: id_GMIPn,id_GMID,id_MZMI,id_MZME,id_GMEPN,id_GMEPD,id_GMEZMI,id_GMED
      type (type_diagnostic_variable_id) :: id_GMIDC,id_GMEDC,id_PN_LLOSS,id_PD_LLOSS,id_ZI_LLOSS,id_ZE_LLOSS,id_fcomm_resp
      type (type_diagnostic_variable_id) :: id_pd_jlim,id_pd_nlim,id_pd_felim,id_pd_silim,id_pd_silim2,id_pn_jlim,id_pn_nlim,id_pn_felim
      type (type_diagnostic_variable_id) :: id_fregen,id_fregensi,id_slowdetflux,id_fscal_part
      type (type_diagnostic_variable_id) :: id_ZI_MES_N,id_ZI_MES_D,id_ZI_MES_C,id_ZI_MESDC,id_ZE_MES_N,id_ZE_MES_D,id_ZE_MES_C,id_ZE_MESDC
      type (type_diagnostic_variable_id) :: id_ZI_EXCR,id_ZI_RESP,id_ZI_GROW,id_ZE_EXCR,id_ZE_RESP,id_ZE_GROW
      type (type_diagnostic_variable_id) :: id_C_PROD,id_C_CONS,id_N_PROD,id_N_CONS,id_foxy_prod,id_foxy_cons,id_foxy_anox
      type (type_horizontal_diagnostic_variable_id) :: id_f_sbenin_c,id_f_sbenin_n,id_f_sbenin_fe
      type (type_bottom_state_variable_id) :: id_ZSEDC,id_ZSEDN,id_ZSEDP,id_ZSEDFE

      ! Parameters
      logical :: jliebig
      real(rk) :: xxi,xaln,xald,xnln,xfln,xnld,xsld,xfld,xvpn,xvpd,xsin0,xnsi0,xuif,xthetam,xthetamd
      real(rk) :: xkmi,xpmipn,xpmid,xkme,xpmed,xpmepn,xpmepd,xpmezmi,zpmed,xgmi,xgme
      real(rk) :: xthetad,xphi,xthetapn,xthetazme,xthetazmi,xthetapd,xbetan,xbetac,xkc,xridg_r0,jq10
      integer :: jmpn,jmpd,jmzmi,jmzme,jphy,jmd
      real(rk) :: xmetapn,xmetapd,xmetazmi,xmetazme,xmpn,xmpd,xmzmi,xmzme,xkphn,xkphd,xkzmi,xkzme
      real(rk) :: xmd,xmdc,xsdiss
      real(rk) :: xfdfrac1,xfdfrac2,xfdfrac3,xrfn
      real(rk) :: xthetanit,xthetarem,xo2min
      real(rk) :: wg,wdep
      integer  :: seafloor

   contains

      procedure :: initialize
      procedure :: do
      procedure :: do_bottom

  end type type_medusa_pelagic

contains

   subroutine initialize(self,configunit)

   class(type_medusa_pelagic),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   !Register parameters
   call self%get_parameter(self%xxi, 'xxi', 'mol N g C-1','C : N conversion factor', default=0.01257_rk)
   call self%get_parameter(self%xaln, 'xaln', 'g C(g chl)-1 (W m-2)-1 d-1','chl-specific initial slope of P-I curve (non-diatoms)', default=15.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xald, 'xald', 'g C(g chl)-1 (W m-2)-1 d-1','chl-specific initial slope of P-I curve (diatoms)', default=11.25_rk,scale_factor=d_per_s)
   call self%get_parameter(self%jliebig, 'jliebig', 'Liebig''s minimum law for nutrient limitation', default=.false.)
   call self%get_parameter(self%xnln, 'xnln', 'mmol N m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.5_rk)
   call self%get_parameter(self%xfln, 'xfln', 'mmol Fe m-3','Fe nutrient uptake half-saturation constant (non-diatoms)', default=0.00033_rk)
   call self%get_parameter(self%xnld, 'xnld', 'mmol N m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.75_rk)
   call self%get_parameter(self%xsld, 'xsld', 'mmol Si m-3','Si nutrient uptake half-saturation constant (diatoms)', default=3.0_rk)
   call self%get_parameter(self%xfld, 'xfld', 'mmol Fe m-3','Fe nutrient uptake half-saturation constant (diatoms)', default=0.00067_rk)
   call self%get_parameter(self%xvpn, 'xvpn', 'd-1','Maximum phytoplankton growth rate (non-diatoms)', default=0.53_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xvpd, 'xvpd', 'd-1','Maximum phytoplankton growth rate (diatoms)', default=0.50_rk,scale_factor=d_per_s)
   call self%get_parameter(self%jphy, 'jphy','-','Temperature regulation (phyto growth): 1-Eppley,2-q10',default=1)
   call self%get_parameter(self%jq10, 'jq10','-','q10 factor for temperature regulation option 2',default=1.5_rk)
   call self%get_parameter(self%xsin0, 'xsin0', 'mol Si mol N-1','minimum diatom Si : N ratio', default=0.2_rk)
   call self%get_parameter(self%xnsi0, 'xnsi0', 'mol N mol Si-1','minimum diatom N : Si ratio', default=0.2_rk)
   call self%get_parameter(self%xuif, 'xuif', '-','hypothetical growth ratio at Inf Si : N ratio', default=1.5_rk)
   call self%get_parameter(self%xthetam, 'xthetam', 'g chl g C-1','maximum Chl : C ratio (non-diatoms)', default=0.05_rk)
   call self%get_parameter(self%xthetamd, 'xthetamd', 'g chl g C-1','maximum Chl : C ratio (diatoms)', default=0.05_rk)
   call self%get_parameter(self%xkmi, 'xkmi', 'mmol N m-3','microzooplankton grazing half-saturation constant', default=0.8_rk)
   call self%get_parameter(self%xpmipn, 'xpmipn', '-','microzooplankton grazing preference on non-diatoms', default=0.75_rk)
   call self%get_parameter(self%xpmid, 'xpmid', '-','microzooplankton grazing preference on detritus', default=0.25_rk)
   call self%get_parameter(self%xkme, 'xkme', 'mmol N m-3','mesozooplankton grazing half-saturation constant', default=0.3_rk)
   call self%get_parameter(self%xpmepn, 'xpmepn', '-','mesozooplankton grazing preference on non-diatoms', default=0.15_rk)
   call self%get_parameter(self%xpmepd, 'xpmepd', '-','mesozooplankton grazing preference on diatoms', default=0.35_rk)
   call self%get_parameter(self%xpmezmi, 'xpmezmi', '-','mesozooplankton grazing preference on microzooplankton', default=0.35_rk)
   call self%get_parameter(self%xpmed, 'xpmed', '-','mesozooplankton grazing preference on detritus', default=0.15_rk)
   call self%get_parameter(self%xgmi, 'xgmi', 'd-1','maximum microzooplankton grazing rate', default=2.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xgme, 'xgme', 'd-1','maximum mesozooplankton grazing rate', default=0.5_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xthetad, 'xthetad', 'mol C mol N-1','detritus C : N ratio', default=6.625_rk)
   call self%get_parameter(self%xphi, 'xphi', '-','zooplankton grazing inefficiency', default=0.20_rk)
   call self%get_parameter(self%xthetapn, 'xthetapn', 'mol C mol N-1','phytoplankton C:N ratio (non-diatoms)', default=6.625_rk)
   call self%get_parameter(self%xthetapd, 'xthetapd', 'mol C mol N-1','phytoplankton C:N ratio (diatoms)', default=6.625_rk)
   call self%get_parameter(self%xbetan, 'xbetan', '-','zooplankton N assimilation efficiency', default=0.77_rk)
   call self%get_parameter(self%xthetazmi, 'xthetazmi', 'mol C mol N-1','microzooplankton C:N ratio', default=5.625_rk)
   call self%get_parameter(self%xbetac, 'xbetac', '-','zooplankton C assimilation efficiency', default=0.64_rk)
   call self%get_parameter(self%xkc, 'xkc', '-','zooplankton net C growth efficiency', default=0.8_rk)
   call self%get_parameter(self%xthetazme, 'xthetazme', 'mol C mol N-1','mesozooplankton C:N ratio', default=5.625_rk)
   call self%get_parameter(self%xmetapn, 'xmetapn', 'd-1','phytoplankton loss rate (non-diatoms)', default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmetapd, 'xmetapd', 'd-1','phytoplankton loss rate (diatoms)', default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmetazmi, 'xmetazmi', 'd-1','microzooplankton loss rate', default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmetazme, 'xmetazme', 'd-1','mesozooplankton loss rate', default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmpn,'xmpn','d-1','phytoplankton maximum loss rate (non-diatoms)', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmpd,'xmpd','d-1','phytoplankton maximum loss rate (diatoms)', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmzmi,'xmzmi','d-1','microzooplankton maximum loss rate', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmzme,'xmzme','d-1','mesozooplankton maximum loss rate', default=0.2_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xkphn,'xkphn','mmol N m-3','phytoplankton loss half-saturation constant (non-diatoms)', default=0.5_rk)
   call self%get_parameter(self%xkphd,'xkphd','mmol N m-3','phytoplankton loss half-saturation constant (diatoms)', default=0.5_rk)
   call self%get_parameter(self%xkzmi,'xkzmi','mmol N m-3','microzooplankton loss half-saturation constant', default=0.5_rk)
   call self%get_parameter(self%xkzme,'xkzme','mmol N m-3','mesozooplankton loss half-saturation constant', default=0.75_rk)
   call self%get_parameter(self%jmpn,'jmpn','-','mortality formulation (non-diatoms): 1-linear, 2-quadratic, 3-hyperbolic, 4-sigmoid', default = 3)
   call self%get_parameter(self%jmpd,'jmpd','-','mortality formulation (diatoms): 1-linear, 2-quadratic, 3-hyperbolic, 4-sigmoid', default = 3)
   call self%get_parameter(self%jmzmi,'jmzmi','-','mortality formulation (non-diatoms): 1-linear, 2-quadratic, 3-hyperbolic, 4-sigmoid', default = 3)
   call self%get_parameter(self%jmzme,'jmzme','-','mortality formulation (non-diatoms): 1-linear, 2-quadratic, 3-hyperbolic, 4-sigmoid', default = 3)
   call self%get_parameter(self%jmd, 'jmd','-','Temperature regulation (detritus remin): 1-Eppley,2-q10',default=1)
   call self%get_parameter(self%xmd,'xmd','d-1','detrital N remineralisation rate', default=0.0158_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xmdc,'xmdc','d-1','detrital C remineralisation rate', default=0.0127_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsdiss,'xsdiss','d-1','diatom frustule dissolution rate', default=0.006_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xfdfrac1,'xfdfrac1','-','fast detritus fraction of diatom losses',default=0.33_rk)
   call self%get_parameter(self%xfdfrac2,'xfdfrac2','-','fast detritus fraction of mesozooplankton losses',default=1._rk)
   call self%get_parameter(self%xfdfrac3,'xfdfrac3','-','fast detritus fraction of mesozooplankton grazing',default=0.8_rk)
   call self%get_parameter(self%xrfn,'xrfn','mmol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=0.03_rk)
   call self%get_parameter(self%xridg_r0,'xridg_r0','-','CaCO3 : POC export rain ratio scalar, Ridgwell et al (2007)',default=0.026_rk)
   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)
   call self%get_parameter(self%wg,'wg','m d-1','detritus sinking rate (<0 for sinking)', default=-2.5_rk)
   call self%get_parameter(self%wdep,'wdep','m d-1','detritus deposition rate', default=2.5_rk,scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_ZCHN,'CHN','mg chl/m**3', 'chlorophyll in non-diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'CHD','mg chl/m**3', 'chlorophyll in diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPHN,'PHN','mmol N/m**3', 'non-diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPHD,'PHD','mmol N/m**3', 'diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPDS,'PDS','mmol Si/m**3', 'diatom phytoplankton (silicon)', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIN,'DIN','mmol N/m**3', 'nitrogen nutrient', minimum=0.0_rk) !no river dilution?
   call self%register_state_variable(self%id_ZFER,'FER','mmol Fe/m**3', 'iron nutrient', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSIL,'SIL','mmol Si/m**3', 'silicic acid', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDET,'DET','mmol N/m**3', 'slow-sinking detritus (N)',minimum=0.0_rk,vertical_movement=self%wg*d_per_s)
   call self%register_state_variable(self%id_ZDTC,'DTC','mmol C/m**3', 'slow-sinking detritus (C)',minimum=0.0_rk,vertical_movement=self%wg*d_per_s)
   call self%register_state_variable(self%id_ZZMI,'ZMI','mmol N/m**3', 'microzooplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZZME,'ZME','mmol N/m**3', 'mesozooplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZALK,'ALK','meq/m**3', 'total alkalinity', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIC,'DiC','mmol C/m**3', 'dissolved inorganic carbon', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZOXY,'OXY','mmol O_2/m**3', 'dissolved oxygen')

   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZPHN)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZPHD)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZDIN)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZDET)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZZMI)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZZME)

   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZPHN, scale_factor=self%xthetapn)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZPHD, scale_factor=self%xthetapd)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZZMI, scale_factor=self%xthetazmi)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZZME, scale_factor=self%xthetazme)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZDTC)
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZDIC)

   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_ZPDS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_ZSIL)

   call self%register_state_dependency(self%id_tempn,'tempn','mmol N/m**3', 'fast-sinking detritus (N)')
   call self%register_state_dependency(self%id_tempc,'tempc','mmol C/m**3', 'fast-sinking detritus (C)')
   call self%register_state_dependency(self%id_tempsi,'tempsi','mmol Si/m**3', 'fast-sinking detritus (Si)')
   call self%register_state_dependency(self%id_tempfe,'tempfe','mmol Fe/m**3', 'fast-sinking detritus (Fe)')
   call self%register_state_dependency(self%id_tempca,'tempca','mmol CaCO3/m**3', 'fast-sinking detritus (CaCO3)')

   call self%get_parameter(self%seafloor,'seafloor','-','seafloor handling: 1-inorganic returns, 2-organic returns, 3-coupled benthic model', default = 3)
   if (self%seafloor .eq. 3) then
         call self%register_state_dependency(self%id_ZSEDC,'BEN_C','mmol C m-2', 'sediment (C)')
         call self%register_state_dependency(self%id_ZSEDN,'BEN_N','mmol N m-2', 'sediment (N)')
         call self%register_state_dependency(self%id_ZSEDP,'BEN_P','mmol P m-2', 'sediment (P)')
         call self%register_state_dependency(self%id_ZSEDFE,'BEN_FE','mmol Fe m-2', 'sediment (Fe)')
         call self%request_coupling_to_model(self%id_ZSEDC, 'BEN', standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_ZSEDN, 'BEN', standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_ZSEDP, 'BEN', standard_variables%total_phosphorus)
   end if
   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)

   call self%register_dependency(self%id_om_cal,'OM_CAL3','-','calcite saturation')
   call self%register_dependency(self%id_xpar,standard_variables%downwelling_photosynthetic_radiative_flux)

   call self%register_diagnostic_variable(self%id_fscal_part,'fscal_part','nmol C cm-2 s-1','carbon in suspended particles')
   call self%register_diagnostic_variable(self%id_prn,'PRN','mmolN/m3/d','Non-diatom primary production')
   call self%register_diagnostic_variable(self%id_prd,'PRD','mmolN/m3/d','Diatom primary production')
   call self%register_diagnostic_variable(self%id_mpn,'MPN','mmolN/m3/d','Non-diatom non-grazing losses')
   call self%register_diagnostic_variable(self%id_mpd,'MPD','mmolN/m3/d','Diatom non-grazing losses')
   call self%register_diagnostic_variable(self%id_OPAL,'OPAL','mmolSi/m3/d','Diatom biogenic opal production')
   call self%register_diagnostic_variable(self%id_OPALDISS,'OPALDISS','mmolSi/m3/d','Diatom biogenic opal dissolution')
   call self%register_diagnostic_variable(self%id_GMIPn,'GMIPn','mmolN/m3/d','Microzoo grazing on non-diatoms')
   call self%register_diagnostic_variable(self%id_GMID,'GMID','mmolN/m3/d','Microzoo grazing on detritus')
   call self%register_diagnostic_variable(self%id_GMEPN,'GMEPN','mmolN/m3/d','Mesozoo grazing on non-diatoms')
   call self%register_diagnostic_variable(self%id_GMEPD,'GMEPD','mmolN/m3/d','Mesozoo grazing on diatoms')
   call self%register_diagnostic_variable(self%id_GMEZMI,'GMEZMI','mmolN/m3/d','Mesozoo grazing on microzoo')
   call self%register_diagnostic_variable(self%id_GMED,'GMED','mmolN/m3/d','Mesozoo grazing on detritus')
   call self%register_diagnostic_variable(self%id_MZMI,'MZMI','mmolN/m3/d','Microzoo non-grazing losses')
   call self%register_diagnostic_variable(self%id_MZME,'MZME','mmolN/m3/d','Mesozoo non-grazing losses')
   call self%register_diagnostic_variable(self%id_detn,'DETN','mmolN/m3/d','Slow detritus creation')
   call self%register_diagnostic_variable(self%id_detc,'DETC','mmolC/m3/d','Slow detritus C creation')
   call self%register_diagnostic_variable(self%id_MDET,'MDET','mmolN/m3/d','Detritus non-grazing losses')
   call self%register_diagnostic_variable(self%id_MDETC,'MDETC','mmolC/m3/d','Detritus non-grazing losses, carbon')
   call self%register_diagnostic_variable(self%id_GMIDC,'GMIDC','mmolC/m3/d','Microzoo grazing on detritus, carbon')
   call self%register_diagnostic_variable(self%id_GMEDC,'GMEDC','mmolC/m3/d','Mesozoo  grazing on detritus, carbon')
   call self%register_diagnostic_variable(self%id_fregen,'regen_slow','mmolN/m3/d','Total slow remin flux N 3D')
   call self%register_diagnostic_variable(self%id_fregensi,'regenSi_slow','mmolSi/m3/d','Total slow remin flux Si 3D')
   call self%register_diagnostic_variable(self%id_slowdetflux,'flux_slow','mmolN/m3/d','Total slow detrital flux 3D')
   call self%register_diagnostic_variable(self%id_pd_jlim,'PD_JLIM','-','Diatom light limitation')
   call self%register_diagnostic_variable(self%id_pd_nlim,'PD_NLIM','-','Diatom N limitation')
   call self%register_diagnostic_variable(self%id_pd_felim,'PD_FELIM','-','Diatom Fe limitation')
   call self%register_diagnostic_variable(self%id_pd_silim,'PD_SILIM','-','Diatom Si limitation')
   call self%register_diagnostic_variable(self%id_pd_silim2,'PD_SILIM2','-','Diatom Si uptake limitation')
   call self%register_diagnostic_variable(self%id_pn_jlim,'PN_JLIM','-','Non-diatom light limitation')
   call self%register_diagnostic_variable(self%id_pn_nlim,'PN_NLIM','-','Non-diatom N limitation')
   call self%register_diagnostic_variable(self%id_pn_felim,'PN_FELIM','-','Non-diatom Fe limitation')
   call self%register_diagnostic_variable(self%id_PN_LLOSS,'PN_LLOSS','mmolN/m3/d','Non-diatom linear losses')
   call self%register_diagnostic_variable(self%id_PD_LLOSS,'PD_LLOSS','mmolN/m3/d','Diatom linear losses')
   call self%register_diagnostic_variable(self%id_ZI_LLOSS,'ZI_LLOSS','mmolN/m3/d','Microzooplankton linear losses')
   call self%register_diagnostic_variable(self%id_ZE_LLOSS,'ZE_LLOSS','mmolN/m3/d','Mesozooplankton linear losses')
   call self%register_diagnostic_variable(self%id_ZI_MES_N,'ZI_MES_N','mmolN/m3/d','Microzoo messy feeding loss to N')
   call self%register_diagnostic_variable(self%id_ZI_MES_D,'ZI_MES_D','mmolN/m3/d','Microzoo messy feeding loss to D')
   call self%register_diagnostic_variable(self%id_ZI_MES_C,'ZI_MES_C','mmolC/m3/d','Microzoo messy feeding loss to C')
   call self%register_diagnostic_variable(self%id_ZI_MESDC,'ZI_MESDC','mmolC/m3/d','Microzoo messy feeding loss to Dc')
   call self%register_diagnostic_variable(self%id_ZE_MES_N,'ZE_MES_N','mmolN/m3/d','Mesozoo messy feeding loss to N')
   call self%register_diagnostic_variable(self%id_ZE_MES_D,'ZE_MES_D','mmolN/m3/d','Mesozoo messy feeding loss to D')
   call self%register_diagnostic_variable(self%id_ZE_MES_C,'ZE_MES_C','mmolC/m3/d','Mesozoo messy feeding loss to C')
   call self%register_diagnostic_variable(self%id_ZE_MESDC,'ZE_MESDC','mmolC/m3/d','Mesozoo messy feeding loss to Dc')
   call self%register_diagnostic_variable(self%id_ZI_EXCR,'ZI_EXCR','mmolN/m3/d','Microzoo excretion')
   call self%register_diagnostic_variable(self%id_ZI_RESP,'ZI_RESP','mmolN/m3/d','Microzoo respiration')
   call self%register_diagnostic_variable(self%id_ZI_GROW,'ZI_GROW','mmolN/m3/d','Microzoo growth')
   call self%register_diagnostic_variable(self%id_ZE_EXCR,'ZE_EXCR','mmolN/m3/d','Mesozoo excretion')
   call self%register_diagnostic_variable(self%id_C_PROD,'C_PROD','mmolC/m3/d','Carbon production')
   call self%register_diagnostic_variable(self%id_C_CONS,'C_CONS','mmolC/m3/d','Carbon consumption')
   call self%register_diagnostic_variable(self%id_N_PROD,'N_PROD','mmolN/m3/d','Nitogen production')
   call self%register_diagnostic_variable(self%id_N_CONS,'N_CONS','mmolN/m3/d','Nitrogen consumption')
   call self%register_diagnostic_variable(self%id_foxy_prod,'O2_PROD','mmolO2/m3/d','Oxygen production')
   call self%register_diagnostic_variable(self%id_foxy_cons,'O2_CONS','mmolO2/m3/d','Oxygen consumption')
   call self%register_diagnostic_variable(self%id_foxy_anox,'O2_ANOX','mmolO2/m3/d','Unrealised oxygen consumption')
   call self%register_diagnostic_variable(self%id_ZE_RESP,'ZE_RESP','mmolN/m3/d','Mesozoo respiration')
   call self%register_diagnostic_variable(self%id_ZE_GROW,'ZE_GROW','mmolN/m3/d','Mesozoo growth')
   call self%register_diagnostic_variable(self%id_fcomm_resp,'COM_RESP','mmolC/m3/d','Community respiration')
   call self%register_diagnostic_variable(self%id_f_sbenin_c,'sbenin_c','mmolC/m2/d','Benthic input slow carbon',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_f_sbenin_n,'sbenin_n','mmolN/m2/d','Benthic input slow nitrogen',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_f_sbenin_fe,'sbenin_fe','mmolFe/m2/d','Benthic input slow iron',source=source_do_bottom)

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   class(type_medusa_pelagic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

! !LOCAL VARIABLES:

    real(rk) :: ZCHN,ZCHD,ZPHN,ZPHD,ZPDS,ZDIN,ZFER,ZSIL,ZDET,ZDTC,ZZMI,ZZME,ZALK,ZDIC,ZOXY,loc_T,par
    real(rk) :: fthetan,fthetad,faln,fald !scaled chl/biomass ratio
    real(rk) :: fnln,ffln ! non-diatom Qn/Qf terms
    real(rk) :: fnld,fsld,ffld ! diatom Qn/Qs/Qf terms
    real(rk) :: fun_T,fun_q10,xvpnT,xvpdT,fchn1,fchn,fjln,fchd1,fchd,fjld
    real(rk) :: fsin,fnsi,fsin1,fnsi1,fnsi2,fprn,fprd,fsld2,frn,frd,fprds
    real(rk) :: fpnlim,fpdlim !nutrient limitation of primary production
    real(rk) :: fmi1,fmi,fgmipn,fgmid,fgmidc,finmi,ficmi,fstarmi,fmith,fmigrow,fmiexcr,fmiresp
    real(rk) :: fme1,fme,fgmepn,fgmepd,fgmepds,fgmezmi,fgmed,fgmedc,finme,ficme,fstarme,fmeth,fmegrow,fmeexcr,fmeresp
    real(rk) :: fdpn2,fdpd2,fdpds2,fdzmi2,fdzme2,fdpn,fdpd,fdpds,fdzmi,fdzme
    real(rk) :: fdd,fddc,fsdiss,fjlim_pn,fjlim_pd
    real(rk) :: fslowc,fslown,fregen,fregensi,fregenc,ftempn,ftempsi,ftempfe,ftempc,fq1,fcaco3,ftempca
    real(rk) :: fn_prod,fn_cons,fs_cons,fs_prod,fc_cons,fc_prod,fa_prod,fa_cons,fo2_ccons,fo2_ncons,fo2_cons,fo2_prod
    real(rk) :: rsmall, om_cal
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values

    _GET_(self%id_ZCHN,ZCHN)
    _GET_(self%id_ZCHD,ZCHD)
    _GET_(self%id_ZPHN,ZPHN)
    _GET_(self%id_ZPHD,ZPHD)
    _GET_(self%id_ZPDS,ZPDS)
    _GET_(self%id_ZDIN,ZDIN)
    _GET_(self%id_ZFER,ZFER)
    _GET_(self%id_ZSIL,ZSIL)
    _GET_(self%id_ZDET,ZDET)
    _GET_(self%id_ZDTC,ZDTC)
    _GET_(self%id_ZZMI,ZZMI)
    _GET_(self%id_ZZME,ZZME)
    _GET_(self%id_ZALK,ZALK)
    _GET_(self%id_ZDIC,ZDIC)
    _GET_(self%id_ZOXY,ZOXY)
    _GET_(self%id_temp,loc_T)
    _GET_(self%id_xpar,par)

   _SET_DIAGNOSTIC_(self%id_slowdetflux,-self%wg * ZDET)

   rsmall = 0.5_rk * EPSILON( 1._rk )

   !PHYTOPLANKTON GROWTH
   !Chlorophyll
   if (ZPHN .gt. rsmall) then
    fthetan = max(tiny(ZCHN),(ZCHN * self%xxi) / (ZPHN + tiny(ZPHN)))
    faln = self%xaln * fthetan
   else
    fthetan = 0._rk
    faln    = 0._rk
   end if

   if (ZPHD .gt. rsmall) then
    fthetad = max(tiny(ZCHD),(ZCHD * self%xxi) / (ZPHD + tiny(ZPHD)))
    fald = self%xald * fthetad
   else
    fthetad = 0._rk
    fald    = 0._rk
   end if

  !Temperature limitation
   fun_T = 1.066_rk**loc_T
   fun_q10 = self%jq10**((loc_T - 0._rk) / 10._rk)
   if (self%jphy .eq. 1) then
     xvpnT = self%xvpn * fun_T
     xvpdT = self%xvpd * fun_T
   elseif (self%jphy .eq. 2) then
     xvpnT = self%xvpn * fun_Q10
     xvpdT = self%xvpd * fun_Q10
   else
     xvpnT = self%xvpn
     xvpdT = self%xvpd
   endif
  !Phytoplankton light limitation
   fchn1 = (xvpnT * xvpnT) + (faln * faln * par * par)
   if (fchn1 .gt. rsmall) then
    fchn = xvpnT / (sqrt(fchn1) + tiny(fchn1))
   else
    fchn = 0._rk
   endif
   fjln = fchn * faln * par !non-diatom J term
   fjlim_pn = fjln / xvpnT

   fchd1 = (xvpdT * xvpdT) + (fald * fald * par * par)
   if (fchd1 .gt. rsmall) then
    fchd = xvpdT / (sqrt(fchd1) + tiny(fchd1))
   else
    fchd = 0._rk
   end if
   fjld = fchd * fald * par !diatom J term
   fjlim_pd = fjld / xvpdT

   ! Phytoplankton nutrient limitation
   !! Non-diatoms (N, Fe)
   fnln = ZDIN / (ZDIN + self%xnln) !non-diatom Qn term
   ffln = ZFER / (ZFER + self%xfln) !non-diatom Qf term
   !! Diatoms (N, Si, Fe)
   fnld = ZDIN / (ZDIN + self%xnld) !diatom Qn term
   fsld = ZSIL / (ZSIL + self%xsld) !diatom Qs term
   ffld = ZFER / (ZFER + self%xfld) !diatom Qf term

   ! Primary production (non-diatoms)
   if (self%jliebig.eqv..false.) then
     fpnlim = fnln * ffln
   else
     fpnlim = min(fnln,ffln)
   end if

   fprn = fjln * fpnlim

   !Primary production (diatoms)
   if (self%jliebig.eqv..false.) then
     fpdlim = fnld * ffld
   else
     fpdlim = min(fnld,ffld)
   end if

  if (ZPHD .gt. rsmall .and. ZPDS .gt. rsmall) then

   fsin = 0._rk
   if ( zphd .gt. rsmall) fsin = ZPDS / ZPHD
   fnsi = 0._rk
   if ( zpds .gt. rsmall) fnsi = ZPHD / ZPDS
   fsin1 = 3.0_rk * self%xsin0   ! = 0.6
   fnsi1 = 1.0_rk / fsin1        ! = 1.667
   fnsi2 = 1.0_rk / self%xsin0   ! = 5.0
   if (fsin .le. self%xsin0) then
      fprd = 0._rk
      fsld2 = 0._rk
   elseif (fsin .lt. fsin1) then
      fprd = self%xuif * ((fsin - self%xsin0) / (fsin + tiny(fsin))) * (fjld * fpdlim)
      fsld2 = self%xuif * ((fsin - self%xsin0) / (fsin + tiny(fsin)))
   elseif (fsin .ge. fsin1) then
      fprd = fjld * fpdlim
      fsld2 = 1.0_rk
   end if

  !Silicon
   if (fsin.lt.fnsi1) then
     fprds = (fjld * fsld)
   elseif (fsin.lt.fnsi2) then
     fprds = self%xuif * ((fnsi - self%xnsi0) / (fnsi + tiny(fnsi))) * (fjld * fsld)
   else
     fprds = 0._rk
   end if
 else
   fsin = 0._rk
   fnsi = 0._rk
   fprd = 0._rk
   fsld2 = 0._rk
   fprds = 0._rk
 end if

   _SET_DIAGNOSTIC_(self%id_pd_jlim, fjlim_pd * ZPHD)
   _SET_DIAGNOSTIC_(self%id_pd_nlim, fnld * ZPHD)
   _SET_DIAGNOSTIC_(self%id_pd_felim, ffld * ZPHD)
   _SET_DIAGNOSTIC_(self%id_pd_silim, fsld2 * ZPHD)
   _SET_DIAGNOSTIC_(self%id_pd_silim2, fsld * ZPHD)

   _SET_DIAGNOSTIC_(self%id_pn_jlim, fjlim_pn * ZPHN)
   _SET_DIAGNOSTIC_(self%id_pn_nlim, fnln * ZPHN)
   _SET_DIAGNOSTIC_(self%id_pn_felim, ffln * ZPHN)

   _SET_DIAGNOSTIC_(self%id_prn, (fprn * ZPHN) / d_per_s)
   _SET_DIAGNOSTIC_(self%id_prd, (fprd * ZPHD) / d_per_s)

  !Chlorophyll production
  frn = (self%xthetam * fchn * fnln * ffln) / (fthetan + tiny(fthetan))
  frd = (self%xthetamd * fchd * fnld * ffld * fsld2) / (fthetad + tiny(fthetad))  

  !ZOOPLANKTON GRAZING
  !Microzooplankton
  fmi1 = (self%xkmi * self%xkmi) + (self%xpmipn * ZPHN * ZPHN) + (self%xpmid * ZDET * ZDET)
  fmi = self%xgmi * ZZMI / fmi1
  fgmipn = fmi * self%xpmipn * ZPHN * ZPHN !grazing on non-diatoms
  fgmid = fmi * self%xpmid * ZDET * ZDET   !grazing on detrital nitrogen
  fgmidc = rsmall
  if (ZDET .gt. rsmall) fgmidc = (ZDTC / (ZDET + tiny(ZDET))) * fgmid !ROAM formulation

  _SET_DIAGNOSTIC_(self%id_GMIPn,fgmipn / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMID,fgmid / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMIDC,fgmidc / d_per_s)

  finmi = (1.0_rk - self%xphi) * (fgmipn + fgmid)
  ficmi = (1.0_rk - self%xphi) * ((self%xthetapn * fgmipn) + fgmidc)
  fstarmi = (self%xbetan * self%xthetazmi) / (self%xbetac * self%xkc) !the ideal food C : N ratio for microzooplankton
  fmith = (ficmi / (finmi + tiny(finmi)))
  if (fmith .ge. fstarmi) then
     fmigrow = self%xbetan * finmi
     fmiexcr = 0.0_rk
  else
     fmigrow = (self%xbetac * self%xkc * ficmi) / self%xthetazmi
     fmiexcr = ficmi * ((self%xbetan / (fmith + tiny(fmith))) - ((self%xbetac * self%xkc) / self%xthetazmi))
  end if
  fmiresp = (self%xbetac * ficmi) - (self%xthetazmi * fmigrow) !Respiration

  !Mesozooplankton
  fme1 = (self%xkme * self%xkme) + (self%xpmepn * ZPHN * ZPHN) + (self%xpmepd * ZPHD * ZPHD) + (self%xpmezmi * ZZMI * ZZMI) + (self%xpmed * ZDET * ZDET)
  fme = self%xgme * ZZME / fme1
  fgmepn = fme * self%xpmepn * ZPHN * ZPHN
  fgmepd = fme * self%xpmepd * ZPHD * ZPHD
  fgmepds = fsin * fgmepd
  fgmezmi = fme * self%xpmezmi * ZZMI * ZZMI
  fgmed = fme * self%xpmed * ZDET * ZDET
  fgmedc = rsmall
  if (ZDET .gt. rsmall) fgmedc = (ZDTC / (ZDET + tiny(ZDET))) * fgmed !ROAM formulation

  _SET_DIAGNOSTIC_(self%id_GMEPN,fgmepn / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMEPD,fgmepd / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMEZMI,fgmezmi / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMED,fgmed / d_per_s)
  _SET_DIAGNOSTIC_(self%id_GMEDC,fgmedc / d_per_s)

  finme = (1.0_rk - self%xphi) * (fgmepn + fgmepd + fgmezmi + fgmed)
  ficme = (1.0_rk - self%xphi) * ((self%xthetapn * fgmepn) + (self%xthetapd * fgmepd) + (self%xthetazmi * fgmezmi) + fgmedc) 
  fstarme = (self%xbetan * self%xthetazme) / (self%xbetac * self%xkc)
  fmeth = (ficme / (finme + tiny(finme)))
  if (fmeth .ge. fstarme) then
     fmegrow = self%xbetan * finme
     fmeexcr = 0.0_rk
  else
     fmegrow = (self%xbetac * self%xkc * ficme) / self%xthetazme
     fmeexcr = ficme * ((self%xbetan / (fmeth + tiny(fmeth))) - ((self%xbetac * self%xkc) / self%xthetazme))
  end if
  fmeresp = (self%xbetac * ficme) - (self%xthetazme * fmegrow)

  _SET_DIAGNOSTIC_(self%id_ZI_MES_N, self%xphi * (fgmipn + fgmid) / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_MES_D, (1._rk - self%xbetan) * finmi / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_MES_C, self%xphi * ((self%xthetapn * fgmipn) + fgmidc) / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_MESDC, (1._rk - self%xbetac) * ficmi / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_MES_N, self%xphi * (fgmepn + fgmepd + fgmezmi + fgmed) / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_MES_D, (1._rk - self%xbetan) * finme / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_MES_C, self%xphi * ((self%xthetapn * fgmepn) + (self%xthetapd * fgmepd) + (self%xthetazmi * fgmezmi) + fgmedc) / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_MESDC, (1._rk - self%xbetac) * ficme / d_per_s)

  _SET_DIAGNOSTIC_(self%id_ZI_EXCR, fmiexcr / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_RESP, fmiresp / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_GROW, fmigrow / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_EXCR, fmeexcr / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_RESP, fmeresp / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_GROW, fmegrow / d_per_s)

  !Plankton metabolic losses
  !Linear loss processes assumed to be metabolic in origin
  fdpn2 = self%xmetapn * ZPHN
  fdpd2 = self%xmetapd * ZPHD
  fdpds2 = self%xmetapd * ZPDS
  fdzmi2 = self%xmetazmi * ZZMI
  fdzme2 = self%xmetazme * ZZME

  _SET_DIAGNOSTIC_(self%id_PN_LLOSS, fdpn2  / d_per_s)
  _SET_DIAGNOSTIC_(self%id_PD_LLOSS, fdpd2 / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZI_LLOSS, fdzmi2 / d_per_s)
  _SET_DIAGNOSTIC_(self%id_ZE_LLOSS, fdzme2 / d_per_s)

  !Plankton mortality losses
  ! non-diatom phytoplankton
  if (self%jmpn == 1) then
      fdpn = self%xmpn * ZPHN                                     !! linear
  elseif (self%jmpn == 2) then
      fdpn = self%xmpn * ZPHN * ZPHN                              !! quadratic
  elseif (self%jmpn == 3) then
      fdpn = self%xmpn * ZPHN * (ZPHN / (self%xkphn + ZPHN))      !! hyperbolic
  elseif (self%jmpn == 4) then
      fdpn = self%xmpn * ZPHN * &                                 !! sigmoid
                 ((ZPHN * ZPHN) / (self%xkphn + (ZPHN * ZPHN)))
  end if
  ! diatom phytoplankton
  if (self%jmpd == 1) then
      fdpd = self%xmpd * ZPHD               !! linear
  elseif (self%jmpd == 2) then
      fdpd = self%xmpd * ZPHD * ZPHD        !! quadratic
  elseif (self%jmpd == 3) then
      fdpd = self%xmpd * ZPHD * &           !! hyperbolic
                  (ZPHD / (self%xkphd + ZPHD))
  elseif (self%jmpd == 4) then
      fdpd = self%xmpd * ZPHD * &           !! sigmoid
                  ((ZPHD * ZPHD) / (self%xkphd + (ZPHD * ZPHD)))
  end if
  fdpds = fdpd * fsin

  _SET_DIAGNOSTIC_(self%id_mpn, fdpn / d_per_s)
  _SET_DIAGNOSTIC_(self%id_mpd, fdpd / d_per_s)

  ! microzooplankton
  if (self%jmzmi == 1) then
     fdzmi = self%xmzmi * ZZMI            !! linear
  elseif (self%jmzmi == 2) then
     fdzmi = self%xmzmi * ZZMI * ZZMI     !! quadratic
  elseif (self%jmzmi == 3) then
     fdzmi = self%xmzmi * ZZMI * &        !! hyperbolic
                  (ZZMI / (self%xkzmi + ZZMI))
  elseif (self%jmzmi == 4) then
     fdzmi = self%xmzmi * ZZMI * &        !! sigmoid
                  ((ZZMI * ZZMI) / (self%xkzmi + (ZZMI * ZZMI)))
  end if
  ! mesozooplankton
  if (self%jmzme == 1) then
     fdzme = self%xmzme * ZZME            !! linear
  elseif (self%jmzme == 2) then
     fdzme = self%xmzme * ZZME * ZZME     !! quadratic
  elseif (self%jmzme == 3) then
     fdzme = self%xmzme * ZZME * &        !! hyperbolic
                  (ZZME / (self%xkzme + ZZME))
  elseif (self%jmzme == 4) then
     fdzme = self%xmzme * ZZME * &        !! sigmoid
                  ((ZZME * ZZME) / (self%xkzme + (ZZME * ZZME)))
  end if

  _SET_DIAGNOSTIC_(self%id_MZMI, fdpn / d_per_s)
  _SET_DIAGNOSTIC_(self%id_MZME, fdpd / d_per_s)

  !Detritus remineralisation (temperature-dependent)

  if (self%jmd == 1) then
    fdd = self%xmd * fun_T * ZDET
    fddc = self%xmdc * fun_T * ZDTC
  elseif (self%jmd == 2) then
    fdd  = self%xmd  * fun_Q10 * ZDET
    fddc = self%xmdc * fun_Q10 * ZDTC
  else
    fdd  = self%xmd  * ZDET
    fddc = self%xmdc * ZDTC
  end if

  _SET_DIAGNOSTIC_(self%id_MDET, fdd / d_per_s)
  _SET_DIAGNOSTIC_(self%id_MDETC, fddc / d_per_s)

  !Original contains accelerated detrital remineralisation in the bottom box

  !Diatom frustule dissolution
  fsdiss = self%xsdiss * ZPDS
  _SET_DIAGNOSTIC_(self%id_OPALDISS, fsdiss / d_per_s)

  !!Aeolian iron deposition and seafloor iron addition - should be dealt with through model inputs.

  !Slow detritus creation
  fslown  = fdpn + fdzmi + ((1._rk - self%xfdfrac1) * fdpd) + ((1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetan) * (finmi + finme))
  fslowc  = (self%xthetapn * fdpn) + (self%xthetazmi * fdzmi) + (self%xthetapd * (1._rk - self%xfdfrac1) * fdpd) + (self%xthetazme * (1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetac) * (ficmi + ficme))

  _SET_DIAGNOSTIC_(self%id_detn, fslown / d_per_s)
  _SET_DIAGNOSTIC_(self%id_detc, fslowc / d_per_s)

  !Nutrient regeneration !Should be saved as diagnostics
  fregen = (( (self%xphi * (fgmipn + fgmid)) +                            &  ! messy feeding
  (self%xphi * (fgmepn + fgmepd + fgmezmi + fgmed)) +                     &  ! messy feeding
  fmiexcr + fmeexcr + fdd +                                               &  ! excretion + D remin.
  fdpn2 + fdpd2 + fdzmi2 + fdzme2))                                          ! linear mortality
 
 _SET_DIAGNOSTIC_(self%id_fregen,fregen / d_per_s)

  !silicon
  fregensi = (( fsdiss + ((1._rk - self%xfdfrac1) * fdpds) +              &  ! dissolution + non-lin. mortality
  ((1._rk - self%xfdfrac3) * fgmepds) +                                   &  ! egestion by zooplankton
  fdpds2))                                                                   ! linear mortality

 _SET_DIAGNOSTIC_(self%id_fregensi,fregensi / d_per_s)

  !carbon
!  fregenc  = (( (self%xphi * ((self%xthetapn * fgmipn) + fgmidc)) +       &  ! messy feeding
!  (self%xphi * ((self%xthetapn * fgmepn) + (self%xthetapd * fgmepd) +     &  ! messy feeding
! (self%xthetazmi * fgmezmi) + fgmedc)) +                                 &  ! messy feeding
!  fmiresp + fmeresp + fddc +                                              &  ! respiration + D remin.
!  (self%xthetapn * fdpn2) + (self%xthetapd * fdpd2) +                     &  ! linear mortality
!  (self%xthetazmi * fdzmi2) + (self%xthetazme * fdzme2)))                    ! linear mortality

  ! Fast-sinking detritus terms
  ! nitrogen:   diatom and mesozooplankton mortality
  ftempn = (self%xfdfrac1 * fdpd)  + (self%xfdfrac2 * fdzme)
  _SET_ODE_(self%id_tempn,ftempn)
  ! silicon:    diatom mortality and grazed diatoms
  ftempsi = (self%xfdfrac1 * fdpds) + (self%xfdfrac3 * fgmepds)
  _SET_ODE_(self%id_tempsi,ftempsi)
  ! iron:       diatom and mesozooplankton mortality
  ftempfe = ((self%xfdfrac1 * fdpd) + (self%xfdfrac2 * fdzme)) * self%xrfn
  _SET_ODE_(self%id_tempfe,ftempfe)
  ! carbon:     diatom and mesozooplankton mortality
  ftempc = (self%xfdfrac1 * self%xthetapd * fdpd) + (self%xfdfrac2 * self%xthetazme * fdzme)
  _SET_ODE_(self%id_tempc,ftempc)

  ! CaCO3: Ridgwell et al. (2007) submodel, uses FULL 3D omega calcite to regulate rain ratio

  _GET_(self%id_om_cal,om_cal)
 ! if (om_cal .ge. 1._rk) then !get f3_omcal!
 !    fq1 = (om_cal - 1._rk)**0.81_rk
 ! else
    fq1 = 0._rk
 ! endif
  fcaco3 = self%xridg_r0 * fq1
  ftempca = ftempc * fcaco3
  _SET_ODE_(self%id_tempca,ftempca)

  !LOCAL SMS TRENDS
  ! chlorophyll
  _SET_ODE_(self%id_ZCHN,((frn * fprn * ZPHN) - fgmipn - fgmepn - fdpn - fdpn2) * (fthetan / self%xxi))
  _SET_ODE_(self%id_ZCHD,((frd * fprd * ZPHD) - fgmepd - fdpd - fdpd2) * (fthetad / self%xxi))

  ! phytoplankton
  _SET_ODE_(self%id_ZPHN,(fprn * ZPHN) - fgmipn - fgmepn - fdpn - fdpn2 )
  _SET_ODE_(self%id_ZPHD,(fprd * ZPHD) - fgmepd - fdpd - fdpd2 )
  _SET_ODE_(self%id_ZPDS,(fprds * ZPDS) - fgmepds - fdpds - fsdiss - fdpds2 )

  _SET_DIAGNOSTIC_(self%id_OPAL, (fprds * ZPDS) / d_per_s)

  ! zooplankton
  _SET_ODE_(self%id_ZZMI, fmigrow - fgmezmi - fdzmi - fdzmi2 )
  _SET_ODE_(self%id_ZZME, fmegrow - fdzme - fdzme2 )

  ! detritus
  _SET_ODE_(self%id_ZDET,fslown-fgmid-fgmed-fdd)

  !dissolved inorganic nitrogen
   fn_cons = - (fprn * ZPHN) - (fprd * ZPHD)                         ! primary production
   fn_prod = + (self%xphi * (fgmipn + fgmid))                     &  ! messy feeding remin.
             + (self%xphi * (fgmepn + fgmepd + fgmezmi + fgmed))  &  ! messy feeding remin.
             + fmiexcr + fmeexcr + fdd                            &  ! excretion and remin.
             + fdpn2 + fdpd2 + fdzmi2 + fdzme2                       ! metab. losses

   _SET_DIAGNOSTIC_(self%id_N_PROD, fn_prod / d_per_s)
   _SET_DIAGNOSTIC_(self%id_N_CONS, fn_cons / d_per_s)
   _SET_ODE_(self%id_ZDIN,fn_prod + fn_cons)

  ! dissolved silicic acid
   fs_cons = - (fprds * ZPDS)                                       ! opal production
   fs_prod = + fsdiss                                             &  ! opal dissolution
             + ((1.0_rk - self%xfdfrac1) * fdpds)                    &  ! mort. loss
             + ((1.0_rk - self%xfdfrac3) * fgmepds)                  &  ! egestion of grazed Si
             + fdpds2                                                ! metab. losses
   _SET_ODE_(self%id_ZSIL,fs_prod + fs_cons)

  ! dissolved iron

   _SET_ODE_(self%id_ZFER, self%xrfn * (fn_prod + fn_cons))


  ! detrital carbon
   _SET_ODE_(self%id_ZDTC, fslowc - fgmidc - fgmedc - fddc)
                 
  ! dissolved inorganic carbon
   fc_cons = - (self%xthetapn * fprn * ZPHN) - (self%xthetapd * fprd * ZPHD)                      ! primary production
   fc_prod = + (self%xthetapn * self%xphi * fgmipn) + (self%xphi * fgmidc)                     &  ! messy feeding remin
             + (self%xthetapn * self%xphi * fgmepn) + (self%xthetapd * self%xphi * fgmepd)     &  ! messy feeding remin
             + (self%xthetazmi * self%xphi * fgmezmi) + (self%xphi * fgmedc)                   &  ! messy feeding remin
             + fmiresp + fmeresp + fddc                                                        &
             + (self%xthetapn * fdpn2)                                                         &  ! resp., remin., losses
             + (self%xthetapd * fdpd2) + (self%xthetazmi * fdzmi2)                             &  ! losses
             + (self%xthetazme * fdzme2)                                                          ! losses

   _SET_DIAGNOSTIC_(self%id_fscal_part, (self%xthetapn * ZPHN + self%xthetapd * ZPHD + self%xthetazmi * ZZMI + self%xthetazme * ZZME + self%xthetad * ZDET) * 0.002_rk)

   _SET_DIAGNOSTIC_(self%id_fcomm_resp, fc_prod / d_per_s) !! community respiration (does not include CaCO3 terms)

   fc_prod = fc_prod - ftempca         ! CaCO3

   _SET_DIAGNOSTIC_(self%id_C_PROD, fc_prod / d_per_s)
   _SET_DIAGNOSTIC_(self%id_C_CONS, fc_cons / d_per_s)
   _SET_ODE_(self%id_ZDIC,fc_prod + fc_cons)

  ! alkalinity

   fa_cons = -2._rk * ftempca                                                    ! CaCO3 production

  _SET_ODE_(self%id_ZALK, fa_cons)

  ! oxygen
   fo2_prod = + (self%xthetanit * fprn * ZPHN)                                            & ! Pn primary production, N
              + (self%xthetanit * fprd * ZPHD)                                            & ! Pd primary production, N
              + (self%xthetarem * self%xthetapn * fprn * ZPHN)                            & ! Pn primary production, C
              + (self%xthetarem * self%xthetapd * fprd * ZPHD)                              ! Pd primary production, C
   fo2_ncons = - (self%xthetanit * self%xphi * fgmipn)                                    & ! Pn messy feeding remin., N
               - (self%xthetanit * self%xphi * fgmid)                                     & ! D  messy feeding remin., N
               - (self%xthetanit * self%xphi * fgmepn)                                    & ! Pn messy feeding remin., N
               - (self%xthetanit * self%xphi * fgmepd)                                    & ! Pd messy feeding remin., N
               - (self%xthetanit * self%xphi * fgmezmi)                                   & ! Zi messy feeding remin., N
               - (self%xthetanit * self%xphi * fgmed)                                     & ! D  messy feeding remin., N
               - (self%xthetanit * fmiexcr)                                               & ! microzoo excretion, N
               - (self%xthetanit * fmeexcr)                                               & ! mesozoo  excretion, N
               - (self%xthetanit * fdd)                                                   & ! slow detritus remin., N
               - (self%xthetanit * fdpn2)                                                 & ! Pn  losses, N
               - (self%xthetanit * fdpd2)                                                 & ! Pd  losses, N
               - (self%xthetanit * fdzmi2)                                                & ! Zmi losses, N
               - (self%xthetanit * fdzme2)                                                  ! Zme losses, N

   fo2_ccons = - (self%xthetarem * self%xthetapn * self%xphi * fgmipn)                    & ! Pn messy feeding remin., C
               - (self%xthetarem * self%xphi * fgmidc)                                    & ! D  messy feeding remin., C
               - (self%xthetarem * self%xthetapn * self%xphi * fgmepn)                    & ! Pn messy feeding remin., C
               - (self%xthetarem * self%xthetapd * self%xphi * fgmepd)                    & ! Pd messy feeding remin., C
               - (self%xthetarem * self%xthetazmi * self%xphi * fgmezmi)                  & ! Zi messy feeding remin., C
               - (self%xthetarem * self%xphi * fgmedc)                                    & ! D  messy feeding remin., C
               - (self%xthetarem * fmiresp)                                               & ! microzoo respiration, C
               - (self%xthetarem * fmeresp)                                               & ! mesozoo  respiration, C
               - (self%xthetarem * fddc)                                                  & ! slow detritus remin., C
               - (self%xthetarem * self%xthetapn * fdpn2)                                 & ! Pn  losses, C
               - (self%xthetarem * self%xthetapd * fdpd2)                                 & ! Pd  losses, C
               - (self%xthetarem * self%xthetazmi * fdzmi2)                               & ! Zmi losses, C
               - (self%xthetarem * self%xthetazme * fdzme2)                                 ! Zme losses, C

   fo2_cons = fo2_ncons + fo2_ccons

   if (ZOXY .lt. self%xo2min) then                     ! deficient O2; production fluxes only
      _SET_ODE_(self%id_ZOXY, fo2_prod )
      _SET_DIAGNOSTIC_(self%id_foxy_anox, fo2_cons / d_per_s)
      _SET_DIAGNOSTIC_(self%id_foxy_prod, fo2_prod / d_per_s)
   else                                                ! sufficient O2; production + consumption fluxes
      _SET_ODE_(self%id_ZOXY, fo2_prod + fo2_cons )
      _SET_DIAGNOSTIC_(self%id_foxy_prod, fo2_prod / d_per_s)
      _SET_DIAGNOSTIC_(self%id_foxy_cons, fo2_cons / d_per_s)
   endif

   _LOOP_END_

   end subroutine do

   subroutine do_bottom(self,_ARGUMENTS_DO_)

   class(type_medusa_pelagic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_

! !LOCAL VARIABLES

    real(rk) :: ZSEDC,ZSEDN,ZSEDFE,ZDET,ZDTC
    real(rk) :: fluxc,fluxn,fluxfe
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _HORIZONTAL_LOOP_BEGIN_

    _GET_(self%id_ZDET,ZDET)
    _GET_(self%id_ZDTC,ZDTC)

     fluxc = self%wdep * ZDTC
     fluxn = self%wdep * ZDET
     fluxfe = self%wdep * ZDET * self%xrfn

    _SET_BOTTOM_EXCHANGE_(self%id_ZDET,-fluxn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDTC,-fluxc)

    _SET_BOTTOM_ODE_(self%id_ZSEDC,  + fluxc)
    _SET_BOTTOM_ODE_(self%id_ZSEDN,  + fluxn)
    _SET_BOTTOM_ODE_(self%id_ZSEDFE, + fluxfe)

     _SET_BOTTOM_ODE_(self%id_ZSEDP, + fluxc/106._rk)

    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_sbenin_c,  + fluxc / d_per_s)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_sbenin_n,  + fluxn / d_per_s)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_sbenin_fe, + fluxfe / d_per_s)

    _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

  end module medusa_pelagic
