#include "fabm_driver.h"

!
!*********************************************************    
!                FABM-MEDUSA pelagic module                      
!*********************************************************                                                                            

module medusa_pelagic

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_pelagic
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZCHN,id_ZCHD,id_ZPHN,id_ZPHD,id_ZPDS,id_ZDIN,id_ZFER,id_ZSIL,id_ZDET,id_ZDTC,id_ZZMI,id_ZZME,id_ZDIC,id_ZALK,id_ZOXY
      type (type_dependency_id)            :: id_temp,id_par,id_depth,id_salt,id_dz
      type (type_diagnostic_variable_id)   :: id_dPAR,id_ftempc
      type (type_horizontal_diagnostic_variable_id) ::  id_fair
      type (type_horizontal_dependency_id)    :: id_apress,id_wnd

      ! Parameters
      logical :: jliebig
      real(rk) :: xxi,xaln,xald,xnln,xfln,xnld,xsld,xfld,xvpn,xvpd,xsin0,xnsi0,xuif,xthetam,xthetamd
      real(rk) :: xkmi,xpmipn,xpmid,xkme,xpmed,xpmepn,xpmepd,xpmezmi,zpmed,xgmi,xgme
      real(rk) :: xthetad,xphi,xthetapn,xthetazme,xthetazmi,xthetapd,xbetan,xbetac,xkc,xridg_r0,jq10
      integer :: jmpn,jmpd,jmzmi,jmzme,jphy,jmd,jiron
      real(rk) :: xmetapn,xmetapd,xmetazmi,xmetazme,xmpn,xmpd,xmzmi,xmzme,xkphn,xkphd,xkzmi,xkzme
      real(rk) :: xmd,xmdc,xsdiss
      real(rk) :: xk_FeL,xLgT,xk_sc_Fe
      real(rk) :: xfdfrac1,xfdfrac2,xfdfrac3,xrfn
      real(rk) :: xthetanit,xthetarem,xo2min
      real(rk) :: wg

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: fast_detritus

  end type

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
   call self%get_parameter(self%xfln, 'xfln', 'mmol Fe m-3','Fe nutrient uptake half-saturation constant (non-diatoms)', default=0.33_rk)
   call self%get_parameter(self%xnld, 'xnld', 'mmol N m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.75_rk)
   call self%get_parameter(self%xsld, 'xsld', 'mmol Si m-3','Si nutrient uptake half-saturation constant (diatoms)', default=3.0_rk)
   call self%get_parameter(self%xfld, 'xfld', 'mmol Fe m-3','Fe nutrient uptake half-saturation constant (diatoms)', default=0.67_rk)
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
   call self%get_parameter(self%xk_FeL,'xk_FeL','-','dissociation constant for (Fe+ligand)',default=100.0_rk)
   call self%get_parameter(self%xLgT,'xLgT','umol m-3','total ligand concentration',default=1.0_rk)
   call self%get_parameter(self%xk_sc_Fe,'xk_sc_Fe','d-1','scavenging rate of "free" Fe',default=0.001_rk,scale_factor=d_per_s)
   call self%get_parameter(self%jiron,'jiron','-','iron scavenging scheme: 1-Dutkiewicz et al. (2005),2-Moore et al. (2004),3-Moore et al. (2008),4-Galbraith et al. (2010)',default=1)
   call self%get_parameter(self%xfdfrac1,'xfdfrac1','-','fast detritus fraction of diatom losses',default=0.33_rk)
   call self%get_parameter(self%xfdfrac2,'xfdfrac2','-','fast detritus fraction of mesozooplankton losses',default=1._rk)
   call self%get_parameter(self%xfdfrac3,'xfdfrac3','-','fast detritus fraction of mesozooplankton grazing',default=0.8_rk)
   call self%get_parameter(self%xrfn,'xrfn','umol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=30._rk)
   call self%get_parameter(self%xridg_r0,'xridg_r0','-','CaCO3 : POC export rain ratio scalar, Ridgwell et al (2007)',default=0.026_rk)
   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)
   call self%get_parameter(self%wg,'wg','m d-1','detritus sinking rate (<0 for sinking)', default=-2.5_rk, scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_ZCHN,'ZCHN','mg chl/m**3', 'chlorophyll in non-diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZCHD','mg chl/m**3', 'chlorophyll in diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPHN,'ZPHN','mmol N/m**3', 'non-diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPHD,'ZPHD','mmol N/m**3', 'diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZPDS,'ZPDS','mmol Si/m**3', 'diatom phytoplankton (silicon)', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIN,'ZDIN','mmol N/m**3', 'nitrogen nutrient', minimum=0.0_rk) !no river dilution?
   call self%register_state_variable(self%id_ZFER,'ZFER','mmol Fe/m**3', 'iron nutrient', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSIL,'ZSIL','mmol Si/m**3', 'silicic acid', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDET,'ZDET','mmol N/m**3', 'slow-sinking detritus (N)', minimum=0.0_rk,vertical_movement=self%wg)
   call self%register_state_variable(self%id_ZDTC,'ZDTC','mmol C/m**3', 'slow-sinking detritus (C)', minimum=0.0_rk,vertical_movement=self%wg)
   call self%register_state_variable(self%id_ZZMI,'ZZMI','mmol N/m**3', 'microzooplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZZME,'ZZME','mmol N/m**3', 'mesozooplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZALK,'ZALK','meq/m**3', 'total alkalinity', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIC,'ZDIC','mmol C/m**3', 'dissolved inorganic carbon', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZOXY,'ZOXY','mmol O_2/m**3', 'dissolved oxygen', minimum=0.0_rk)
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W m-2',       'photosynthetically active radiation', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fair,'fair','mmol O_2/m^2/d','Air-sea flux of oxygen')
   call self%register_diagnostic_variable(self%id_ftempc,'ftempc','mmol C/m^3/d','fast detritus carbon')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_apress, standard_variables%surface_air_pressure)
   call self%register_dependency(self%id_wnd,standard_variables%wind_speed)
   call self%register_dependency(self%id_depth, standard_variables%depth)
   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   class(type_medusa_pelagic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

! !LOCAL VARIABLES:

    real(rk) :: ZCHN,ZCHD,ZPHN,ZPHD,ZPDS,ZDIN,ZFER,ZSIL,ZDET,ZDTC,ZZMI,ZZME,ZALK,ZDIC,ZOXY,loc_T,par,depth
    real(rk) :: fthetan,fthetad,faln,fald !scaled chl/biomass ratio
    real(rk) :: fnln,ffln ! non-diatom Qn/Qf terms
    real(rk) :: fnld,fsld,ffld ! diatom Qn/Qs/Qf terms
    real(rk) :: fun_T,fun_q10,xvpnT,xvpdT,fchn1,fchn,fjln,fchd1,fchd,fjld
    real(rk) :: fsin,fnsi,fsin1,fnsi1,fnsi2,fprn,fprd,fsld2,frn,frd,fprds
    real(rk) :: fpnlim,fpdlim !nutrient limitation of primary production
    real(rk) :: fmi1,fmi,fgmipn,fgmid,fgmidc,finmi,ficmi,fstarmi,fmith,fmigrow,fmiexcr,fmiresp
    real(rk) :: fme1,fme,fgmepn,fgmepd,fgmepds,fgmezmi,fgmed,fgmedc,finme,ficme,fstarme,fmeth,fmegrow,fmeexcr,fmeresp
    real(rk) :: fdpn2,fdpd2,fdpds2,fdzmi2,fdzme2,fdpn,fdpd,fdpds,fdzmi,fdzme
    real(rk) :: fdd,fddc,fsdiss
    real(rk) :: xFeT,xb_coef_tmp,xb2M4ac,xLgF,xFel,xFeF,xFree,ffescav,xmaxFeF,fdeltaFe
    real(rk) :: fslowc,fslown,fregen,fregensi,fregenc,ftempn,ftempsi,ftempfe,ftempc,fq1,fcaco3,ftempca
    real(rk) :: fn_prod,fn_cons,fs_cons,fs_prod,fc_cons,fc_prod,fa_prod,fa_cons,fo2_ccons,fo2_ncons,fo2_cons,fo2_prod

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
    _GET_(self%id_par,par) !check PAR // what about self-shading?
    _GET_(self%id_depth,depth)

   !PHYTOPLANKTON GROWTH
   !Chlorophyll
   fthetan = max(tiny(ZCHN),(ZCHN * self%xxi) / (ZPHN + tiny(ZPHN)))
   faln = self%xaln * fthetan
   fthetad = max(tiny(ZCHD),(ZCHD * self%xxi) / (ZPHD + tiny(ZPHD)))
   fald = self%xald * fthetad

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
   fchn = xvpnT / (sqrt(fchn1) + tiny(fchn1))
   fjln = fchn * faln * par !non-diatom J term

   fchd1 = (xvpdT * xvpdT) + (fald * fald * par * par)
   fchd = xvpdT / (sqrt(fchd1) + tiny(fchd1))
   fjld = fchd * fald * par !diatom J term

   ! Phytoplankton nutrient limitation
   !! Non-diatoms (N, Fe)
   fnln = ZDIN / (ZDIN + self%xnln) !non-diatom Qn term
   ffln = 1. !ZFER / (ZFER + self%xfln) !non-diatom Qf term
   !! Diatoms (N, Si, Fe)
   fnld = ZDIN / (ZDIN + self%xnld) !diatom Qn term
   fsld = ZSIL / (ZSIL + self%xsld) !diatom Qs term
   ffld = 1. !ZFER / (ZFER + self%xfld) !diatom Qf term

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

   fsin = ZPDS / ZPHD
   fnsi = ZPHD / ZPDS
   fsin1 = 3.0_rk * self%xsin0 !! = 0.6
   fnsi1 = 1.0_rk / fsin1      !! = 1.667
   fnsi2 = 1.0_rk / self%xsin0 !! = 5.0
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

  !Chlorophyll production
  frn = (self%xthetam * fchn * fnln * ffln) / (fthetan + tiny(fthetan))
  frd = (self%xthetamd * fchd * fnld * ffld * fsld2) / (fthetad + tiny(fthetad))  

  !ZOOPLANKTON GRAZING
  !Microzooplankton
  fmi1 = (self%xkmi * self%xkmi) + (self%xpmipn * ZPHN * ZPHN) + (self%xpmid * ZDET * ZDET)
  fmi = self%xgmi * ZZMI / fmi1
  fgmipn = fmi * self%xpmipn * ZPHN * ZPHN !grazing on non-diatoms
  fgmid = fmi * self%xpmid * ZDET * ZDET   !grazing on detrital nitrogen
  fgmidc = (ZDTC / (ZDET + tiny(ZDET))) * fgmid !ROAM formulation
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
  fgmedc = (ZDTC / (ZDET + tiny(ZDET))) * fgmed !ROAM formulation
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

  !Plankton metabolic losses
  !Linear loss processes assumed to be metabolic in origin
  fdpn2 = self%xmetapn * ZPHN
  fdpd2 = self%xmetapd * ZPHD
  fdpds2 = self%xmetapd * ZPDS
  fdzmi2 = self%xmetazmi * ZZMI
  fdzme2 = self%xmetazme * ZZME

  !Plankton mortality losses
  ! non-diatom phytoplankton
  if (self%jmpn == 1) fdpn = self%xmpn * ZPHN                                     !! linear
  if (self%jmpn == 2) fdpn = self%xmpn * ZPHN * ZPHN                              !! quadratic
  if (self%jmpn == 3) fdpn = self%xmpn * ZPHN * (ZPHN / (self%xkphn + ZPHN))      !! hyperbolic
  if (self%jmpn == 4) fdpn = self%xmpn * ZPHN * &                                 !! sigmoid
                 ((ZPHN * ZPHN) / (self%xkphn + (ZPHN * ZPHN)))
  ! diatom phytoplankton
  if (self%jmpd == 1) fdpd = self%xmpd * ZPHD               !! linear
  if (self%jmpd == 2) fdpd = self%xmpd * ZPHD * ZPHD        !! quadratic
  if (self%jmpd == 3) fdpd = self%xmpd * ZPHD * &           !! hyperbolic
                  (ZPHD / (self%xkphd + ZPHD))
  if (self%jmpd == 4) fdpd = self%xmpd * ZPHD * &           !! sigmoid
                  ((ZPHD * ZPHD) / (self%xkphd + (ZPHD * ZPHD)))
          fdpds = fdpd * fsin
  ! microzooplankton
  if (self%jmzmi == 1) fdzmi = self%xmzmi * ZZMI            !! linear
  if (self%jmzmi == 2) fdzmi = self%xmzmi * ZZMI * ZZMI     !! quadratic
  if (self%jmzmi == 3) fdzmi = self%xmzmi * ZZMI * &        !! hyperbolic
                  (ZZMI / (self%xkzmi + ZZMI))
  if (self%jmzmi == 4) fdzmi = self%xmzmi * ZZMI * &        !! sigmoid
                  ((ZZMI * ZZMI) / (self%xkzmi + (ZZMI * ZZMI)))
  ! mesozooplankton
  if (self%jmzme == 1) fdzme = self%xmzme * ZZME            !! linear
  if (self%jmzme == 2) fdzme = self%xmzme * ZZME * ZZME     !! quadratic
  if (self%jmzme == 3) fdzme = self%xmzme * ZZME * &        !! hyperbolic
                  (ZZME / (self%xkzme + ZZME))
  if (self%jmzme == 4) fdzme = self%xmzme * ZZME * &        !! sigmoid
                  ((ZZME * ZZME) / (self%xkzme + (ZZME * ZZME)))

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

  !Original contains accelerated detrital remineralisation in the bottom box (how do I let model know it is a bottom box?)

  !Diatom frustule dissolution
  fsdiss = self%xsdiss * ZPDS

  !IRON CHEMISTRY AND FRACTIONATION
  xFeT = ZFER * 1.e3_rk !total iron concentration (mmolFe/m3 -> umolFe/m3)
  xb_coef_tmp = self%xk_FeL * (self%xLgT - xFeT) - 1.0_rk !this is F1 in Yool et al (2013)
  xb2M4ac = max(((xb_coef_tmp * xb_coef_tmp) + (4.0_rk * self%xk_FeL * self%xLgT)), 0._rk)
  xLgF = 0.5_rk * (xb_coef_tmp + (xb2M4ac**0.5_rk)) / self%xk_FeL ! "free" ligand concentration
  xFeL = self%xLgT - xLgF ! ligand-bound iron concentration
  xFeF = (xFeT - xFeL) * 1.e-3_rk! "free" iron concentration (and convert to mmolFe/m3)
  xFree = xFeF / ZFER
  ! Scavenging of iron
  if (self%jiron == 1) then
     ffescav = self%xk_sc_Fe * xFeF
     xmaxFeF = min((xFeF * 1.e3_rk), 0.3_rk)        ! = umol/m3
     fdeltaFe = (xFeT - (xFeL + xmaxFeF)) * 1.e-3   ! = mmol/m3
     ffescav     = ffescav + fdeltaFe               ! = mmol/m3/d

     if ((depth .gt. 1000._rk) .and. (xFeT .lt. 0.5_rk)) then
        ffescav = 0._rk
     endif
  elseif (self%jiron == 2) then
  elseif (self%jiron == 3) then
  elseif (self%jiron == 4) then
  else
     ffescav = 0._rk
  end if


  !!Aeolian iron deposition and seafloor iron addition - should be dealt with through model inputs.

  !Slow detritus creation
  fslown  = fdpn + fdzmi + ((1._rk - self%xfdfrac1) * fdpd) + ((1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetan) * (finmi + finme))
  fslowc  = (self%xthetapn * fdpn) + (self%xthetazmi * fdzmi) + (self%xthetapd * (1._rk - self%xfdfrac1) * fdpd) + (self%xthetazme * (1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetac) * (ficmi + ficme))

  !Nutrient regeneration !Not necessary at all, should be saved as diagnostic of 
  !nitrogen
  fregen = (( (self%xphi * (fgmipn + fgmid)) +                            &  ! messy feeding
  (self%xphi * (fgmepn + fgmepd + fgmezmi + fgmed)) +                     &  ! messy feeding
  fmiexcr + fmeexcr + fdd +                                               &  ! excretion + D remin.
  fdpn2 + fdpd2 + fdzmi2 + fdzme2))                                          ! linear mortality
  !silicon
  fregensi = (( fsdiss + ((1._rk - self%xfdfrac1) * fdpds) +              &  ! dissolution + non-lin. mortality
  ((1._rk - self%xfdfrac3) * fgmepds) +                                   &  ! egestion by zooplankton
  fdpds2))                                                                   ! linear mortality
  !carbon
  fregenc  = (( (self%xphi * ((self%xthetapn * fgmipn) + fgmidc)) +       &  ! messy feeding
  (self%xphi * ((self%xthetapn * fgmepn) + (self%xthetapd * fgmepd) +     &  ! messy feeding
  (self%xthetazmi * fgmezmi) + fgmedc)) +                                 &  ! messy feeding
  fmiresp + fmeresp + fddc +                                              &  ! respiration + D remin.
  (self%xthetapn * fdpn2) + (self%xthetapd * fdpd2) +                     &  ! linear mortality
  (self%xthetazmi * fdzmi2) + (self%xthetazme * fdzme2)))                    ! linear mortality

  ! Fast-sinking detritus terms
  ! nitrogen:   diatom and mesozooplankton mortality
  ftempn = (self%xfdfrac1 * fdpd)  + (self%xfdfrac2 * fdzme)
  ! silicon:    diatom mortality and grazed diatoms
  ftempsi = (self%xfdfrac1 * fdpds) + (self%xfdfrac3 * fgmepds)
  ! iron:       diatom and mesozooplankton mortality
  ftempfe = ((self%xfdfrac1 * fdpd) + (self%xfdfrac2 * fdzme)) * self%xrfn
  ! carbon:     diatom and mesozooplankton mortality
  ftempc = (self%xfdfrac1 * self%xthetapd * fdpd) + (self%xfdfrac2 * self%xthetazme * fdzme)
  _SET_DIAGNOSTIC_(self%id_ftempc,ftempc)

  ! CaCO3: Ridgwell et al. (2007) submodel, uses FULL 3D omega calcite to regulate rain ratio
  !if (f3_omcal .ge. 1._rk) then !get f3_omcal!
  !   fq1 = (f3_omcal - 1._rk)**0.81
  !else
  !   fq1 = 0.
  !endif
  !   fcaco3 = self%xridg_r0 * fq1
  !ftempca = ftempc * fcaco3

  !Fast-sinking detritus magic...

  !LOCAL SMS TRENDS
  ! chlorophyll
  _SET_ODE_(self%id_ZCHN,((frn * fprn * ZPHN) - fgmipn - fgmepn - fdpn - fdpn2) * (fthetan / self%xxi))
  _SET_ODE_(self%id_ZCHD,((frd * fprd * ZPHD) - fgmepd - fdpd - fdpd2) * (fthetad / self%xxi))

  ! phytoplankton
  _SET_ODE_(self%id_ZPHN,(fprn * ZPHN) - fgmipn - fgmepn - fdpn - fdpn2 )
  _SET_ODE_(self%id_ZPHD,(fprd * ZPHD) - fgmepd - fdpd - fdpd2 )
  _SET_ODE_(self%id_ZPDS,(fprds * ZPDS) - fgmepds - fdpds - fsdiss - fdpds2 )

  ! zooplankton
  _SET_ODE_(self%id_ZZMI, fmigrow - fgmezmi - fdzmi - fdzmi2 )
  _SET_ODE_(self%id_ZZME, fmegrow - fdzme - fdzme2 )

  ! detritus
  _SET_ODE_(self%id_ZDET,fdpn + ((1._rk - self%xfdfrac1) * fdpd) + fdzmi + ((1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetan) * (finmi + finme))-fgmid-fgmed-fdd)
        
         !+ ffast2slown                                        ! seafloor fast->slow

  !dissolved inorganic nitrogen
   fn_cons = - (fprn * ZPHN) - (fprd * ZPHD)                         ! primary production
   fn_prod = + (self%xphi * (fgmipn + fgmid))                     &  ! messy feeding remin.
             + (self%xphi * (fgmepn + fgmepd + fgmezmi + fgmed))  &  ! messy feeding remin.
             + fmiexcr + fmeexcr + fdd                            &  ! excretion and remin.
            !+ freminn (fast-detritus contribution!!!)            &      
             + fdpn2 + fdpd2 + fdzmi2 + fdzme2                       ! metab. losses

   _SET_ODE_(self%id_ZDIN,fn_prod + fn_cons)

  ! dissolved silicic acid
   fs_cons = - (fprds * ZPDS)                                       ! opal production
   fs_prod = + fsdiss                                             &  ! opal dissolution
             + ((1.0 - self%xfdfrac1) * fdpds)                    &  ! mort. loss
             + ((1.0 - self%xfdfrac3) * fgmepds)                  &  ! egestion of grazed Si
             !+ freminsi (fast-detritus contribution!!!)          &  ! fast dissolution
             + fdpds2                                                ! metab. losses
   _SET_ODE_(self%id_ZSIL,fs_prod + fs_cons)

  ! dissolved iron
   _SET_ODE_(self%id_ZFER, (self%xrfn * (fn_prod-fn_cons)) - ffescav)
                                               !+ ffetop     &
                                               !+ ffebot     &
                                                 

  ! detrital carbon
   _SET_ODE_(self%id_ZDTC, (self%xthetapn * fdpn) + ((1._rk - self%xfdfrac1) * (self%xthetapd * fdpd)) + (self%xthetazmi * fdzmi) + ((1._rk - self%xfdfrac2) * (self%xthetazme * fdzme)) + ((1._rk - self%xbetac) * (ficmi + ficme))- fgmidc - fgmedc - fddc)
                 !+ ffast2slowc                             ! seafloor fast->slow

  ! dissolved inorganic carbon
   fc_cons = - (self%xthetapn * fprn * ZPHN) - (self%xthetapd * fprd * ZPHD)                      ! primary production
   fc_prod = + (self%xthetapn * self%xphi * fgmipn) + (self%xphi * fgmidc)                     &  ! messy feeding remin
             + (self%xthetapn * self%xphi * fgmepn) + (self%xthetapd * self%xphi * fgmepd)     &  ! messy feeding remin
             + (self%xthetazmi * self%xphi * fgmezmi) + (self%xphi * fgmedc)                   &  ! messy feeding remin
             + fmiresp + fmeresp + fddc &
             !+ freminc 
             + (self%xthetapn * fdpn2)                                                         &  ! resp., remin., losses
             + (self%xthetapd * fdpd2) + (self%xthetazmi * fdzmi2)                              &  ! losses
             + (self%xthetazme * fdzme2)                                                          ! losses
               
  ! fc_prod = fc_prod - ftempca + freminca         ! CaCO3
                 
   _SET_ODE_(self%id_ZDIC,fc_prod + fc_cons)

  ! alkalinity
  ! fa_prod =  2._rk * freminca                                                   ! CaCO3 dissolution
  ! fa_cons = -2._rk * ftempca                                                    ! CaCO3 production

  ! _SET_ODE_(self%id_ZALK, fa_prod + fa_cons)

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
               !- (self%xthetanit * freminn)                                               & ! fast detritus remin., N
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
               !- (self%xthetarem * freminc)                                               & ! fast detritus remin., C
               - (self%xthetarem * self%xthetapn * fdpn2)                                 & ! Pn  losses, C
               - (self%xthetarem * self%xthetapd * fdpd2)                                 & ! Pd  losses, C
               - (self%xthetarem * self%xthetazmi * fdzmi2)                               & ! Zmi losses, C
               - (self%xthetarem * self%xthetazme * fdzme2)                                 ! Zme losses, C

   fo2_cons = fo2_ncons + fo2_ccons

   if (zoxy .lt. self%xo2min) then                     ! deficient O2; production fluxes only
      _SET_ODE_(self%id_ZOXY, fo2_prod )
   else                                                ! sufficient O2; production + consumption fluxes
      _SET_ODE_(self%id_ZOXY, fo2_prod + fo2_cons )
   endif

   !+ air-sea fluxes...
  _SET_DIAGNOSTIC_(self%id_dPAR,par)

   _LOOP_END_

   end subroutine do

   subroutine fast_detritus(self,_ARGUMENTS_VERTICAL_)
   class(type_medusa_pelagic),intent(in) :: self

   _DECLARE_ARGUMENTS_VERTICAL_  

   real(rk) :: dz,fq0,fq1,fq2,fq3,fq4,fq5,fq6,fq7,fq8,fprotf
   real(rk) :: xmassc = 12.011_rk
   real(rk) :: xmassca = 100.086_rk
   real(rk) :: xmasssi = 60.084_rk
   real(rk) :: xprotca = 0.07_rk
   real(rk) :: xprotsi = 0.026_rk
   real(rk) :: xfastc = 188._rk
   real(rk) :: freminc,freminn,freminfe,freminsi,freminca
   real(rk) :: ffastc=0._rk,ffastn=0._rk,ffastca=0._rk,ffastsi=0._rk
   real(rk) :: ftempc
   

   _VERTICAL_LOOP_BEGIN_

   _GET_(self%id_dz,dz)
   _GET_(self%id_ftempc,ftempc)

   !Carbon
   fq0      = ffastc      !! how much organic C enters this box        (mol)
   fq1      = (fq0 * xmassc)     !! how much it weighs                        (mass)
   fq2      = (ffastca * xmassca)        !! how much CaCO3 enters this box            (mass)
   fq3      = (ffastsi * xmasssi)        !! how much opal enters this box            (mass)
   fq4      = (fq2 * xprotca) + (fq3 * xprotsi) !! total protected organic C                 (mass)

   !! this next term is calculated for C but used for N and Fe as well
   !! it needs to be protected in case ALL C is protected

   if (fq4.lt.fq1) then
     fprotf   = (fq4 / (fq1 + tiny(fq1)))      !! protected fraction of total organic C     (non-dim)
   else
     fprotf   = 1._rk                            !! all organic C is protected                (non-dim)
   endif

   fq5      = (1._rk - fprotf)                    !! unprotected fraction of total organic C   (non-dim)
   fq6      = (fq0 * fq5)                       !! how much organic C is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected C leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))            !! how much total C leaves this box          (mol)
   freminc  = (fq0 - fq8) / dz            !! C remineralisation in this box            (mol)
   ffastc = fq8
                           
   !Nitrogen
   fq0      = ffastn   !! how much organic N enters this box        (mol)
   fq5      = (1._rk - fprotf)   !! unprotected fraction of total organic N   (non-dim)
   fq6      = (fq0 * fq5)      !! how much organic N is unprotected         (mol)
   fq7      = (fq6 * exp(-(dz / xfastc)))     !! how much unprotected N leaves this box    (mol)
   fq8      = (fq7 + (fq0 * fprotf))            !! how much total N leaves this box          (mol)
   freminn  = (fq0 - fq8) / dz                !! N remineralisation in this box            (mol)
   ffastn = fq8 
   
    ffastc  = ffastc + ftempc

   _VERTICAL_LOOP_END_

   end subroutine fast_detritus

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

   class(type_medusa_pelagic),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: pt,ps,o2,o2_sato,o2_schmidt,kwo2,o2_sat,pp0,kw660,wnd,o2flux,o2sat
   real(rk) :: a0 = 2.00907_rk
   real(rk) :: a1 = 3.22014_rk
   real(rk) :: a2 = 4.05010_rk
   real(rk) :: a3 = 4.94457_rk
   real(rk) :: a4 = -2.56847E-1_rk
   real(rk) :: a5 = 3.88767_rk
   real(rk) :: b0 = -6.24523E-3_rk
   real(rk) :: b1 = -7.37614E-3_rk
   real(rk) :: b2 = -1.03410E-2_rk
   real(rk) :: b3 = -8.17083E-3_rk
   real(rk) :: c0 = -4.88682E-7_rk
   real(rk) :: tt,tk,ts,ts2,ts3,ts4,ts5
   real(rk) :: ans1, ans2
  ! Wanninkhof (2014) coefficients
   real(rk) :: as0 = 1920.4_rk
   real(rk) :: as1 = -135.6_rk
   real(rk) :: as2 = 5.2121_rk
   real(rk) :: as3 = -0.10939_rk
   real(rk) :: as4 = 0.00093777_rk
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk) :: a(7)
   real(rk) :: b(7)
   real(rk) :: tmp_k
   data a(1) / 0.166_rk /  ! Liss & Merlivat (1986)    [approximated]
   data a(2) / 0.3_rk /    ! Wanninkhof (1992)         [sans enhancement]
   data a(3) / 0.23_rk /   ! Nightingale et al. (2000) [good]
   data a(4) / 0.23_rk /   ! Nightingale et al. (2000) [better]
   data a(5) / 0.222_rk /  ! Nightingale et al. (2000) [best]
   data a(6) / 0.337_rk /  ! OCMIP-2                   [sans variability]
   data a(7) / 0.251_rk /  ! Wanninkhof (2014)         [assumes 6h avg winds]
   data b(1) / 0.133_rk /
   data b(2) / 0.0_rk /
   data b(3) / 0.0_rk /
   data b(4) / 0.1_rk /
   data b(5) / 0.333_rk /
   data b(6) / 0.0_rk /
   data b(7) / 0.0_rk /


   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,pt)
   _GET_(self%id_salt,ps)
   _GET_(self%id_ZOXY,o2)
   _GET_HORIZONTAL_(self%id_apress,pp0)
   _GET_HORIZONTAL_(self%id_wnd,wnd)

      o2 = o2/1000._rk

! Calculate gas transfer velocity (cm/h)

      tmp_k = (a(7) * wnd**2) + (b(7) * wnd)

! Convert tmp_k from cm/h to m/s
      kw660 = tmp_k / (3600._rk * 100._rk)

 !note: air-sea fluxes must be corrected for sea ice

      tt   = 298.15_rk - pt
      tk   = 273.15_rk + pt
      ts   = log(tt / tk)
      ts2  = ts**2_rk
      ts3  = ts**3_rk
      ts4  = ts**4_rk
      ts5  = ts**5_rk

      ans1 = a0 + a1*ts + a2*ts2 + a3*ts3 + a4*ts4 + a5*ts5  &
             + ps*(b0 + b1*ts + b2*ts2 + b3*ts3)             &
             + c0*(ps*ps)

      ans2 = exp(ans1)

!  Convert from ml/l to mol/m3
   o2_sato = (ans2 / 22391.6_rk) * 1000.0_rk

   o2_schmidt = as0 + pt*(as1 + pt*(as2 + pt*(as3 + pt*as4)))
   kwo2 = kw660 * (660._rk / o2_schmidt)**0.5_rk
   o2sat = o2_sato * pp0 / 101325._rk
   !print*,pp0
   
   o2flux = kwo2 * (o2sat - o2)
   o2flux = o2flux *1000.

   _SET_SURFACE_EXCHANGE_(self%id_ZOXY, o2flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fair, o2flux)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

  end module medusa_pelagic
