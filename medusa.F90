#include "fabm_driver.h"

!
!********************************************    
!                FABM-MEDUSA                       
!********************************************                                                                            

module fabm_medusa

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZCHN,id_ZCHD,id_ZPHN,id_ZPHD,id_ZPDS,id_ZDIN,id_ZFER,id_ZSIL,id_ZDET,id_ZDTC,id_ZZMI,id_ZZME
      type (type_dependency_id)            :: id_temp,id_par
  !    type (type_diagnostic_variable_id)   ::

      ! Parameters
      logical :: jliebig
      real(rk) :: xxi,xaln,xald,xnln,xfln,xnld,xsld,xfld,xvpn,xvpd,xsin0,xnsi0,xuif,xthetam,xthetamd
      real(rk) :: xkmi,xpmipn,xpmid,xkme,xpmepn,xpmepd,xpmezmi,zpmed,xgmi,xgme,xthetad,xphi,xthetapn,xthetazme,xthetazmi,xthetapd
      real(rk) :: xmetapn,xmetapd,xmetazmi,xmetazme,xmpn,xmpd,xmzmi,xmzme,xkphn,xkphd,xkzmi,xkzme
      real(rk) :: xmd,xmdc,xsdiss
      real(rk) :: xk_FeL,xLgT,xk_sc_Fe
      real(rk) :: xfdfrac1,xfdfrac2,xfdfrac3,xrfn
   contains
      procedure :: initialize
      procedure :: do

  end type

contains

   subroutine initialize(self,configunit)
! SELF%DT???
   class(type_medusa),intent(inout),target :: self
   integer,               intent(in)           :: configunit   
   !Register parameters
   call self%get_parameter(self%xxi, 'xxi', 'molN gC-1','C : N conversion factor', default=0.01257_rk)
   call self%get_parameter(self%xaln, 'xaln', 'gC(g chl)-1 (W m-2)-1 d-1','chl-specific initial slope of P-I curve (non-diatoms)', default=15.0_rk)
   call self%get_parameter(self%xald, 'xald', 'gC(g chl)-1 (W m-2)-1 d-1','chl-specific initial slope of P-I curve (diatoms)', default=11.25_rk)
   call self%get_parameter(self%jliebig, 'jliebig', 'Liebig''s minimum law for nutrient limitation', default=.false.)
   call self%get_parameter(self%xnln, 'xnln', 'mmolN m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.5_rk)
   call self%get_parameter(self%xfln, 'xfln', 'mmolFe m-3','Fe nutrient uptake half-saturation constant (non-diatoms)', default=0.33_rk)
   call self%get_parameter(self%xnld, 'xnld', 'mmolN m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.75_rk)
   call self%get_parameter(self%xsld, 'xsld', 'mmolSi m-3','Si nutrient uptake half-saturation constant (diatoms)', default=3.0_rk)
   call self%get_parameter(self%xfld, 'xfld', 'mmolFe m-3','Fe nutrient uptake half-saturation constant (diatoms)', default=0.67_rk)
   call self%get_parameter(self%xvpn, 'xvpn', 'd-1','Maximum phytoplankton growth rate (non-diatoms)', default=0.53_rk)
   call self%get_parameter(self%xvpd, 'xvpd', 'd-1','Maximum phytoplankton growth rate (diatoms)', default=0.50_rk)
   call self%get_parameter(self%xsin0, 'xsin0', 'molSi molN-1','minimum diatom Si : N ratio', default=0.2_rk)
   call self%get_parameter(self%xnsi0, 'xnsi0', 'molN molSi-1','minimum diatom N : Si ratio', default=0.2_rk)
   call self%get_parameter(self%xuif, 'xuif', '-','hypothetical growth ratio at Inf Si : N ratio', default=1.5_rk)
   call self%get_parameter(self%xthetam, 'xthetam', 'g chl gC-1','maximum Chl : C ratio (non-diatoms)', default=0.05_rk)
   call self%get_parameter(self%xthetamd, 'xthetamd', 'g chl gC-1','maximum Chl : C ratio (diatoms)', default=0.05_rk)
   call self%get_parameter(self%xkmi, 'xkmi', 'mmolN m-3','microzooplankton grazing half-saturation constant', default=0.8_rk)
   call self%get_parameter(self%xpmipn, 'xpmipn', '-','microzooplankton grazing preference on non-diatoms', default=0.75_rk)
   call self%get_parameter(self%xpmid, 'xpmid', '-','microzooplankton grazing preference on detritus', default=0.25_rk)
   call self%get_parameter(self%xkme, 'xkme', 'mmolN m-3','mesozooplankton grazing half-saturation constant', default=0.3_rk)
   call self%get_parameter(self%xpmepn, 'xpmepn', '-','mesozooplankton grazing preference on non-diatoms', default=0.15_rk)
   call self%get_parameter(self%xpmepd, 'xpmepd', '-','mesozooplankton grazing preference on diatoms', default=0.35_rk)
   call self%get_parameter(self%xpmezmi, 'xpmezmi', '-','mesozooplankton grazing preference on microzooplankton', default=0.35_rk)
   call self%get_parameter(self%xpmed, 'xpmed', '-','mesozooplankton grazing preference on detritus', default=0.15_rk)
   call self%get_parameter(self%xgmi, 'xgmi', 'd-1','maximum microzooplankton grazing rate', default=2.0_rk)
   call self%get_parameter(self%xgme, 'xgme', 'd-1','maximum mesozooplankton grazing rate', default=0.5_rk)
   call self%get_parameter(self%xthetad, 'xthetad', 'molC molN-1','detritus C : N ratio', default=6.625_rk)
   call self%get_parameter(self%xphi, 'xphi', '-','zooplankton grazing inefficiency', default=0.20_rk)
   call self%get_parameter(self%xthetapn, 'xthetapn', 'molC molN-1','phytoplankton C:N ratio (non-diatoms)', default=6.625_rk)
   call self%get_parameter(self%xthetapd, 'xthetapd', 'molC molN-1','phytoplankton C:N ratio (diatoms)', default=6.625_rk)
   call self%get_parameter(self%xbetan, 'xbetan', '-','zooplankton N assimilation efficiency', default=0.77_rk)
   call self%get_parameter(self%xthetazmi, 'xthetazmi', 'molC molN-1','microzooplankton C:N ratio', default=5.625_rk)
   call self%get_parameter(self%xbetac, 'xbetac', '-','zooplankton C assimilation efficiency', default=0.64_rk)
   call self%get_parameter(self%xkc, 'xkc', '-','zooplankton net C growth efficiency', default=0.8_rk)
   call self%get_parameter(self%xthetazme, 'xthetazme', 'molC molN-1','mesozooplankton C:N ratio', default=5.625_rk)
   call self%get_parameter(self%xmetapn, 'xmetapn', 'd-1','phytoplankton loss rate (non-diatoms)', default=0.02_rk)
   call self%get_parameter(self%xmetapd, 'xmetapd', 'd-1','phytoplankton loss rate (diatoms)', default=0.02_rk)
   call self%get_parameter(self%xmetazmi, 'xmetazmi', 'd-1','microzooplankton loss rate', default=0.02_rk)
   call self%get_parameter(self%xmetazme, 'xmetazme', 'd-1','mesozooplankton loss rate', default=0.02_rk)
   call self%get_parameter(self%xmpn,'xmpn','d-1','phytoplankton maximum loss rate (non-diatoms)', default=0.1_rk)
   call self%get_parameter(self%xmpd,'xmpd','d-1','phytoplankton maximum loss rate (diatoms)', default=0.1_rk)
   call self%get_parameter(self%xmzmi,'xmzmi','d-1','microzooplankton maximum loss rate', default=0.1_rk)
   call self%get_parameter(self%xmzme,'xmzme','d-1','mesozooplankton maximum loss rate', default=0.2_rk)
   call self%get_parameter(self%xkphn,'xkphn','mmolN m-3','phytoplankton los half-saturation constant (non-diatoms)', default=0.5_rk)
   call self%get_parameter(self%xkphd,'xkphd','mmolN m-3','phytoplankton los half-saturation constant (diatoms)', default=0.5_rk)
   call self%get_parameter(self%xkzmi,'xkzmi','mmolN m-3','microzooplankton loss half-saturation constant', default=0.5_rk)
   call self%get_parameter(self%xkzme,'xkzme','mmolN m-3','mesozooplankton loss half-saturation constant', default=0.75_rk)
   call self%get_parameter(self%xmd,'xmd','d-1','detrital N remineralisation rate', default=0.0158_rk)
   call self%get_parameter(self%xmdc,'xmdc','d-1','detrital C remineralisation rate', default=0.0127_rk)
   call self%get_parameter(self%xsdiss,'xsdiss','d-1','diatom frustule dissolution rate', default=0.006_rk)
   call self%get_parameter(self%xk_FeL,'xk_FeL','-','dissociation constant for (Fe+ligand)',default=100.0_rk)
   call self%get_parameter(self%xLgT,'xLgT','umol m-3','total ligand concentration',default=1.0_rk)
   call self%get_parameter(self%xk_sc_Fe,'xk_sc_Fe','d-1','scavenging rate of "free" Fe',default=0.001_rk)
   call self%get_parameter(self%xfdfrac1,'xfdfrac1','-','fast detritus fraction of diatom losses',default=0.33_rk)
   call self%get_parameter(self%xfdfrac2,'xfdfrac2','-','fast detritus fraction of mesozooplankton losses',default=1._rk)
   call self%get_parameter(self%xfdfrac3,'xfdfrac3','-','fast detritus fraction of mesozooplankton grazing',default=0.8_rk)
   call self%get_parameter(self%xrfn,'xrfn','umolFe molN-1 m','phytoplankton Fe : N uptake ratio',default=30._rk)
   ! Register state variables
   call self%register_state_variable(self%id_ZCHN,'ZCHN','mg chl/m**3', 'chlorophyll in non-diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZCHD','mg chl/m**3', 'chlorophyll in diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZPHN','mmolN/m**3', 'non-diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZPHD','mmolN/m**3', 'diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZPDS','mmolSi/m**3', 'diatom phytoplankton (silicon)', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIN,'ZDIN','mmolN/m**3', 'nitrogen nutrient', minimum=0.0_rk) !no river dilution?
   call self%register_state_variable(self%id_ZFER,'ZFER','mmolFe/m**3', 'iron nutrient', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSIL,'ZSIL','mmolSi/m**3', 'silicic acid', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDET,'ZDET','mmolN/m**3', 'slow-sinking detritus (N)', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDTC,'ZDTC','mmolC/m**3', 'slow-sinking detritus (C)', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZZMI,'ZZMI','mmolN/m**3', 'microzooplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZZME,'ZZME','mmolN/m**3', 'mesozooplankton', minimum=0.0_rk)
   ! Register diagnostic variables
     
   ! Register environmental dependencies

   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   class(type_medusa), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

! !LOCAL VARIABLES:

    real(rk) :: ZCHN,ZCHD,ZPHN,ZPHD,ZPDS,ZDIN,ZFER,ZSIL,ZDET,ZDTC,ZZMI,ZZME,loc_T,par
    real(rk) :: fthetan,fthetad,faln,fald !scaled chl/biomass ratio
    real(rk) :: fnln,ffln ! non-diatom Qn/Qf terms
    real(rk) :: fnld,fsld,ffld ! diatom Qn/Qs/Qf terms
    real(rk) :: fun_T,xvpnT,xvpdT,fchn1,fchn,fjln,fchd1,fchd,fjld
    real(rk) :: fsin,fnsi,fsin1,fnsi1,fnsi2,fprn,fprd,fsld2,frn,frd,fprds
    real(rk) :: fpnlim,fpdlim !nutrient limitation of primary production
    real(rk) :: fmi1,fmi,fgmipn,fgmid,fgmidc,finmi,ficmi,fstarmi,fmith,fmigrow,fmiexcr,fmiresp
    real(rk) :: fme1,fme,fgmepn,fgmepd,fgmepds,fgmezmi,fgmed,fgmedc,finme,ficme,fstarme,fmeth,fmegrow,fmeexcr,fmeresp
    real(rk) :: fdpn2,fdpd2,fdpds2,fdzmi2,fdzme2,fdpn,fdpd,fdzmi,fdzme
    real(rk) :: fdd,fddc,fsdiss
    real(rk) :: xFeT,xb_coef_tmp,xb2M4ac,xLgF,xFel,xFeF,xFree,ffescav
    real(rk) :: fslown,fregen,fregensi,fregenc,ftempn,ftempsi,ftempfe,ftempc

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
    _GET_(self%id_temp,loc_T)
    _GET_(self%id_par,par) !check PAR // what about self-shading?

   !PHYTOPLANKTON GROWTH
   !Chlorophyll
   fthetan = (ZCHN * self%xxi) / ZPHN
   faln = self%xaln * fthetan
   fthetad = (ZCHD * self%xxi) / ZPHD
   fald = self%xald * fthetad

  !Temperature limitation
   fun_T = 1.066_rk**loc_T

   xvpnT = self%xvpn * fun_T
   xvpdT = self%xvpd * fun_T

  !Phytoplankton light limitation
   fchn1 = (xvpnT * xvpnT) + (faln * faln * par * par)
   fchn = xvpnT / sqrt(fchn1)
   fjln = fchn * faln * par !non-diatom J term

   fchd1 = (xvpdT * xvpdT) + (fald * fald * par * par)
   fchd = xvpdT / sqrt(fchd1)
   fjld = fchd * fald * par !diatom J term

   ! Phytoplankton nutrient limitation
   !! Non-diatoms (N, Fe)
   fnln = ZDIN / (ZDIN + self%xnln) !non-diatom Qn term
   ffln = ZFER / (ZFER + self%xfln) !non-diatom Qf term
   !! Diatoms (N, Si, Fe)
   fnld = ZDIN / (ZDIN + self%xnld) !diatom Qn term
   fsld = ZSIL / (ZSIL + self%xsld) !diatom Qs term
   ffld = ZFER / (ZFER + self%xfld) !diatom Qf term

   ! Primary production (non-diatoms)
   if (self%njliebig=.false.) then
     fpnlim = fnln * ffln
   else
     fpnlim = min(fnln,ffln)
   end if

   fprn = fjln * fpnlim

   !Primary production (diatoms)
   if (self%nliebig=.false.) then
     fpdlim = fnld * fsld
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
      fprd = self%xuif * ((fsin - self%xsin0) / fsin) * (fjld * fpdlim)
      fsld2 = self%xuif * ((fsin - self%xsin0) / fsin)
   elseif (fsin .ge. fsin1) then
      fprd = fjld * fpdlim
      fsld2 = 1.0_rk
   end if
  !Silicon

   if (fsin.lt.fnsi1) then
     fprds = (fjld * fsld)
   elseif (fsin.lt.fnsi2) then
     fprds = self%xuif * ((fnsi - self%xnsi0) / fnsi) * (fjld * fsld)
   else
     fprds = 0._rk
   end if

  !Chlorophyll production
  frn = (self%xthetam * fchn * fnln * ffln) / fthetan
  frd = (self%xthetamd * fchd * fnld * ffld * fsld2) / fthetad  

  !ZOOPLANKTON GRAZING
  !Microzooplankton
  fmi1 = (self%xkmi * self%xkmi) + (self%xpmipn * ZPHN * ZPHN) + (self%xpmid * ZDET * ZDET)
  fmi = self%xgmi * ZZMI / fmi1
  fgmipn = fmi * self%xpmipn * ZPHN * ZPHN !grazing on non-diatoms
  fgmid = fmi * self%xpmid * ZDET * ZDET   !grazing on detrital nitrogen
  fgmidc = self%xthetad * fgmid !grazing on detrital carbon: non-ROAM formulation (see switch in original code to be implemented)
  finmi = (1.0_rk - self%xphi) * (fgmipn + fgmid)
  ficmi = (1.0_rk - self%xphi) * ((self%xthetapn * fgmipn) + fgmidc)
  fstarmi = (self%xbetan * self%xthetazmi) / (self%xbetac * self%xkc) !the ideal food C : N ratio for microzooplankton
  fmith = ficmi / finmi
  if (fmith .ge. fstarmi) then
     fmigrowth = self%xbetan * finmi
     fmiexcr = 0.0_rk
  else
     fmigrow = (self%xbetac * self%xkc * ficmi) / self%xthetazmi
     fmiexcr = ficmi * ((self%xbetan / fmith) - ((self%xbetac * self%xkc) / self%xthetazmi))
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
  fgmedc = self%xthetad * fgmed !grazing on detrital carbon: non-ROAM formulation (see switch in original code to be implemented)
  finme = (1.0_rk - self%xphi) * (fgmepn + fgmepd + fgmezmi + fgmed)
  ficme = (1.0_rk - self%xphi) * ((self%xthetapn * fgmepn) + (self%xthetapd * fgmepd) + (self%xthetazmi * fgmezmi) + fgmedc) 
  fstarme = (self%xbetan * self%xthetazme) / (self%xbetac * self%xkc)
  fmeth = ficme / finme
  if (fmeth .ge. fstarme) then
     fmegrow = self%xbetan * finme
     fmeexcr = 0.0_rk
  else
     fmegrow = (self%xbetac * self%xkc * ficme) / self%xthetazme
     fmeexcr = ficme * ((self%xbetan / fmeth) - ((self%xbetac * self%xkc) / self%xthetazme))
  end if
  fmeresp = (self%xbetac * ficme) - (self%xthetazme * fmegrow)

  !Plankton metabolic losses
  !Linear loss processes assumed to be metabolic in origin
  fdpn2 = self%xmetapn * ZPHN
  fdpd2 = self%xmetapd * ZPHD
  fdpds2 = self%xmetapd * ZPDS
  fdzmi2 = self%xmetazmi * ZZMI
  fdzme2 = self%xmetazme * ZZME

  !Plankton mortality losses !NB: currently hyperbolic mortality term only (18/08/2017)
  fdpn = self%xmpn * ZPHN * (ZPHN / (self%xkphn + ZPHN)) !non-diatom phytoplankton
  fdpd = self%xmpd * ZPHD * (ZPHD / (self%xkphd + ZPHD)) !diatom phytoplankton
  fdzmi = self%xmzmi * ZZMI * (ZZMI / (self%xkzmi + ZZMI)) !microzooplankton
  fdzme = self%xmzme * ZZME * (ZZME / (self%xkzme + ZZME)) !mesozooplankton

  !Detritus remineralisation (temperature-dependent) !another option is Q10-based
  fdd = self%xmd * fun_T * ZDET
  fddc = self%xmdc * fun_T * ZDTC

  !Original contains accelerated detrital remineralisation in the bottom box (how do I let model know it is a bottom box?)

  !Diatom frustule dissolution
  fsdiss = self%xsdiss * ZPDS

  !IRON CHEMISTRY AND FRACTIONATION
  xFeT = ZFER * 1.e3_rk !total iron concentration (mmolFe/m3 -> umolFe/m3)
  xb_coef_tmp = self%xk_FeL * (self%xLgT - xFeT) - 1.0_rk !this is F1 in Yool et al (2013)
  xb2M4ac = max(((xb_coef_tmp * xb_coef_tmp) + (4.0_rk * self%xk_FeL * self%xLgT)), 0._rk)
  xLgF = 0.5_rk * (xb_coef_tmp + (xb2M4ac**0.5_rk)) / self%xk_FeL ! "free" ligand concentration
  xFeL = xLgT - xLgF ! ligand-bound iron concentration
  xFeF = (xFeT - xFeL) * 1.e-3_rk! "free" iron concentration (and convert to mmolFe/m3)
  xFree = xFeF / ZFER
  ! Scavenging of iron (option 1 from original code)
  ffescav = self%xk_sc_Fe * xFeF
  !!Further (optional) implicit "scavenging" by Mick Follows, caps concentration of total Fe, and not included here yet...
  !!Aeolian iron deposition and seafloor iron addition - should be dealt with through model inputs.
  ! Stop scavenging for inron at depths greater than 1000 m - how to do this?
 
  !Slow detritus creation
  fslown  = fdpn + fdzmi + ((1._rk - self%xfdfrac1) * fdpd) + ((1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetan) * (finmi + finme))
  fslowc  = (self%xthetapn * fdpn) + (self%xthetazmi * fdzmi) + (self%xthetapd * (1._rk - self%xfdfrac1) * fdpd) + (self%xthetazme * (1._rk - self%xfdfrac2) * fdzme) + ((1._rk - self%xbetac) * (ficmi + ficme))

  !Nutrient regeneration
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

  !_SET_ODE_(self%id_..,)

   _LOOP_END_

   end subroutine do


  end module fabm_medusa

