#include "fabm_driver.h"
module medusa_carbonate

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_medusa_carbonate
      type (type_state_variable_id)     :: id_O3c,id_TA,id_bioalk
      type (type_dependency_id)         :: id_ETW, id_X1X, id_dens, id_pres
      type (type_horizontal_dependency_id) :: id_wnd,id_PCO2A
      type (type_diagnostic_variable_id) :: id_ph,id_pco2,id_CarbA,id_BiCarb,id_Carb,id_TA_diag
      type (type_diagnostic_variable_id) :: id_Om_cal,id_Om_arg

   contains
     procedure :: initialize
     procedure :: do
     procedure :: do_surface  
   end type

   public :: CO2_dynamics,co2dyn,CaCO3_Saturation

contains
    
    subroutine initialize(self,configunit)

      class (type_ersem_carbonate), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

     call self%register_state_variable(self%id_O3c,'c','mmol C/m^3','total dissolved inorganic carbon', 2200._rk,minimum=0._rk)
     call self%register_state_variable(self%id_TA,'TA','umol/kg','total alkalinity',2300._rk,minimum=1.e-4_rk)

     call self%register_diagnostic_variable(self%id_ph,    'pH',    '-',      'pH',standard_variable=standard_variables%ph_reported_on_total_scale)
     call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  '1e-6',    'partial pressure of CO2')
     call self%register_diagnostic_variable(self%id_CarbA, 'CarbA', 'mmol/m^3','carbonic acid concentration')
     call self%register_diagnostic_variable(self%id_BiCarb,'BiCarb','mmol/m^3','bicarbonate concentration')
     call self%register_diagnostic_variable(self%id_Carb,  'Carb',  'mmol/m^3','carbonate concentration')
     call self%register_diagnostic_variable(self%id_Om_cal,'Om_cal','-','calcite saturation')
     call self%register_diagnostic_variable(self%id_Om_arg,'Om_arg','-','aragonite saturation')

     call self%register_dependency(self%id_ETW, standard_variables%temperature)
     call self%register_dependency(self%id_X1X, standard_variables%practical_salinity)
     call self%register_dependency(self%id_dens,standard_variables%density)
     call self%register_dependency(self%id_pres,standard_variables%pressure)
     call self%register_dependency(self%id_wnd,  standard_variables%wind_speed)
     call self%register_dependency(self%id_PCO2A,standard_variables%mole_fraction_of_carbon_dioxide_in_air)

    end subroutine

    subroutine do(self,_ARGUMENTS_DO_)
    class (type_ersem_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_


      _LOOP_BEGIN_
         _GET_(self%id_O3C,O3C)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_(self%id_dens,density)
         _GET_(self%id_pres,pres)

   call CO2_dynamics( Temp, Sal, Depth, DIC, ALK, pCO2a, &                  ! inputs 
      pco2w, ph, h2co3, hco3, co3, henry, om_cal, om_arg, TDIC, TALK, &     ! outputs 
      dcf, iters) 

      _LOOP_END_
   end subroutine

   subroutine CO2_dynamics( T, S, Z, DIC, TALK, pco2a, pco2w, ph, h2co3, bicarb, &
                            carb, henry, om_cal, om_arg, TCO2, TA, dcf, iters ) 

      real(rk),intent( in )    :: T,S,Z
      real(rk),intent( in )    :: DIC        ! total dissolved inorganic carbon (mmol.m-3) 
      real(rk),intent( in )    :: TALK       ! total alkalinity (meq.m-3) 
      real(rk),intent( in )    :: pco2a      ! atmospheric pCO2 
      real(rk),intent( inout ) :: pco2w      ! seawater pCO2 
      real(rk),intent( inout ) :: ph         ! seawater pH 
      real(rk),intent( inout ) :: h2co3      ! seawater carbonic acid (H2CO3) 
      real(rk),intent( inout ) :: bicarb     ! seawater bicarbonate ion (HCO3) 
      real(rk),intent( inout ) :: carb       ! seawater carbonate ion (CO3) 
      real(rk),intent( inout ) :: henry      ! Henry constant 
      real(rk),intent( inout ) :: om_cal     ! Omega calcite 
      real(rk),intent( inout ) :: om_arg     ! Omega aragonite 
      real(rk),intent( inout ) :: TCO2       ! total dissolved inorganic carbon (mol.kg-1) 
      real(rk),intent( inout ) :: TA         ! total alkalinity (eq.kg-1) 
      real(rk),intent( inout ) :: dcf        ! density conversion factor 
      integer,intent( inout ) :: iters      ! # iterations to convergence 
      real(rk) :: a, b, c 
      real(rk) :: ca, bc, cb 
      real(rk) :: pco2water, fairco2 

          a   =  8.24493e-1_rk - 4.0899e-3_rk*T +  7.6438e-5_rk*T**2 - 8.2467e-7_rk*T**3 + 5.3875e-9_rk*T**4 
          b   = -5.72466e-3_rk + 1.0227e-4_rk*T - 1.6546e-6_rk*T**2  
          c   = 4.8314e-4_rk
          dcf = (999.842594_rk + 6.793952e-2_rk*T- 9.095290e-3_rk*T**2 + 1.001685e-4_rk*T**3 & 
                - 1.120083e-6_rk*T**4 + 6.536332e-9_rk*T**5+a*S+b*S**1.5_rk+c*S**2)/1.0e3_rk 

          TA    = TALK / (1.0e6_rk*dcf) 
          TCO2  = DIC  / (1.0e6_rk*dcf) 

! Call the parent routine for the carbonate system 

      call CO2DYN ( TCO2, TA, T, S, pco2a, &     ! inputs 
          pco2water, pH, HENRY, ca, bc, cb, iters )  ! outputs 

          pco2w  = pco2water * (1.0e6_rk)  ! partial pressure of co2 in water  
          TA     = TA * (1.0e6_rk)         ! total alkalinity (umol/kg) 
          h2co3  = ca * (1.0e6_rk*dcf)     ! carbonic acid concentration (mmol/m3) 
          bicarb = bc * (1.0e6_rk*dcf)     ! bicarbonate ion concentration (mmol/m3) 
          carb   = cb * (1.0e6_rk*dcf)     ! carbonate ion concentration (mmol/m3) 
          TCO2   = TCO2 * (1.0e6_rk)       ! total C or DIC in units of umol/kg 

! Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states

      call CaCO3_Saturation ( T, S, Z, cb, &  ! inputs
          om_cal, om_arg )                    ! outputs

         _SET_DIAGNOSTIC_(self%id_Om_cal,Om_cal)
         _SET_DIAGNOSTIC_(self%id_Om_arg,Om_arg)

   end subroutine CO2_dynamics
!----------------------------------------------------------------------
!----------------------------------------------------------------------

  subroutine CO2SET(P,T,S,AKVAL,NKVAL,IC) 

      INTEGER             :: MAXK, MAXCON, NKVAL, ICON, IC, IK 
      PARAMETER(MAXK=4,MAXCON=3) 
      REAL(rk), DIMENSION(MAXK)        :: A, B, C 
      REAL(rk), DIMENSION(MAXK,MAXCON) :: A0, A1, A2 
      REAL(rk), DIMENSION(MAXK,MAXCON) :: B0, B1, B2 
      REAL(rk), DIMENSION(NKVAL)       :: AKVAL 
      REAL(rk)            :: P,T,S,VAL,TK, delta, kappa, Rgas 
      REAL(rk)            :: dlogTK, S2, sqrtS, S15, k1, k2, kb 

      DATA Rgas/83.131_rk/ 
      DATA A/-167.8108_rk, 290.9097_rk, 207.6548_rk, 148.0248_rk/ 
      DATA B/9345.17_rk, -14554.21_rk, -11843.79_rk, -8966.9_rk/ 
      DATA C/23.3585_rk, -45.0575_rk, -33.6485_rk, -24.4344_rk/ 
      DATA (A0(1,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (A0(2,ICON),ICON=1,MAXCON) /0.0221_rk, 0.5709_rk, -45.8076_rk/ 
      DATA (A0(3,ICON),ICON=1,MAXCON) /0.9805_rk, 1.4853_rk, -39.5492_rk/ 
      DATA (A0(4,ICON),ICON=1,MAXCON) /0.0473_rk, 0.5998_rk, 0.5998_rk/ 
      DATA (A1(1,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (A1(2,ICON),ICON=1,MAXCON) /34.02_rk, -84.25_rk, 1935.07_rk/ 
      DATA (A1(3,ICON),ICON=1,MAXCON) /-92.65_rk, -192.69_rk, 1590.14_rk/ 
      DATA (A1(4,ICON),ICON=1,MAXCON) /49.10_rk, -75.25_rk, -75.25_rk/ 
      DATA (A2(1,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (A2(2,ICON),ICON=1,MAXCON) /2*0.0_rk,6.9513_rk/ 
      DATA (A2(3,ICON),ICON=1,MAXCON) /2*0.0_rk,6.1523_rk/ 
      DATA (A2(4,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B0(1,ICON),ICON=1,MAXCON) /3*0.023517_rk/ 
      DATA (B0(2,ICON),ICON=1,MAXCON) /0.0_rk,-0.01632_rk,-0.01566_rk/ 
      DATA (B0(3,ICON),ICON=1,MAXCON) /-0.03294_rk,-0.05058_rk,-0.04997_rk/ 
      DATA (B0(4,ICON),ICON=1,MAXCON) /0.0_rk, -0.01767_rk, -0.01767_rk/ 
      DATA (B1(1,ICON),ICON=1,MAXCON) /3*-2.3656e-4_rk/ 
      DATA (B1(2,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B1(3,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B1(4,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B2(1,ICON),ICON=1,MAXCON) /3*4.7036e-7_rk/ 
      DATA (B2(2,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B2(3,ICON),ICON=1,MAXCON) /3*0.0_rk/ 
      DATA (B2(4,ICON),ICON=1,MAXCON) /3*0.0_rk/ 

      TK=T+273.15_rk 

      DO 100 IK=1,NKVAL 
        VAL=A(IK) + B(IK)/TK + C(IK)*LOG(TK) 
        VAL=VAL + (A0(IK,IC) + A1(IK,IC)/TK + A2(IK,IC)*LOG(TK))*SQRT(S) 
        VAL=VAL + (B0(IK,IC) + B1(IK,IC)*TK + B2(IK,IC)*TK*TK)*S 
        AKVAL(IK)=EXP(VAL)
100    CONTINUE 

      IF (IC .EQ. 3) THEN 

        dlogTK = log(TK)
        S2 = S*S
        sqrtS = sqrt(S)
        S15 = S**1.5_rk

        k1=10.**(-1.*(3670.7/TK - 62.008 + 9.7944*dlogTK - &
     &          0.0118 * S + 0.000116*S2))

        k2=10.**(-1.*(1394.7/TK + 4.777 - &
     &          0.0184*S + 0.000118*S2))

        delta=-25.5+0.1271*T
        kappa=(-3.08+0.0877*T)/1000.0
        k1=k1*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))

        delta=-15.82-0.0219*T
        kappa=(1.13-0.1475*T)/1000.0
        k2=k2*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))

        kb=exp((-8966.90 - 2890.53*sqrtS - 77.942*S + &
     &          1.728*S15 - 0.0996*S2)/TK + &
     &          (148.0248 + 137.1942*sqrtS + 1.62142*S) + &
     &          (-24.4344 - 25.085*sqrtS - 0.2474*S) * &
     &          dlogTK + 0.053105*sqrtS*TK)

        delta=-29.48+0.1622*T-0.002608*T**2.0 
        kappa=-2.84/1000.0 
        kb=kb*exp((-delta+0.5*kappa*P)*P/(Rgas*TK)) 

        delta=-25.60+0.2324*T-0.0036246*T**2 
        kappa=(-5.13+0.0794*T)/1000.0 

      AKVAL(1)=AKVAL(1)*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))
      AKVAL(2) = k1
      AKVAL(3) = k2
      AKVAL(4) = kb 

      end if
   end subroutine
!----------------------------------------------------------------------
!----------------------------------------------------------------------
   subroutine CaCO3_Saturation (Tc, S, D, CO3, Om_cal, Om_arg)
         real(rk),intent(in)  :: Tc, S, D, CO3
         real(rk),intent(out) :: Om_cal, Om_arg

        real(rk) Tk, Ca
        real(rk) logKspc, Kspc
        real(rk) logKspa, Kspa
        real(rk) tmp1, tmp2, tmp3
        real(rk) dV, dK, P, R

        real(rk),parameter :: Kelvin = 273.15_rk

        Tk = Tc + Kelvin
        Ca = 0.01028_rk    ! Currently oceanic mean value at S=25, needs refining)
        Ca = 0.010279_rk * (S / 35.0_rk)  ! Ca varies with salinity (cf. Feeley et al., 2004; Yool et al., 2010)
        R = 83.131_rk      !(cm3.bar.mol-1.K-1)
        P = D/10._rk    !pressure in bars

! calculate K for calcite
        tmp1 = -171.9065_rk - (0.077993_rk*Tk) + (2839.319_rk/Tk) + 71.595_rk*log10(Tk) 
        tmp2 = + (-0.77712_rk + (0.0028426_rk*Tk) + (178.34_rk/Tk))*SQRT(S)
        tmp3 = - (0.07711_rk*S) + (0.0041249_rk*(S**1.5_rk))
        logKspc = tmp1 + tmp2 + tmp3
        Kspc = 10._rk**logKspc

! correction for pressure for calcite
        dV = -48.76_rk + 0.5304_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspc = Kspc*exp(tmp1)
        logKspc = log10(Kspc)

        tmp1 = -171.945_rk - 0.077993_rk*Tk + 2903.293_rk / Tk + 71.595_rk* log10(Tk)
        tmp2 = + (-0.068393_rk + 0.0017276_rk*Tk + 88.135_rk/Tk)*SQRT(S)
        tmp3 = - 0.10018_rk*S + 0.0059415_rk*S**1.5_rk
        logKspa = tmp1 + tmp2 + tmp3
        Kspa = 10._rk**logKspa

! correction for pressure for aragonite
        dV = -46._rk + 0.530_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspa = Kspa*exp(tmp1)
        logKspa = log10(Kspa)

! calculate saturation states
        Om_cal = (CO3 * Ca) / Kspc
        Om_arg = (CO3 * Ca) / Kspa

      end subroutine CaCO3_SATURATION

end module
