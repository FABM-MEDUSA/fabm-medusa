#include "fabm_driver.h"
module medusa_carbonate

!******************************************************************************************
!                       MEDUSA carbonate chemistry
!******************************************************************************************
!
! Calculate the carbonate chemistry for the whole ocean
!
! Global simulations imply calculations on the first timestep and every month subsequently;
! the resulting 3D field of omega calcite is used to determine the depth of the CCD
! This is achieved by implemeting FABM capability for scheduling calls to subroutines
!
! Air-sea exchange of carbon dioxide is also calculated in this module
!------------------------------------------------------------------------------------------
!
   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   type,extends(type_base_model),public :: type_medusa_carbonate
      type (type_dependency_id)            :: id_ZALK
      type (type_state_variable_id)        :: id_ZDIC
      type (type_dependency_id)            :: id_temp,id_salt,id_dens,id_pres,id_depth,id_dz
      type (type_horizontal_dependency_id) :: id_PCO2A,id_kw660,id_fr_i
      type (type_diagnostic_variable_id)   :: id_ph,id_pco2,id_CarbA,id_BiCarb,id_Carb,id_TA_diag
      type (type_diagnostic_variable_id)   :: id_Om_cal,id_Om_arg
      type (type_horizontal_diagnostic_variable_id) ::  id_fairco2,id_pco2s,id_ATM_PCO2

   contains
     procedure :: initialize
     procedure :: do
     procedure :: do_surface

   end type

   public :: CO2_dynamics,co2dyn,polyco,CaCO3_Saturation

contains
    
    subroutine initialize(self,configunit)

      class (type_medusa_carbonate), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

     call self%register_state_dependency(self%id_ZDIC,'DiC','mmol C/m^3','dissolved inorganic carbon')
     call self%register_dependency(self%id_ZALK,'ALK','meq/m**3','total alkalinity')

     call self%register_diagnostic_variable(self%id_ph,    'PH3',    '-',      'Ocean pH 3D',standard_variable=standard_variables%ph_reported_on_total_scale)
     call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  '1e-6',    'partial pressure of CO2')
     call self%register_diagnostic_variable(self%id_CarbA, 'CarbA', 'mmol/m^3','carbonic acid concentration')
     call self%register_diagnostic_variable(self%id_BiCarb,'BiCarb','mmol/m^3','bicarbonate concentration')
     call self%register_diagnostic_variable(self%id_Carb,  'Carb',  'mmol/m^3','carbonate concentration')
     call self%register_diagnostic_variable(self%id_Om_cal,'OM_CAL3','-','Omega calcite 3D')
     call self%register_diagnostic_variable(self%id_Om_arg,'OM_ARG','-','aragonite saturation')
     self%id_ph%link%target%prefill = prefill_previous_value
     self%id_pco2%link%target%prefill = prefill_previous_value
     self%id_CarbA%link%target%prefill = prefill_previous_value
     self%id_BiCarb%link%target%prefill = prefill_previous_value
     self%id_Carb%link%target%prefill = prefill_previous_value
     self%id_Om_cal%link%target%prefill = prefill_previous_value
     self%id_Om_arg%link%target%prefill = prefill_previous_value

     call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
     call self%register_dependency(self%id_temp, standard_variables%temperature)
     call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
     call self%register_dependency(self%id_PCO2A,standard_variables%mole_fraction_of_carbon_dioxide_in_air)
     call self%register_dependency(self%id_depth, standard_variables%depth)
     call self%register_dependency(self%id_dens,standard_variables%density)
     call self%register_dependency(self%id_pres,standard_variables%pressure)
     call self%register_horizontal_dependency(self%id_kw660, 'KW660', 'm/s', 'gas transfer velocity')
     call self%register_dependency(self%id_fr_i,standard_variables%ice_area_fraction)
     call self%register_diagnostic_variable(self%id_fairco2,'CO2FLUX','mmol C/m^2/d','Air-sea CO2 flux')
     call self%register_diagnostic_variable(self%id_pco2s,'OCN_PCO2','-','Surface ocean pCO2')
     call self%register_diagnostic_variable(self%id_ATM_PCO2,'ATM_PCO2','ppmv','Atmospheric pCO2')
    end subroutine

    subroutine do(self,_ARGUMENTS_DO_)
    class (type_medusa_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ZDIC,ZALK,temp,salt,density,pres,depth
      real(rk) :: Om_cal,Om_arg,dcf,henry,TDIC,TALK
      real(rk) :: ph,pco2a,pco2w,h2co3,hco3,co3,k0co2
      integer :: iters

      _LOOP_BEGIN_
         _GET_(self%id_ZDIC,ZDIC)
         _GET_(self%id_ZALK,ZALK)
         _GET_(self%id_temp,temp)
         _GET_(self%id_salt,salt)
         _GET_HORIZONTAL_(self%id_pco2a,pco2a)
         _GET_(self%id_dens,density)

         _GET_(self%id_pres,pres)
         _GET_(self%id_depth,depth)

    call CO2_dynamics(temp,salt,depth,ZDIC,ZALK,pCO2a,pco2w,ph,h2co3,hco3,co3,henry,om_cal,om_arg,TDIC,TALK,dcf,iters)

    _SET_DIAGNOSTIC_(self%id_ph,pH)
    _SET_DIAGNOSTIC_(self%id_pco2,pco2w)
    _SET_DIAGNOSTIC_(self%id_CarbA, h2co3)
    _SET_DIAGNOSTIC_(self%id_Bicarb,hco3)
    _SET_DIAGNOSTIC_(self%id_Carb,  co3)
    _SET_DIAGNOSTIC_(self%id_Om_cal,Om_cal)
    _SET_DIAGNOSTIC_(self%id_Om_arg,Om_arg)

      _LOOP_END_
   end subroutine

   subroutine CO2_dynamics( T, S, Pr, DIC, TALK, pco2a, pco2w, ph, h2co3, bicarb, &
                            carb, henry, om_cal, om_arg, TCO2, TA, dcf, iters ) 

      real(rk),intent( in )    :: T,S,Pr
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
      integer,intent( inout ) :: iters       ! iterations to convergence
      real(rk) :: a, b, c 
      real(rk) :: ca, bc, cb 
      real(rk) :: pco2water, fairco2

          a   =  8.24493e-1_rk - 4.0899e-3_rk*T +  7.6438e-5_rk*T**2 - 8.2467e-7_rk*T**3 + 5.3875e-9_rk*T**4 
          b   = -5.72466e-3_rk + 1.0227e-4_rk*T - 1.6546e-6_rk*T**2  
          c   = 4.8314e-4_rk
          dcf = (999.842594_rk + 6.793952e-2_rk*T- 9.095290e-3_rk*T**2 + 1.001685e-4_rk*T**3 & 
                - 1.120083e-6_rk*T**4 + 6.536332e-9_rk*T**5+a*S+b*S**1.5_rk+c*S**2)

          TA    = TALK / (1.0e3_rk*dcf) 
          TCO2  = DIC  / (1.0e3_rk*dcf) 

! Call the parent routine for the carbonate system 

      call CO2DYN ( TCO2, TA, T, S, pco2a, &     ! inputs 
          pco2water, pH, HENRY, ca, bc, cb, iters )  ! outputs 

          pco2w  = pco2water * (1.0e6_rk)  ! partial pressure of co2 in water  
          TA     = TA * (1.0e6_rk)         ! total alkalinity (umol/kg) 
          h2co3  = ca * (1.0e3_rk*dcf)     ! carbonic acid concentration (mmol/m3) 
          bicarb = bc * (1.0e3_rk*dcf)     ! bicarbonate ion concentration (mmol/m3) 
          carb   = cb * (1.0e3_rk*dcf)     ! carbonate ion concentration (mmol/m3) 
          TCO2   = TCO2 * (1.0e6_rk)       ! total C or DIC in units of umol/kg 

! Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states

      call CaCO3_Saturation ( T, S, Pr, cb, &  ! inputs
          om_cal, om_arg )                    ! outputs

   end subroutine CO2_dynamics

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine CO2dyn ( TCO2, TA, T, S, PCO2A, PCO2, PH, HENRY, ca, bc, cb, iters ) 
!     This subroutine acts as an interface to the Haltafall iteration, setting options etc. 

      REAL(rk), INTENT( in )    :: TCO2, TA, T, S, PCO2A 
      REAL(rk), INTENT( inout ) :: PCO2, PH, HENRY, ca, bc, cb 
      INTEGER,  INTENT( inout ) :: iters
      REAL(rk) :: PRSS 
      INTEGER  :: MCONC, MKVAL, ICONST, ICALC 
      PARAMETER ( MCONC = 9,MKVAL = 4 ) 
      REAL(rk), DIMENSION(MKVAL) :: AKVAL 
      REAL(rk), DIMENSION(MCONC) :: CONCS 
      ICONST = 6 
      PRSS = 1._rk
      CONCS(1) = TCO2 
      CONCS(2) = TA 
      ICALC = 1 
      CALL POLYCO(PRSS,T,S,CONCS,MCONC,AKVAL,MKVAL,ICALC,ICONST,iters) 
      if(iters.eq.25) then
         CONCS(3) = PCO2A * 1E-06_rk ! atmospheric value
         CONCS(4) = 8.1_rk           ! global average pH
         CONCS(5) = TCO2 * 0.005_rk  ! some "standard" values plucked
         CONCS(6) = TCO2 * 0.865_rk  ! from pages 5-6 of Zeebe & Wolf-
         CONCS(7) = TCO2 * 0.130_rk  ! Gladow (2001)
         AKVAL(1) = 0.5E-01_rk       ! trial-and-error value
      endif

      PCO2  = CONCS(3) 
      PH    = CONCS(4) 
      ca    = CONCS(5) 
      bc    = CONCS(6) 
      cb    = CONCS(7) 
      HENRY = AKVAL(1)

      end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine POLYCO(PD,TD,SD,CONCS,NCONC,AKVAL,NKVAL,ICALC,ICONST,iters) 
! MASTER SUBROUTINE FOR CALCULATION OF THE CO2 SYSTEM THERMODYNAMICS 

! EXPLANATION OF POLYCO PARAMETERS 

!       AXY (16/08/11): Yuri change for non-surface carbonate chemistry  
!       P - PRESSURE IN ATMOSPHERES (P<>1 CODED only for ICONST=3 and ICONST=6) 
!       T - TEMPERATURE IN DEG.C 
!       S - SALINITY IN PPT 
!      CONCS(1) - TOTAL C (MOL/KG) 
!      CONCS(2) - TOTAL ALKALINITY (MOL/KG) 
!      CONCS(3) - PCO2 (ATM) 
!      CONCS(4) - PH 
!      CONCS(5) - {H2CO3} (MOL/KG) 
!      CONCS(6) - {HCO3} (MOL/KG) 
!      CONCS(7) - {CO3} (MOL/KG) 
!      CONCS(8) - CARBONATE ALKALINITY  ) FOR ICONST = 4,5,6 
!      CONCS(9) - BORATE ALKALINITY     )       ONLY 
!      NCONC - SIZE OF CONCS ARRAY (7 FOR ICONST=1,2,3; 9 FOR ICONST 
!      AKVAL(1) - KP (HENRY'S LAW CONSTANT) (MOL/KG/ATM) 
!      AKVAL(2) - K1C (H2CO3 DISSOCIATION) (MOL/KG) 
!      AKVAL(3) - K2C (HCO3 DISSOCIATION) (MOL/KG) 
!      AKVAL(4) - KB (B(OH)3 DISSOCIATION) (MOL/KG)  FOR ICONST=4,5,6 
!      NKVAL - SIZE OF AKVAL ARRAY (3 FOR ICONST=1,2,3; 4 FOR ICONST= 
!      ICALC - SELECTION OF THE TWO INPUT PARAMETERS: 
!      ICALC = 1  TOTAL C AND ALKALINITY 
!      ICALC = 2  TOTAL C AND PCO2 
!      ICALC = 3  TOTAL C AND PH 
!      ICALC = 4  ALKALINITY AND PCO2 
!      ICALC = 5  ALKALINITY AND PH 
!      ICALC = 6  PCO2 AND PH 
!      ICALC = 7  CALCULATE CONSTANTS AKVAL ONLY 
!      ICONST - SELECTION OF PH SCALE AND COMPONENTS: 
!      ICONST = 1  NBS PH SCALE 
!      ICONST = 2  HANSSON'S SCALE (SWS WITHOUT FLUORIDE) 
!      ICONST = 3  SWS PH SCALE 
!      ICONST = 4  AS 1 BUT INCLUDING BORATE IN THE CALCULATION 
!      ICONST = 5  AS 2 BUT INCLUDING BORATE IN THE CALCULATION 
!      ICONST = 6  AS 3 BUT INCLUDING BORATE IN THE CALCULATION 
!  NOTE: FOR ICONST=1,2,3 CONCS(2) REPRESENTS CARBONATE ALKALINITY SINC 
!        BORATE IS NOT INCLUDED IN THE CALCULATION. FOR ICONST=4,5,6 CO 
!        REPRESENTS TOTAL ALKALINITY (CARBONATE + BORATE), THE COMPONEN 
!        WHICH ARE GIVEN IN CONCS(8) AND CONCS(9) 

      REAL(rk) :: PMIN, PMAX, SMIN, SMAX, TMIN, TMAX, &  
     &     PD, TD, SD, P, T, S, BTOT 
      INTEGER  :: MINJC, MAXJC, MINJK, MAXJK, MINCAL, MAXCAL, MINCON,  & 
     &     MAXCON, NCONC, NKVAL, ICALC, ICONST, IC, iters 
      LOGICAL  :: BORON 

      PARAMETER(MINJC=7,MAXJC=9,MINJK=3,MAXJK=4) 
      PARAMETER(MINCAL=1,MAXCAL=7,MINCON=1,MAXCON=6) 
      PARAMETER(PMIN=0.99999_rk,PMAX=1.00001_rk,SMIN=0.0_rk,  & 
     &  SMAX=45.0_rk,TMIN=-3.0_rk,TMAX=40.0_rk) 
      REAL(rk), DIMENSION(NKVAL) :: AKVAL 
      REAL(rk), DIMENSION(NCONC) :: CONCS 

      P = PD
      S = SD
      T = TD

      IF(ICALC.LT.MINCAL.OR.ICALC.GT.MAXCAL) STOP 'POLYCO - ICALC OUT OR RANGE' 
      IF(ICONST.LT.MINCON.OR.ICONST.GT.MAXCON) STOP 'POLYCO - ICONST OUT OF RANGE' 

      BORON=(ICONST.GT.3)
      IF(BORON) THEN 

        IC=ICONST-3 
        BTOT=0.0004128_rk*S/35.0_rk 

        IF(NCONC.NE.MAXJC) STOP 'POLYCO - WRONG NCONC VALUE' 
        IF(NKVAL.NE.MAXJK) STOP 'POLYCO - WRONG NKVAL VALUE' 

      ELSE 

        IC=ICONST 

        IF(NCONC.NE.MINJC) STOP 'POLYCO - WRONG NCONC VALUE' 
        IF(NKVAL.NE.MINJK) STOP 'POLYCO - WRONG NKVAL VALUE' 

      ENDIF 

      CALL CO2SET(P,T,S,AKVAL,NKVAL,IC) 

      IF(ICALC.LT.MAXCAL)  & 
     & CALL CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT,iters) 

      end subroutine

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

        k1=10**(-1*(3670.7_rk/TK - 62.008_rk + 9.7944_rk*dlogTK - &
     &          0.0118_rk * S + 0.000116_rk*S2))

        k2=10**(-1*(1394.7_rk/TK + 4.777_rk - &
     &          0.0184_rk*S + 0.000118_rk*S2))

        delta=-25.5_rk+0.1271_rk*T
        kappa=(-3.08_rk+0.0877_rk*T)/1000._rk
        k1=k1*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

        delta=-15.82_rk-0.0219_rk*T
        kappa=(1.13_rk-0.1475_rk*T)/1000._rk
        k2=k2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

        kb=exp((-8966.90_rk - 2890.53_rk*sqrtS - 77.942_rk*S + &
     &          1.728_rk*S15 - 0.0996_rk*S2)/TK + &
     &          (148.0248_rk + 137.1942_rk*sqrtS + 1.62142_rk*S) + &
     &          (-24.4344_rk - 25.085_rk*sqrtS - 0.2474_rk*S) * &
     &          dlogTK + 0.053105_rk*sqrtS*TK)

        delta=-29.48_rk+0.1622_rk*T-0.002608_rk*T**2.0_rk
        kappa=-2.84_rk/1000._rk
        kb=kb*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

        delta=-25.60_rk+0.2324_rk*T-0.0036246_rk*T**2
        kappa=(-5.13_rk+0.0794_rk*T)/1000._rk

      AKVAL(1)=AKVAL(1)*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))
      AKVAL(2) = k1
      AKVAL(3) = k2
      AKVAL(4) = kb 

      end if
   end subroutine
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT,iters) 
! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI 
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975) 
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE. 

      INTEGER             :: NCONC, NKVAL, ICALC, II, KARL, LQ, iters 
      INTEGER             :: COUNTER,C_CHECK,C_SW,III 

      REAL(rk)            :: CTOT,ALK,PCO2,PH,H2CO3,HCO3,CO3,ALKC 
      REAL(rk)            :: ALKB,AKP,AK1C,AK2C,AKB,BTOT 
      REAL(rk)            :: AKR,AHPLUS 
      REAL(rk)            :: PROD,tol1,tol2,tol3,tol4,steg,fak 
      REAL(rk)            :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z 

      REAL(rk), DIMENSION(NCONC) :: CONCS 
      REAL(rk), DIMENSION(NKVAL) :: AKVAL 
      REAL(rk), DIMENSION(9)     :: CONCS2 
      REAL(rk), DIMENSION(4)     :: AKVAL2 

      real(rk),DIMENSION(:) :: old_CONCS(NCONC) 

      EQUIVALENCE (CTOT  , CONCS2(1)), (ALK   , CONCS2(2)),  & 
     &            (PCO2  , CONCS2(3)), (PH    , CONCS2(4)),  & 
     &            (H2CO3 , CONCS2(5)), (HCO3  , CONCS2(6)),  & 
     &            (CO3   , CONCS2(7)), (ALKC  , CONCS2(8)),  & 
     &            (ALKB  , CONCS2(9)),                       & 
     &            (AKP   , AKVAL2(1)), (AK1C  , AKVAL2(2)),  & 
     &            (AK2C  , AKVAL2(3)), (AKB   , AKVAL2(4)) 

      LOGICAL             :: BORON,DONE 


!  DERIVING PH REQUIRES FOLLOWING LOOP TO CONVERGE. 
!  THIS SUBROUTINE RELIES ON CONVERGENCE.  IF THE ENVIRONMENTAL 
!  CONDITIONS DO NOT ALLOW FOR CONVERGENCE (IN 3D MODEL THIS IS 
!  LIKELY TO OCCUR NEAR LOW SALINITY REGIONS) THE MODEL WILL  
!  BE STUCK IN THE LOOP.  TO AVOID THIS A CONVERGENCE CONDITION 
!  IS PUT IN PLACE TO SET A FLAGG OF -99 IN THE PH VAR FOR NON CONVEGENCE.  
!  THE MODEL IS THEN ALLOWED TO CONTINUE. 'COUNTER, C_SW,C_CHECK' ARE  
!  THE LOCAL VARS USED. 
! C_SW = condition of convergence 0=yes, 1= no 
! COUNTER = number of iterations 
! C_CHECK = maximum number of iterations 

! SET COUNTER AND SWITCH TO ZERO AND OFF 

      COUNTER=0 
      C_SW=0 

! FROM EXPERIENCE IF THE ITERATIONS IN THE FOLLOWING DO LOOP 
! EXCEEDS 15 CONVERGENCE WILL NOT OCCUR.  THE OVERHEAD OF 25 ITERATIONS  
! IS OK FOR SMALL DOMAINS WITH 1/10 AND 1/15 DEG RESOLUTION. 
! I RECOMMEND A LOWER VALUE OF 15 FOR HIGHER RESOLUTION OR LARGER DOMAINS. 

      C_CHECK=25 

      DO 100 II=1,NCONC
        CONCS2(II)=CONCS(II) 

! IF CONVERGENCE IS NOT ACHIEVED THE LOCAL ARRAY CONCS MUST BE STORED TO
! ALLOW THE MODEL TO CONTINUE. THEREFORE .... 

! UPDATE OLD_CONCS 
        old_CONCS(II)=CONCS(II) 

100    CONTINUE 

      DO 110 II=1,NKVAL 
        AKVAL2(II)=AKVAL(II) 

110    CONTINUE 

      AKR = AK1C/AK2C 
      AHPLUS=10.0_rk**(-PH) 
      PROD=AKR*AKP*PCO2 

      IF(BORON) THEN 
        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN 

!         *** ALK, BTOT AND CTOT OR PCO2 FIXED *** 
!         *** ITERATIVE CALCULATION NECESSARY HERE 
!         SET INITIAL GUESSES AND TOLERANCE 

          H2CO3=PCO2*AKP 
          CO3=ALK/10.0_rk 
          AHPLUS=1.0e-8_rk 
          ALKB=BTOT 
          TOL1=ALK/1.0e5_rk 
          TOL2=H2CO3/1.0e5_rk 
          TOL3=CTOT/1.0e5_rk 
          TOL4=BTOT/1.0e5_rk 

!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS 

          KARL=1 
          STEG=2.0_rk 
          FAK=1.0_rk 
          STEGBY=0.4_rk 
10        DONE=.TRUE. 

! SET COUNTER UPDATE. FLAG 10 IS THE POINT OF RETURN FOR  
! THE CONVERGENCE CONDITION 

          COUNTER=COUNTER+1 
          iters  = COUNTER 

! CHECK IF CONVERGENCE HAS OCCURED IN THE NUMBER OF  
! ACCEPTABLE ITTERATIONS.  SET C_SW TO LET MODEL KNOW
! WHAT TO DO AT THE END OF THE SUBROUTINE 

          if(counter.ge.c_check)then
             c_sw=1
             GOTO 123

          endif

          IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***

            Y=AHPLUS*AHPLUS*CO3/(AK1C*AK2C)
            IF(ABS(Y-H2CO3).GT.TOL2) THEN

              CO3=CO3*H2CO3/Y
              DONE=.FALSE.

            ENDIF

          ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED *** 

            Y=CO3*(1.0_rk+AHPLUS/AK2C+AHPLUS*AHPLUS/(AK1C*AK2C))
            IF(ABS(Y-CTOT).GT.TOL3) THEN

              CO3=CO3*CTOT/Y
              DONE=.FALSE.

            ENDIF

          ENDIF

          Y=ALKB*(1.0_rk+AHPLUS/AKB)
          IF(ABS(Y-BTOT).GT.TOL4) THEN

            ALKB=ALKB*BTOT/Y
            DONE=.FALSE.

          ENDIF


! Alkalinity is equivalent to -(total H+), so the sign of W is opposite 
! to that normally used 

          Y=CO3*(2.0_rk+AHPLUS/AK2C)+ALKB 

          IF(ABS(Y-ALK).GT.TOL1) THEN 

            DONE=.FALSE. 
            X=LOG(AHPLUS) 
 

            IF ( Y-ALK .GE. 0.0_rk ) THEN 

            W=1.0_rk 

            ELSE 

            W=-1.0_rk 

            ENDIF 

            IF(W.GE.0.0_rk) THEN 

              X1=X 
              Y1=Y 

            ELSE 

              X2=X 
              Y2=Y 

            ENDIF 

            LQ=KARL 

            IF(LQ.EQ.1) THEN 

              KARL=2*NINT(W) 

            ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0.) THEN 

              FAK=0.5_rk
              KARL=3 

            ENDIF 

            IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN 

              W=(X2-X1)/(Y2-Y1) 
              X=X1+W*(ALK-Y1) 

            ELSE 

              STEG=STEG*FAK 
              X=X+STEG*W 

            ENDIF

            AHPLUS=EXP(X)

          ENDIF

          IF(.NOT.DONE) GOTO 10

          HCO3=CO3*AHPLUS/AK2C

          IF(ICALC.EQ.4) THEN

            CTOT=H2CO3+HCO3+CO3

          ELSEIF(ICALC.EQ.1) THEN

            H2CO3=HCO3*AHPLUS/AK1C
            PCO2=H2CO3/AKP

          ENDIF

          PH=-LOG10(AHPLUS)
          ALKC=ALK-ALKB

        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***

          Y=SQRT(PROD*(PROD-4.0_rk*AKP*PCO2+4.0_rk*CTOT))

          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0_rk
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2.0_rk*CO3
          AHPLUS=AK1C*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB)
          ALK=ALKC+ALKB

        ELSEIF(ICALC.EQ.3) THEN 
!         *** CTOT, PH AND BTOT FIXED *** 

          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C) 
          CO3=FACTOR*AK1C*AK2C 
          HCO3=FACTOR*AK1C*AHPLUS 
          H2CO3=FACTOR*AHPLUS*AHPLUS 
          PCO2=H2CO3/AKP 
          ALKC=HCO3+2.0_rk*CO3 
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB) 
          ALK=ALKC+ALKB 

        ELSEIF(ICALC.EQ.5) THEN 
!         *** ALK, PH AND BTOT FIXED *** 

          ALKB=BTOT/(1.0_rk+AHPLUS/AKB) 
          ALKC=ALK-ALKB 
          HCO3=ALKC/(1.0_rk+2.0_rk*AK2C/AHPLUS) 
          CO3=HCO3*AK2C/AHPLUS 
          H2CO3=HCO3*AHPLUS/AK1C 
          PCO2=H2CO3/AKP 
          CTOT=H2CO3+HCO3+CO3 

        ELSEIF(ICALC.EQ.6) THEN 
!         *** PCO2, PH AND BTOT FIXED *** 

          ALKB=BTOT/(1.0_rk+AHPLUS/AKB) 
          H2CO3=PCO2*AKP 
          HCO3=H2CO3*AK1C/AHPLUS 
          CO3=HCO3*AK2C/AHPLUS 
          CTOT=H2CO3+HCO3+CO3 
          ALKC=HCO3+2.0_rk*CO3 
          ALK=ALKC+ALKB 

        ENDIF 

      ELSE 

        IF(ICALC.EQ.1) THEN 
!         *** CTOT AND ALK FIXED ***
 
          TERM=4.0_rk*ALK+CTOT*AKR-ALK*AKR 
          Z=SQRT(TERM*TERM+4.0_rk*(AKR-4.0_rk)*ALK*ALK) 
          CO3=(ALK*AKR-CTOT*AKR-4.0_rk*ALK+Z)/(2.0_rk*(AKR-4.0_rk)) 
          HCO3=(CTOT*AKR-Z)/(AKR-4.0_rk) 
          H2CO3=CTOT-ALK+CO3 
          PCO2=H2CO3/AKP 
          PH=-LOG10(AK1C*H2CO3/HCO3)

        ELSEIF(ICALC.EQ.2) THEN

!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4.0_rk*AKP*PCO2+4.0_rk*CTOT)) 
          H2CO3=PCO2*AKP 
          HCO3=(Y-PROD)/2.0_rk 
          CO3=CTOT-H2CO3-HCO3 
          ALK=HCO3+2.0_rk*CO3 
          PH=-LOG10(AK1C*H2CO3/HCO3) 

        ELSEIF(ICALC.EQ.3) THEN 
!         *** CTOT AND PH FIXED *** 

          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C) 
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP 
          ALK=HCO3+2.0_rk*CO3 

        ELSEIF(ICALC.EQ.4) THEN 
!         *** ALK AND PCO2 FIXED ***

          TERM=SQRT((8.0_rk*ALK+PROD)*PROD) 
          CO3=ALK/2.0_rk+PROD/8.0_rk-TERM/8.0_rk 
          HCO3=-PROD/4.0_rk+TERM/4.0_rk 
          CTOT=CO3+HCO3+H2CO3 
          PH=-LOG10(AK1C*H2CO3/HCO3) 

        ELSEIF(ICALC.EQ.5) THEN 
!         *** ALK AND PH FIXED *** 

          HCO3=ALK/(1.0_rk+2.0_rk*AK2C/AHPLUS) 
          CO3=HCO3*AK2C/AHPLUS 
          H2CO3=HCO3*AHPLUS/AK1C 
          PCO2=H2CO3/AKP 
          CTOT=H2CO3+HCO3+CO3 

        ELSEIF(ICALC.EQ.6) THEN 

!         *** PCO2 AND PH FIXED *** 

          H2CO3=PCO2*AKP 
          HCO3=H2CO3*AK1C/AHPLUS 
          CO3=HCO3*AK2C/AHPLUS 
          CTOT=H2CO3+HCO3+CO3 
          ALK=HCO3+2.0_rk*CO3 

        ENDIF 

      ENDIF 

 
      DO 120 II=1,NCONC 

        CONCS(II)=CONCS2(II) 

120    CONTINUE 


! C_SW IS SET AT 0 TO START 
! THEN IF NON CONVERGENCE C_SW SET TO 1 

123   IF(C_SW.EQ.1)THEN 
! IF NON CONVERGENCE, THE MODEL REQUIRES CONCS TO CONTAIN USABLE VALUES. 
! BEST OFFER BEING THE OLD CONCS VALUES WHEN CONVERGENCE HAS BEEN  
! ACHIEVED 

        DO II=1,NCONC
           CONCS(II)=OLD_CONCS(II) 
        END DO 

! SPECIFIC CARBONATE VALUES TO PUSH CODE ON THROUGH THE 
! NON CONVERGENCE CONDITIONS 
! PCO2W = 0 SO CO2 UPTAKE WILL BE ENCOURAGED 
!       CONCS(3)=O3C(III)*0.005_fp8/1.e6_fp8 

! RESET SWITCH FOR NEXT CALL TO THIS SUBROUTINE 

        C_SW=0 

      ENDIF 

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
      if (D .GT. 0._rk) then
        dV = -48.76_rk + 0.5304_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspc = Kspc*exp(tmp1)
        logKspc = log10(Kspc)
      end if

        tmp1 = -171.945_rk - 0.077993_rk*Tk + 2903.293_rk / Tk + 71.595_rk* log10(Tk)
        tmp2 = + (-0.068393_rk + 0.0017276_rk*Tk + 88.135_rk/Tk)*SQRT(S)
        tmp3 = - 0.10018_rk*S + 0.0059415_rk*S**1.5_rk
        logKspa = tmp1 + tmp2 + tmp3
        Kspa = 10._rk**logKspa

! correction for pressure for aragonite
      if (D .GT. 0._rk) then
        dV = -46._rk + 0.5304_rk*Tc
        dK = -11.76_rk/1.e3_rk + (0.3692_rk/1.e3_rk) * Tc
        tmp1 = -(dV/(R*Tk))*P + (0.5_rk*dK/(R*Tk))*P*P
        Kspa = Kspa*exp(tmp1)
        logKspa = log10(Kspa)
      end if

! calculate saturation states
        Om_cal = (CO3 * Ca) / Kspc
        Om_arg = (CO3 * Ca) / Kspa

   end subroutine CaCO3_SATURATION

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

   class(type_medusa_carbonate),intent(in) :: self
   !**********************************
   ! Air-sea exchange of CO2
   !**********************************
   _DECLARE_ARGUMENTS_DO_SURFACE_

    real(rk) :: ZALK,ZDIC,temp,salt,sc,fwind,flux,kw660,fr_i
    real(rk) :: henry, pco2a, pco2, dcf, a, b, c, ca, bc, cb, ph, TA, TCO2
    integer :: iters

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_HORIZONTAL_(self%id_kw660,kw660)
   _GET_HORIZONTAL_(self%id_fr_i,fr_i)
   _GET_(self%id_ZALK,ZALK)
   _GET_(self%id_ZDIC,ZDIC)
   _GET_HORIZONTAL_(self%id_pco2a,pco2a)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ATM_PCO2,pco2a)  

          a   =  8.24493e-1_rk - 4.0899e-3_rk*temp +  7.6438e-5_rk*temp**2 - 8.2467e-7_rk*temp**3 + 5.3875e-9_rk*temp**4 
          b   = -5.72466e-3_rk + 1.0227e-4_rk*temp - 1.6546e-6_rk*temp**2  
          c   = 4.8314e-4_rk
          dcf = (999.842594_rk + 6.793952e-2_rk*temp- 9.095290e-3_rk*temp**2 + 1.001685e-4_rk*temp**3 & 
                - 1.120083e-6_rk*temp**4 + 6.536332e-9_rk*temp**5+a*salt+b*salt**1.5_rk+c*salt**2) 

          TA    = ZALK  / (1.0e3_rk*dcf) 
          TCO2  = ZDIC  / (1.0e3_rk*dcf)

    call CO2dyn ( TCO2, TA, temp, salt, pco2a, PCO2, PH, HENRY, ca, bc, cb, iters )

    sc    = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.0432190_rk*temp**3._rk
    fwind = kw660 * (sc/660._rk)**(-0.5_rk)
    fwind = fwind * 24._rk/100._rk

    ! Calculate air-sea flux, correct for sea-ice
    flux = fwind * henry * ( pco2a - pco2 * 1.0e6_rk) * dcf / 1000._rk
    flux = (1._rk - fr_i) * flux
    _SET_SURFACE_EXCHANGE_(self%id_ZDIC,flux /86400._rk)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fairco2,flux)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pco2s,PCO2 * 1.0e6_rk)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module
