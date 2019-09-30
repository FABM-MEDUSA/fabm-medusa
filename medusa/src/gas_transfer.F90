#include "fabm_driver.h"

!
!********************************************************************
!            FABM-MEDUSA gas transfer velocity
!********************************************************************
!     Original author: Andrew Yool (15/10/04, revised 12/06/15)
! 
! This subroutine uses near-surface wind speed to calculate gas
! transfer velocity for use in CO2 and O2 exchange calculations.
!
! Note that the parameterisation of Wanninkhof quoted here is a
! truncation of the original equation.  It excludes a chemical
! enhancement function (based on temperature), although such
! temperature dependence is reported negligible by Etcheto & 
! Merlivat (1988).
!
! Note also that in calculating scalar wind, the variance of the
! wind over the period of a timestep is ignored.  Some authors,
! for instance OCMIP-2, favour including some reference to the
! variability of wind.  However, their wind fields are averaged
! over relatively long time periods, and so this issue may be
! safely (!) ignored here.
! 
!     eqn       choice of parameterisation:
!               1 = Liss & Merlivat (1986)    [approximated]
!               2 = Wanninkhof (1992)         [sans enhancement]
!               3 = Nightingale et al. (2000) [good]
!               4 = Nightingale et al. (2000) [better]
!               5 = Nightingale et al. (2000) [best]
!               6 = OCMIP-2                   [sans variability]
!               7 = Wanninkhof (2014)         [assumes 6h avg winds]
!
! Wanninkhof, R. (2014). Relationship between wind speed and gas 
! exchange over the ocean revisited. LIMNOLOGY AND OCEANOGRAPHY-METHODS
! 12, 351-362, doi:10.4319/lom.2014.12.351

module gas_transfer

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_gas_transfer
        type (type_horizontal_dependency_id) :: id_wnd
        type (type_horizontal_diagnostic_variable_id)   :: id_kw660,id_wnd_diag
        integer :: eqn

   contains

     procedure :: initialize
     procedure :: do_surface

   end type

contains
    
    subroutine initialize(self,configunit)

     class (type_gas_transfer), intent(inout), target :: self
     integer,                      intent(in)            :: configunit

     call self%get_parameter(self%eqn,'eqn','-','choice of gas transfer coefficients', default = 7)
     call self%register_dependency(self%id_wnd,  standard_variables%wind_speed)

     call self%register_diagnostic_variable(self%id_kw660,'KW660','m/s','Piston velocity', output=output_time_step_averaged, source=source_do_surface)
     call self%register_diagnostic_variable(self%id_wnd_diag,'WIND','m/s','Surface scalar wind',source=source_do_surface)
    end subroutine

    subroutine do_surface(self,_ARGUMENTS_DO_)
    class (type_gas_transfer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: wnd,kw660
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

           _GET_HORIZONTAL_(self%id_wnd,wnd)
           _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wnd_diag,wnd)

! Calculate gas transfer velocity (cm/h)
           tmp_k = (a(self%eqn) * wnd**2._rk) + (b(self%eqn) * wnd)

! Convert tmp_k from cm/h to m/s
           kw660 = tmp_k / (3600._rk * 100._rk)

  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kw660,kw660)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module gas_transfer
