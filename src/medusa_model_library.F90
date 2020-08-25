module medusa_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use medusa_light
   use medusa_pelagic
   use medusa_fast_detritus
   use medusa_oxygen
   use medusa_carbonate
   use medusa_benthic
   use gas_transfer
   use medusa_iron_scav
   use medusa_ccd

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: medusa_model_factory

contains

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('pelagic');       allocate(type_medusa_pelagic::model)
         case ('fast_detritus'); allocate(type_medusa_fast_detritus::model)
         case ('oxygen');        allocate(type_medusa_oxygen::model)
         case ('carbonate');     allocate(type_medusa_carbonate::model)
         case ('benthic');       allocate(type_medusa_benthic::model)
         case ('gas_transfer');  allocate(type_gas_transfer::model)
         case ('light');         allocate(type_medusa_light::model)
         case ('iron_scav');     allocate(type_medusa_iron_scav::model)
         case ('ccd');           allocate(type_medusa_ccd::model)
        ! Add new models here
      end select
   end subroutine create

end module
