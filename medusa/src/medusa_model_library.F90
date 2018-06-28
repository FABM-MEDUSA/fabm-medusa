module medusa_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use medusa_light
   use medusa_pelagic
   use medusa_fast_detritus
   use medusa_oxygen
   use medusa_carbonate
   use medusa_benthic
   use gas_transfer
   use medusa_iron_scav

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: medusa_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('medusa_pelagic');                          allocate(type_medusa_pelagic::model)
         case ('medusa_fast_detritus');                    allocate(type_medusa_fast_detritus::model)
         case ('medusa_oxygen');                           allocate(type_medusa_oxygen::model)
         case ('medusa_carbonate');                        allocate(type_medusa_carbonate::model)
         case ('medusa_benthic');                          allocate(type_medusa_benthic::model)
         case ('gas_transfer');                            allocate(type_gas_transfer::model)
         case ('medusa_light');                            allocate(type_medusa_light::model)
         case ('medusa_iron_scav');                        allocate(type_medusa_iron_scav::model)
        ! Add new models here
      end select
   end subroutine create

end module
