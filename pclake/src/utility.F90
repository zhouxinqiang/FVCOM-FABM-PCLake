#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pclake_utility
!
! !INTERFACE:
   module pclake_utility
!
!  Originally implemented by Fenjuan Hu (fenjuan@bios.au.dk)
!
! !USES:
   use fabm_types
   use fabm_expressions
   implicit none
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public uFunTmAbio
   public uFunTmBio
   public uFunTmVeg
   public DayLength
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter                 :: cTmRef=20.0_rk
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   pure real(rk) function uFunTmAbio(uTm,cTheta)

! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cTheta
!
! !REVISION HISTORY:
!  Original author(s): Fenjuan Hu
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmAbio = cTheta**(uTm-cTmRef)
   end function uFunTmAbio
!EOC

!-----------------------------------------------------------------------
!BOP
!
!
! !INTERFACE:
   pure real(rk) function uFunTmBio(uTm,cSigTm,cTmOpt)
!
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cSigTm,cTmOpt
!
! !REVISION HISTORY:
!  original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmBio = exp(-0.5_rk/cSigTm**2 *((uTm-cTmOpt)**2- (cTmRef -cTmOpt)**2))
   end function uFunTmBio
!EOC

!-----------------------------------------------------------------------
!BOP
!
!
! !INTERFACE:
   pure real(rk) function uFunTmVeg(uTm,cQ10Veg)
!

! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cQ10Veg
!
! !REVISION HISTORY:
!  Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmVeg = ((cQ10Veg )** (0.1_rk * (uTm - cTmRef)))
   end function uFunTmVeg
!EOC

!-----------------------------------------------------------------------
!BOP
!
!
! !INTERFACE:
   pure real(rk) function DayLength(latitude,day)
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: latitude
   integer, intent(in)                 :: day
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !Local variables:
   real(rk),parameter :: Pi=3.14159265358979_rk
   real(rk),parameter :: deg2rad=Pi/180._rk
   real(rk),parameter :: ecliptic=23.439_rk*deg2rad
   real(rk)           :: m 
!EOP
!-----------------------------------------------------------------------
!BOC
!  http://www.gandraxa.com/length_of_day.xml
   m = 1._rk-tan(deg2rad*latitude)*tan(ecliptic*cos(0.0172_rk*Day))
!   DayLength = acos(1._rk-m)/Pi*24._rk
   DayLength = acos(1._rk-m)/Pi ! fraction of day with sun shine
   end function DayLength
!EOC

!-----------------------------------------------------------------------

   end module pclake_utility

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
