MODULE tide_mod
   !!======================================================================
   !!                       ***  MODULE  tide_mod  ***
   !! Compute nodal modulations corrections and pulsations
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   !JT USE oce             ! ocean dynamics and tracers
   USE lib_mpp         ! distributed memory computing library

   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE daymod         ! calendar
   USE in_out_manager ! I/O units
   USE ioipsl  , ONLY :   ymds2ju      ! for calendar


   IMPLICIT NONE
   PRIVATE

   PUBLIC   tide_harmo       ! called by tideini and diaharm modules
   PUBLIC   tide_init_Wave   ! called by tideini and diaharm modules
   PUBLIC   tide_init_calendar_options   ! called by tideini and diaharm modules

   ! davbyr: increase maximum number of harmonics from 19 to 34
   INTEGER, PUBLIC, PARAMETER ::   jpmax_harmo = 34   !: maximum number of harmonic

   TYPE, PUBLIC ::    tide
      CHARACTER(LEN=4) ::   cname_tide
      REAL(wp)         ::   equitide
      INTEGER          ::   nutide
      INTEGER          ::   nt, ns, nh, np, np1, shift
      INTEGER          ::   nksi, nnu0, nnu1, nnu2, R
      INTEGER          ::   nformula
   END TYPE tide

   TYPE(tide), PUBLIC, DIMENSION(jpmax_harmo) ::   Wave   !:

   REAL(wp) ::   sh_T, sh_s, sh_h, sh_p, sh_p1             ! astronomic angles
   REAL(wp) ::   sh_xi, sh_nu, sh_nuprim, sh_nusec, sh_R   !
   REAL(wp) ::   sh_I, sh_x1ra, sh_N                       !


   !JT origin angles
   REAL(wp) ::   sh_T_o, sh_s_o, sh_h_o, sh_p_o, sh_p1_o, sh_N_o                ! astronomic angles
   !REAL(wp) ::   sh_xi_o, sh_nu_o, sh_nuprim_o, sh_nusec_o, sh_R_o   !
   !REAL(wp) ::   sh_I_o, sh_x1ra_o                    !
   !JT origin angles

   !!JT
   INTEGER(KIND=8)  ::  days_since_origin
   LOGICAL  ::   ln_astro_verbose 
   !LOGICAL  ::   ln_tide_360_cal 
   !LOGICAL  ::   ln_tide_drift_time_cont_manual
   LOGICAL  ::   ln_tide_drift                  ! Do we want to run with "drifting" tides? (Namelist)
   LOGICAL  ::   ln_tide_compress               ! Do we want to run with "compressed" tides? (Namelist)
   INTEGER  ::   nn_tide_orig_yr,nn_tide_orig_mn,nn_tide_orig_dy            !JT

   !!JT

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tide_init_Wave
#     include "tide.h90"
   END SUBROUTINE tide_init_Wave


   SUBROUTINE tide_init_calendar_options

      INTEGER                              ::   ios


      ln_tide_drift = .FALSE.
      ln_tide_compress = .FALSE.

      NAMELIST/nam_tides360/ ln_tide_drift,ln_tide_compress,ln_astro_verbose,&
        & nn_tide_orig_yr,nn_tide_orig_mn,nn_tide_orig_dy

      ! read in Namelist. 
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
      READ   ( numnam_ref, nam_tides360, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_tides360 in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_tides360, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 ) CALL ctl_nam ( ios , 'nam_tides360 in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_tides360 )


      IF( lwp ) THEN
        WRITE(numout,*) " "
        WRITE(numout,*) "tide_harmo: nam_tides360 - 360 day tides "
        WRITE(numout,*) "~~~~~~~~~~~~~~~~~~~~~"
        WRITE(numout,*) "       tides360: allow tides to drift through year: ln_tide_drift = ",ln_tide_drift
        WRITE(numout,*) "       tides360: Compress tides, so around a 360 day year: ln_tide_compress = ",ln_tide_compress
        WRITE(numout,*) "       tides360:           USE ln_tide_compress  WITH CARE. INCOMPLETE."
        WRITE(numout,*) "       tides360: Increase output verbosity: ln_astro_verbose = ",ln_astro_verbose
        !WRITE(numout,*) "       tides360: Calculate time between origin and gregorian and 360 manually: ln_tide_drift_time_cont_manual = ",ln_tide_drift_time_cont_manual
        WRITE(numout,*) "       tides360: 360 day origin date year: nn_tide_orig_yr = ",nn_tide_orig_yr
        WRITE(numout,*) "       tides360: 360 day origin date month: nn_tide_orig_mn = ",nn_tide_orig_mn
        WRITE(numout,*) "       tides360: 360 day origin date day: nn_tide_orig_dy = ",nn_tide_orig_dy
        WRITE(numout,*) " "
      ENDIF

 
      IF( nleapy == 30 ) THEN
          IF ( ln_tide_drift .AND. ln_tide_compress ) THEN
              CALL ctl_stop( 'tide_harmo: nam_tides360: if 360 day calendar ln_tide_drift and ln_tide_compress cannot be true' )
          ENDIF
          

          IF ( ln_tide_drift   ) THEN
              WRITE(numout,*) "       tides360: Tides continuous so equinoctal tides drift through the year,"
              WRITE(numout,*) "                 as the S2-K2 beating occurs 5 days later every year."
          ENDIF

          IF ( ln_tide_compress   ) THEN
              WRITE(numout,*) "       tides360: The Tropical Year (and so some tidal periods) are compressed,"
              WRITE(numout,*) "                 so the tides repeat with an annual cycle, so the "
              WRITE(numout,*) "                 the S2-K2 beating is fixed relative to the calendar, but the "
              WRITE(numout,*) "                 M2 period varies slightly."
              WRITE(numout,*) "                 Use with care, as this requires more work."
          ENDIF

          IF ( ( .NOT. ln_tide_drift  ) .AND. ( .NOT. ln_tide_compress ) ) THEN
              WRITE(numout,*) "       tides360: Use the default NEMO tide code, where the tides are reset "
              WRITE(numout,*) "                 at the beginning of each month, leading to a slight discontinuity"
              WRITE(numout,*) "                 in the tides, and making tidal analysis difficult."
          ENDIF

      ELSE        
          WRITE(numout,*) "       tides360: Gregorian calendar so using standard tides"
      ENDIF

      IF ( ln_tide_compress   )  CALL astronomic_angle_origin

   END SUBROUTINE tide_init_calendar_options


   SUBROUTINE tide_harmo( pomega, pvt, put , pcor, ktide ,kc)

      !! Externally called by sbctide.F90/sbc_tide
      !! Externally named: omega_tide, v0tide, utide, ftide, ntide, nb_harmo
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide            ! Indice of tidal constituents
      INTEGER                , INTENT(in ) ::   kc               ! Total number of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pomega           ! pulsation in radians/s
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pvt, put, pcor   !
      !!----------------------------------------------------------------------
      !

!      INTEGER                              ::   ios


!      ln_tide_drift = .FALSE.
!      ln_tide_compress = .FALSE.

!      NAMELIST/nam_tides360/ ln_tide_drift,ln_tide_compress,ln_astro_verbose,&
!        & nn_tide_orig_yr,nn_tide_orig_mn,nn_tide_orig_dy

!      ! read in Namelist. 
!      !!----------------------------------------------------------------------
!      !
!      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
!      READ   ( numnam_ref, nam_tides360, IOSTAT=ios, ERR= 901 )
!901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_tides360 in reference namelist' )

!      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
!      READ  ( numnam_cfg, nam_tides360, IOSTAT = ios, ERR = 902 )
!902   IF( ios > 0 ) CALL ctl_nam ( ios , 'nam_tides360 in configuration namelist' )
!      IF(lwm) WRITE ( numond, nam_tides360 )


!      IF( lwp ) THEN
!        WRITE(numout,*) " "
!        WRITE(numout,*) "tide_harmo: nam_tides360 - 360 day tides "
!        WRITE(numout,*) "~~~~~~~~~~~~~~~~~~~~~"
!        WRITE(numout,*) "       tides360: allow tides to drift through year: ln_tide_drift = ",ln_tide_drift
!        WRITE(numout,*) "       tides360: Compress tides, so around a 360 day year: ln_tide_compress = ",ln_tide_compress
!        WRITE(numout,*) "       tides360:           USE ln_tide_compress  WITH CARE. INCOMPLETE."
!        WRITE(numout,*) "       tides360: Increase output verbosity: ln_astro_verbose = ",ln_astro_verbose
!        !WRITE(numout,*) "       tides360: Calculate time between origin and gregorian and 360 manually: ln_tide_drift_time_cont_manual = ",ln_tide_drift_time_cont_manual
!        WRITE(numout,*) "       tides360: 360 day origin date year: nn_tide_orig_yr = ",nn_tide_orig_yr
!        WRITE(numout,*) "       tides360: 360 day origin date month: nn_tide_orig_mn = ",nn_tide_orig_mn
!        WRITE(numout,*) "       tides360: 360 day origin date day: nn_tide_orig_dy = ",nn_tide_orig_dy
!        WRITE(numout,*) " "
!      ENDIF

! 
!      IF( nleapy == 30 ) THEN
!          IF ( ln_tide_drift .AND. ln_tide_compress ) THEN
!              CALL ctl_stop( 'tide_harmo: nam_tides360: if 360 day calendar ln_tide_drift and ln_tide_compress cannot be true' )
!          ENDIF
!          

!          IF ( ln_tide_drift   ) THEN
!              WRITE(numout,*) "       tides360: Tides continuous so equinoctal tides drift through the year,"
!              WRITE(numout,*) "                 as the S2-K2 beating occurs 5 days later every year."
!          ENDIF

!          IF ( ln_tide_compress   ) THEN
!              WRITE(numout,*) "       tides360: The Tropical Year (and so some tidal periods) are compressed,"
!              WRITE(numout,*) "                 so the tides repeat with an annual cycle, so the "
!              WRITE(numout,*) "                 the S2-K2 beating is fixed relative to the calendar, but the "
!              WRITE(numout,*) "                 M2 period varies slightly."
!              WRITE(numout,*) "                 Use with care, as this requires more work."
!          ENDIF

!          IF ( ( .NOT. ln_tide_drift  ) .AND. ( .NOT. ln_tide_compress ) ) THEN
!              WRITE(numout,*) "       tides360: Use the default NEMO tide code, where the tides are reset "
!              WRITE(numout,*) "                 at the beginning of each month, leading to a slight discontinuity"
!              WRITE(numout,*) "                 in the tides, and making tidal analysis difficult."
!          ENDIF

!      ELSE        
!          WRITE(numout,*) "       tides360: Gregorian calendar so using standard tides"
!      ENDIF
    
      CALL astronomic_angle
      CALL tide_pulse( pomega, ktide ,kc )
      CALL tide_vuf  ( pvt, put, pcor, ktide ,kc )

      !
   END SUBROUTINE tide_harmo


   SUBROUTINE astronomic_angle
      !!----------------------------------------------------------------------
      !!  tj is time elapsed since 1st January 1900, 0 hour, counted in julian
      !!  century (e.g. time in days divide by 36525)
      !!----------------------------------------------------------------------
      REAL(wp) ::   cosI, p, q, t2, t4, sin2I, s2, tgI2, P1, sh_tgn2, at1, at2
      REAL(wp) ::   zqy , zsy, zday, zdj, zhfrac

      
        ! JT
        ! Tides are added as boundary conditions, and as tidal potential. 
        !
        ! For the Boundaries, the complex tide amplitudes are give for each point. 
        ! This gives the amplitude and the phase for each consititent. 
        ! The tidal frequency is calculated from Wave in tide.h90 via tide_pulse below.
        ! A the start(ish) of everyday, astronomic_angle is called via tide_harmo 
        ! from SBC/sbctide.F90.
        ! The astronomic_angle specifies the location of the moon and sun etc at the given 
        ! model time. these are used to update the tidal phase. 
        ! 
        ! In the bdytides.F90 the function bdy_dta_tides calls
        ! tide_init_elevation and tide_init_velocities (also in bdytides.F90) 
        ! every day. This uses the values from astro angles to update the phase of the 
        ! tidal constiuents as read in from the boundary files. 
        ! 
        ! The tidal potential in (re)initialised every day in sbstide.F90 in the function
        ! tide_init_potential. This uses the values from astro angles:
        ! (v0tide + utide) and produces amp_pot and phi_pot.
        ! These are then used in SBC/updtide.F90 (every timestep?) to set pot_astro.
        !
        ! Both SBC/sbctide.F90 and bdy_dta_tides calculate zoff+z_arg which is the number of seconds since the beginning of the day.
        ! the tidal phases are then corrected for this reset with the v0tide parameter, calucate by tide_vuf.
        ! nodal correction is much smaller, with ftide (which affects the amplitude), and utide (which affects the phase).
        ! 
        ! 
        ! As the phase of the tidal constituents for both the boundaries and the tidal potential
        ! are adjusted by the astronomic_angle, we can adapt this one module to adapt the tides 
        ! for 360 day calendars.
        ! 
        ! There are different approaches to tides in a 360 day calendar.
        !   1) (current), the tides are effectively reset to the first of the month.
        !       therefore skip 31st's and repeat 29th and 30th of Feb
        !       its happy with extra days of the month (doesn't crash for 30th Feb)
        !       Tide is anchored to correct part of the year, but extra/missing days
        !       are unrealistic, add noise to the system, and make least square tidal analysis fail
        !
        !   2) Start the tides at the begining of the run and then let run continuously.
        !       The tides drift throughout the year, so the equinox's are not at the correct part of the year.
        !       
        !       This is the approach set up below
        !
        !       2b) Adapt the equations to use decimal years (they sort of do, as they use day of year) 
        !           This would make the counting forward and backward from the origin easier 
        !           (the final step (going from DOY to mon and yr) would be removed)
        !
        !   4) Adapt the equations that affect the location of the moon and tides. 
        !       This very likely to be a hugely complex job, that would affect the amphidromic systems,
        !       As you're likely to need to change many/all of the tidal constants. this is then likely
        !       to change the tidal frequencies, and so the the tidal wave speed, and hence the amphidromes, 
        !       co-tide and co-phase lines. 
        !
        !       This approach is not followed
        !
        !
        ! To make the tide continueous for 360 and 365.25 day calendars,
        ! firstly, I make temporary working yr/mn/day integers.
        ! for the gregorian calendar these are simply set to nyear, nmonth and nday.
        ! 
        ! For a 360 day calendar, I count the days from 1900/1/1 to the current day 
        ! according to the the 360 day calendar.
        ! I then count forward that many days according to the gregorian calendar.
        ! therefore every 30 day month of the 360d model run, the tides move forward 30 days.



    ! Questions:
    ! Are the better ways of adding offets to dates in Nemo/Fortran? i.e. python timedelta from datetime
    !       is there a leap year function in NEMO/Fortran?
    ! Not sure if the algorithm is very stable for different origins.
    !   nleap is corrected for 1900 not being a leap year. Needs updated for a different origin year
    !   Does it work if its not starting on a leap year?
    !   Does it work if its after the 28th Feb?
    !   Does it work if the origin is after the start date of the run?
    !   When adjusting the DOY and Y for the number of leap years, what happens if its more than 365 leap years?
    !       
    ! h mean solar Longitude and s mean lunar Longitude and are functions of zhfrac, zday and zsy, 
    !   but the coeffiencents are not 1/86400:1:365
    !   zday is the DOY corrected for the number of leap years since 1900. 
    !   So can run from 20 to 385. this wouldn't matter if the coefficients were 1/86400:1:365
    !   should zsy and zday be updated so zday is between e.g. 1 and 365??
    !
    ! What are the impacts on the NWS if the tide drifts?
    ! What is the impact on the NWS if the tide repeats/skips days?
    !   can this make the model go unstable?
    
   

    ! New variables defined for new code
      INTEGER  ::   yr_org,mn_org,dy_org            !JT
      REAL(wp) ::   sec_grg                         !JT
      INTEGER  ::   yr_grg,mn_grg,dy_grg            !JT
      INTEGER  ::   yr_360,mn_360,dy_360            !JT
      INTEGER  ::   yr_wrk,mn_wrk,dy_wrk            !JT

    

      LOGICAL  ::   ln_tide_drift_time_cont    ! Do we correct for a continueous time


      ! INTEGER(KIND=8)  ::  days_since_origin  ! added to module
      INTEGER  ::  init_yr, day_in_init_yr,nleap,init_doy
      INTEGER  ::  init_doy_inc_l,yg_is_leap_mod,doy_grg
      INTEGER,DIMENSION(12) ::  idayt, idays
      !INTEGER  ::    ji

      !REAL(wp) ::   fjulday_org       !: current julian day 
      ! REAL(wp) ::   days_since_origin_ymds2ju
      INTEGER(KIND=8) ::   days_since_origin_ymds2ju_int



      REAL(wp) ::   jul_org_greg,jul_org_360,jul_pres_360
         
      DATA idayt/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      
      
      ! Currently hardcode the verbosity and the of the code
      ! how to I read the calendar type
      !ln_tide_360_cal = .TRUE.
      !IF ( nleapy == 30)  ln_tide_360_cal = .TRUE.



      !! Nameslist values


      IF (ln_astro_verbose .AND. lwp) THEN
        WRITE(numout,*) 'astro '
        WRITE(numout,*) 'astro -------------------------------------------------'
      ENDIF
      

      ln_tide_drift_time_cont = .FALSE. ! the same the original code

      IF( nleapy == 30 ) THEN
        IF ( ln_tide_drift ) THEN               
            ln_tide_drift_time_cont = .TRUE.
        ENDIF
        IF ( ln_tide_compress ) THEN               
            ! ##################################################################
            ! ##################################################################
            ! ##################################################################
            ! #####
            ! #####  For the 360 day tide constituents, 
            ! #####  We only use days_since_origin for v0tide in tide_vuf.
            ! #####  
            ! #####  To use the correct tide nodal correction (utide)
            ! #####    (which is a small ajustment)
            ! #####  use keep that linked to the gregorian dates.
            ! #####
            ! #####
            ! #####
            ! #####  Therefore, we keep yr_wrk, mn_wrk, dy_wrk to equal
            ! #####    nyear, nmonth, nday
            ! #####
            ! ##################################################################
            ! ##################################################################
            ! ##################################################################

            ln_tide_drift_time_cont = .FALSE.

            ! ##################################################################
            ! ##################################################################
            ! ##################################################################
            ! #####
            ! #####  NEMO4.0.4
            ! #####  BUT!!! need to calc days_since_origin, so need to 
            ! #####  set ln_tide_drift_time_cont too true, then reset 
            ! #####  yr_wrk, mn_wrk, dy_wrk to equal nyear, nmonth, nday
            ! #####
            ! #####
            ! ##################################################################
            ! ##################################################################
            ! ##################################################################
            ln_tide_drift_time_cont = .TRUE.
        ENDIF

      ELSE
        ln_tide_drift_time_cont = .FALSE.
      ENDIF


      IF (ln_astro_verbose .AND. lwp) THEN
        WRITE(numout,*) 'astro ln_tide_drift_time_cont = ',ln_tide_drift_time_cont
      ENDIF

      !IF( ln_tide_360_cal ) THEN
      !IF( nleapy == 30 ) THEN
      IF( ln_tide_drift_time_cont ) THEN


        ! clear and set current dates. 
        yr_360 = nyear
        mn_360 = nmonth
        dy_360 = nday


        yr_grg = 0
        mn_grg = 0
        dy_grg = 0

        yr_wrk = 0
        mn_wrk = 0
        dy_wrk = 0

        ! Set the origin in the name list

        yr_org = nn_tide_orig_yr
        mn_org = nn_tide_orig_mn
        dy_org = nn_tide_orig_dy

        
        !IF (ln_tide_drift_time_cont_manual) THEN



!            IF (ln_astro_verbose .AND. lwp) THEN
!                WRITE(numout,*) 'astro: yr_360,yr_org,((yr_360-yr_org)*360)', yr_360,yr_org,((yr_360-yr_org)*360)
!                WRITE(numout,*) 'astro: mn_360,mn_org,((mn_360-mn_org)*30)', mn_360,mn_org,((mn_360-mn_org)*30)
!                WRITE(numout,*) 'astro: dy_360,dy_org,(dy_360-dy_org)', dy_360,dy_org,(dy_360-dy_org)
!            ENDIF
!            
!            ! how many days from 1900 in the 360 day calendar
!            days_since_origin = ((yr_360-yr_org)*360) + ((mn_360-mn_org)*30) + (dy_360-dy_org)
!
!            ! first guess of what year this would be for the same numbers of days from 1/1/1900 in a gregorian calendar
!            init_yr = yr_org + days_since_origin/365
!
!            ! was the initial estimated year a leap year? how many days in this year?
!            day_in_init_yr = 365
!            if (MOD(init_yr,4) == 0) day_in_init_yr = 366
!
!
!
!            !CALL ymds2ju_JT (yr_org, mn_org, dy_org, 0.0, fjulday_org,360.)
!
!            !IF (ln_astro_verbose) THEN
!            !  IF(lwp) THEN
!            !    WRITE(numout,*) 'astro: ymds2ju_JT yr_org, mn_org, dy_org,fjulday_org', yr_org, mn_org, dy_org,fjulday_org
!            !  ENDIF
!            !ENDIF
!
!
!            !CALL ymds2ju( yr_org, mn_org, dy_org, 0.0, fjulday_org )  ! we assume that we start run at 00:00
!            !IF( ABS(fjulday_org - REAL(NINT(fjulday_org),wp)) < 0.1 / rday )   fjulday_org = REAL(NINT(fjulday_org),wp)   ! avoid truncation error
!            !fjulday_org = fjulday_org + 1.                             ! move back to the day at nit000 (and not at nit000 - 1)
!
!            !days_since_origin_ymds2ju_int = AINT(fjulday - fjulday_org)
!
!            IF (ln_astro_verbose .AND. lwp) THEN             
!                WRITE(numout,*) 'astro: days_since_origin,init_yr,day_in_init_yr', days_since_origin,init_yr,day_in_init_yr
!                !WRITE(numout,*) 'astro: fjulday_org', fjulday_org
!                !WRITE(numout,*) 'astro: fjulday', fjulday
!                !WRITE(numout,*) 'astro: fjulday - fjulday_org', fjulday - fjulday_org
!                !WRITE(numout,*) 'astro: days_since_origin_ymds2ju_int', days_since_origin_ymds2ju_int
!            ENDIF
!
!
!            ! how many leap years since the origin. 
!            nleap = (yr_360-yr_org)/4 - 1 !1900 is not a leap year
!            
!            ! initial estimate of the day of year
!            init_doy = MOD(days_since_origin,365)
!            
!            ! correct the initial estimate for the DOY for the number of leap days since the origin
!            init_doy_inc_l = init_doy - nleap
!
!
!            IF (ln_astro_verbose .AND. lwp) THEN
!                WRITE(numout,*) 'astro: nleap,init_doy,init_doy_inc_l',nleap,init_doy,init_doy_inc_l
!            ENDIF
!            
!
!            ! The number of leap days could pull the  DOY before 0.
!            ! in which case decrement the year, and reset the DOY.
!            ! of the origin is 365 leap years ago, and initial DOY could be adjusted by more than one year..
!            ! Unlikely to be a prob, but need to remember if planning very long control runs. Need to think about this.
!
!            IF (init_doy_inc_l .LT. 0) THEN
!                init_doy_inc_l = init_doy_inc_l+365
!                init_yr = init_yr - 1 
!                IF (MOD(init_yr, 4) == 0 ) THEN
!                    init_doy_inc_l = init_doy_inc_l + 1
!                ENDIF
!            ENDIF
!
!            
!            ! This gives the year and the day of year in the gregorian calendar
!            yr_grg = init_yr    
!            doy_grg = init_doy_inc_l
!            yg_is_leap_mod = MOD(yr_grg, 4)
!
!            IF (ln_astro_verbose .AND. lwp) THEN
!                WRITE(numout,*) 'astro: yr_grg,doy_grg,yg_is_leap_mod',yr_grg,doy_grg,yg_is_leap_mod
!            ENDIF
!
!
!            ! Convert from day of year to month and day in the gregorian calendar.
!            !   dayjul code adapted
!            !   this perhaps should be a function, but not sure how to write one
!            !   there may be this code functionality elsewhere in NEMO
!            !!----------------------------------------------------------------------
!            
!
!            ! what is the DOY of the first day of the month for each month.
!            !   correct for leap years.
!            
!            idays(1) = 0.
!            idays(2) = 31.
!            inc = 0.
!            IF( yg_is_leap_mod == 0.)   inc = 1.
!
!            DO ji = 3, 12
!                idays(ji)=idayt(ji)+inc
!            END DO
!        
!            ! cycle through the months.
!            !   if the DOY is greater than the DOY of the first Day of Month
!            !       Note the month. Calculate day of month by subtraction.
!            !   Once beyond the correct month, the if statement won't be true, so wont calculate.
!
!            DO ji = 1, 12
!                IF ( doy_grg .GE. idays(ji) )  THEN
!                    mn_grg = ji
!                    dy_grg = doy_grg-idays(ji) +1
!                ENDIF
!            END DO
!
!
!
!
!
!            IF(ln_astro_verbose .AND. lwp) THEN
!                WRITE(numout,*) 'astro: mn_grg,dy_grg',mn_grg,dy_grg
!                WRITE(numout,*) ' '
!                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : yr_360,mn_360,dy_360,yr_grg,mn_grg,dy_grg,doy_grg =',yr_360,mn_360,dy_360,yr_grg,mn_grg,dy_grg,doy_grg
!
!                WRITE(numout,*) ' '
!            ENDIF
!            
!
!            
!            IF (ln_astro_verbose .AND. lwp)  WRITE(numout,*) 'tide_mod_astro_ang_meth_1,',yr_grg, mn_grg, dy_grg


        !ELSE ! ln_tide_drift_time_cont_manual
            
            
            ! number of days since 15th October 1582, for namelist origin, in both calendars, and for current model day.
            
            CALL ymds2ju_JT( yr_org,mn_org,dy_org, 0. ,jul_org_greg,365.24 )
            CALL ymds2ju_JT( yr_org,mn_org,dy_org, 0. ,jul_org_360,360. )
            CALL ymds2ju_JT( yr_360,mn_360,dy_360, 0. ,jul_pres_360,360. )

            ! Calculate the days since the origin: days_since_origin_ymds2ju_int
            ! How many days between the current day, and the origin, in the 360 day calendar.
            days_since_origin_ymds2ju_int = jul_pres_360 - jul_org_360

            IF (ln_astro_verbose .AND. lwp) THEN
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : jul_org_360,jul_pres_360,jul_pres_360 - jul_org_360 =',jul_org_360,jul_pres_360,jul_pres_360 - jul_org_360
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : days_since_origin_ymds2ju_int, days_since_origin_ymds2ju_int mod 360 =',days_since_origin_ymds2ju_int,MOD( days_since_origin_ymds2ju_int ,360 )
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : yr_org,mn_org,dy_org, jul_org_greg =',yr_org,mn_org,dy_org, jul_org_greg
            ENDIF

            !add days_since_origin_ymds2ju_int days to the origin in the gregorian calendar.
            CALL ju2ymds_JT( days_since_origin_ymds2ju_int + jul_org_greg, yr_grg, mn_grg, dy_grg, sec_grg,365.24 )

            IF (ln_astro_verbose .AND. lwp) THEN
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : yr_grg, mn_grg, dy_grg =',yr_grg, mn_grg, dy_grg
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : yr_360, mn_360, dy_360 =',yr_360, mn_360, dy_360
                WRITE(numout,*) 'tide_mod_astro_ang 360_corr : yr_org, mn_org, dy_org =',yr_org, mn_org, dy_org
            ENDIF



            
            IF (ln_astro_verbose .AND. lwp)  WRITE(numout,*) 'tide_mod_astro_ang_meth_2,',yr_grg, mn_grg, dy_grg

        !ENDIF !ln_tide_drift_time_cont_manual

        ! for 360 calendars, work with the pseudo gregorian dates
        yr_wrk = yr_grg
        mn_wrk = mn_grg
        dy_wrk = dy_grg

        days_since_origin = days_since_origin_ymds2ju_int

        
        !IF (ln_tide_compress) THEN        
        !    yr_wrk = nyear
        !    mn_wrk = nmonth
        !    dy_wrk = nday
        !ENDIF

      ELSE

        ! for gregorian calendars, work with the model gregorian dates
        yr_wrk = nyear
        mn_wrk = nmonth
        dy_wrk = nday

      ENDIF

      ! continue with original code, using working year, month and day.

      !
      zqy = AINT( (yr_wrk-1901.)/4. )        ! leap years since 1901
      zsy = yr_wrk - 1900.                   ! years since 1900
      !
      zdj  = dayjul( yr_wrk, mn_wrk, dy_wrk )  ! day number of year
      zday = zdj + zqy - 1.                 ! day number of year + No of leap yrs
                                            ! i.e. what would doy if every year = 365 day??
      !
      zhfrac = nsec_day / 3600.             ! The seconds of the day/3600

       
      !
      !----------------------------------------------------------------------
      !  Sh_n Longitude of ascending lunar node
      !----------------------------------------------------------------------
      sh_N=(259.1560564-19.328185764*zsy-.0529539336*zday-.0022064139*zhfrac)*rad
      !----------------------------------------------------------------------
      ! T mean solar angle (Greenwhich time)
      !----------------------------------------------------------------------
      sh_T=(180.+zhfrac*(360./24.))*rad
      !----------------------------------------------------------------------
      ! h mean solar Longitude
      !----------------------------------------------------------------------
      sh_h=(280.1895014-.238724988*zsy+.9856473288*zday+.0410686387*zhfrac)*rad
      !----------------------------------------------------------------------
      ! s mean lunar Longitude
      !----------------------------------------------------------------------
      sh_s=(277.0256206+129.38482032*zsy+13.176396768*zday+.549016532*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p1 Longitude of solar perigee
      !----------------------------------------------------------------------
      sh_p1=(281.2208569+.01717836*zsy+.000047064*zday+.000001961*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p Longitude of lunar perigee
      !----------------------------------------------------------------------
      sh_p=(334.3837214+40.66246584*zsy+.111404016*zday+.004641834*zhfrac)*rad



      IF(ln_astro_verbose .AND. lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'tide_mod_astro_ang : yr_wrk,mn_wrk,dy_wrk=',yr_wrk,mn_wrk,dy_wrk
          WRITE(numout,*) 'tide_mod_astro_ang : nyear, nmonth, nday,nsec_day=',nyear, nmonth, nday,nsec_day
          WRITE(numout,*) 'tide_mod_astro_ang : sh_N,sh_T,sh_h,sh_s,sh_p1,sh_p=', sh_N,sh_T,sh_h,sh_s,sh_p1,sh_p
          WRITE(numout,*) 'tide_mod_astro_ang : zsy,zday,zhfrac,rad=', zsy,zday,zhfrac,rad
          WRITE(numout,*) 'tide_mod_astro_ang : zqy ,zdj,yr_wrk, mn_wrk, dy_wrk =', zqy ,zdj,yr_wrk, mn_wrk, dy_wrk
          WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF



      sh_N = MOD( sh_N ,2*rpi )
      sh_s = MOD( sh_s ,2*rpi )
      sh_h = MOD( sh_h, 2*rpi )
      sh_p = MOD( sh_p, 2*rpi )
      sh_p1= MOD( sh_p1,2*rpi )

      cosI = 0.913694997 -0.035692561 *cos(sh_N)

      sh_I = ACOS( cosI )

      sin2I   = sin(sh_I)
      sh_tgn2 = tan(sh_N/2.0)

      at1=atan(1.01883*sh_tgn2)
      at2=atan(0.64412*sh_tgn2)

      sh_xi=-at1-at2+sh_N

      IF( sh_N > rpi )   sh_xi=sh_xi-2.0*rpi

      sh_nu = at1 - at2

      !----------------------------------------------------------------------
      ! For constituents l2 k1 k2
      !----------------------------------------------------------------------

      tgI2 = tan(sh_I/2.0)
      P1   = sh_p-sh_xi

      t2 = tgI2*tgI2
      t4 = t2*t2
      sh_x1ra = sqrt( 1.0-12.0*t2*cos(2.0*P1)+36.0*t4 )

      p = sin(2.0*P1)
      q = 1.0/(6.0*t2)-cos(2.0*P1)
      sh_R = atan(p/q)

      p = sin(2.0*sh_I)*sin(sh_nu)
      q = sin(2.0*sh_I)*cos(sh_nu)+0.3347
      sh_nuprim = atan(p/q)

      s2 = sin(sh_I)*sin(sh_I)
      p  = s2*sin(2.0*sh_nu)
      q  = s2*cos(2.0*sh_nu)+0.0727
      sh_nusec = 0.5*atan(p/q)
      !
   END SUBROUTINE astronomic_angle






   SUBROUTINE astronomic_angle_origin
      !!----------------------------------------------------------------------
      !!  tj is time elapsed since 1st January 1900, 0 hour, counted in julian
      !!  century (e.g. time in days divide by 36525)
      !!----------------------------------------------------------------------
      REAL(wp) ::   cosI, p, q, t2, t4, sin2I, s2, tgI2, P1, sh_tgn2, at1, at2
      REAL(wp) ::   zqy , zsy, zday, zdj, zhfrac

      
    
   

    ! New variables defined for new code
      INTEGER  ::   yr_wrk,mn_wrk,dy_wrk            !JT

    

      ! for gregorian calendars, work with the model gregorian dates
      yr_wrk = nn_tide_orig_yr
      mn_wrk = nn_tide_orig_mn
      dy_wrk = nn_tide_orig_dy


      !
      zqy = AINT( (yr_wrk-1901.)/4. )        ! leap years since 1901
      zsy = yr_wrk - 1900.                   ! years since 1900
      !
      zdj  = dayjul( yr_wrk, mn_wrk, dy_wrk )  ! day number of year
      zday = zdj + zqy - 1.                 ! day number of year + No of leap yrs
                                            ! i.e. what would doy if every year = 365 day??
      !
      zhfrac = nsec_day / 3600.             ! The seconds of the day/3600

       
      !
      !----------------------------------------------------------------------
      !  Sh_n Longitude of ascending lunar node
      !----------------------------------------------------------------------
      sh_N_o=(259.1560564-19.328185764*zsy-.0529539336*zday-.0022064139*zhfrac)*rad
      !----------------------------------------------------------------------
      ! T mean solar angle (Greenwhich time)
      !----------------------------------------------------------------------
      sh_T_o=(180.+zhfrac*(360./24.))*rad
      !----------------------------------------------------------------------
      ! h mean solar Longitude
      !----------------------------------------------------------------------
      sh_h_o=(280.1895014-.238724988*zsy+.9856473288*zday+.0410686387*zhfrac)*rad
      !----------------------------------------------------------------------
      ! s mean lunar Longitude
      !----------------------------------------------------------------------
      sh_s_o=(277.0256206+129.38482032*zsy+13.176396768*zday+.549016532*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p1 Longitude of solar perigee
      !----------------------------------------------------------------------
      sh_p1_o=(281.2208569+.01717836*zsy+.000047064*zday+.000001961*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p Longitude of lunar perigee
      !----------------------------------------------------------------------
      sh_p_o=(334.3837214+40.66246584*zsy+.111404016*zday+.004641834*zhfrac)*rad



      IF(ln_astro_verbose .AND. lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'tide_mod_astro_ang_orig,yr_wrk,mn_wrk,dy_wrk,nsec_day,=',yr_wrk,mn_wrk,dy_wrk,nsec_day
          WRITE(numout,*) 'tide_mod_astro_ang_orig,sh_N_o,sh_T_o,sh_h_o,sh_s_o,sh_p1_o,sh_p_o,', sh_N_o,sh_T_o,sh_h_o,sh_s_o,sh_p1_o,sh_p_o
          WRITE(numout,*) 'tide_mod_astro_ang_orig,zqy ,zdj,zsy,zday,zhfrac,rad,', zqy ,zdj,zsy,zday,zhfrac,rad

          WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF



      sh_N_o = MOD( sh_N_o ,2*rpi )
      sh_s_o = MOD( sh_s_o ,2*rpi )
      sh_h_o = MOD( sh_h_o, 2*rpi )
      sh_p_o = MOD( sh_p_o, 2*rpi )
      sh_p1_o= MOD( sh_p1_o,2*rpi )


!      cosI = 0.913694997 -0.035692561 *cos(sh_N_o)
!
!
!
!!      REAL(wp) ::   cosI, p, q, t2, t4, sin2I, s2, tgI2, P1, sh_tgn2, at1, at2
!!      REAL(wp) ::   zqy , zsy, zday, zdj, zhfrac
!
!
!      sh_I_o = ACOS( cosI )
!
!      sin2I   = sin(sh_I_o)
!      sh_tgn2 = tan(sh_N_o/2.0)
!
!      at1=atan(1.01883*sh_tgn2)
!      at2=atan(0.64412*sh_tgn2)
!
!      sh_xi_o=-at1-at2+sh_N
!
!      IF( sh_N_o > rpi )   sh_xi_o=sh_xi_o-2.0*rpi
!
!      sh_nu_o = at1 - at2
!
!      !----------------------------------------------------------------------
!      ! For constituents l2 k1 k2
!      !----------------------------------------------------------------------
!
!      tgI2 = tan(sh_I_o/2.0)
!      P1   = sh_p_o-sh_xi_o
!
!      t2 = tgI2*tgI2
!      t4 = t2*t2
!      sh_x1ra_o = sqrt( 1.0-12.0*t2*cos(2.0*P1)+36.0*t4 )
!
!      p = sin(2.0*P1)
!      q = 1.0/(6.0*t2)-cos(2.0*P1)
!      sh_R = atan(p/q)
!
!      p = sin(2.0*sh_I)*sin(sh_nu)
!      q = sin(2.0*sh_I)*cos(sh_nu)+0.3347
!      sh_nuprim_o = atan(p/q)
!
!      s2 = sin(sh_I_o)*sin(sh_I_o)
!      p  = s2*sin(2.0*sh_nu_o)
!      q  = s2*cos(2.0*sh_nu_o)+0.0727
!      sh_nusec_o = 0.5*atan(p/q)



      !
   END SUBROUTINE astronomic_angle_origin















   SUBROUTINE tide_pulse( pomega, ktide ,kc )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_pulse  ***
      !!                      
      !! ** Purpose : Compute tidal frequencies
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in ) ::   kc       ! Total number of tidal constituents
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide    ! Indice of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pomega   ! pulsation in radians/s
      !
      INTEGER  ::   jh
      REAL(wp) ::   zscale
      REAL(wp) ::   zomega_T =  13149000.0_wp                ! Mean Solar Day  ! degrees/century
      REAL(wp) ::   zomega_s =    481267.892_wp              ! Sidereal Month  ! degrees/century
      REAL(wp) ::   zomega_h  !=     36000.76892_wp           ! Tropical Year   ! degrees/century
      REAL(wp) ::   zomega_p =      4069.0322056_wp          ! Moons Perigee   ! degrees/century
      REAL(wp) ::   zomega_n =      1934.1423972_wp          ! Regression of Lunar Nodes  ! degrees/century
      REAL(wp) ::   zomega_p1=         1.719175_wp           ! Perihelion      ! degrees/century
      !!----------------------------------------------------------------------

      zomega_h = 36000.76892_wp                  ! Tropical Year   ! degrees/century
      IF (( nleapy == 30 ) .AND. ln_tide_compress)      zomega_h = 36525.0_wp   ! 360 day Tropical Year   ! degrees/century (360 deg/360 days= 1deg/day. cent = 36525

      !
      zscale =  rad / ( 36525._wp * 86400._wp ) ! Convert to radians per second. 
      !
      DO jh = 1, kc
         pomega(jh) = (  zomega_T * Wave( ktide(jh) )%nT   &
            &          + zomega_s * Wave( ktide(jh) )%ns   &
            &          + zomega_h * Wave( ktide(jh) )%nh   &
            &          + zomega_p * Wave( ktide(jh) )%np   &
            &          + zomega_p1* Wave( ktide(jh) )%np1  ) * zscale
      END DO
    
      IF (ln_astro_verbose .AND. lwp) THEN

          WRITE(numout,*) 'astro tide_pulse nleapy:',nleapy
          WRITE(numout,*) 'astro tide_pulse zomega_h:',zomega_h
          WRITE(numout,*) 'astro tide_pulse if zomega_h = 36000.76892 for 365.24 day year'
          WRITE(numout,*) 'astro tide_pulse if zomega_h = 36525.00000 for 360.00 day year'

 
          DO jh = 1, kc
              WRITE(numout,*) 'astro tide_pulse const, pomega, period(hr):',Wave(ktide(jh))%cname_tide, pomega(jh),2*rpi/(3600.0_wp*pomega(jh))
          END DO

      ENDIF
      !
   END SUBROUTINE tide_pulse


   SUBROUTINE tide_vuf( pvt, put, pcor, ktide ,kc )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_vuf  ***
      !!                      
      !! ** Purpose : Compute nodal modulation corrections
      !!
      !! ** Outputs : vt: Phase of tidal potential relative to Greenwich (radians)
      !!              ut: Phase correction u due to nodal motion (radians)
      !!              ft: Nodal correction factor
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in ) ::   kc               ! Total number of tidal constituents
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide            ! Indice of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pvt, put, pcor   !

      !
      INTEGER ::   jh   ! dummy loop index
      !!----------------------------------------------------------------------

      !! JT for compress
      REAL(wp)                             :: hours_since_origin
      REAL(wp), DIMENSION(kc) ::   pomega           ! pulsation in radians/s
      REAL(wp), DIMENSION(kc) ::   freq_per_day,   v0linearslope,v0linearintercept             ! pulsation in radians/s  !offset,cycle_reset,freq,per_hr
      !! JT for compress



      IF ( ln_tide_compress) THEN

        CALL tide_pulse( pomega, ktide ,kc )

        DO jh = 1, kc
          !per_hr(jh) = (2*rpi/pomega(jh))/3600.
          !freq(jh) = (2*rpi/per_hr(jh))
          !freq_per_day(jh) = freq(jh)*24
          freq_per_day(jh) = pomega(jh) * 86400.0_wp
          !cycle_reset(jh) = mod(hours_since_origin*freq(jh),2.*rpi)
          v0linearslope(jh) =   - mod (   (-freq_per_day(jh)), (2*rpi)  )
          IF(ln_astro_verbose .AND. lwp) WRITE(numout,*) 'astro tide_vuf 1:',jh,kc,ktide(jh),v0linearslope(jh),freq_per_day(jh), pomega(jh),(2*rpi/pomega(jh))/3600.! * 86400.0_wp,freq(jh)*24,per_hr(jh),freq(jh)
        ENDDO


!        !offset(1) = 0.10789890_wp
!        !offset(2) = 1.10897897_wp
!        !offset(3) = 2.11005903_wp
!        !offset(4) = 0.00000000_wp
!        !offset(5) = 3.47632710_wp
!        !offset(6) = 0.16751976_wp
!        !offset(7) = -0.05503165_wp
!        !offset(8) = 0.94604842_wp
!        !offset(9) = 6.10534877_wp
!        !offset(10) = 0.21579780_wp
!        !offset(11) = 0.00000000_wp
!        !offset(12) = 0.00000000_wp
!        !offset(13) = 0.00000000_wp
!        !offset(14) = 0.00000000_wp
!        !offset(15) = 3.14159265_wp
!        !offset(16) = 0.21833313_wp
!        !offset(17) = 5.50043837_wp
!        !offset(18) = 2.24841149_wp
!        !offset(19) = 0.01800173_wp

!        !v0linearintercept(1) = 0.11044027_wp
!        !v0linearintercept(2) = 1.11152799_wp
!        !v0linearintercept(3) = 2.11261570_wp
!        !v0linearintercept(4) = 0.00000000_wp
!        !v0linearintercept(5) = 3.49727335_wp
!        !v0linearintercept(6) = 0.17784035_wp
!        !v0linearintercept(7) = 6.21578523_wp
!        !v0linearintercept(8) = 0.93368764_wp
!        !v0linearintercept(9) = 6.10534496_wp
!        !v0linearintercept(10) = 0.22088055_wp
!        !v0linearintercept(11) = 0.00000000_wp
!        !v0linearintercept(12) = 0.00000000_wp
!        !v0linearintercept(13) = 0.00000000_wp
!        !v0linearintercept(14) = 0.00000000_wp
!        !v0linearintercept(15) = 3.14159265_wp

!        !v0linearintercept(1) = v0linearintercept(1) - 0.000000_wp
!        !v0linearintercept(2) = v0linearintercept(2) - 0.000000_wp
!        !v0linearintercept(3) = v0linearintercept(3) - 0_wp
!        !v0linearintercept(4) = v0linearintercept(4) - 0.165795_wp
!        !v0linearintercept(5) = v0linearintercept(5) + 2.821252_wp
!        !v0linearintercept(6) = v0linearintercept(6) + 0.479504_wp
!        !v0linearintercept(7) = v0linearintercept(7) - 2.175621_wp
!        !v0linearintercept(8) = v0linearintercept(8) + 1.900267_wp
!        !v0linearintercept(9) = v0linearintercept(9) + 0.107633_wp
!        !v0linearintercept(10) = v0linearintercept(10) - 0.000000_wp
!        !v0linearintercept(11) = v0linearintercept(11) - 0.000000_wp
!        !v0linearintercept(12) = v0linearintercept(12) - 0.225730_wp
!        !v0linearintercept(13) = v0linearintercept(13) - 0.238641_wp
!        !v0linearintercept(14) = v0linearintercept(14) - 3.005851_wp
!        !v0linearintercept(15) = v0linearintercept(15) - 0.000000_wp

!        !v0linearintercept(1) =   0.11044026999999999_wp
!        !v0linearintercept(2) =   1.11152798999999990_wp
!        !v0linearintercept(3) =   2.11261570000000010_wp
!        !v0linearintercept(4) =  -0.16579500000000000_wp
!        !v0linearintercept(5) =   6.31852534999999980_wp
!        !v0linearintercept(6) =   0.65734435000000002_wp
!        !v0linearintercept(7) =   4.04016423000000020_wp
!        !v0linearintercept(8) =   2.83395464000000000_wp
!        !v0linearintercept(9) =   6.21297795999999990_wp
!        !v0linearintercept(10) =  0.22088055000000001_wp
!        !v0linearintercept(11) =  0.00000000000000000_wp
!        !v0linearintercept(12) = -0.22572999999999999_wp
!        !v0linearintercept(13) = -0.23864099999999999_wp
!        !v0linearintercept(14) = -3.00585099999999980_wp
!        !v0linearintercept(15) =  3.14159265000000020_wp

!        v0linearintercept( 1) =   0.2208805500_wp   -  (rpi* 68.0_wp/180.0_wp) !   M2  1
!        v0linearintercept( 2) =   3.1186126191_wp  !   N2  2
!        v0linearintercept( 3) =   0.9305155436_wp  !  2N2  3
!        v0linearintercept( 4) =   0.0194858941_wp  !   S2  4
!        v0linearintercept( 5) =  -2.5213114949_wp  !   K2  5
!        v0linearintercept( 6) =   6.5970532125_wp  !   K1  6
!        v0linearintercept( 7) =   1.1115279900_wp  !   O1  7
!        v0linearintercept( 8) =   0.1104402700_wp  !   Q1  8
!        !     v0linearintercept( 9) =   4.2269096542_wp  !   P1  9
!        !v0linearintercept( 9) =  -2.0351042402_wp  !   P1  9  compress3
!        !v0linearintercept( 9) =  -2.0351042402_wp  - 2.6179938779914944 !   P1  9  compress4

!        v0linearintercept( 9) =   rpi* 345.0_wp/180.0_wp -  (rpi* 140.0_wp/180.0_wp) !   P1  9  compress4

!        v0linearintercept(10) =   3.1415926500_wp  !   M4 10
!        v0linearintercept(11) =   0.0000000000_wp  !   Mf 11
!        v0linearintercept(12) =   0.0000000000_wp  !   Mm 12
!        v0linearintercept(13) =   0.0000000000_wp  ! Msqm 13
!        v0linearintercept(14) =   0.0000000000_wp  !  Mtm 14
!        v0linearintercept(15) =  -0.0230244122_wp  !   S1 15
!        v0linearintercept(16) =   4.2565208698_wp  !  MU2 16
!        v0linearintercept(17) =   6.5001767059_wp  !  NU2 17
!        v0linearintercept(18) =   0.0000000000_wp    -  (rpi* 113.0_wp/180.0_wp) !   L2 18
!        v0linearintercept(19) =   0.0092971808_wp  !   T2 19  + rpi/2.

!        !v0linearintercept(1) = v0linearintercept(1) - 0.034975_wp    ! M2
!        !v0linearintercept(2) = v0linearintercept(2) - 0.030244_wp    ! N2
!        !v0linearintercept(3) = v0linearintercept(3) - 0.036046_wp    ! 2N2
!        !v0linearintercept(4) = v0linearintercept(4) + 0.002092_wp    ! S2
!        !v0linearintercept(5) = v0linearintercept(5) - 0.273826_wp    ! K2
!        !v0linearintercept(6) = v0linearintercept(6) - 0.144677_wp    ! K1
!        !v0linearintercept(7) = v0linearintercept(7) + 0.031938_wp    ! O1
!        !v0linearintercept(8) = v0linearintercept(8) - 0.812030_wp    ! Q1
!        !v0linearintercept(9) = v0linearintercept(9) + 2.109118_wp    ! P1
!        !v0linearintercept(10) = v0linearintercept(10) + 0.070021_wp    ! M4
!        !v0linearintercept(11) = v0linearintercept(11) - 0.000000_wp    ! Mf
!        !v0linearintercept(12) = v0linearintercept(12) - 0.000000_wp    ! Mm
!        !v0linearintercept(13) = v0linearintercept(13) - 0.000000_wp    ! Msqm
!        !v0linearintercept(14) = v0linearintercept(14) - 0.000000_wp    ! Mtm
!        !v0linearintercept(15) = v0linearintercept(15) - 0.035676_wp    ! S1
!        !v0linearintercept(16) = v0linearintercept(16) + 0.007598_wp    ! MU2
!        !v0linearintercept(17) = v0linearintercept(17) - 0.043060_wp    ! NU2
!        !v0linearintercept(18) = v0linearintercept(18) + 0.023561_wp    ! L2
!        !v0linearintercept(19) = v0linearintercept(19) + 0.025624_wp    ! T2

!        v0linearintercept(1) = v0linearintercept(1) - (rpi*2.003909_wp/180.0_wp)    ! M2
!        v0linearintercept(2) = v0linearintercept(2) - (rpi*1.732874_wp/180.0_wp)    ! N2
!        v0linearintercept(3) = v0linearintercept(3) - (rpi*2.065265_wp/180.0_wp)    ! 2N2
!        v0linearintercept(4) = v0linearintercept(4) + (rpi*0.119842_wp/180.0_wp)    ! S2
!        v0linearintercept(5) = v0linearintercept(5) - (rpi*15.689068_wp/180.0_wp)    ! K2
!        v0linearintercept(6) = v0linearintercept(6) - (rpi*8.289390_wp/180.0_wp)    ! K1
!        v0linearintercept(7) = v0linearintercept(7) + (rpi*1.829931_wp/180.0_wp)    ! O1
!        v0linearintercept(8) = v0linearintercept(8) - (rpi*46.525902_wp/180.0_wp)    ! Q1
!        v0linearintercept(9) = v0linearintercept(9) + (rpi*120.843575_wp/180.0_wp)    ! P1
!        v0linearintercept(10) = v0linearintercept(10) + (rpi*4.011896_wp/180.0_wp)    ! M4
!        v0linearintercept(11) = v0linearintercept(11) - (rpi*0.000000_wp/180.0_wp)    ! Mf
!        v0linearintercept(12) = v0linearintercept(12) - (rpi*0.000000_wp/180.0_wp)    ! Mm
!        v0linearintercept(13) = v0linearintercept(13) - (rpi*0.000000_wp/180.0_wp)    ! Msqm
!        v0linearintercept(14) = v0linearintercept(14) - (rpi*0.000000_wp/180.0_wp)    ! Mtm
!        v0linearintercept(15) = v0linearintercept(15) - (rpi*2.044069_wp/180.0_wp)    ! S1
!        v0linearintercept(16) = v0linearintercept(16) + (rpi*0.435315_wp/180.0_wp)    ! MU2
!        v0linearintercept(17) = v0linearintercept(17) - (rpi*2.467160_wp/180.0_wp)    ! NU2
!        v0linearintercept(18) = v0linearintercept(18) + (rpi*1.349939_wp/180.0_wp)    ! L2
!        v0linearintercept(19) = v0linearintercept(19) + (rpi*1.468170_wp/180.0_wp)    ! T2


!        ! wave data.

!        !Wave( 1) = tide(  'M2'     , 0.242297 ,    2   ,  2 , -2 ,  2 ,  0 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
!        !Wave( 2) = tide(  'N2'     , 0.046313 ,    2   ,  2 , -3 ,  2 ,  1 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
!        !Wave( 3) = tide( '2N2'     , 0.006184 ,    2   ,  2 , -4 ,  2 ,  2 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
!        !Wave( 4) = tide(  'S2'     , 0.113572 ,    2   ,  2 ,  0 ,  0 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,     0   )
!        !Wave( 5) = tide(  'K2'     , 0.030875 ,    2   ,  2 ,  0 ,  2 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   , -2   , 0 ,   235   )
!        !!              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
!        !Wave( 6) = tide(  'K1'     , 0.142408 ,    1   ,  1 ,  0 ,  1 ,  0 ,  0  ,  -90  ,  0   ,  0   , -1   ,  0   , 0 ,   227   )
!        !Wave( 7) = tide(  'O1'     , 0.101266 ,    1   ,  1 , -2 ,  1 ,  0 ,  0  ,  +90  ,  2   , -1   ,  0   ,  0   , 0 ,    75   )
!        !Wave( 8) = tide(  'Q1'     , 0.019387 ,    1   ,  1 , -3 ,  1 ,  1 ,  0  ,  +90  ,  2   , -1   ,  0   ,  0   , 0 ,    75   )
!        !Wave( 9) = tide(  'P1'     , 0.047129 ,    1   ,  1 ,  0 , -1 ,  0 ,  0  ,  +90  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )
!        !!              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
!        !Wave(10) = tide(  'M4'     , 0.000000 ,    4   ,  4 , -4 ,  4 ,  0 ,  0  ,    0  ,  4   , -4   ,  0   ,  0   , 0 ,    1    )
!        !!              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
!        !Wave(11) = tide(  'Mf'     , 0.042017 ,    0   ,  0 ,  2 ,  0 ,  0 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
!        !Wave(12) = tide(  'Mm'     , 0.022191 ,    0   ,  0 ,  1 ,  0 , -1 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,   73    )
!        !Wave(13) = tide(  'Msqm'   , 0.000667 ,    0   ,  0 ,  4 , -2 ,  0 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
!        !Wave(14) = tide(  'Mtm'    , 0.008049 ,    0   ,  0 ,  3 ,  0 , -1 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
!        !!              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
!        !Wave(15) = tide(  'S1'     , 0.000000 ,    1   ,  1 ,  0 ,  0 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )   
!        !Wave(16) = tide(  'MU2'    , 0.005841 ,    2   ,  2 , -4 ,  4 ,  0 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,   78    )
!        !Wave(17) = tide(  'NU2'    , 0.009094 ,    2   ,  2 , -3 ,  4 , -1 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,   78    ) 
!        !Wave(18) = tide(  'L2'     , 0.006694 ,    2   ,  2 , -1 ,  2 , -1 ,  0  , +180  ,  2   , -2   ,  0   ,  0   , 0 ,  215    )
!        !Wave(19) = tide(  'T2'     , 0.006614 ,    2   ,  2 ,  0 , -1 ,  0 ,  1  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )

!        !name list
!        !  clname(1)='Q1'
!        !  clname(2)='O1'
!        !  clname(3)='P1'
!        !  clname(4)='S1'
!        !  clname(5)='K1'
!        !  clname(6)='2N2'
!        !  clname(7)='MU2'
!        !  clname(8)='N2'
!        !  clname(9)='NU2'
!        !  clname(10)='M2'
!        !  clname(11)='L2'
!        !  clname(12)='T2'
!        !  clname(13)='S2'
!        !  clname(14)='K2'
!        !  clname(15)='M4'

!        ! ktide 8,7,9,15

!        !ktide = 
!        !8
!        !7
!        !9
!        !15
!        !6
!        !3
!        !16
!        !2
!        !17
!        !1
!        !18
!        !19
!        !4
!        !5
!        !10












!        !NEMO4

!!        clname(1)='Q1'
!!        clname(2)='O1'
!!        clname(3)='P1'
!!        clname(4)='S1'
!!        clname(5)='K1'
!!        clname(6)='2N2'
!!        clname(7)='MU2'
!!        clname(8)='N2'
!!        clname(9)='NU2'
!!        clname(10)='M2'
!!        clname(11)='L2'
!!        clname(12)='T2'
!!        clname(13)='S2'
!!        clname(14)='K2'
!!        clname(15)='M4'
!!        ktide = [10,9,11,12,8,23,21,15,22,14,18,19,16,17,28]


!        v0linearintercept( 1) =   0.1104402700_wp  !   Q1  8
!        v0linearintercept( 2) =   1.1115279900_wp  !   O1  7
!        v0linearintercept( 3) =   rpi* 345.0_wp/180.0_wp -  (rpi* 140.0_wp/180.0_wp) !   P1  9  compress4
!        v0linearintercept( 4) =  -0.0230244122_wp  !   S1 15
!        v0linearintercept( 5) =   6.5970532125_wp  !   K1  6
!        v0linearintercept( 6) =   0.9305155436_wp  !  2N2  3
!        v0linearintercept( 7) =   4.2565208698_wp  !  MU2 16
!        v0linearintercept( 8) =   3.1186126191_wp  !   N2  2
!        v0linearintercept( 9) =   6.5001767059_wp  !  NU2 17
!        v0linearintercept(10) =   0.2208805500_wp   -  (rpi* 68.0_wp/180.0_wp) !   M2  1
!        v0linearintercept(11) =   0.0000000000_wp    -  (rpi* 113.0_wp/180.0_wp) !   L2 18
!        v0linearintercept(12) =   0.0092971808_wp  !   T2 19  + rpi/2.
!        v0linearintercept(13) =   0.0194858941_wp  !   S2  4
!        v0linearintercept(14) =  -2.5213114949_wp  !   K2  5
!        v0linearintercept(15) =   3.1415926500_wp  !   M4 10



!        v0linearintercept( 1) = v0linearintercept( 1) - (rpi*46.525902_wp/180.0_wp)   ! Q1
!        v0linearintercept( 2) = v0linearintercept( 2) + (rpi*1.829931_wp/180.0_wp)    ! O1
!        v0linearintercept( 3) = v0linearintercept( 3) + (rpi*120.843575_wp/180.0_wp)  ! P1
!        v0linearintercept( 4) = v0linearintercept( 4) - (rpi*2.044069_wp/180.0_wp)    ! S1
!        v0linearintercept( 5) = v0linearintercept( 5) - (rpi*8.289390_wp/180.0_wp)    ! K1
!        v0linearintercept( 6) = v0linearintercept( 6) - (rpi*2.065265_wp/180.0_wp)    ! 2N2
!        v0linearintercept( 7) = v0linearintercept( 7) + (rpi*0.435315_wp/180.0_wp)    ! MU2
!        v0linearintercept( 8) = v0linearintercept( 8) - (rpi*1.732874_wp/180.0_wp)    ! N2
!        v0linearintercept( 9) = v0linearintercept( 9) - (rpi*2.467160_wp/180.0_wp)    ! NU2
!        v0linearintercept(10) = v0linearintercept(10) - (rpi*2.003909_wp/180.0_wp)    ! M2
!        v0linearintercept(11) = v0linearintercept(11) + (rpi*1.349939_wp/180.0_wp)    ! L2
!        v0linearintercept(12) = v0linearintercept(12) + (rpi*1.468170_wp/180.0_wp)    ! T2
!        v0linearintercept(13) = v0linearintercept(13) + (rpi*0.119842_wp/180.0_wp)    ! S2
!        v0linearintercept(14) = v0linearintercept(14) - (rpi*15.689068_wp/180.0_wp)   ! K2
!        v0linearintercept(14) = v0linearintercept(15) + (rpi*4.011896_wp/180.0_wp)    ! M4






        DO jh = 1, kc
          IF(ln_astro_verbose .AND. lwp) WRITE(numout,*) 'astro tide_vuf 2:',jh,days_since_origin,v0linearslope(jh),v0linearintercept(ktide(jh))!,cycle_reset(jh)! ,offset(jh) 
        ENDDO
      ENDIF



      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !! JT
      !!!  phi_tide(ib)=phi_tide(ib)+v0tide(itide)+utide(itide)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------

      DO jh = 1, kc
         !  Phase of the tidal potential relative to the Greenwhich 
         !  meridian (e.g. the position of the fictuous celestial body). Units are radian:
         !  Linear with time

         IF ( ln_tide_compress ) THEN

            v0linearintercept(jh) = sh_T_o * Wave( ktide(jh) )%nT    &
                &    + sh_s_o * Wave( ktide(jh) )%ns    &
                &    + sh_h_o * Wave( ktide(jh) )%nh    &
                &    + sh_p_o * Wave( ktide(jh) )%np    &
                &    + sh_p1_o* Wave( ktide(jh) )%np1   &
                &    +          Wave( ktide(jh) )%shift * rad



             !JT pvt(jh) = mod(  ( (v0linearslope(jh)*days_since_origin) + v0linearintercept(  ktide(jh)  ) ),  2*rpi)-(2*rpi)
             pvt(jh) = mod(  ( (v0linearslope(jh)*days_since_origin) + v0linearintercept(  jh  ) ),  2*rpi)-(2*rpi)
         ELSE
             pvt(jh) = sh_T * Wave( ktide(jh) )%nT    &
                &    + sh_s * Wave( ktide(jh) )%ns    &
                &    + sh_h * Wave( ktide(jh) )%nh    &
                &    + sh_p * Wave( ktide(jh) )%np    &
                &    + sh_p1* Wave( ktide(jh) )%np1   &
                &    +        Wave( ktide(jh) )%shift * rad
         ENDIF

         
         
         !
         !  Phase correction u due to nodal motion. Units are radian:
         !  Cyclical with time. Much smaller terms than pvt. 
         put(jh) = sh_xi     * Wave( ktide(jh) )%nksi   &
            &    + sh_nu     * Wave( ktide(jh) )%nnu0   &
            &    + sh_nuprim * Wave( ktide(jh) )%nnu1   &
            &    + sh_nusec  * Wave( ktide(jh) )%nnu2   &
            &    + sh_R      * Wave( ktide(jh) )%R





         !JT
         pvt(jh) = mod( pvt(jh),  2*rpi)
         put(jh) = mod( put(jh),  2*rpi)
         !JT


         !  Nodal correction factor:
         pcor(jh) = nodal_factort( Wave( ktide(jh) )%nformula )


      END DO


      IF(ln_astro_verbose .AND. lwp) THEN
          DO jh = 1, kc
              WRITE(numout,*) 'astro tide_vuf 3,',jh,Wave(ktide(jh))%cname_tide, pomega(jh),2*rpi/(3600.0_wp*pomega(jh)),pvt(jh), put(jh), pcor(jh), mod(pvt(jh),2*rpi), mod(put(jh),2*rpi), 2*rpi
          END DO
      ENDIF



      !
   END SUBROUTINE tide_vuf


   RECURSIVE FUNCTION nodal_factort( kformula ) RESULT( zf )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kformula
      !
      REAL(wp) :: zf
      REAL(wp) :: zs, zf1, zf2
      !!----------------------------------------------------------------------
      !
      SELECT CASE( kformula )
      !
      CASE( 0 )                  !==  formule 0, solar waves
         zf = 1.0
         !
      CASE( 1 )                  !==  formule 1, compound waves (78 x 78)
         zf=nodal_factort(78)
         zf = zf * zf
         !
      CASE ( 2 )                 !==  formule 2, compound waves (78 x 0)  ===  (78) 
       zf1= nodal_factort(78)
       zf = nodal_factort( 0)
       zf = zf1 * zf
       !
      CASE ( 4 )                 !==  formule 4,  compound waves (78 x 235) 
         zf1 = nodal_factort( 78)
         zf  = nodal_factort(235)
         zf  = zf1 * zf
         !
      CASE ( 5 )                 !==  formule 5,  compound waves (78 *78 x 235)
         zf1 = nodal_factort( 78)
         zf  = nodal_factort(235)
         zf  = zf * zf1 * zf1
         !
      CASE ( 6 )                 !==  formule 6,  compound waves (78 *78 x 0)
         zf1 = nodal_factort(78)
         zf  = nodal_factort( 0)
         zf  = zf * zf1 * zf1 
         !
      CASE( 7 )                  !==  formule 7, compound waves (75 x 75)
         zf = nodal_factort(75)
         zf = zf * zf
         !
      CASE( 8 )                  !==  formule 8,  compound waves (78 x 0 x 235)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(  0)
         zf2 = nodal_factort(235)
         zf  = zf * zf1 * zf2
         !
      CASE( 9 )                  !==  formule 9,  compound waves (78 x 0 x 227)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(  0)
         zf2 = nodal_factort(227)
         zf  = zf * zf1 * zf2
         !
      CASE( 10 )                 !==  formule 10,  compound waves (78 x 227)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(227)
         zf  = zf * zf1
         !
      CASE( 11 )                 !==  formule 11,  compound waves (75 x 0)
!!gm bug???? zf 2 fois !
         zf = nodal_factort(75)
         zf1 = nodal_factort( 0)
         zf = zf * zf1
         !
      CASE( 12 )                 !==  formule 12,  compound waves (78 x 78 x 78 x 0) 
         zf1 = nodal_factort(78)
         zf  = nodal_factort( 0)
         zf  = zf * zf1 * zf1 * zf1
         !
      CASE( 13 )                 !==  formule 13, compound waves (78 x 75)
         zf1 = nodal_factort(78)
         zf  = nodal_factort(75)
         zf  = zf * zf1
         !
      CASE( 14 )                 !==  formule 14, compound waves (235 x 0)  ===  (235)
         zf  = nodal_factort(235)
         zf1 = nodal_factort(  0)
         zf  = zf * zf1
         !
      CASE( 15 )                 !==  formule 15, compound waves (235 x 75) 
         zf  = nodal_factort(235)
         zf1 = nodal_factort( 75)
         zf  = zf * zf1
         !
      CASE( 16 )                 !==  formule 16, compound waves (78 x 0 x 0)  ===  (78)
         zf  = nodal_factort(78)
         zf1 = nodal_factort( 0)
         zf  = zf * zf1 * zf1
         !
      CASE( 17 )                 !==  formule 17,  compound waves (227 x 0) 
         zf1 = nodal_factort(227)
         zf  = nodal_factort(  0)
         zf  = zf * zf1
         !
      CASE( 18 )                 !==  formule 18,  compound waves (78 x 78 x 78 )
         zf1 = nodal_factort(78)
         zf  = zf1 * zf1 * zf1
         !
      CASE( 19 )                 !==  formule 19, compound waves (78 x 0 x 0 x 0)  ===  (78)
!!gm bug2 ==>>>   here identical to formule 16,  a third multiplication by zf1 is missing
         zf  = nodal_factort(78)
         zf1 = nodal_factort( 0)
         zf = zf * zf1 * zf1
         !
         
      !--- davbyr 11/2017
      CASE( 20 )                 !==  formule 20,  compound waves ( 78 x 78 x 78 x 78 )
         zf1 = nodal_factort(78)
         zf  = zf1 * zf1 * zf1 * zf1
      !--- END davbyr
      CASE( 73 )                 !==  formule 73
         zs = sin(sh_I)
         zf = (2./3.-zs*zs)/0.5021
         !
      CASE( 74 )                 !==  formule 74
         zs = sin(sh_I)
         zf = zs * zs / 0.1578
         !
      CASE( 75 )                 !==  formule 75
         zs = cos(sh_I/2)
         zf = sin(sh_I) * zs * zs / 0.3800
         !
      CASE( 76 )                 !==  formule 76
         zf = sin(2*sh_I) / 0.7214
         !
      CASE( 77 )                 !==  formule 77
         zs = sin(sh_I/2)
         zf = sin(sh_I) * zs * zs / 0.0164
         !
      CASE( 78 )                 !==  formule 78
         zs = cos(sh_I/2)
         zf = zs * zs * zs * zs / 0.9154
         !
      CASE( 79 )                 !==  formule 79
         zs = sin(sh_I)
         zf = zs * zs / 0.1565
         !
      CASE( 144 )                !==  formule 144
         zs = sin(sh_I/2)
         zf = ( 1-10*zs*zs+15*zs*zs*zs*zs ) * cos(sh_I/2) / 0.5873
         !
      CASE( 149 )                !==  formule 149
         zs = cos(sh_I/2)
         zf = zs*zs*zs*zs*zs*zs / 0.8758
         !
      CASE( 215 )                !==  formule 215
         zs = cos(sh_I/2)
         zf = zs*zs*zs*zs / 0.9154 * sh_x1ra
         !
      CASE( 227 )                !==  formule 227 
         zs = sin(2*sh_I)
         zf = sqrt( 0.8965*zs*zs+0.6001*zs*cos (sh_nu)+0.1006 )
         !
      CASE ( 235 )               !==  formule 235 
         zs = sin(sh_I)
         zf = sqrt( 19.0444*zs*zs*zs*zs + 2.7702*zs*zs*cos(2*sh_nu) + .0981 )
         !
      END SELECT
      !
   END FUNCTION nodal_factort


   FUNCTION dayjul( kyr, kmonth, kday )
      !!----------------------------------------------------------------------
      !!  *** THIS ROUTINE COMPUTES THE JULIAN DAY (AS A REAL VARIABLE)
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in) ::   kyr, kmonth, kday
      !
      INTEGER,DIMENSION(12) ::  idayt, idays
      INTEGER  ::   inc, ji
      REAL(wp) ::   dayjul, zyq
      !
      DATA idayt/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      !!----------------------------------------------------------------------
      !
      idays(1) = 0.
      idays(2) = 31.
      inc = 0.
      zyq = MOD( kyr-1900. , 4. )
      IF( zyq == 0.)   inc = 1.
      DO ji = 3, 12
         idays(ji)=idayt(ji)+inc
      END DO
      dayjul = idays(kmonth) + kday
      !
   END FUNCTION dayjul






    SUBROUTINE ju2ymds_JT (julian,year,month,day,sec,one_year)
    !---------------------------------------------------------------------
    !  IMPLICIT NONE
    !-
    !  REAL,INTENT(IN) :: julian
    !-
    !  INTEGER,INTENT(OUT) :: year,month,day
    !  REAL,INTENT(OUT)    :: sec
    !-
    !  INTEGER :: julian_day
    !  REAL    :: julian_sec
    !---------------------------------------------------------------------
    !  julian_day = INT(julian)
    !  julian_sec = (julian-julian_day)*one_day
    !-
    !  CALL ju2ymds_internal(julian_day,julian_sec,year,month,day,sec)
    !---------------------
    !END SUBROUTINE ju2ymds
    !-
    !===
    !-
    !SUBROUTINE ju2ymds_internal (julian_day,julian_sec,year,month,day,sec)
    !---------------------------------------------------------------------
    !- This subroutine computes from the julian day the year,
    !- month, day and seconds
    !-
    !- In 1968 in a letter to the editor of Communications of the ACM
    !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    !- and Thomas C. Van Flandern presented such an algorithm.
    !-
    !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
    !-
    !- In the case of the Gregorian calendar we have chosen to use
    !- the Lilian day numbers. This is the day counter which starts
    !- on the 15th October 1582. This is the day at which Pope
    !- Gregory XIII introduced the Gregorian calendar.
    !- Compared to the true Julian calendar, which starts some 7980
    !- years ago, the Lilian days are smaler and are dealt with easily
    !- on 32 bit machines. With the true Julian days you can only the
    !- fraction of the day in the real part to a precision of a 1/4 of
    !- a day with 32 bits.
    !---------------------------------------------------------------------
      IMPLICIT NONE


    !-
      
      REAL,INTENT(IN)    :: julian,one_year
    !-
      INTEGER,INTENT(OUT) :: year,month,day
      REAL,INTENT(OUT)    :: sec
    !-
      INTEGER :: l,n,i,jd,j,d,m,y,ml
      INTEGER :: add_day
      REAL :: eps_day

      REAL,PARAMETER :: one_day = 86400.0
    !---------------------------------------------------------------------

      INTEGER :: julian_day
      REAL    :: julian_sec
      INTEGER :: mon_len(12)
    !---------------------------------------------------------------------
    
    IF  ( (one_year > 365.0).AND.(one_year < 366.0) )  THEN
      mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    ELSE IF ( ABS(one_year-365.0) <= EPSILON(one_year) ) THEN
      mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    ELSE IF ( ABS(one_year-366.0) <= EPSILON(one_year) ) THEN
      mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    ELSE IF ( ABS(one_year-360.0) <= EPSILON(one_year) ) THEN
      mon_len(:)=(/30,30,30,30,30,30,30,30,30,30,30,30/)
    ENDIF



      julian_day = INT(julian)
      julian_sec = (julian-julian_day)*one_day

      eps_day = SPACING(one_day)
    !  lock_one_year = .TRUE.
    !-
      jd = julian_day
      sec = julian_sec
      IF (sec > (one_day-eps_day)) THEN
        add_day = INT(sec/one_day)
        sec = sec-add_day*one_day
        jd = jd+add_day
      ENDIF
      IF (sec < -eps_day) THEN
        sec = sec+one_day
        jd = jd-1
      ENDIF
    !-
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
    !-- Gregorian
        jd = jd+2299160
    !-
        l = jd+68569
        n = (4*l)/146097
        l = l-(146097*n+3)/4
        i = (4000*(l+1))/1461001
        l = l-(1461*i)/4+31
        j = (80*l)/2447
        d = l-(2447*j)/80
        l = j/11
        m = j+2-(12*l)
        y = 100*(n-49)+i+l
      ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year)) &
     &         .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
    !-- No leap or All leap
        !if ( ABS(one_year-365.0) <= EPSILON(one_year) ) mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
        !if ( ABS(one_year-366.0) <= EPSILON(one_year) ) mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
        y = jd/NINT(one_year)
        l = jd-y*NINT(one_year)
        m = 1
        ml = 0
        DO WHILE (ml+mon_len(m) <= l)
          ml = ml+mon_len(m)
          m = m+1
        ENDDO
        d = l-ml+1
      ELSE
    !-- others
        ml = NINT(one_year/12.)
        y = jd/NINT(one_year)
        l = jd-y*NINT(one_year)
        m = (l/ml)+1
        d = l-(m-1)*ml+1
      ENDIF
    !-
      day = d
      month = m
      year = y
    !------------------------------


    END SUBROUTINE ju2ymds_JT








    !-
    SUBROUTINE ymds2ju_JT (year,month,day,sec,julian,one_year)
    !---------------------------------------------------------------------
    !  IMPLICIT NONE
    !-
    !  INTEGER,INTENT(IN) :: year,month,day
    !  REAL,INTENT(IN)    :: sec
    !-
    !  REAL,INTENT(OUT) :: julian
    !-
    !  INTEGER :: julian_day
    !  REAL    :: julian_sec
    !---------------------------------------------------------------------
    !  CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
    !-
    !---------------------
    !END SUBROUTINE ymds2ju
    !-
    !===
    !-
    !SUBROUTINE ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
    !---------------------------------------------------------------------
    !- Converts year, month, day and seconds into a julian day
    !-
    !- In 1968 in a letter to the editor of Communications of the ACM
    !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    !- and Thomas C. Van Flandern presented such an algorithm.
    !-
    !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
    !-
    !- In the case of the Gregorian calendar we have chosen to use
    !- the Lilian day numbers. This is the day counter which starts
    !- on the 15th October 1582.
    !- This is the day at which Pope Gregory XIII introduced the
    !- Gregorian calendar.
    !- Compared to the true Julian calendar, which starts some
    !- 7980 years ago, the Lilian days are smaler and are dealt with
    !- easily on 32 bit machines. With the true Julian days you can only
    !- the fraction of the day in the real part to a precision of
    !- a 1/4 of a day with 32 bits.
    !---------------------------------------------------------------------
      IMPLICIT NONE
    !-
      INTEGER,INTENT(IN) :: year,month,day
      REAL,INTENT(IN)    :: sec
      REAL,INTENT(IN)    :: one_year
    !-
    !  INTEGER,INTENT(OUT) :: julian_day
    !  REAL,INTENT(OUT)    :: julian_sec
      REAL,INTENT(OUT) :: julian
      INTEGER          :: julian_day
      REAL             :: julian_sec
    !-
      INTEGER :: jd,m,y,d,ml
      REAL,PARAMETER :: one_day = 86400.0


          INTEGER :: mon_len(12)  
      
        IF  ( (one_year > 365.0).AND.(one_year < 366.0) )  THEN
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
        ELSE IF ( ABS(one_year-365.0) <= EPSILON(one_year) ) THEN
          mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
        ELSE IF ( ABS(one_year-366.0) <= EPSILON(one_year) ) THEN
          mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
        ELSE IF ( ABS(one_year-360.0) <= EPSILON(one_year) ) THEN
          mon_len(:)=(/30,30,30,30,30,30,30,30,30,30,30,30/)
        ENDIF


    !---------------------------------------------------------------------
     ! lock_one_year = .TRUE.
    !-


      m = month
      y = year
      d = day


        !---------------------------------------------------------------------


    !-
    !- We deduce the calendar from the length of the year as it
    !- is faster than an INDEX on the calendar variable.
    !-
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
    !-- "Gregorian"
        jd = (1461*(y+4800+INT((m-14)/12)))/4 &
     &      +(367*(m-2-12*(INT((m-14)/12))))/12 &
     &      -(3*((y+4900+INT((m-14)/12))/100))/4 &
     &      +d-32075
        jd = jd-2299160
      ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year))  &
     &         .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
    !-- "No leap" or "All leap"
        ml = SUM(mon_len(1:m-1))
        jd = y*NINT(one_year)+ml+(d-1)
      ELSE
    !-- Calendar with regular month
        ml = NINT(one_year/12.)
        jd = y*NINT(one_year)+(m-1)*ml+(d-1)
      ENDIF
    !-
      julian_day = jd
      julian_sec = sec

      julian = julian_day+julian_sec/one_day

    !------------------------------
    END SUBROUTINE ymds2ju_JT








   !!======================================================================
END MODULE tide_mod
