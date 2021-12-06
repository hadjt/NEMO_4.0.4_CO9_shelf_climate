MODULE diaharm_fast 
   !!======================================================================
   !!                       ***  MODULE  example  ***
   !! Ocean physics:  On line harmonic analyser
   !!                 
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!                   :                Calculate harmonic analysis
   !!----------------------------------------------------------------------
   !!   harm_ana        :
   !!   harm_ana_init   :
   !!   NB: 2017-12 : add 3D harmonic analysis of velocities
   !!                 integration of Maria Luneva's development
   !!   'key_3Ddiaharm'
   !!----------------------------------------------------------------------

   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal bdy forcing
   USE daymod          ! calendar
   USE tideini
   USE restart
   USE ioipsl, ONLY : ju2ymds    ! for calendar
   !
   !
   USE timing          ! preformance summary
   USE zdf_oce

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC dia_harm_fast                                      ! routine called in step.F90 module
   LOGICAL, PUBLIC, PARAMETER :: lk_diaharm_fast  = .TRUE.   ! to be run or not
   LOGICAL, PUBLIC :: lk_diaharm_2D   ! = .TRUE.   ! to run 2d
   LOGICAL, PUBLIC :: lk_diaharm_3D   ! = .TRUE.   ! to run 3d

   !! * Module variables
   INTEGER, PARAMETER ::  nharm_max  = jpmax_harmo  ! max number of harmonics to be analysed 
   INTEGER, PARAMETER ::  nhm_max    = 2*nharm_max+1 
   INTEGER, PARAMETER ::  nvab       = 2 ! number of 3D variables
   INTEGER            ::  nharm
   INTEGER            ::  nhm 
   INTEGER ::                 & !!! ** toto namelist (namtoto) **
      nflag  =  1                ! default value of nflag 
   REAL(wp), DIMENSION(nharm_max) ::                & 
      om_tide                     ! tidal frequencies ( rads/sec)
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:)   ::                & 
      bzz,c,x    ! work arrays
   REAL(wp) :: cca,ssa,zm,bt,dd_cumul
!
   REAL(wp), PUBLIC ::   fjulday_startharm       !: Julian Day since start of harmonic analysis
   REAL(wp), PUBLIC, ALLOCATABLE,DIMENSION(:) :: anau, anav, anaf   ! nodel/phase corrections used by diaharmana
   REAL(WP), ALLOCATABLE,SAVE,DIMENSION(:,:)   :: cc,a
!
   INTEGER ::  nvar_2d, nvar_3d    !: number of 2d and 3d variables to analyse
   INTEGER, ALLOCATABLE,DIMENSION(:) :: m_posi_2d, m_posi_3d

!  Name of variables used in the restart
   CHARACTER( LEN = 10 ), DIMENSION(5), PARAMETER :: m_varName2d = (/'ssh','u2d','v2d','ubfr','vbfr'/)
   CHARACTER( LEN = 10 ), DIMENSION(4), PARAMETER :: m_varName3d = (/'rho','u3d','v3d','w3d'/)
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:  ) :: g_cosamp2D, g_sinamp2D, g_cumul_var2D
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:) :: g_cosamp3D, g_sinamp3D, g_cumul_var3D  
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:)       :: g_out2D,h_out2D  ! arrays for output
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:)     :: g_out3D,h_out3D  ! arrays for 3D output
!
!  NAMELIST
   LOGICAL, PUBLIC :: ln_diaharm_store           !: =T  Stores data for harmonic Analysis
   LOGICAL, PUBLIC :: ln_diaharm_compute         !: =T  Compute harmonic Analysis
   LOGICAL, PUBLIC :: ln_diaharm_read_restart   !: =T  Read restart from a previous run 
   !JT
   LOGICAL, PUBLIC :: ln_diaharm_multiyear   !: =T  Read restart from a previous run 
   INTEGER, PUBLIC :: nn_diaharm_multiyear   !: =T  Read restart from a previous run 
   LOGICAL, PUBLIC :: ln_diaharm_update_nodal_daily   !: =T  update the nodes every day
   LOGICAL, PUBLIC :: ln_diaharm_fast
   LOGICAL, PUBLIC :: ln_diaharm_postproc_vel


   !JT
   LOGICAL, PUBLIC :: ln_ana_ssh, ln_ana_uvbar, ln_ana_bfric, ln_ana_rho, ln_ana_uv3d, ln_ana_w3d
   INTEGER ::   nb_ana        ! Number of harmonics to analyse
   CHARACTER (LEN=4), DIMENSION(jpmax_harmo) ::   tname   ! Names of tidal constituents ('M2', 'K1',...)
   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_all ! INDEX within the full set of constituents (tide.h90)
   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_sub ! INDEX within the subset of constituents pass in input

   !! * Substitutions

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! or LIM 2.0 , UCL-LOCEAN-IPSL (2005)
   !! or  TOP 1.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/module_example,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_harm_fast( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Harmonic analyser
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!        !  02-08  (Author names)  brief description of modifications
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * arguments
      INTEGER, INTENT( in  ) ::   &  
         kt                          ! describe it!!!

      !! * local declarations
      INTEGER  :: ji, jk, jj          ! dummy loop arguments
      INTEGER  :: jh, i1, i2, jgrid
      INTEGER  :: j2d, j3d
      REAL(WP) :: sec2start,sec2start_old
      CHARACTER (len=40) :: tmp_name
      !!--------------------------------------------------------------------

      !JT IF( nn_timing == 1 )   CALL timing_start( 'dia_harm_fast' )
      IF( ln_timing )   CALL timing_start( 'dia_harm_fast' )
      IF( kt == nit000   )   CALL harm_ana_init(kt)    ! Initialization (first time-step only)


      IF ( ln_diaharm_update_nodal_daily ) THEN
         !IF (MOD(kt,nint(86400./rdt)) == 0) THEN


         IF( nsec_day == NINT(0.5_wp * rdt) .OR. kt == nit000 ) THEN      
            DO jh = 1, nb_ana
               !JT anau(jh) = 3.141579*utide ( ntide_sub(jh) )/180.
               !JT anav(jh) = 3.141579*v0tide( ntide_sub(jh) )/180.
               anau(jh) = utide ( ntide_sub(jh) )
               anav(jh) = v0tide( ntide_sub(jh) )
               anaf(jh) = ftide ( ntide_sub(jh) )
            END DO

            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'harm_ana : update nodes?',ln_diaharm_update_nodal_daily
            IF(lwp) WRITE(numout,*) 'harm_ana : date, time',ndastp, nhour, nminute

         ENDIF
      ENDIF

     IF ( ln_diaharm_fast .and. ln_diaharm_store .and. ( lk_diaharm_2D .or. lk_diaharm_3D) ) THEN

          ! this bit done every time step
          nhm=2*nb_ana+1
          c(1) = 1.0

          sec2start_old = nint( (fjulday-fjulday_startharm)*86400._wp ) 
          sec2start = nsec_day - NINT(0.5_wp * rdt)

          sec2start = adatrj * 86400._wp

          
          IF(lwp) WRITE(numout,*) 'diaharm_fast: sec2start = ',nint( (fjulday-fjulday_startharm)*86400._wp ),nsec_day - NINT(0.5_wp * rdt),adatrj * 86400._wp

          IF( iom_use('tide_t') ) CALL iom_put( 'tide_t', sec2start )




          !IF(lwp) WRITE(numout,*) "ztime NEW", kt, sec2start, fjulday_startharm

          DO jh=1,nb_ana
             c(2*jh  ) = anaf(jh)*cos( sec2start*om_tide(jh) + anau(jh) + anav(jh) )
             c(2*jh+1) = anaf(jh)*sin( sec2start*om_tide(jh) + anau(jh) + anav(jh) )

             c(2*jh  ) = anaf(jh)*cos( sec2start*om_tide(jh) + anau(jh) )
             c(2*jh+1) = anaf(jh)*sin( sec2start*om_tide(jh) + anau(jh) )


             
             !IF(lwp) WRITE(numout,*) 'diaharm_fast: analwave,',kt,tname(jh),sec2start,sec2start/3600.,sec2start_old,sec2start_old/3600,c(2*jh),c(2*jh+1),om_tide(jh),anau(jh),anav(jh)
          ENDDO 

          !IF(lwp) WRITE(numout,*) "c init", c, "c end", sec2start, om_tide(1), anau(1), anav(1),"end nodal"


          ! CUMULATE
          DO ji=1,jpi         ! loop lon
             DO jj=1,jpj      ! loop lat
                DO jh=1,nhm   ! loop harmonic

                   DO j2d=1,nvar_2d
                      IF ( m_posi_2d(j2d) .eq. 1 ) dd_cumul = c(jh) * sshn(ji,jj) * ssmask (ji,jj)             ! analysis elevation
                      IF ( m_posi_2d(j2d) .eq. 2 ) dd_cumul = c(jh) * un_b(ji,jj) * ssumask(ji,jj)             ! analysis depth average velocities 
                      IF ( m_posi_2d(j2d) .eq. 3 ) dd_cumul = c(jh) * vn_b(ji,jj) * ssvmask(ji,jj)
                      !JT IF ( m_posi_2d(j2d) .eq. 4 ) dd_cumul = c(jh) * bfrua(ji,jj) * un(ji,jj,mbku(ji,jj)) * ssumask(ji,jj) ! analysis bottom friction
                      !JT IF ( m_posi_2d(j2d) .eq. 5 ) dd_cumul = c(jh) * bfrva(ji,jj) * vn(ji,jj,mbkv(ji,jj)) * ssvmask(ji,jj)
                      g_cumul_var2D(jh,ji,jj,j2d) = g_cumul_var2D(jh,ji,jj,j2d) + dd_cumul
                   ENDDO

                   DO j3d=1,nvar_3d
                      DO jk=1,jpkm1
                         IF ( m_posi_3d(j3d) .eq. 1 ) dd_cumul = c(jh) *  rhd(ji,jj,jk)               * tmask(ji,jj,jk)   
                         IF ( m_posi_3d(j3d) .eq. 2 ) dd_cumul = c(jh) * ( un(ji,jj,jk)-un_b(ji,jj) ) * umask(ji,jj,jk) 
                         IF ( m_posi_3d(j3d) .eq. 3 ) dd_cumul = c(jh) * ( vn(ji,jj,jk)-vn_b(ji,jj) ) * vmask(ji,jj,jk)
                         IF ( m_posi_3d(j3d) .eq. 4 ) dd_cumul = c(jh) *   wn(ji,jj,jk)               * wmask(ji,jj,jk)
                         g_cumul_var3D(jh,ji,jj,jk,j3d) = g_cumul_var3D(jh,ji,jj,jk,j3d) + dd_cumul
                      ENDDO
                   ENDDO

                ENDDO     ! end loop harmonic
             ENDDO        ! end loop lat
          ENDDO           ! end loop lon

          ! Compute nodal factor cumulative cross-product
          DO i1=1,nhm
             DO i2=1,nhm
                cc(i1,i2)=cc(i1,i2)+c(i1)*c(i2)
             ENDDO
          ENDDO

          ! Output RESTART
          IF( kt == nitrst ) THEN
             CALL harm_rst_write(kt) ! Dump out data for a restarted run 
          ENDIF

          ! At End of run
          IF ( kt ==  nitend ) THEN

             IF(lwp) WRITE(numout,*)
             IF(lwp) WRITE(numout,*) 'harm_ana : harmonic analysis of tides at end of run'
             IF(lwp) WRITE(numout,*) '~~~~~~~~~'

             IF( ln_diaharm_compute ) THEN

                 ! INITIALISE TABLE TO 0
                 IF ( nvar_2d .gt. 0 ) THEN
                    g_cosamp2D = 0.0_wp
                    g_sinamp2D = 0.0_wp
                 ENDIF
                 IF ( nvar_3d .gt. 0 ) THEN
                    g_cosamp3D = 0.0_wp
                    g_sinamp3D = 0.0_wp
                 ENDIF

                 ! FIRST OUTPUT 2D VARIABLES
                 DO jgrid=1,nvar_2d    ! loop number of 2d variables (ssh, U2d, V2d, UVfric) to analyse harmonically
                    DO ji=1,jpi        ! loop lon
                       DO jj=1,jpj     ! loop lat
                          bt = 1.0_wp; bzz(:) = 0.0_wp
                          DO jh=1,nhm  ! loop harmonic
                             bzz(jh) = g_cumul_var2D(jh,ji,jj,jgrid)
                             bt = bt*bzz(jh)
                          ENDDO
                          ! Copy back original cumulated nodal factor
                          a(:,:) = cc(:,:)
    !                     now do gaussian elimination of the system
    !                     a * x = b
    !                     the matrix x is (a0,a1,b1,a2,b2 ...)
    !                     the matrix a and rhs b solved here for x
                          x=0.0_wp
                          IF(bt.ne.0.) THEN
                            CALL gelim( a, bzz, x, nhm )
    !                       Backup output in variables
                            DO jh=1,nb_ana
                               g_cosamp2D(jh,ji,jj,jgrid) = x(jh*2  )
                               g_sinamp2D(jh,ji,jj,jgrid) = x(jh*2+1)
                            ENDDO
                            g_cosamp2D   ( 0,ji,jj,jgrid) = x(1)
                            g_sinamp2D   ( 0,ji,jj,jgrid) = 0.0_wp
                          ENDIF     ! bt.ne.0.
                       ENDDO        ! jj
                    ENDDO           ! ji
                 ENDDO              ! jgrid

                 ! SECOND OUTPUT 3D VARIABLES
                 DO jgrid=1,nvar_3d     ! loop number of 3d variables rho, U, V, W
                    DO jk=1,jpkm1       ! loop over vertical level
                       DO ji=1,jpi      ! loop over lon
                          DO jj=1,jpj   ! loop over lat
                             bt = 1.0_wp; bzz(:) = 0.0_wp
                             DO jh=1,nhm
                                bzz(jh) = g_cumul_var3D(jh,ji,jj,jk,jgrid)
                                bt = bt*bzz(jh)
                             ENDDO
                             ! Copy back original cumulated nodal factor
                             a(:,:) = cc(:,:)                      
    !                        now do gaussian elimination of the system
    !                        a * x = b
    !                        the matrix x is (a0,a1,b1,a2,b2 ...)
    !                        the matrix a and rhs b solved here for x
                             x=0.0_wp
                             IF(bt.ne.0.) THEN
                               CALL gelim( a, bzz, x, nhm )
    !                          Backup output in variables
                               DO jh=1,nb_ana
                                  g_cosamp3D(jh,ji,jj,jk,jgrid) = x(jh*2  )
                                  g_sinamp3D(jh,ji,jj,jk,jgrid) = x(jh*2+1)
                               ENDDO
                               g_cosamp3D   ( 0,ji,jj,jk,jgrid) = x(1)
                               g_sinamp3D   ( 0,ji,jj,jk,jgrid) = 0.0_wp
                            ENDIF     ! bt.ne.0.
                          ENDDO       ! jj
                       ENDDO          ! ji
                    ENDDO             ! jk
                 ENDDO                ! jgrid

                 CALL harm_ana_out     ! output analysis (last time step)

             ELSE    ! ln_harmana_compute = False 
                 IF(lwp) WRITE(numout,*) " Skipping Computing harmonics at last step"

             ENDIF   ! ln_harmana_compute 
          ENDIF      ! kt ==  nitend

     ENDIF

      !JT IF( nn_timing == 1 )   CALL timing_stop( 'dia_harm_fast' )
      IF( ln_timing  )   CALL timing_stop( 'dia_harm_fast' )

   END SUBROUTINE dia_harm_fast 

   SUBROUTINE harm_ana_init( kt )   !JT
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
      !! * local declarations

      INTEGER, INTENT(in) ::   kt     ! ocean time-step  !JT
      !!
      INTEGER ::   ji, jk, jh  ! dummy loop indices
      INTEGER ::   ios                  ! Local integer output status for namelist read
      INTEGER ::   k2d, k3d             ! dummy number of analysis


      

      NAMELIST/nam_diaharm_fast/ ln_diaharm_fast, ln_diaharm_store, ln_diaharm_compute, ln_diaharm_read_restart, ln_ana_ssh, ln_ana_uvbar, ln_ana_bfric, ln_ana_rho, ln_ana_uv3d, ln_ana_w3d, &
               & tname,ln_diaharm_multiyear,nn_diaharm_multiyear,ln_diaharm_update_nodal_daily,ln_diaharm_postproc_vel
      !!----------------------------------------------------------------------
      !JT
      ln_diaharm_fast = .FALSE.
      ln_diaharm_multiyear = .FALSE.
      nn_diaharm_multiyear = 20
      !JT
      lk_diaharm_2D    = .TRUE.   ! to run 2d
      lk_diaharm_3D    = .TRUE.   ! to run 3d

      ln_diaharm_store = .TRUE.

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'harm_init : initialization of harmonic analysis of tides'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

      ! GET NAMELIST DETAILS
      REWIND( numnam_ref )              ! Namelist nam_diaharm_fast in reference namelist : Tidal harmonic analysis
      READ  ( numnam_ref, nam_diaharm_fast, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm_fast in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist nam_diaharm_fast in configuration namelist : Tidal harmonic analysis
      READ  ( numnam_cfg, nam_diaharm_fast, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 ) CALL ctl_nam ( ios , 'nam_diaharm_fast in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_diaharm_fast )


      !
      IF(lwp) THEN
         WRITE(numout,*) 'Tidal diagnostics_fast'
         WRITE(numout,*) '   Fast Harmonic analysis ?:         ln_diaharm_fast= ', ln_diaharm_fast
         WRITE(numout,*) '   Store output in restart?:         ln_diaharm_store= ', ln_diaharm_store
         WRITE(numout,*) '   Compute analysis?:         ln_diaharm_compute= ', ln_diaharm_compute
         WRITE(numout,*) '   Read in restart? : ln_diaharm_read_restart = ', ln_diaharm_read_restart
         WRITE(numout,*) '   SSH harmonic analysis: ln_ana_ssh = ', ln_ana_ssh
         WRITE(numout,*) '   Barotropic Velocities harmonic analysis: ln_ana_uvbar = ', ln_ana_uvbar
         WRITE(numout,*) '   Bed Friction for harmonic analysis (not implemented): ln_ana_bfric = ', ln_ana_bfric
         WRITE(numout,*) '   Density harmonic analysis: ln_ana_rho = ', ln_ana_rho
         WRITE(numout,*) '   3D velocities harmonic analysis: ln_ana_uv3d = ', ln_ana_uv3d
         WRITE(numout,*) '   Vertical Velocities harmonic analysis: ln_ana_w3d = ', ln_ana_w3d
         WRITE(numout,*) '   Names of harmonics: tname = ', tname
         WRITE(numout,*) '   Max number of harmonics: jpmax_harmo = ', jpmax_harmo
         WRITE(numout,*) '   Number of Harmonics: nb_harmo = ', nb_harmo
         WRITE(numout,*) '   Multi-year harmonic analysis: ln_diaharm_multiyear = ', ln_diaharm_multiyear
         WRITE(numout,*) '   Multi-year harmonic analysis - number of years: nn_diaharm_multiyear = ', nn_diaharm_multiyear
         WRITE(numout,*) '   Multi-year harmonic analysis - number of years: ln_diaharm_update_nodal_daily = ', ln_diaharm_update_nodal_daily
         WRITE(numout,*) '   Number of Harmonics: nyear, nmonth = ', nyear, nmonth
         WRITE(numout,*) '   Post-process velocity stats: ln_diaharm_postproc_vel = ', ln_diaharm_postproc_vel

      ENDIF
      ! JT


      IF ( ln_diaharm_multiyear ) THEN
          ln_diaharm_store = .True.
          IF(lwp) WRITE(numout,*) '   Multi-year harmonic analysis ', nyear,nn_diaharm_multiyear,nmonth
          IF ((mod(nyear,nn_diaharm_multiyear) == 0) .AND. ( nmonth == 1 )) THEN ! Jan, year = 1980,2000,2020,2040, restart tidal calculation
              ln_diaharm_read_restart = .FALSE.
              IF(lwp) WRITE(numout,*) '   Read in restart? : ln_diaharm_read_restart = ', ln_diaharm_read_restart
          ELSE
              ln_diaharm_read_restart = .TRUE.
              IF(lwp) WRITE(numout,*) '   Read in restart? : ln_diaharm_read_restart = ', ln_diaharm_read_restart
          ENDIF



          IF(lwp) WRITE(numout,*) '   Multi-year harmonic analysis ', nyear,nn_diaharm_multiyear,nmonth
          IF ((mod(nyear,nn_diaharm_multiyear) == (nn_diaharm_multiyear - 1)) .AND. ( nmonth == 12 )) THEN ! Dec year = 1999,2019,2039,2040, restart tidal calculation
              ln_diaharm_compute = .TRUE.
              IF(lwp) WRITE(numout,*) '   Compute analysis?:         ln_diaharm_compute= ', ln_diaharm_compute
          ELSE
              ln_diaharm_compute = .FALSE.
              IF(lwp) WRITE(numout,*) '   Compute analysis?:         ln_diaharm_compute= ', ln_diaharm_compute
          ENDIF

      ENDIF
      IF ( kt < 10 ) THEN
        ln_diaharm_read_restart = .FALSE.
        IF(lwp) WRITE(numout,*) '   kt = ',kt
        IF(lwp) WRITE(numout,*) '   kt < 10, so setting ln_diaharm_read_restart to .FALSE.'
        IF(lwp) WRITE(numout,*) '   Read in restart? : ln_diaharm_read_restart = ', ln_diaharm_read_restart        
      ENDIF

      ! JT


      ! GET NUMBER OF HARMONIC TO ANALYSE - from diaharm.F90
      nb_ana = 0
      DO jk=1,jpmax_harmo
         DO ji=1,nb_harmo
            IF(TRIM(tname(jk)) == TRIM(Wave( ntide(ji) )%cname_tide) ) THEN
               nb_ana=nb_ana+1
            ENDIF
         END DO
      END DO
      !
      IF(lwp) THEN
         WRITE(numout,*) '        Namelist nam_diaharm_fast'
         WRITE(numout,*) '        nb_ana    = ', nb_ana
         CALL flush(numout)
      ENDIF
      !
      IF (nb_ana > nharm_max) THEN
        IF(lwp) WRITE(numout,*) ' E R R O R harm_ana : nb_ana must be lower than nharm_max, stop'
        IF(lwp) WRITE(numout,*) ' nharm_max = ', nharm_max
        nstop = nstop + 1
      ENDIF


      ALLOCATE(ntide_all(nb_ana))
      ALLOCATE(ntide_sub(nb_ana))

      DO jk=1,nb_ana
       DO ji=1,nb_harmo
          IF (TRIM(tname(jk)) .eq. Wave( ntide(ji) )%cname_tide ) THEN
             ntide_sub(jk) = ji
             ntide_all(jk) = ntide(ji)
             EXIT
          END IF
       END DO
      END DO

     ALLOCATE( anau(nb_ana) )
     ALLOCATE( anav(nb_ana) )
     ALLOCATE( anaf(nb_ana) )


    IF( ln_diaharm_fast ) THEN

          ! SEARCH HOW MANY VARIABLES 2D AND 3D TO PROCESS
          nvar_2d = 0; nvar_3d = 0
          IF ( ln_ana_ssh   ) nvar_2d = nvar_2d + 1       ! analysis elevation
          IF ( ln_ana_uvbar ) nvar_2d = nvar_2d + 2       ! analysis depth-averaged velocity
          IF ( ln_ana_bfric ) nvar_2d = nvar_2d + 2       ! analysis bottom friction 
                
          IF ( ln_ana_rho   ) nvar_3d = nvar_3d + 1       ! analysis density
          IF ( ln_ana_uv3d  ) nvar_3d = nvar_3d + 2       ! analysis 3D horizontal velocities
          IF ( ln_ana_w3d   ) nvar_3d = nvar_3d + 1       ! analysis 3D vertical velocity

          ! CHECK IF SOMETHING TO RUN
          IF ( nvar_2d .eq. 0 ) lk_diaharm_2D = .FALSE.   ! no 2d to run
          IF ( nvar_3d .eq. 0 ) lk_diaharm_3D = .FALSE.   ! no 3d to run
    !      IF ( nvar_2d .gt. 0 .and. nvar_3d .gt. 0 ) lk_diaharm_fast = .FALSE.
    !      IF ( .NOT. ln_diaharm_store ) lk_diaharm_fast = .FALSE.

          IF ( ln_diaharm_store .and. ( lk_diaharm_2D .or. lk_diaharm_3D) ) THEN

             ! DO ALLOCATIONS
             IF ( lk_diaharm_2D ) THEN
                ALLOCATE( g_cumul_var2D(nb_ana*2+1,jpi,jpj,    nvar_2d) )
                ALLOCATE( g_cosamp2D( 0:nb_ana*2+1,jpi,jpj,    nvar_2d) )
                ALLOCATE( g_sinamp2D( 0:nb_ana*2+1,jpi,jpj,    nvar_2d) )
                ALLOCATE( g_out2D (jpi,jpj) )
                ALLOCATE( h_out2D (jpi,jpj) )
                ALLOCATE( m_posi_2d( nvar_2d ) ); m_posi_2d(:)=0
             ENDIF
     
             IF ( lk_diaharm_3D ) THEN
                ALLOCATE( g_cumul_var3D(nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
                ALLOCATE( g_cosamp3D( 0:nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
                ALLOCATE( g_sinamp3D( 0:nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
                ALLOCATE( g_out3D (jpi,jpj,jpk) )
                ALLOCATE( h_out3D (jpi,jpj,jpk) )
                ALLOCATE( m_posi_3d( nvar_3d ) ); m_posi_3d(:)=0
             ENDIF

             ALLOCATE( cc(nb_ana*2+1,nb_ana*2+1) )
             ALLOCATE( a (nb_ana*2+1,nb_ana*2+1) )
             ALLOCATE( bzz(nb_ana*2+1) )
             ALLOCATE( x  (nb_ana*2+1) )
             ALLOCATE( c  (nb_ana*2+1) )
             ! END ALLOCATE 

             ! STORE INDEX OF WHAT TO PRODUCE DEPENDING ON ACTIVATED LOGICAL
             ! MAKES THINGS EASIER AND FASTER LATER
             ! !!! UGLY !!!
             jh = 1; k2d = 0; 
             IF ( ln_ana_ssh   ) THEN
                k2d = k2d + 1; m_posi_2d(k2d) = jh
                IF(lwp) WRITE(numout,*) "   - ssh harmonic analysis activated (ln_ana_ssh)"
             ENDIF
             jh = jh + 1
             IF ( ln_ana_uvbar ) THEN
                k2d = k2d + 1; m_posi_2d(k2d) = jh
                jh  = jh  + 1 
                k2d = k2d + 1; m_posi_2d(k2d) = jh
                IF(lwp) WRITE(numout,*) "   - barotropic currents harmonic analysis activated (ln_ana_uvbar)"
             ELSE
                jh  = jh  + 1
             ENDIF
             jh = jh + 1
             IF ( ln_ana_bfric ) THEN
                k2d = k2d + 1; m_posi_2d(k2d) = jh
                jh  = jh  + 1; 
                k2d = k2d + 1; m_posi_2d(k2d) = jh
                IF(lwp) WRITE(numout,*) "   - bottom friction harmonic analysis activated (ln_ana_vbfr)"
             ELSE
                jh  = jh  + 1
             ENDIF

             ! and for 3D
             jh = 1; k3d = 0; 
             IF ( ln_ana_rho  ) THEN
                k3d = k3d + 1; m_posi_3d(k3d) = jh
                IF(lwp) WRITE(numout,*) "   - 3D density harmonic analysis activated (ln_ana_rho)"
             ENDIF
             jh = jh + 1
             IF ( ln_ana_uv3d )  THEN
                k3d = k3d + 1; m_posi_3d(k3d) = jh
                jh  = jh  + 1 
                k3d = k3d + 1; m_posi_3d(k3d) = jh
                IF(lwp) WRITE(numout,*) "   - 3D horizontal currents harmonic analysis activated (ln_ana_uv3d)"
             ELSE
                jh  = jh  + 1
             ENDIF
             jh = jh + 1
             IF ( ln_ana_w3d ) THEN
                k3d = k3d + 1; m_posi_3d(k3d) = jh
                IF(lwp) WRITE(numout,*) "   - 3D vertical currents harmonic analysis activated (ln_ana_w3d)"
             ENDIF

             ! SELECT AND STORE FREQUENCIES
             IF(lwp)    WRITE(numout,*) 'Analysed frequency  : ',nb_ana ,'Frequency '
             DO jh=1,nb_ana
                om_tide(jh) = omega_tide( ntide_sub(jh) ) 
                IF(lwp) WRITE(numout,*) '        - ',tname(jh),' ',om_tide(jh), (2*rpi/3600.)/om_tide(jh),"hr"
             ENDDO

             ! READ RESTART IF 
             IF ( ln_diaharm_read_restart ) THEN
                IF (lwp) WRITE(numout,*) "Reading previous harmonic data from previous run. kt = ",kt
                ! Need to read in  bssh bz, cc anau anav and anaf 
                call harm_rst_read  ! This reads in from the previous day
                                    ! Currrently the data in in assci format
             ELSE 

                IF (lwp) WRITE(numout,*) "Starting harmonic analysis from Fresh. kt = ",kt 
     
                IF ( lk_diaharm_2D ) g_cumul_var2D(:,:,:,:  ) = 0.0_wp
                IF ( lk_diaharm_3D ) g_cumul_var3D(:,:,:,:,:) = 0.0_wp
                cc           = 0.0_wp
                a    (:,:)   = 0.0_wp ! NB
                bzz  (:)     = 0.0_wp
                x    (:)     = 0.0_wp
                c    (:)     = 0.0_wp
                anau (:)     = 0.0_wp
                anav (:)     = 0.0_wp
                anaf (:)     = 0.0_wp

                DO jh = 1, nb_ana
                   anau(jh) = utide ( ntide_sub(jh) )
                   anav(jh) = v0tide( ntide_sub(jh) )
                   anaf(jh) = ftide ( ntide_sub(jh) )
                END DO

                fjulday_startharm=fjulday !Set this at very start and store
                !JT this is a mistake - only works on daily cycles, should use fjulnsec_dayday

                IF (lwp) THEN
                   WRITE(numout,*) '--------------------------'
                   WRITE(numout,*) '   - Output anaf for check'
                   WRITE(numout,*) 'ANA F', anaf
                   WRITE(numout,*) 'ANA U', anau
                   WRITE(numout,*) 'ANA V', anav
                   WRITE(numout,*) 'fjulday',fjulday
                   WRITE(numout,*) 'fjulday_startharm',fjulday_startharm
                   WRITE(numout,*) 'nsec_day',nsec_day
                   WRITE(numout,*) 'kt',kt
                   WRITE(numout,*) '--------------------------'
                ENDIF

             ENDIF

          ELSE

             IF (lwp) WRITE(numout,*) "No variable setup for harmonic analysis"

          ENDIF
      ENDIF

   END SUBROUTINE harm_ana_init
!
   SUBROUTINE gelim (a,b,x,n)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Guassian elimination
      !!
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
        implicit none
!
        integer  :: n
        REAL(WP) :: b(nb_ana*2+1), a(nb_ana*2+1,nb_ana*2+1)
        REAL(WP) :: x(nb_ana*2+1)
        INTEGER  :: row,col,prow,pivrow,rrow
        REAL(WP) :: atemp
        REAL(WP) :: pivot
        REAL(WP) :: m

        do row=1,n-1
           pivrow=row
           pivot=a(row,n-row+1)
           do prow=row+1,n
              if (abs(a(prow,n-row+1)).gt.abs(pivot)  ) then
                 pivot=a(prow,n-row+1)
                 pivrow=prow
              endif
           enddo
!	swap row and prow
           if ( pivrow .ne. row ) then
              atemp=b(pivrow)
              b(pivrow)=b(row)
              b(row)=atemp
              do col=1,n
                 atemp=a(pivrow,col)
                 a(pivrow,col)=a(row,col)
                 a(row,col)=atemp
              enddo
           endif

           do rrow=row+1,n
              if (a(row,row).ne.0) then
   
                 m=-a(rrow,n-row+1)/a(row,n-row+1)
                 do col=1,n
                    a(rrow,col)=m*a(row,col)+a(rrow,col)
                 enddo
                 b(rrow)=m*b(row)+b(rrow)
              endif
           enddo
        enddo
!	back substitution now

        x(1)=b(n)/a(n,1)
        do row=n-1,1,-1
           x(n-row+1)=b(row)
           do col=1,(n-row)
              x(n-row+1)=(x(n-row+1)-a(row,col)*x(col)) 
           enddo

           x(n-row+1)=(x(n-row+1)/a(row,(n-row)+1))
        enddo

        return
   END SUBROUTINE gelim

   SUBROUTINE harm_ana_out
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
        USE dianam          ! build name of file (routine)
 
      !! * local declarations
      INTEGER :: ji, jj, jk, jgrid, jh    ! dummy loop indices
!      INTEGER :: nh_T
!      INTEGER :: nid_harm
!      CHARACTER (len=40) :: clhstnamt, clop1, clop2 ! temporary names 
!      CHARACTER (len=40) :: clhstnamu, clhstnamv    ! temporary names 
      CHARACTER (len=40) :: suffix
      CHARACTER (len=40) :: tmp_name
!      REAL(wp) :: zsto1, zsto2, zout, zmax, zjulian, zdt, zmdi  ! temporary scalars

      REAL(wp), ALLOCATABLE,DIMENSION(:,:,:)       :: amp_u2d,phi_u2d, amp_v2d,phi_v2d  ! arrays for output

      REAL(wp)   :: tmp_u_amp ,tmp_v_amp ,tmp_u_phi ,tmp_v_phi
      REAL(wp)   :: a_u, b_u, a_v, b_v, twodelta, delta, alpha2, alpha, qmin, qmax, ecc,thetamax, thetamin
      REAL(wp)   :: Qc, Qac, gc,gac, Phi_Ua, dir_Ua, polarity
      REAL(wp)   :: tmpreal

      REAL(wp), ALLOCATABLE,DIMENSION(:,:)         :: tmp_u_amp_mat,tmp_v_amp_mat,tmp_u_phi_mat,tmp_v_phi_mat
!      REAL(wp), ALLOCATABLE,DIMENSION(:,:)         :: a_u_mat,b_u_mat,a_v_mat,b_v_mat,qmax_mat,qmin_mat,ecc_mat
!      REAL(wp), ALLOCATABLE,DIMENSION(:,:)         :: thetamax_mat,thetamin_mat,Qc_mat,Qac_mat,gc_mat,gac_mat
!      REAL(wp), ALLOCATABLE,DIMENSION(:,:)         :: Phi_Ua_mat,dir_Ua_mat,polarity_mat



      IF (ln_diaharm_postproc_vel .AND. ln_ana_uvbar)  THEN
         ALLOCATE( amp_u2d(jh,jpi,jpj),amp_v2d(jh,jpi,jpj),phi_u2d(jh,jpi,jpj),phi_v2d(jh,jpi,jpj) )


         ALLOCATE(tmp_u_amp_mat(jpi,jpj),tmp_v_amp_mat(jpi,jpj),tmp_u_phi_mat(jpi,jpj),tmp_v_phi_mat(jpi,jpj))
!         ALLOCATE(a_u_mat(jpi,jpj),b_u_mat(jpi,jpj),a_v_mat(jpi,jpj),b_v_mat(jpi,jpj))
!         ALLOCATE(qmax_mat(jpi,jpj),qmin_mat(jpi,jpj),ecc_mat(jpi,jpj))
!         ALLOCATE(thetamax_mat(jpi,jpj),thetamin_mat(jpi,jpj),Qc_mat(jpi,jpj),Qac_mat(jpi,jpj))
!         ALLOCATE(gc_mat(jpi,jpj),gac_mat(jpi,jpj),Phi_Ua_mat(jpi,jpj),dir_Ua_mat(jpi,jpj),polarity_mat(jpi,jpj))

      endif

      do jgrid=1,nvar_2d
          do jh=1,nb_ana
             h_out2D = 0.0
             g_out2D = 0.0
             do jj=1,nlcj
                do ji=1,nlci
                   cca=g_cosamp2D(jh,ji,jj,jgrid)
                   ssa=g_sinamp2D(jh,ji,jj,jgrid)
                   h_out2D(ji,jj)=sqrt(cca**2+ssa**2)
                   IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN 
                      g_out2D(ji,jj)= 0.0_wp
                   ELSE
                      g_out2D(ji,jj)=(180.0/rpi)*atan2(ssa,cca)       
                   ENDIF 
                   IF (h_out2D(ji,jj).ne.0) THEN
                       h_out2D(ji,jj)=h_out2D(ji,jj)/anaf(jh)
                   ENDIF
                   IF (g_out2D(ji,jj).ne.0) THEN  !Correct and take modulus
                       !JT 
                       !JT  g_out2D(ji,jj) = g_out2D(ji,jj) + MOD( (anau(jh)+anav(jh))/rad , 360.0)
                       !JT 
                       g_out2D(ji,jj) = g_out2D(ji,jj) + MOD( (anau(jh))/rad , 360.0)
                       if (g_out2D(ji,jj).gt.360.0) then
                           g_out2D(ji,jj)=g_out2D(ji,jj)-360.0
                       else if (g_out2D(ji,jj).lt.0.0) then
                           g_out2D(ji,jj)=g_out2D(ji,jj)+360.0
                       endif
                   ENDIF
                enddo
             enddo
             !
             ! NETCDF OUTPUT
             suffix = TRIM( m_varName2d( m_posi_2d(jgrid) ) )
             IF(lwp) WRITE(numout,*) "harm_ana_out", suffix

             tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'amp_'//TRIM(suffix)
             IF( iom_use(TRIM(tmp_name)) )  THEN
                IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name),'; shape = ', SHAPE(h_out2D)
                IF(lwp) WRITE(numout,*) "harm_ana_out names", tmp_name,tname(jh),' ',om_tide(jh), (2*rpi/3600.)/om_tide(jh),"hr"
                CALL iom_put( TRIM(tmp_name), h_out2D(:,:) )
             ELSE
                IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
             ENDIF

             tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'pha_'//TRIM(suffix)
             IF( iom_use(TRIM(tmp_name)) )  THEN
                IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name),'; shape = ', SHAPE(g_out2D)
                CALL iom_put( TRIM(tmp_name), g_out2D(:,:) )
             ELSE
                IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
             ENDIF



             IF (ln_diaharm_postproc_vel .AND. ln_ana_uvbar)  THEN

               !IF (m_posi_2d(jgrid) == 2) THEN
               IF (TRIM(suffix) == TRIM('u2d')) THEN
                  if (lwp)  WRITE(numout,*) "harm_ana_out ln_diaharm_postproc_vel: "//TRIM(Wave(ntide_all(jh))%cname_tide)//' u2d  '//TRIM(suffix)
                  amp_u2d(jh,:,:) = h_out2D(:,:)
                  phi_u2d(jh,:,:) = rpi*g_out2D(:,:)/180.0
               ENDIF

               !IF (m_posi_2d(jgrid) == 3) THEN
               IF (TRIM(suffix) == TRIM('v2d')) THEN
                  if (lwp)  WRITE(numout,*) "harm_ana_out ln_diaharm_postproc_vel: "//TRIM(Wave(ntide_all(jh))%cname_tide)//' v2d  '//TRIM(suffix)
                  amp_v2d(jh,:,:) = h_out2D(:,:)
                  phi_v2d(jh,:,:) = rpi*g_out2D(:,:)/180.0
               ENDIF
             ENDIF

             CALL FLUSH(numout)


          enddo

         suffix = TRIM( m_varName2d( m_posi_2d(jgrid) ) )
         tmp_name='TA_'//TRIM(suffix)//'_off'
         IF( iom_use(TRIM(tmp_name)) )  THEN
            IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
            CALL iom_put( TRIM(tmp_name), g_cosamp2D( 0,:,:,jgrid))
         ELSE
            IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
         ENDIF

         CALL FLUSH(numout)

      enddo
!
! DO THE SAME FOR 3D VARIABLES
!
      do jgrid=1,nvar_3d
          do jh=1,nb_ana
             h_out3D = 0.0
             g_out3D = 0.0
             DO jk=1,jpkm1
                do jj=1,nlcj
                   do ji=1,nlci
                      cca=g_cosamp3D(jh,ji,jj,jk,jgrid)
                      ssa=g_sinamp3D(jh,ji,jj,jk,jgrid)
                      h_out3D(ji,jj,jk)=sqrt(cca**2+ssa**2)
                      IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN
                         g_out3D(ji,jj,jk) = 0.0_wp
                      ELSE
                         g_out3D(ji,jj,jk) = (180.0/rpi)*atan2(ssa,cca)
                      ENDIF
                      IF (h_out3D(ji,jj,jk).ne.0) THEN
                          h_out3D(ji,jj,jk) = h_out3D(ji,jj,jk)/anaf(jh)
                      ENDIF
                      IF (g_out3D(ji,jj,jk).ne.0) THEN  !Correct and take modulus
                          !JT                        
                          !JT g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk) + MOD( (anau(jh)+anav(jh))/rad , 360.0)
                          !JT 
                          g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk) + MOD( (anau(jh))/rad , 360.0)
                          if      (g_out3D(ji,jj,jk).gt.360.0) then
                                   g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk)-360.0
                          else if (g_out3D(ji,jj,jk).lt.0.0) then
                                   g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk)+360.0
                          endif
                      ENDIF
                   enddo    ! ji
                enddo       ! jj
             ENDDO          ! jk
             !
             ! NETCDF OUTPUT
             suffix = TRIM( m_varName3d( m_posi_3d(jgrid) ) )
             IF(lwp) WRITE(numout,*) "harm_ana_out", suffix

             tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'amp_'//TRIM(suffix)
             IF( iom_use(TRIM(tmp_name)) )  THEN
                IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name),'; shape = ', SHAPE(h_out3D)
                CALL iom_put( TRIM(tmp_name), h_out3D(:,:,:) )
             ELSE
                IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
             ENDIF

             tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'pha_'//TRIM(suffix)
             IF( iom_use(TRIM(tmp_name)) )  THEN
                IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name),'; shape = ', SHAPE(g_out3D)
                CALL iom_put(tmp_name, g_out3D(:,:,:) )
             ELSE
                IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
             ENDIF

          enddo             ! jh 

         suffix = TRIM( m_varName3d( m_posi_3d(jgrid) ) )
         tmp_name='TA_'//TRIM(suffix)//'_off'
         IF( iom_use(TRIM(tmp_name)) )  THEN
            IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
            CALL iom_put( TRIM(tmp_name), g_cosamp3D( 0,:,:,:,jgrid))
         ELSE
            IF(lwp) WRITE(numout,*) "harm_ana_out: not requested: ",TRIM(tmp_name)
         ENDIF

      enddo                 ! jgrid

     CALL FLUSH(numout)

      IF (ln_diaharm_postproc_vel .AND. ln_ana_uvbar)  THEN
         IF(lwp) WRITE(numout,*) "diaharm_fast: Postprocess barotropic velocity tidal parameters"
         CALL FLUSH(numout)
         DO jh=1,nb_ana


            tmp_u_amp_mat(:,:) = 0.
            tmp_v_amp_mat(:,:) = 0.
            tmp_u_phi_mat(:,:) = 0.
            tmp_v_phi_mat(:,:) = 0.

!            a_u_mat(:,:) = 0.
!            b_u_mat(:,:) = 0.
!            a_v_mat(:,:) = 0.
!            b_v_mat(:,:) = 0.

!            qmax_mat(:,:) = 0.
!            qmin_mat(:,:) = 0.

!            ecc_mat(:,:) = 0
!            thetamax_mat(:,:) =0.
!            thetamin_mat(:,:) = 0.

!            Qc_mat(:,:) = 0.
!            Qac_mat(:,:) = 0.
!            gc_mat(:,:) = 0.
!            gac_mat(:,:) = 0.

!            Phi_Ua_mat(:,:) = 0.
!            dir_Ua_mat(:,:) = 0.
!            polarity_mat(:,:) = 0.


!             DO jj = 2, nlcj - 1
!                DO ji = 2, nlci - 1

!             do jj=2,nlcj
!                do ji=2,nlci
                    !IF ((ssumask(ji,jj) + ssumask(ji-1,jj)) == 0 ) CYCLE
                    !IF ((ssvmask(ji,jj) + ssvmask(ji,jj-1)) == 0 ) CYCLE

!                    IF ( ((ssumask(ji,jj) + ssumask(ji-1,jj)) > 0 ) .AND. ((ssvmask(ji,jj) + ssvmask(ji,jj-1)) > 0 ) ) THEN
!                        tmp_u_amp = ((amp_u2d(jh,ji,jj)*ssumask(ji,jj)) + (amp_u2d(jh,ji-1,jj)*ssumask(ji-1,jj)))/(ssumask(ji,jj) + ssumask(ji-1,jj))
!                        tmp_v_amp = ((amp_v2d(jh,ji,jj)*ssvmask(ji,jj)) + (amp_v2d(jh,ji,jj-1)*ssvmask(ji,jj-1)))/(ssvmask(ji,jj) + ssvmask(ji,jj-1))
!                        ! WORK ON THE WRAP AROUND
!                        tmp_u_phi = ((phi_u2d(jh,ji,jj)*ssumask(ji,jj)) + (phi_u2d(jh,ji-1,jj)*ssumask(ji-1,jj)))/(ssumask(ji,jj) + ssumask(ji-1,jj))
!                        tmp_v_phi = ((phi_v2d(jh,ji,jj)*ssvmask(ji,jj)) + (phi_v2d(jh,ji,jj-1)*ssvmask(ji,jj-1)))/(ssvmask(ji,jj) + ssvmask(ji,jj-1))

             do jj=1,nlcj
                do ji=1,nlci

!                        tmp_u_amp = ((amp_u2d(jh,ji,jj)) + (amp_u2d(jh,ji-1,jj)))/(2.)
!                        tmp_v_amp = ((amp_v2d(jh,ji,jj)) + (amp_v2d(jh,ji,jj-1)))/(2.)
!                        ! WORK ON THE WRAP AROUND
!                        tmp_u_phi = ((phi_u2d(jh,ji,jj)) + (phi_u2d(jh,ji-1,jj)))/(2.)
!                        tmp_v_phi = ((phi_v2d(jh,ji,jj)) + (phi_v2d(jh,ji,jj-1)))/(2.)



                        tmp_u_amp = (amp_u2d(jh,ji,jj)) 
                        tmp_v_amp = (amp_v2d(jh,ji,jj)) 
                        ! WORK ON THE WRAP AROUND
                        tmp_u_phi = (phi_u2d(jh,ji,jj)) 
                        tmp_v_phi = (phi_v2d(jh,ji,jj)) 



!                        a_u = tmp_U_amp * cos(tmp_U_phi)
!                        b_u = tmp_U_amp * sin(tmp_U_phi)
!                        a_v = tmp_V_amp * cos(tmp_V_phi)
!                        b_v = tmp_V_amp * sin(tmp_V_phi)

!                        twodelta =  atan2( (tmp_V_amp**2  * sin( 2*(tmp_U_phi - tmp_V_phi)  ) ) , (   tmp_U_amp**2   +   tmp_V_amp**2  * cos( 2*(tmp_U_phi - tmp_V_phi)  )     ) )
!                        delta = twodelta/2.

!                        !alpha2 = sqrt( tmp_U_amp**4 + tmp_V_amp**4 + 2*tmp_U_amp**2*tmp_V_amp**2*cos(2*(tmp_U_phi - tmp_V_phi))  )

!                        tmpreal = tmp_U_amp**4 + tmp_V_amp**4 + 2*tmp_U_amp**2*tmp_V_amp**2*cos(2*(tmp_U_phi - tmp_V_phi)) 
!                        if (tmpreal < 0) CYCLE
!                        alpha2 = sqrt( tmp_U_amp**4 + tmp_V_amp**4 + 2*tmp_U_amp**2*tmp_V_amp**2*cos(2*(tmp_U_phi - tmp_V_phi))  )
!                        if (alpha2 < 0) CYCLE
!                        alpha= sqrt( alpha2 )


!                        !major and minor axis of the ellipse
!                        !qmax = sqrt( (tmp_U_amp**2 + tmp_V_amp**2 + alpha**2)/2 )
!                        !tmpreal =  (tmp_U_amp**2 + tmp_V_amp**2 - alpha**2)/2
!                        !qmin = 0
!                        !if (tmpreal > 0) qmin = sqrt( (tmp_U_amp**2 + tmp_V_amp**2 - alpha**2)/2 )   ! but always positive.

!                        tmpreal =  (tmp_U_amp**2 + tmp_V_amp**2 - alpha**2)/2
!                        if (tmpreal < 0) CYCLE
!                        qmin = sqrt( (tmp_U_amp**2 + tmp_V_amp**2 - alpha**2)/2 )   ! but always positive.

!                        !eccentricity of ellipse
!                        ecc = (qmax - qmin)/(qmax + qmin)
!                        ! Angle of major and minor ellipse
!                        thetamax = atan2((  tmp_V_amp * cos((tmp_U_phi - tmp_V_phi) - delta)   ) , ( tmp_U_amp * cos( delta) )  )
!                        thetamin = thetamax + rpi/2.



!                        ! Rotary current components: Pugh A3.10
!                        ! Clockwise (c) and anticlockwise (ac) rotating rotate_wind_vectors
!                        ! so   Qc = clockwise     = anticyclonic = negative
!                        ! and Qac = anticlockwise = cyclonic     = negative

!                        tmpreal = tmp_U_amp**2 + tmp_V_amp**2 - (2*tmp_U_amp*tmp_V_amp*sin( tmp_V_phi - tmp_U_phi))
!                        if (tmpreal < 0) CYCLE
!                        Qc  = 0.5*sqrt( tmp_U_amp**2 + tmp_V_amp**2 - (2*tmp_U_amp*tmp_V_amp*sin( tmp_V_phi - tmp_U_phi))  )

!                        tmpreal = tmp_U_amp**2 + tmp_V_amp**2 + (2*tmp_U_amp*tmp_V_amp*sin( tmp_V_phi - tmp_U_phi)) 
!                        if (tmpreal < 0) CYCLE
!                        Qac = 0.5*sqrt( tmp_U_amp**2 + tmp_V_amp**2 + (2*tmp_U_amp*tmp_V_amp*sin( tmp_V_phi - tmp_U_phi))  )


!                        gc  = atan2(  (  (  tmp_U_amp*sin( tmp_U_phi ) ) +  (tmp_V_amp*cos( tmp_V_phi)  ) )  ,  (  (tmp_U_amp*cos( tmp_U_phi ))  -  (tmp_V_amp*sin( tmp_V_phi ))  )  )
!                        gac = atan2(  (  ( -tmp_U_amp*sin( tmp_U_phi ) ) +  (tmp_V_amp*cos( tmp_V_phi)  ) )  ,  (  (tmp_U_amp*cos( tmp_U_phi ))  +  (tmp_V_amp*sin( tmp_V_phi ))  )  )

!                        !Pugh A3.2
!                        Phi_Ua = -0.5*(gac - gc)
!                        dir_Ua = 0.5*(gac + gc)  ! positive from x axis
!                        polarity = (Qac - Qc)/qmax



                        tmp_u_amp_mat(ji,jj) = tmp_u_amp
                        tmp_v_amp_mat(ji,jj) = tmp_v_amp
                        tmp_u_phi_mat(ji,jj) = tmp_u_phi
                        tmp_v_phi_mat(ji,jj) = tmp_v_phi


!                        a_u_mat(ji,jj) = a_u
!                        b_u_mat(ji,jj) = b_u
!                        a_v_mat(ji,jj) = a_v
!                        b_v_mat(ji,jj) = b_v

!                        qmax_mat(ji,jj) = qmax
!                        qmin_mat(ji,jj) = qmin

!                        ecc_mat(ji,jj) = ecc
!                        thetamax_mat(ji,jj) = thetamax
!                        thetamin_mat(ji,jj) = thetamin

!                        Qc_mat(ji,jj) = Qc
!                        Qac_mat(ji,jj) = Qac
!                        gc_mat(ji,jj) = gc
!                        gac_mat(ji,jj) = gac

!                        Phi_Ua_mat(ji,jj) = Phi_Ua
!                        dir_Ua_mat(ji,jj) = dir_Ua
!                        polarity_mat(ji,jj) = polarity

!                    ENDIF
                END DO
             END DO


!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_u_amp_t_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!               IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!               CALL iom_put( TRIM(tmp_name), tmp_u_amp_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_v_amp_t_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), tmp_v_amp_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_u_phi_t_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), tmp_u_phi_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_v_phi_t_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), tmp_v_phi_mat(:,:))
!            ENDIF



!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_a_u_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!               IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!               CALL iom_put( TRIM(tmp_name), a_u_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_a_v_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), a_v_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_b_u_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), b_u_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_b_v_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), b_v_mat(:,:))
!            ENDIF

!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_qmax_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), qmax_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_qmin_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), qmin_mat(:,:))
!            ENDIF

!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_ecc_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), ecc_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_thetamax_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), thetamax_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_thetamin_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), thetamin_mat(:,:))
!            ENDIF

!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_Qc_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), Qc_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_Qac_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), Qac_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_gc_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), gc_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_gac_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), gac_mat(:,:))
!            ENDIF


!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_Phi_Ua_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), Phi_Ua_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_dir_Ua_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), dir_Ua_mat(:,:))
!            ENDIF
!            tmp_name=TRIM(Wave(ntide_all(jh))%cname_tide)//'_polarity_uvbar'
!            IF( iom_use(TRIM(tmp_name)) ) THEN
!              IF(lwp) WRITE(numout,*) "harm_ana_out: iom_put: ",TRIM(tmp_name)
!              CALL iom_put( TRIM(tmp_name), polarity_mat(:,:))
!            ENDIF

            tmp_u_amp_mat(:,:) = 0.
            tmp_v_amp_mat(:,:) = 0.
            tmp_u_phi_mat(:,:) = 0.
            tmp_v_phi_mat(:,:) = 0.

!            a_u_mat(:,:) = 0.
!            b_u_mat(:,:) = 0.
!            a_v_mat(:,:) = 0.
!            b_v_mat(:,:) = 0.

!            qmax_mat(:,:) = 0.
!            qmin_mat(:,:) = 0.

!            ecc_mat(:,:) = 0
!            thetamax_mat(:,:) =0.
!            thetamin_mat(:,:) = 0.

!            Qc_mat(:,:) = 0.
!            Qac_mat(:,:) = 0.
!            gc_mat(:,:) = 0.
!            gac_mat(:,:) = 0.

!            Phi_Ua_mat(:,:) = 0.
!            dir_Ua_mat(:,:) = 0.
!            polarity_mat(:,:) = 0.


         END DO
         IF(lwp) WRITE(numout,*) "diaharm_fast: Finshed postprocessing 2d velocity tidal parameters"
      ENDIF

     CALL FLUSH(numout)

      IF (ln_diaharm_postproc_vel .AND. ln_ana_uv3d)  THEN
           IF(lwp) WRITE(numout,*) "diaharm_fast: Postprocess 3d velocity tidal parameters"
      ENDIF


     CALL FLUSH(numout)

! to output tidal parameters, u and v on t grid
!
!                                  !==  standard Cd  ==!
!         DO jj = 2, jpjm1
!            DO ji = 2, jpim1
!               imk = k_mk(ji,jj)    ! ocean bottom level at t-points
!               zut = un(ji,jj,imk) + un(ji-1,jj,imk)     ! 2 x velocity at t-point
!               zvt = vn(ji,jj,imk) + vn(ji,jj-1,imk)
!               !                                                           ! here pCd0 = mask*boost * drag
!               pCdU(ji,jj) = - pCd0(ji,jj) * SQRT(  0.25 * ( zut*zut + zvt*zvt ) + pke0  )
!            END DO
!         END DO



      IF (ln_diaharm_postproc_vel .AND. ln_ana_uvbar)  THEN

         DEALLOCATE(amp_u2d, amp_v2d, phi_u2d, phi_v2d )


         DEALLOCATE(tmp_u_amp_mat, tmp_v_amp_mat, tmp_u_phi_mat, tmp_v_phi_mat )
!         DEALLOCATE(a_u_mat, b_u_mat, a_v_mat, b_v_mat, qmax_mat, qmin_mat, ecc_mat )
!         DEALLOCATE(thetamax_mat, thetamin_mat, Qc_mat, Qac_mat, gc_mat, gac_mat )
!         DEALLOCATE(Phi_Ua_mat, dir_Ua_mat, polarity_mat )

      endif

      IF(lwp) WRITE(numout,*) "diaharm_fast: Deallocated 2d velocity tidal parameters"

      CALL FLUSH(numout)
!
   END SUBROUTINE harm_ana_out
!
   SUBROUTINE harm_rst_write(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To write out cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   restart files will be dated by default
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      INTEGER             ::   jh, j2d, j3d
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name

      !restart file
      DO j2d=1,nvar_2d
         CALL iom_rstput( kt, nitrst, numrow, 'Mean_'//TRIM(m_varName2d( m_posi_2d(j2d) )), g_cumul_var2D( 1, :, :, j2d ) )
         DO jh=1,nb_ana
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_cos', g_cumul_var2D( jh*2  , :, :, j2d ) )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_sin', g_cumul_var2D( jh*2+1, :, :, j2d ) )
         ENDDO
      ENDDO

      DO j3d=1,nvar_3d
         !JT CALL iom_rstput( kt, nitrst, numrow, 'Mean_'//TRIM(m_varName2d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         CALL iom_rstput( kt, nitrst, numrow, 'Mean_'//TRIM(m_varName3d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         DO jh=1,nb_ana
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_cos', g_cumul_var3D( jh*2  , :, :, :, j3d ) )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_sin', g_cumul_var3D( jh*2+1, :, :, :, j3d ) )
         ENDDO
      ENDDO

      IF(lwp) THEN
        IF( kt > 999999999 ) THEN ; WRITE(clkt, *       ) kt
        ELSE                      ; WRITE(clkt, '(i8.8)') kt
        ENDIF
        clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
        clpath = TRIM(cn_ocerst_outdir)
        IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
        IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for writing: ',TRIM(clpath)//clname

        WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
        OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
        WRITE(66) cc
        WRITE(66) anau
        WRITE(66) anav
        WRITE(66) anaf
        WRITE(66) fjulday_startharm
        CLOSE(66)
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   harm_rst_write: DONE '
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) anau
        WRITE(numout,*) anav
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_write

   SUBROUTINE harm_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To read in  cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name
      INTEGER             ::   jh, j2d, j3d

      IF( nit000 > 999999999 ) THEN ; WRITE(clkt, *       ) nit000-1
      ELSE                      ; WRITE(clkt, '(i8.8)') nit000-1
      ENDIF
      clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
      clpath = TRIM(cn_ocerst_outdir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'

      IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for reading: ',TRIM(clpath)//clname

      DO j2d=1,nvar_2d
         CALL iom_get( numror,jpdom_autoglo, 'Mean_'//TRIM(m_varName2d( m_posi_2d(j2d) )), g_cumul_var2D( 1, :, :, j2d ) )
         IF(lwp) WRITE(numout,*) "2D", j2d, m_posi_2d(j2d), m_varName2d( m_posi_2d(j2d) )
         DO jh=1,nb_ana
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_cos', g_cumul_var2D( jh*2  , :, :, j2d ) )
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_sin', g_cumul_var2D( jh*2+1, :, :, j2d ) )
         ENDDO
      ENDDO

      DO j3d=1,nvar_3d
         !JT  CALL iom_get( numror,jpdom_autoglo, 'Mean_'//TRIM(m_varName2d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         CALL iom_get( numror,jpdom_autoglo, 'Mean_'//TRIM(m_varName3d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         IF(lwp) WRITE(numout,*) "3D", j3d,  m_posi_3d(j3d), m_varName3d( m_posi_3d(j3d) )

         DO jh=1,nb_ana
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_cos', g_cumul_var3D( jh*2  , :, :, :, j3d ) )
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_sin', g_cumul_var3D( jh*2+1, :, :, :, j3d ) )
         ENDDO
      ENDDO

      WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
      OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
      READ(66) cc
      READ(66) anau
      READ(66) anav
      READ(66) anaf
      READ(66) fjulday_startharm
      CLOSE(66)

      IF(lwp) THEN
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   Checking anaf is correct'
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_read

   !!======================================================================

END MODULE diaharm_fast 
