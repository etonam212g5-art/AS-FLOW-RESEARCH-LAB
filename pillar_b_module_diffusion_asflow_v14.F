MODULE module_diffusion_asflow_v14

! ====================================================================
! AS-FLOW V14 — SEPARATION COMPLETE Theta/omega
! Ref: AS_FLOW_THEORY.md (EDP-J, EDP-Theta)
!
! V14 vs V13-A:
!   - Theta est une variable prognostique separee (3 composantes)
!   - Beltrami COMPLET applique sur Theta UNIQUEMENT (pas sur u,v,w)
!   - EDP de J: S*J - k*J*(J-1) + nu*Lap(J)  (meme nu)
!   - Biot-Savart diagnostique BS(Theta,0) et BS(omega,0) aux pas output
! ====================================================================

CONTAINS

! ====================================================================
! SOUS-ROUTINE 1 : INITIALISATION (appelee une fois a itimestep=1)
! Theta_0 = curl(u,v,w) / max(J, 1.0)
! ====================================================================
   SUBROUTINE asflow_v14_init(u, v, w, j_asflow, &
                              theta_x, theta_y, theta_z, &
                              dx, dy, rdz, &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte

      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: u, v, w
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: rdz
      REAL, INTENT(IN) :: dx, dy

      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           j_asflow

      REAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           theta_x, theta_y, theta_z

      INTEGER :: i, j, k
      REAL :: rdx, rdy, j_val
      REAL :: curl_x, curl_y, curl_z
      REAL :: dudy, dudz, dvdx, dvdz, dwdx, dwdy

      rdx = 1.0 / dx
      rdy = 1.0 / dy

      ! Initialiser a zero
      theta_x = 0.0
      theta_y = 0.0
      theta_z = 0.0

      DO j = MAX(jts, jds+1), MIN(jte, jde-1)
         DO k = MAX(kts, kds+1), MIN(kte, kde-1)
            DO i = MAX(its, ids+1), MIN(ite, ide-1)
               
               ! J démarre strictement à 1.0 partout
               j_asflow(i,k,j) = 1.0

               ! Curl(u,v,w) par differences centrees
               dwdy = (w(i,k,j+1) - w(i,k,j-1)) * 0.5 * rdy
               dvdz = (v(i,k+1,j) - v(i,k-1,j)) * 0.5 * rdz(i,k,j)
               curl_x = dwdy - dvdz

               dudz = (u(i,k+1,j) - u(i,k-1,j)) * 0.5 * rdz(i,k,j)
               dwdx = (w(i+1,k,j) - w(i-1,k,j)) * 0.5 * rdx
               curl_y = dudz - dwdx

               dvdx = (v(i+1,k,j) - v(i-1,k,j)) * 0.5 * rdx
               dudy = (u(i,k,j+1) - u(i,k,j-1)) * 0.5 * rdy
               curl_z = dvdx - dudy

               ! Theta = curl / max(J, 1.0)
               j_val = MAX(j_asflow(i,k,j), 1.0)
               theta_x(i,k,j) = curl_x / j_val
               theta_y(i,k,j) = curl_y / j_val
               theta_z(i,k,j) = curl_z / j_val

            END DO
         END DO
      END DO

   END SUBROUTINE asflow_v14_init

! ====================================================================
! SOUS-ROUTINE 2 : CALCUL PRINCIPAL (appelee a chaque pas de temps)
! - EDP de J : source stretching + Kevlar + diffusion (khdif)
! - Calcul du Frein de Kevlar : ajoute directement a ru_tendf
! - Diagnostic pur de Theta : Theta = curl(u) / J
! ====================================================================
   SUBROUTINE asflow_v14_calc( &
        u, v, w, j_asflow, mu, &
        j_tend, &
        theta_x, theta_y, theta_z, &
        ru_tendf, rv_tendf, rw_tendf, &
        s_mag_diag, &
        khdif, &
        dx, dy, rdz, rdzw, &
        dt, &
        ids, ide, jds, jde, kds, kde, &
        ims, ime, jms, jme, kms, kme, &
        its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte

      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           u, v, w, j_asflow
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: rdz, rdzw
      REAL, INTENT(IN), DIMENSION(ims:ime, jms:jme) :: mu
      REAL, INTENT(IN) :: dx, dy, khdif, dt

      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           j_tend, ru_tendf, rv_tendf, rw_tendf

      REAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           theta_x, theta_y, theta_z, s_mag_diag

      ! Local variables
      INTEGER :: i, j, k
      REAL :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
      REAL :: S_11, S_22, S_33, S_12, S_13, S_23
      REAL :: S_mag, S_mag_eff
      REAL :: j_val, mu_val, rdx, rdy
      REAL :: diff_J, kevlar_factor
      REAL :: curl_x, curl_y, curl_z

      ! AS-Flow constants (AS_FLOW_THEORY.md, Section VIII)
      REAL, PARAMETER :: k_elastic = 0.05
      REAL, PARAMETER :: asflow_c_stretch = 0.25

      rdx = 1.0 / dx
      rdy = 1.0 / dy

      theta_x = 0.0
      theta_y = 0.0
      theta_z = 0.0
      s_mag_diag = 0.0

      DO j = MAX(jts, jds+1), MIN(jte, jde-1)
         DO k = MAX(kts, kds+1), MIN(kte, kde-1)
            DO i = MAX(its, ids+1), MIN(ite, ide-1)

               ! ============================================
               ! ETAPE 1 : Tenseur de deformation |S_ij| et Curl
               ! ============================================
               dudx = (u(i+1,k,j) - u(i-1,k,j)) * 0.5 * rdx
               dudy = (u(i,k,j+1) - u(i,k,j-1)) * 0.5 * rdy
               dudz = (u(i,k+1,j) - u(i,k-1,j)) * 0.5 * rdz(i,k,j)

               dvdx = (v(i+1,k,j) - v(i-1,k,j)) * 0.5 * rdx
               dvdy = (v(i,k,j+1) - v(i,k,j-1)) * 0.5 * rdy
               dvdz = (v(i,k+1,j) - v(i,k-1,j)) * 0.5 * rdz(i,k,j)

               dwdx = (w(i+1,k,j) - w(i-1,k,j)) * 0.5 * rdx
               dwdy = (w(i,k,j+1) - w(i,k,j-1)) * 0.5 * rdy
               dwdz = (w(i,k+1,j) - w(i,k-1,j)) * 0.5 * rdz(i,k,j)

               S_11 = dudx
               S_22 = dvdy
               S_33 = dwdz
               S_12 = 0.5 * (dudy + dvdx)
               S_13 = 0.5 * (dudz + dwdx)
               S_23 = 0.5 * (dvdz + dwdy)

               S_mag = SQRT(2.0 * (S_11**2 + S_22**2 + S_33**2 &
                          + 2.0 * S_12**2 &
                          + 2.0 * S_13**2 &
                          + 2.0 * S_23**2))
               S_mag_eff = S_mag * asflow_c_stretch

               ! Vorticite (Curl) pour diagnostic
               curl_x = dwdy - dvdz
               curl_y = dudz - dwdx
               curl_z = dvdx - dudy

               ! ============================================
               ! ETAPE 2 : EDP de J
               ! dJ/dt = C*S_mag*J - k_e*J*(J-1) + khdif*Lap(J)
               ! ============================================
               j_val = MAX(j_asflow(i,k,j), 1.0)
               mu_val = mu(i,j)

               diff_J = khdif * ( &
                    (j_asflow(i+1,k,j) - 2.*j_asflow(i,k,j) &
                     + j_asflow(i-1,k,j)) * rdx**2 + &
                    (j_asflow(i,k,j+1) - 2.*j_asflow(i,k,j) &
                     + j_asflow(i,k,j-1)) * rdy**2 )

               j_tend(i,k,j) = j_tend(i,k,j) + mu_val * ( &
                    (j_val * S_mag_eff - k_elastic * j_val * (j_val - 1.0)) / &
                    (1.0 + dt * k_elastic * j_val) + &
                    diff_J )

               ! ============================================
               ! ETAPE 3 : FREIN DE KEVLAR (-k_e * (J-1) * u)
               ! Les tendances WRF sont couplees avec mu(i,j)
               ! ============================================
               
               ! Le frein ne s'applique que si J > 1 (cap de stabilite numerique a 0.05/s)
               kevlar_factor = -MIN(k_elastic * MAX(0.0, j_val - 1.0), 0.05)
               
               ru_tendf(i,k,j) = ru_tendf(i,k,j) + mu_val * kevlar_factor * u(i,k,j)
               rv_tendf(i,k,j) = rv_tendf(i,k,j) + mu_val * kevlar_factor * v(i,k,j)
               rw_tendf(i,k,j) = rw_tendf(i,k,j) + mu_val * kevlar_factor * w(i,k,j)

               ! ============================================
               ! ETAPE 4 : Diagnostic pur de Theta
               ! Theta = Curl(u) / J
               ! ============================================
               theta_x(i,k,j) = curl_x / j_val
               theta_y(i,k,j) = curl_y / j_val
               theta_z(i,k,j) = curl_z / j_val
               s_mag_diag(i,k,j) = S_mag

            END DO
         END DO
      END DO

   END SUBROUTINE asflow_v14_calc

! ====================================================================
! SOUS-ROUTINE 3 : BIOT-SAVART DIAGNOSTIQUE (aux pas d'ecriture)
! Resout nabla^2(psi) = -source par SOR, puis u = curl(psi)
! Appelee pour BS(Theta) et BS(J*Theta)
! ====================================================================
   SUBROUTINE asflow_v14_biot_savart( &
        theta_x, theta_y, theta_z, j_asflow, &
        u_theta_diag, v_theta_diag, &
        u_omega_diag, v_omega_diag, &
        dx, dy, rdz, &
        ids, ide, jds, jde, kds, kde, &
        ims, ime, jms, jme, kms, kme, &
        its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte

      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           theta_x, theta_y, theta_z, j_asflow
      REAL, INTENT(IN) :: dx, dy
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: rdz

      REAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
           u_theta_diag, v_theta_diag, u_omega_diag, v_omega_diag

      ! Locales pour le solveur de Poisson
      REAL, ALLOCATABLE :: psi_x(:,:,:), psi_y(:,:,:), psi_z(:,:,:)
      REAL, ALLOCATABLE :: rhs_x(:,:,:), rhs_y(:,:,:), rhs_z(:,:,:)
      REAL :: rdx2, rdy2, rdz2_local
      REAL :: coeff, psi_new, j_val
      INTEGER :: i, j, k, iter
      INTEGER :: i_s, i_e, j_s, j_e, k_s, k_e

      INTEGER, PARAMETER :: max_iter = 200
      REAL, PARAMETER :: omega_sor = 1.6

      i_s = MAX(its, ids+1) ; i_e = MIN(ite, ide-1)
      j_s = MAX(jts, jds+1) ; j_e = MIN(jte, jde-1)
      k_s = MAX(kts, kds+1) ; k_e = MIN(kte, kde-1)

      ALLOCATE(psi_x(ims:ime,kms:kme,jms:jme))
      ALLOCATE(psi_y(ims:ime,kms:kme,jms:jme))
      ALLOCATE(psi_z(ims:ime,kms:kme,jms:jme))
      ALLOCATE(rhs_x(ims:ime,kms:kme,jms:jme))
      ALLOCATE(rhs_y(ims:ime,kms:kme,jms:jme))
      ALLOCATE(rhs_z(ims:ime,kms:kme,jms:jme))

      rdx2 = 1.0 / (dx * dx)
      rdy2 = 1.0 / (dy * dy)

      ! Initialize outputs
      u_theta_diag = 0.0
      v_theta_diag = 0.0
      u_omega_diag = 0.0
      v_omega_diag = 0.0

      ! ==============================
      ! PARTIE A : BS(Theta) -> u_theta
      ! ==============================
      psi_x = 0.0 ; psi_y = 0.0 ; psi_z = 0.0
      rhs_x = 0.0 ; rhs_y = 0.0 ; rhs_z = 0.0

      DO j = j_s, j_e
        DO k = k_s, k_e
          DO i = i_s, i_e
            rhs_x(i,k,j) = -theta_x(i,k,j)
            rhs_y(i,k,j) = -theta_y(i,k,j)
            rhs_z(i,k,j) = -theta_z(i,k,j)
          END DO
        END DO
      END DO

      DO iter = 1, max_iter
        DO j = j_s, j_e
          DO k = k_s, k_e
            DO i = i_s, i_e
              rdz2_local = rdz(i,k,j) * rdz(i,k,j)
              coeff = 1.0 / (2.0*(rdx2 + rdy2 + rdz2_local))

              psi_new = coeff * ( &
                   rdx2*(psi_x(i+1,k,j)+psi_x(i-1,k,j)) + &
                   rdy2*(psi_x(i,k,j+1)+psi_x(i,k,j-1)) + &
                   rdz2_local*(psi_x(i,k+1,j)+psi_x(i,k-1,j)) &
                   - rhs_x(i,k,j) )
              psi_x(i,k,j) = psi_x(i,k,j) + &
                   omega_sor*(psi_new - psi_x(i,k,j))

              psi_new = coeff * ( &
                   rdx2*(psi_y(i+1,k,j)+psi_y(i-1,k,j)) + &
                   rdy2*(psi_y(i,k,j+1)+psi_y(i,k,j-1)) + &
                   rdz2_local*(psi_y(i,k+1,j)+psi_y(i,k-1,j)) &
                   - rhs_y(i,k,j) )
              psi_y(i,k,j) = psi_y(i,k,j) + &
                   omega_sor*(psi_new - psi_y(i,k,j))

              psi_new = coeff * ( &
                   rdx2*(psi_z(i+1,k,j)+psi_z(i-1,k,j)) + &
                   rdy2*(psi_z(i,k,j+1)+psi_z(i,k,j-1)) + &
                   rdz2_local*(psi_z(i,k+1,j)+psi_z(i,k-1,j)) &
                   - rhs_z(i,k,j) )
              psi_z(i,k,j) = psi_z(i,k,j) + &
                   omega_sor*(psi_new - psi_z(i,k,j))
            END DO
          END DO
        END DO
      END DO

      ! u_theta = curl(psi)
      DO j = j_s, j_e
        DO k = k_s, k_e
          DO i = i_s, i_e
            u_theta_diag(i,k,j) = &
                 (psi_z(i,k,j+1)-psi_z(i,k,j-1))*0.5/dy - &
                 (psi_y(i,k+1,j)-psi_y(i,k-1,j))*0.5*rdz(i,k,j)
            v_theta_diag(i,k,j) = &
                 (psi_x(i,k+1,j)-psi_x(i,k-1,j))*0.5*rdz(i,k,j) - &
                 (psi_z(i+1,k,j)-psi_z(i-1,k,j))*0.5/dx
          END DO
        END DO
      END DO

      ! ==============================
      ! PARTIE B : BS(J*Theta) -> u_omega
      ! ==============================
      psi_x = 0.0 ; psi_y = 0.0 ; psi_z = 0.0

      DO j = j_s, j_e
        DO k = k_s, k_e
          DO i = i_s, i_e
            j_val = MAX(j_asflow(i,k,j), 1.0)
            rhs_x(i,k,j) = -j_val * theta_x(i,k,j)
            rhs_y(i,k,j) = -j_val * theta_y(i,k,j)
            rhs_z(i,k,j) = -j_val * theta_z(i,k,j)
          END DO
        END DO
      END DO

      DO iter = 1, max_iter
        DO j = j_s, j_e
          DO k = k_s, k_e
            DO i = i_s, i_e
              rdz2_local = rdz(i,k,j) * rdz(i,k,j)
              coeff = 1.0 / (2.0*(rdx2 + rdy2 + rdz2_local))

              psi_new = coeff * ( &
                   rdx2*(psi_x(i+1,k,j)+psi_x(i-1,k,j)) + &
                   rdy2*(psi_x(i,k,j+1)+psi_x(i,k,j-1)) + &
                   rdz2_local*(psi_x(i,k+1,j)+psi_x(i,k-1,j)) &
                   - rhs_x(i,k,j) )
              psi_x(i,k,j) = psi_x(i,k,j) + &
                   omega_sor*(psi_new - psi_x(i,k,j))

              psi_new = coeff * ( &
                   rdx2*(psi_y(i+1,k,j)+psi_y(i-1,k,j)) + &
                   rdy2*(psi_y(i,k,j+1)+psi_y(i,k,j-1)) + &
                   rdz2_local*(psi_y(i,k+1,j)+psi_y(i,k-1,j)) &
                   - rhs_y(i,k,j) )
              psi_y(i,k,j) = psi_y(i,k,j) + &
                   omega_sor*(psi_new - psi_y(i,k,j))

              psi_new = coeff * ( &
                   rdx2*(psi_z(i+1,k,j)+psi_z(i-1,k,j)) + &
                   rdy2*(psi_z(i,k,j+1)+psi_z(i,k,j-1)) + &
                   rdz2_local*(psi_z(i,k+1,j)+psi_z(i,k-1,j)) &
                   - rhs_z(i,k,j) )
              psi_z(i,k,j) = psi_z(i,k,j) + &
                   omega_sor*(psi_new - psi_z(i,k,j))
            END DO
          END DO
        END DO
      END DO

      DO j = j_s, j_e
        DO k = k_s, k_e
          DO i = i_s, i_e
            u_omega_diag(i,k,j) = &
                 (psi_z(i,k,j+1)-psi_z(i,k,j-1))*0.5/dy - &
                 (psi_y(i,k+1,j)-psi_y(i,k-1,j))*0.5*rdz(i,k,j)
            v_omega_diag(i,k,j) = &
                 (psi_x(i,k+1,j)-psi_x(i,k-1,j))*0.5*rdz(i,k,j) - &
                 (psi_z(i+1,k,j)-psi_z(i-1,k,j))*0.5/dx
          END DO
        END DO
      END DO

      DEALLOCATE(psi_x, psi_y, psi_z, rhs_x, rhs_y, rhs_z)

   END SUBROUTINE asflow_v14_biot_savart

END MODULE module_diffusion_asflow_v14
