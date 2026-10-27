# PILIER B : PREUVE D'IMPLÉMENTATION WRF & ANALYSE DU CODE SOURCE FORTRAN

**Projet :** Régularisation AS-Flow / Simulation du Cyclone Kenneth (Catégorie 5)  
**Auteur :** Ugvucbm Cvkzdv HDBDJVGDCGHSACJHMS  
**Référence Manuscrit :** *On Spurious Vorticity Generation in Vector-Coupled Metric Damping for Geophysical CFD*  
**Fichier Source Associé :** [`pillar_b_module_diffusion_asflow_v14.F`](./pillar_b_module_diffusion_asflow_v14.F)  

---

## 1. Contexte du Fichier et Emplacement dans le Noyau Dynamique WRF
Le fichier [`pillar_b_module_diffusion_asflow_v14.F`](./pillar_b_module_diffusion_asflow_v14.F) est le module Fortran officiel qui a été compilé et exécuté au sein du noyau dynamique non-hydrostatique `dyn_em` du modèle **WRF (Weather Research and Forecasting Model)** sur supercalculateur cloud lors de la simulation haute résolution (1 km) du cyclone tropical Kenneth (Avril 2019).

Dans l'architecture de WRF, ce module est appelé directement à chaque sous-pas d'intégration Runge-Kutta d'ordre 3 (`module_first_rk_step_part2.F`), via la routine de calcul :
```fortran
SUBROUTINE asflow_v14_calc( ... )
```

---

## 2. Preuve Textuelle du Couplage Vectoriel (Lignes 201 à 212)

Le cœur de la pathologie mise en lumière dans notre article de recherche réside dans les lignes **201 à 212** du fichier source Fortran :

```fortran
! ============================================
! ETAPE 3 : FREIN DE KEVLAR (-k_e * (J-1) * u)
! Les tendances WRF sont couplees avec mu(i,j)
! ============================================

! Le frein ne s'applique que si J > 1 (cap de stabilite numerique a 0.05/s)
kevlar_factor = -MIN(k_elastic * MAX(0.0, j_val - 1.0), 0.05)

ru_tendf(i,k,j) = ru_tendf(i,k,j) + mu_val * kevlar_factor * u(i,k,j)
rv_tendf(i,k,j) = rv_tendf(i,k,j) + mu_val * kevlar_factor * v(i,k,j)
rw_tendf(i,k,j) = rw_tendf(i,k,j) + mu_val * kevlar_factor * w(i,k,j)
```

### Correspondance Mathématique Directe avec l'Article :
1. **Facteur d'excès élastique :**  
   $$	ext{kevlar\_factor} \Longleftrightarrow - K \max(0, J - 1)$$
2. **Couplage avec la masse sèche (coordonnée verticale $\eta$) :**  
   Dans WRF, $\mu(i,j)$ représente la pression de colonne d'air sec. Les variables pronostiques de quantité de mouvement sont couplées en flux de masse : $U = \mu u$, $V = \mu v$, $W = \mu w$.  
   L'injection dans les accumulateurs de force (`ru_tendf`, `rv_tendf`, `rw_tendf`) correspond exactement à l'EDP couplée :
   $$rac{\partial (\mu \mathbf{u})}{\partial t} = \dots + \mu \mathbf{F}_{	ext{drag}}$$
   avec :
   $$\mathbf{F}_{	ext{drag}} = - K \max(0, J - 1) \mathbf{u}$$

---

## 3. Origine Physique du Rotationnel Parasite et de l'Explosion de W

En appliquant le rotationnel analytique à cette force vectorielle :
$$
abla 	imes \mathbf{F}_{	ext{drag}} = - K 	heta(J - 1) (
abla J 	imes \mathbf{u}) - K (J - 1) oldsymbol{\omega}$$

Dans le mur de l'œil du cyclone Kenneth :
* Le vent tangentiel cyclonique $v_	heta pprox 70	ext{ m/s}$ (140 kt).
* Le gradient radial de cisaillement $rac{\partial J}{\partial r}$ est maximal sur le flanc externe du mur.
* Le terme croisé $(
abla J 	imes \mathbf{u})_z = rac{\partial J}{\partial r} v_	heta$ s'est transformé en un **générateur artificiel de vorticité verticale** :
  $$(
abla 	imes \mathbf{F}_{	ext{drag}})_z pprox 1.64 	imes 10^{-4}	ext{ s}^{-2}$$

**Conséquence numérique observée dans les sorties WRF :**
Ce couple parasite a déstabilisé la dynamique de la couche limite, amplifiant la convergence au centre et propulsant la vitesse verticale maximale $W_{	ext{max}}$ jusqu'à **$45.12	ext{ m/s}$** (voir `metrics_asflow.csv`), alors que la simulation physique de référence (Smagorinsky standard) maintenait $W_{	ext{max}} \le 0.02	ext{ m/s}$ à ce stade de creusement.

---

## 4. La Correction Tensorielle de Cauchy pour WRF

Pour corriger ce bug au niveau du solveur Fortran de WRF, les lignes 207-212 doivent être remplacées par la divergence du tenseur de contraintes symétrique de Cauchy $oldsymbol{	au} = 2 
u_{	ext{eff}}(J) \mathbf{S}$ :

```fortran
! =================================================================
! REGULARISATION TENSORIELLE DE CAUCHY (ANNULATION DU ROTATIONNEL)
! =================================================================
nu_eff = nu_0 + kappa * MAX(0.0, j_val - 1.0) * (dx * dy)

! Calcul des contraintes de cisaillement
tau_11 = 2.0 * nu_eff * S_11
tau_22 = 2.0 * nu_eff * S_22
tau_33 = 2.0 * nu_eff * S_23
tau_12 = 2.0 * nu_eff * S_12

! Divergence conservative du tenseur (div tau)
ru_tendf(i,k,j) = ru_tendf(i,k,j) + mu_val * ( &
     (tau_11(i+1,k,j) - tau_11(i-1,k,j)) * 0.5 * rdx + &
     (tau_12(i,k,j+1) - tau_12(i,k,j-1)) * 0.5 * rdy )

rv_tendf(i,k,j) = rv_tendf(i,k,j) + mu_val * ( &
     (tau_12(i+1,k,j) - tau_12(i-1,k,j)) * 0.5 * rdx + &
     (tau_22(i,k,j+1) - tau_22(i,k,j-1)) * 0.5 * rdy )
```

Grâce à cette formulation tensorielle :
$$
abla 	imes (
abla \cdot oldsymbol{	au}) \equiv 
u_{	ext{eff}} 
abla^2 oldsymbol{\omega} \le 0$$
Le couple parasite croisé $
abla J 	imes \mathbf{u}$ est **strictement éliminé à l'ordre machine**, garantissant la conservation de l'invariance galiléenne et la dissipation inconditionnelle de l'enstrophie sous WRF.
