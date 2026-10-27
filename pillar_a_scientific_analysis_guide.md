# GUIDE D'ANALYSE SCIENTIFIQUE ET PHYSIQUE DES FIGURES (PILIER A)

**Projet :** Régularisation AS-Flow / Simulation du Cyclone Kenneth (Catégorie 5)  
**Auteur :** Romeo DAGLAN  
**Référence Manuscrit :** *On Spurious Vorticity Generation in Vector-Coupled Metric Damping for Geophysical CFD*  
**Script Associé :** [`pillar_a_spurious_vorticity_2d.py`](./pillar_a_spurious_vorticity_2d.py)  

---

## 1. Contexte Physique & Matière Première Utilisée

Le benchmark 2D autonome a été calibré avec les **paramètres physiques réels du Cyclone tropical Kenneth** (Avril 2019, sud-ouest de l'Océan Indien), le cyclone le plus intense jamais observé dans le canal du Mozambique :
* **Vitesse maximale des vents ($V_{\max}$) :** $70.0\text{ m/s}$ ($\approx 140\text{ nœuds}$ ou $252\text{ km/h}$).
* **Rayon des vents maximaux ($R_{\max}$) :** $25.0\text{ km}$ (le rayon du mur de l'œil).
* **Domaine de calcul :** $200 \times 200\text{ km}$ centré sur l'œil, avec une résolution $\Delta x = \Delta y = 1\text{ km}$ (identique au domaine emboîté fin D03 de WRF).
* **Constantes AS-Flow du run Cloud :** $K_{\text{elastic}} = 0.05\text{ s}^{-1}$ et $K_{\text{drag}} = 0.05\text{ s}^{-1}$ (extraites du module Fortran `module_diffusion_asflow_v14.F`).

---

## 2. Analyse Panneau par Panneau de la Figure 1 (`figure1_radial_torque_generation.png`)

La Figure 1 présente la coupe radiale 1D (du centre $r=0$ jusqu'à $r=100\text{ km}$) illustrant la chaîne causale menant à l'émergence du couple parasite.

```text
[ (a) Cinématique du Cyclone Kenneth (v_theta & S_mag) ]
                           ↓
[ (b) Émergence de la Métrique J (Seuil J > 1 au mur)  ]
                           ↓
[ (c) Injection Spécieuse de Vorticité vs Tenseur      ]
```

### Panneau (a) : Cinématique du Vent Cyclonique et Taux de Cisaillement
* **Courbe bleue pleine ($v_\theta$ en $\text{m/s}$) :**  
  Représente le profil de vitesse azimutale du cyclone :
  * **Pour $r \le 25\text{ km}$ (dans l'œil) :** Le vent augmente linéairement de $0$ à $70\text{ m/s}$ ($v_\theta = V_{\max} \cdot r / R_{\max}$). Le fluide est en **rotation solide**.
  * **Au mur de l'œil ($r = 25\text{ km}$) :** Pic absolu des vents à $70\text{ m/s}$ ($140\text{ kt}$).
  * **Pour $r > 25\text{ km}$ (flanc externe) :** Décroissance cyclonique en $v_\theta \propto r^{-0.6}$.
* **Courbe violette pointillée ($S_{\text{mag}}$ en $10^{-3}\text{ s}^{-1}$) :**  
  Taux de déformation $S_{\text{mag}} = \sqrt{2 S_{ij} S_{ij}}$ :
  * Dans l'œil ($r < 25\text{ km}$), $S_{\text{mag}} = 0$ car la rotation en bloc ne déforme pas le fluide.
  * Exactement au franchissement du mur de l'œil ($r = 25\text{ km}$), le cisaillement radial explose pour atteindre un pic violent à **$4.26 \times 10^{-3}\text{ s}^{-1}$**.

### Panneau (b) : Émergence de la Variable Métrique $J(r)$
* **Courbe orange ($J(r)$) :**  
  Gouvernée par l'équilibre de transport-relaxation d'AS-Flow : $J = 1 + S_{\text{mag}} / K_{\text{elastic}}$.
* **Ligne grise pointillée ($J = 1.0$) :**  
  Représente le seuil critique d'activation du frein élastique.
* **Interprétation :**  
  Dans l'œil, $J = 1.0$ (aucun freinage). Mais au mur de l'œil ($r \ge 25\text{ km}$), $J$ grimpe brutalement pour culminer à **$1.085$**, puis décroît lentement vers l'extérieur.  
  Cette variation spatiale crée un **gradient spatial aigu $\nabla J = \frac{\partial J}{\partial r} \mathbf{e}_r$** dirigé radialement.

### Panneau (c) : Le Terme Croisé Parasite vs la Solution Tensorielle
* **Courbe rouge épaisse (Terme croisé $-K(\nabla J \times \mathbf{u})_z$) :**  
  Sur le flanc externe du mur de l'œil, le gradient radial est négatif ($\partial_r J < 0$) et le vent est cyclonique ($v_\theta > 0$).  
  Par conséquent :
  $$-K \theta(J-1) (\nabla J \times \mathbf{u})_z = -K \left(\frac{\partial J}{\partial r}\right) v_\theta > 0 \quad \text{(STRICTEMENT POSITIF !)}$$
  Il culmine à **$+1.64 \times 10^{-4}\text{ s}^{-2}$**.
* **Courbe rouge en tirets (Rotationnel total $(\nabla \times \mathbf{F}_{\text{drag}})_z$) :**  
  Inclut le frottement stable $-K(J-1)\omega_z$. Même avec ce frottement, le couple parasite domine largement :
  $$\frac{\text{Couple Parasite}}{\text{Amortissement Physique}} = \mathbf{24.93}$$
  Le modèle injecte 25 fois plus de vorticité artificielle qu'il n'en dissipe !
* **Courbe verte (Régularisation Tensorielle de Cauchy $(\nabla \times \mathbf{F}_{\text{tensor}})_z$) :**  
  Représente la formulation tensorielle symétrique $\boldsymbol{\tau} = 2 \nu(J)\mathbf{S}$.  
  Elle reste parfaitement collée à zéro ($\le 10^{-8}\text{ s}^{-2}$). Le couple parasite est **strictement annihilé**, garantissant un amortissement stable.

---

## 3. Analyse de la Figure 2 : Cartes 2D Comparatives (`figure2_2d_vorticity_comparison.png`)

La Figure 2 présente la projection horizontale dans le plan $(X, Y)$ du cyclone :
* **Panneau de gauche (Modèle Vectoriel Intuitif) :**  
  Montre une couronne d'injection de vorticité parasite (zone rouge/bleue) calquée exactement sur le cercle en tirets noirs du mur de l'œil ($R_{\max} = 25\text{ km}$).  
  Cette couronne agit comme une turbine artificielle qui injecte continuellement de l'énergie de rotation dans le vortex.
* **Panneau de droite (Régularisation Tensorielle de Cauchy) :**  
  La couronne parasite a totalement disparu. Le champ est neutre et purement dissipatif sur l'ensemble du domaine de $200 \times 200\text{ km}$.

---

## 4. Analyse de la Figure 3 : Évolution Temporelle et Spin-Up (`figure3_temporal_enstrophy_evolution.png`)

Ce graphique montre ce qui se passe lorsque l'on intègre l'équation de la vorticité au cours du temps $\partial_t \omega_z = (\nabla \times \mathbf{F})_z$ :
* **Courbe rouge (Modèle Vectoriel) :**  
  En 3 minutes de dynamique, la vorticité au mur de l'œil explose de $1.0 \times 10^{-3}\text{ s}^{-1}$ à près de $3.7 \times 10^{-3}\text{ s}^{-1}$. C'est la signature exacte du **spin-up parasite** qui a provoqué l'explosion des vitesses verticales ($W = 45.12\text{ m/s}$) dans les simulations supercalculateur de Kenneth.
* **Courbe noire en tirets (Référence sans couplage) :** Décroissance visqueuse lente.
* **Courbe verte (Régularisation Tensorielle de Cauchy) :** La vorticité reste parfaitement maîtrisée et stable, confirmant le respect du théorème de dissipation de l'enstrophie de Beale-Kato-Majda.
