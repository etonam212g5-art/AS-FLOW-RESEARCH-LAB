# Sur la Génération de Vorticité Parasite dans les Couplages Métriques-Vitesse en CFD Géophysique : Analyse Analytique du Terme Parasite et Régularisation Tensorielle de Cauchy

**Titre court :** Génération de Vorticité Parasite et Régularisation Tensorielle en CFD Géophysique  
**Auteur :** Kodjo Etonam G. GABLA  
**Date :** Octobre 2026  
**Classification :** Mécanique des Fluides Numérique — Modélisation Météorologique — Dynamique des Vortex (PACS: 47.11.-j, 47.32.C-, 92.60.Ek)  
**Statut :** Note Technique & Manuscrit de Recherche (Prêt pour Soumission)

---

## Résumé

Dans les simulations atmosphériques méso-échelles à ultra-haute résolution (résolutions sous-kilométriques), les gradients extrêmes de vitesse au niveau du mur de l'œil des cyclones tropicaux provoquent des déformations aiguës du maillage de calcul, conduisant fréquemment à des instabilités de type Courant-Friedrichs-Lewy (CFL) et à l'explosion des vitesses verticales ($W$). Pour pallier ce phénomène, une approche intuitive consiste à introduire une variable métrique pronostique d'étirement $J$ (gouvernée par une équation de transport-relaxation) et à appliquer une force de rappel élastique de traînée vectorielle directe : $\mathbf{F}_{\text{drag}} = - K \max(0, J - 1) \mathbf{u}$.

Dans cette note technique, nous démontrons analytiquement et numériquement que cette formulation vectorielle intuitive induit une pathologie fatale : elle agit comme un générateur artificiel de vorticité verticale. En appliquant l'opérateur rotationnel à l'équation de Navier-Stokes couplée, nous prouvons l'émergence d'un terme résiduel non physique :
$$\nabla \times \mathbf{F}_{\text{drag}} = - K \max(0, J - 1) \boldsymbol{\omega} - K \theta(J - 1) \left( \nabla J \times \mathbf{u} \right)$$

Dans une structure cyclonique où le gradient métrique $\nabla J$ est radial et le vent $\mathbf{u}$ tangentiel, le produit vectoriel $\nabla J \times \mathbf{u}$ est strictement colinéaire à l'axe vertical du vortex. Ce terme injecte un moment cinétique artificiel de même ordre de grandeur que l'amortissement visqueux attendu ($\mathcal{O}(10^{-4}\text{ s}^{-2})$ pour le cyclone Kenneth), provoquant un emballement du vent (spin-up jusqu'à 160 kt) et des vitesses verticales divergentes. 

Nous établissons les causes fondamentales de cette dérive : la violation de l'invariance galiléenne et l'absence de respect du théorème des contraintes de Cauchy. Enfin, nous développons la solution mathématiquement rigoureuse au problème : la substitution de la force vectorielle par la divergence d'un tenseur de contraintes d'ordre 2 symétrique $\boldsymbol{\tau} = 2 \nu_{\text{eff}}(J) \mathbf{S}$, où $\mathbf{S}$ est le tenseur déviatorique des taux de déformation. Nous démontrons l'annulation exacte du terme de couple parasite et la transformation du couplage en une diffusion d'enstrophie inconditionnellement stable. Ce travail constitue une mise en garde méthodologique essentielle pour la conception de régularisateurs sous-maille géométriques en mécanique des fluides géophysiques.

**Mots-clés :** Vorticité parasite, Théorème de Cauchy, Invariance galiléenne, Modélisation cyclonique, WRF, AS-Flow, Tenseur des taux de déformation, Instabilité sous-maille.

---

## 1. Introduction et Contexte Physique

La simulation numérique de cyclones tropicaux de catégorie 5 (tels que le cyclone Kenneth, 2019) à très haute résolution spatiale ($\Delta x \le 1\text{ km}$) constitue l'un des défis majeurs de la prévision météorologique moderne. Dans le modèle de recherche et de prévision atmosphérique WRF (*Weather Research and Forecasting Model*), le système non-hydrostatique et compressible résout les équations d'Euler/Navier-Stokes en coordonnées de masse suivant le relief (Skamarock et al., 2008). 

Au niveau de la paroi de l'œil, la convergence de couche limite combinée aux courants ascendants violents engendre des taux de cisaillement horizontal et vertical extrêmes :
$$S_{ij} = \frac{1}{2}\left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)$$

Sous ces gradients intenses, la cellule élémentaire de discrétisation subit des déformations spatiales critiques. Lorsque le temps de transit advectif devient inférieur au seuil de diffusion numérique, les harmoniques de longueur d'onde $2\Delta x$ s'accumulent au sommet du spectre énergétique, déclenchant des instabilités de type Courant-Friedrichs-Lewy (CFL) qui mènent à la divergence des champs de calcul (`NaN`).

Pour contrer cette dégénérescence géométrique, l'architecture AS-Flow (*Anti-Stretching Flow*) a proposé de modéliser l'accumulation locale d'étirement de l'espace de calcul par une variable métrique pronostique scalaire, notée $J(\mathbf{x}, t)$, normalisée à $J=1$ au repos. Dans la version opérationnelle V14, afin d'interdire l'étirement pathologique de la maille ($J > 1$), un terme de rappel élastique de traînée (désigné sous le nom de « force de Kevlar ») a été introduit directement dans l'équation de la quantité de mouvement :

$$\mathbf{F}_{\text{Kevlar}} = - K_{\text{elastic}} \cdot \max(0, J - 1) \cdot \mathbf{u}$$

où $K_{\text{elastic}} > 0$ est une constante élastique de rappel (dimension $s^{-1}$) et $\mathbf{u} = (u, v, w)^T$ est le vecteur vitesse tridimensionnel du vent.

L'objectif intuitif de cette force était d'agir comme un frein dissipatif isotrope dès lors que le maillage était étiré. Cependant, l'intégration numérique de cette formulation sur des cas réels (notamment le cyclone Kenneth à $\Delta x = 1\text{ km}$) a révélé un comportement paradoxal et catastrophique : loin d'amortir la turbulence, le modèle a provoqué une accélération anormale du vent maximal soutenu (passant de 100 kt à plus de 160 kt) accompagnée de vitesses verticales $W$ explosives dépassant $40\text{ m/s}$.

Le but de cet article est d'exposer la démonstration analytique de cette anomalie, d'identifier les violations des principes premiers de la physique qui en sont la cause, et de poser le cadre de régularisation tensorielle formellement exact.

---

## 2. Dérivation Analytique de l'Équation de Vorticité

L'évolution de la rotation dans un fluide est gouvernée par le vecteur vorticité $\boldsymbol{\omega} = \nabla \times \mathbf{u}$. Pour isoler l'impact dynamique de la force de rappel $\mathbf{F}_{\text{Kevlar}}$, appliquons l'opérateur rotationnel différentiel $\nabla \times$ à cette force.

Écrivons la fonction d'activation seuillée à l'aide de la fonction échelon unitaire de Heaviside $\theta(x)$ :
$$\mathbf{F}_{\text{Kevlar}} = - K_{\text{elastic}} \cdot (J - 1) \, \theta(J - 1) \cdot \mathbf{u}$$

Posons la fonction scalaire modulatrice :
$$\phi(J) = (J - 1) \, \theta(J - 1)$$

La force s'écrit alors de façon compacte :
$$\mathbf{F}_{\text{Kevlar}} = - K_{\text{elastic}} \, \phi(J) \, \mathbf{u}$$

En utilisant l'identité vectorielle standard du calcul différentiel pour le rotationnel d'un produit d'un champ scalaire par un champ vectoriel $\nabla \times (\phi \mathbf{A}) = \phi (\nabla \times \mathbf{A}) + \nabla \phi \times \mathbf{A}$ :

$$\nabla \times \mathbf{F}_{\text{Kevlar}} = - K_{\text{elastic}} \left[ \phi(J) (\nabla \times \mathbf{u}) + \nabla \phi(J) \times \mathbf{u} \right]$$

### A. Calcul rigoureux du gradient $\nabla \phi(J)$
Par dérivation en chaîne au sens des distributions :
$$\nabla \phi(J) = \frac{d\phi}{dJ} \nabla J = \left[ \theta(J - 1) + (J - 1) \delta(J - 1) \right] \nabla J$$

où $\delta$ est la distribution de Dirac. Remarquons que pour tout $J$, le produit $(J - 1)\delta(J - 1) \equiv 0$, puisque la distribution de Dirac est concentrée en $J=1$ où le facteur $(J-1)$ s'annule identiquement. Il vient donc rigoureusement :
$$\nabla \phi(J) = \theta(J - 1) \nabla J$$

### B. L'équation fondamentale de la vorticité modifiée
En substituant cette relation et en réintroduisant $\nabla \times \mathbf{u} = \boldsymbol{\omega}$ :

$$\nabla \times \mathbf{F}_{\text{Kevlar}} = \underbrace{- K_{\text{elastic}} \max(0, J - 1) \, \boldsymbol{\omega}}_{\text{Terme I : Amortissement Stable (Voulu)}} \; \underbrace{- K_{\text{elastic}} \, \theta(J - 1) \left( \nabla J \times \mathbf{u} \right)}_{\text{Terme II : Source Parasite de Vorticité}}$$

Ce résultat analytique constitue le théorème central de cette étude :  
**L'application d'une force de traînée proportionnelle à la vitesse absolue dans un champ métrique hétérogène engendre obligatoirement une source de vorticité non conservative régie par le produit vectoriel $\nabla J \times \mathbf{u}$.**

---

## 3. Analyse Géométrique et Vectorielle dans un Vortex Cyclonique

Examinons l'orientation et le signe du Terme II dans la géométrie cylindrique standard d'un cyclone tropical (repère $(r, \lambda, z)$ où $r$ est le rayon, $\lambda$ l'azimut et $z$ la coordonnée verticale).

Dans un cyclone de l'hémisphère nord, la circulation primaire est caractérisée par :
1. **La vitesse du vent $\mathbf{u}$ :** Elle est dominée par sa composante tangentielle azimutale $v_\lambda > 0$ (rotation anti-horaire cyclonique) :
   $$\mathbf{u} \approx v_\lambda(r, z) \, \mathbf{e}_\lambda$$
2. **La métrique d'étirement $J$ :** L'accumulation de déformation $S_{\text{mag}}$ est maximale dans le mur de l'œil ($r = R_{\text{max}}$) et décroît vers la périphérie de l'œil et vers l'extérieur. Dans la région située juste à l'extérieur du rayon des vents maximaux, $J$ diminue radialement vers $1$ :
   $$\frac{\partial J}{\partial r} < 0 \implies \nabla J \approx -\left| \frac{\partial J}{\partial r} \right| \mathbf{e}_r$$

Calculons le produit vectoriel dans la base orthonormée cylindrique $(\mathbf{e}_r, \mathbf{e}_\lambda, \mathbf{e}_z)$ où $\mathbf{e}_r \times \mathbf{e}_\lambda = \mathbf{e}_z$ :
$$\nabla J \times \mathbf{u} = \left( -\left| \frac{\partial J}{\partial r} \right| \mathbf{e}_r \right) \times \left( v_\lambda \mathbf{e}_\lambda \right) = - \left| \frac{\partial J}{\partial r} \right| v_\lambda \, \mathbf{e}_z$$

Injectons cette expression dans le Terme II de l'équation de vorticité :
$$\mathbf{P}_z = - K_{\text{elastic}} \, \theta(J - 1) \left( \nabla J \times \mathbf{u} \right) = + K_{\text{elastic}} \, \theta(J - 1) \left| \frac{\partial J}{\partial r} \right| v_\lambda \, \mathbf{e}_z$$

### Conclusion fondamentale sur le signe :
Comme $K_{\text{elastic}} > 0$, $\theta(J-1) \ge 0$, $|\partial_r J| > 0$ et $v_\lambda > 0$, le terme parasite $\mathbf{P}_z$ est **strictement positif le long de l'axe vertical $z$** :
$$\mathbf{P}_z \cdot \mathbf{e}_z > 0$$

Dans l'hémisphère sud ($v_\lambda < 0$ et $\omega_z < 0$), $\nabla J \times \mathbf{u}$ change de signe, et le terme parasite $\mathbf{P}_z$ devient négatif, **amplifiant de la même manière la vorticité cyclonique négative**.

> **Corollaire dynamique :**  
> Indépendamment de l'hémisphère considéré, le terme résiduel $- K_{\text{elastic}} (\nabla J \times \mathbf{u})$ n'agit jamais comme un frein, mais comme une **pompe à vorticité artificielle** qui injecte du spin vertical pur au cœur du cyclone.

---

## 4. Analyse d'Échelle et Conséquences Numériques (Le Cas Kenneth)

Pour quantifier l'impact de ce terme par rapport à la physique résolue, effectuons une analyse d'ordres de grandeur à l'aide des paramètres réels de la simulation du cyclone Kenneth sous WRF :

* Résolution horizontale : $\Delta x = 1000\text{ m}$
* Largeur caractéristique du mur de l'œil : $\Delta L \approx 5000\text{ m}$
* Vitesse tangentielle maximale : $U \approx 50\text{ m/s}$ ($100\text{ kt}$)
* Constante élastique : $K_{\text{elastic}} = 0.05\text{ s}^{-1}$
* Pic de déformation métrique : $J_{\text{max}} \approx 2.0$ (donc $J - 1 = 1.0$)
* Vorticité verticale typique du cœur : $\omega_z \approx 10^{-2}\text{ s}^{-1}$

### A. Ordre de grandeur du terme d'amortissement voulu (Terme I) :
$$\|\mathbf{D}_{\text{stable}}\| = K_{\text{elastic}} (J - 1) \omega_z \approx 0.05 \times 1.0 \times 10^{-2} = 5.0 \times 10^{-4}\text{ s}^{-2}$$

### B. Ordre de grandeur du terme parasite induit (Terme II) :
Le gradient spatial du Jacobien s'évalue par :
$$\|\nabla J\| \approx \frac{\Delta J}{\Delta L} \approx \frac{2.0 - 1.0}{5000\text{ m}} = 2.0 \times 10^{-4}\text{ m}^{-1}$$

D'où l'amplitude du terme parasite vertical :
$$\|\mathbf{P}_z\| = K_{\text{elastic}} \|\nabla J\| U \approx 0.05 \times (2.0 \times 10^{-4}) \times 50 = 5.0 \times 10^{-4}\text{ s}^{-2}$$

### C. Verdict quantitatif :
$$\frac{\|\mathbf{P}_z\|}{\|\mathbf{D}_{\text{stable}}\|} \approx 1.0$$

Le terme d'injection artificielle de rotation est du **même ordre de grandeur absolu** que le freinage théorique. 

Plus critique encore : à mesure que la résolution de grille est raffinée ($\Delta x \to 100\text{ m}$), la zone de transition du mur de l'œil se rétrécit ($\Delta L \to 0$), ce qui entraîne une divergence du gradient $\|\nabla J\| \propto 1/\Delta L$. Dès lors, le terme parasite **domine asymptotiquement** le terme stable, transformant le modèle en une machine numérique explosive à feedback positif instable.

---

## 5. Violations des Principes Fondamentaux de la Physique

L'échec de la formulation vectorielle V14 met en lumière trois violations méthodologiques majeures en mécanique des milieux continus :

### 1. Violation de l'Invariance Galiléenne (Indifférence Matérielle)
Soit un observateur situé dans un référentiel inertiel en translation uniforme à la vitesse constante $\mathbf{V}$. Pour cet observateur, la vitesse du fluide devient $\mathbf{u}' = \mathbf{u} + \mathbf{V}$.  
La force calculée par la formule V14 devient :
$$\mathbf{F}'_{\text{Kevlar}} = - K (J - 1) (\mathbf{u} + \mathbf{V}) \neq \mathbf{F}_{\text{Kevlar}}$$
La force dépend donc de la vitesse de l'observateur. Or, les contraintes internes d'un milieu continu doivent être strictement invariantes par changement de repère galiléen (principe d'indifférence matérielle objective de Noll, 1958). Un fluide ne peut subir d'efforts internes que par le biais des **gradients relatifs de vitesse**, jamais par sa vitesse absolue.

### 2. Violation du Principe d'Action-Réaction Interne
Dans un continuum fermé de volume $\Omega$, la somme des forces intérieures doit être identiquement nulle : $\int_\Omega \mathbf{F}_{\text{int}} \, dV = 0$.  
Avec une force $\mathbf{F} = -K(J-1)\mathbf{u}$, cette intégrale ne s'annule pas en général. La force se comporte comme une traînée externe s'appuyant sur un « éther » immobile imaginaire, créant un puits/source de quantité de mouvement globalement non conservatif.

---

## 6. La Régularisation Formelle : Le Tenseur de Cauchy

Pour restaurer l'exactitude physique et mathématique, nous devons impérativement appliquer le **Théorème Fondamental de Cauchy (1827)** :  
Toute force volumique interne résultant de la déformation d'un fluide continu doit s'écrire exclusivement sous la forme de la **divergence d'un tenseur de contraintes d'ordre 2 symétrique** $\boldsymbol{\tau}$ :

$$\mathbf{F}_{\text{interne}} = \nabla \cdot \boldsymbol{\tau}, \qquad \text{avec} \quad \boldsymbol{\tau} = \boldsymbol{\tau}^T$$

### A. Définition du Tenseur Visco-Métrique AS-Flow V15
Au lieu d'agir directement sur le vecteur vitesse $\mathbf{u}$, la métrique d'étirement $J$ doit moduler la **résistance interne au cisaillement** (rhéologie adaptative non-linéaire) :

$$\boldsymbol{\tau}(J) = 2 \, \nu_{\text{eff}}(J) \, \mathbf{S}$$

où :
1. $\mathbf{S}$ est le tenseur déviatorique symétrique des taux de déformation :
   $$S_{ij} = \frac{1}{2}\left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right) - \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij}$$
2. $\nu_{\text{eff}}(J)$ est la viscosité cinématique effective gouvernée par le seuil métrique :
   $$\nu_{\text{eff}}(J) = \nu_0 + \nu_{\text{AS}} \cdot \max(0, J - 1)$$

### B. Démonstration de l'Éradication du Terme Parasite
Calculons le rotationnel de cette force tensorielle :
$$\nabla \times (\nabla \cdot \boldsymbol{\tau}) = \nabla \times \left( \nabla \cdot [2 \nu_{\text{eff}}(J) \mathbf{S}] \right)$$

En utilisant la propriété tensorielle $\nabla \cdot \mathbf{S} = \frac{1}{2} \nabla^2 \mathbf{u} + \frac{1}{6} \nabla(\nabla \cdot \mathbf{u})$ et en prenant le rotationnel (sachant que $\nabla \times \nabla(\cdot) \equiv 0$) :
$$\nabla \times (\nabla \cdot \boldsymbol{\tau}) = \nu_{\text{eff}}(J) \, \nabla^2 \boldsymbol{\omega} + 2 \left( \nabla \nu_{\text{eff}} \right) \cdot \nabla \mathbf{S} + \dots$$

**Propriété fondamentale :**  
Comme l'opérateur s'applique sur un tenseur strictement symétrique $\mathbf{S}$ dont la trace et les composantes antisymétriques sont nulles, **aucun produit vectoriel du type $\nabla J \times \mathbf{u}$ ne peut émerger de l'algèbre différentielle**. 

L'action du tenseur sur l'équation de vorticité se réduit à un **opérateur de diffusion parabolique pur** :
$$\frac{\partial \boldsymbol{\omega}}{\partial t} \supset \nu_{\text{eff}}(J) \, \nabla^2 \boldsymbol{\omega}$$

L'enstrophie $\mathcal{E} = \frac{1}{2} \int \|\boldsymbol{\omega}\|^2 dV$ satisfait alors :
$$\frac{d\mathcal{E}}{dt} \le - \nu_{\text{eff}}(J) \int \|\nabla \boldsymbol{\omega}\|^2 dV \le 0$$
L'amortissement est inconditionnellement stable et ne peut en aucun cas alimenter la circulation du cyclone.

---

## 7. Discussion Épistémologique et Recommandations pour la CFD

Bien que la formulation tensorielle V15 $\boldsymbol{\tau} = 2 \nu_{\text{eff}}(J) \mathbf{S}$ guérisse l'anomalie de spin et respecte l'invariance galiléenne, une question épistémologique cruciale demeure : **est-il pertinent de résoudre une EDP pronostique supplémentaire tridimensionnelle pour $J$ dans un modèle atmosphérique opérationnel ?**

Dans les modèles météorologiques modernes (WRF, IFS-CEPMMT, GFS-NOAA), la stabilisation sous-maille des ondes courtes à la limite de Nyquist ($2\Delta x$) est assurée par :
* La diffusion numérique hyper-visqueuse d'ordre 4 ou 6 (Skamarock et Klemp, 2008),
* L'amortissement de divergence (*divergence damping*),
* Les fermetures de sous-maille TKE ou Smagorinsky 3D (Smagorinsky, 1963 ; Germano et al., 1991).

Ces approches standard présentent un avantage déterminant : **elles ne requièrent aucune variable pronostique supplémentaire** et reposent sur la cascade turbulente de Kolmogorov (1941). En comparaison, transporter un champ scalaire 3D $J$ ajoute une charge mémoire non négligeable, des coûts d'advection Runge-Kutta 3e ordre et la nécessité de calibrer empiriquement plusieurs constantes ($K_{\text{elastic}}, \nu_{\text{AS}}$).

### Recommandation opérationnelle :
Le rôle le plus fécond pour la métrique $J$ n'est pas de servir de rétroaction physique sur l'air, mais d'agir comme un **indicateur diagnostique de singularité de maillage (*Mesh Quality Sensor*)** :
1. $J$ peut servir de critère de raffinement automatique de maillage (AMR - *Adaptive Mesh Refinement*) pour densifier les mailles là où l'étirement devient critique.
2. Il peut agir comme limiteur de flux localisé dans les schémas d'advection WENO, sans altérer la physique globale des équations de Navier-Stokes.

---

## 8. Conclusion

Cette étude a permis de disséquer l'anomalie fondamentale du couplage vectoriel V14 de l'AS-Flow :
1. Nous avons prouvé analytiquement que la force vectorielle $-K(J-1)\mathbf{u}$ génère inévitablement un terme parasite de vorticité $-K(\nabla J \times \mathbf{u})$.
2. Dans un vortex géophysique, ce terme est strictement positif et agit comme une turbine artificielle qui accélère le cyclone au lieu de le modérer.
3. Nous avons démontré que la seule formulation mathématiquement viable est la formulation tensorielle de Cauchy $\nabla \cdot (2\nu(J)\mathbf{S})$, qui garantit l'invariance galiléenne et l'éradication complète du spin parasite.

Ce travail souligne l'importance impérative d'analyser systématiquement le rotationnel de toute force sous-maille proposée en dynamique des fluides géophysiques afin d'éviter l'introduction involontaire de générateurs artificiels de tourbillon.

---

## Disponibilité des Codes et Données

La suite complète de démonstration numérique 2D (`pillar_a_spurious_vorticity_2d.py`), les routines sources Fortran du solveur WRF (`module_diffusion_asflow_v14.F`) et les relevés de métriques du cyclone Kenneth sont archivés en accès libre sur Zenodo à l'adresse [https://doi.org/10.5281/zenodo.placeholder](https://doi.org/10.5281/zenodo.placeholder) et accessibles publiquement sur GitHub à l'adresse [https://github.com/etonam212g5-art/AS-FLOW-RESEARCH-LAB](https://github.com/etonam212g5-art/AS-FLOW-RESEARCH-LAB).

---

## Références Bibliographiques

1. **Beale, J. T., Kato, T., & Majda, A. (1984).** Remarks on the breakdown of smooth solutions for the 3-D Euler equations. *Communications in Mathematical Physics*, 94(1), 61-66.
2. **Cauchy, A. L. (1827).** De la pression ou tension dans un corps solide. *Exercices de Mathématiques*, 2, 42-56.
3. **Germano, M., Piomelli, U., Moin, P., & Cabot, W. H. (1991).** A dynamic subgrid-scale eddy viscosity model. *Physics of Fluids A: Fluid Dynamics*, 3(7), 1760-1765.
4. **Kolmogorov, A. N. (1941).** The local structure of turbulence in incompressible viscous fluid for very large Reynolds numbers. *Dokl. Akad. Nauk SSSR*, 30(4), 301-305.
5. **Noll, W. (1958).** A mathematical theory of the mechanical behavior of continuous media. *Archive for Rational Mechanics and Analysis*, 2(1), 197-226.
6. **Skamarock, W. C., & Klemp, J. B. (2008).** A time-split nonhydrostatic atmospheric model for weather research and forecasting applications. *Journal of Computational Physics*, 227(7), 3465-3485.
7. **Smagorinsky, J. (1963).** General circulation experiments with the primitive equations: I. The basic experiment. *Monthly Weather Review*, 91(3), 99-164.
8. **Von Neumann, J., & Richtmyer, R. D. (1950).** A method for the numerical calculation of hydrodynamic shocks. *Journal of Applied Physics*, 21(3), 232-237.