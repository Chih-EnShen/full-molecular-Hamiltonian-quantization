# Catalog: Lagrangians and Hamiltonians at Each Coordinate Stage

This file lists the classical **Lagrangian in velocity form** $L(q,\dot q)$ and the classical **Hamiltonian in momentum form** $H(q,p)$ at the three coordinate stages of the paper:

1. lab-frame Cartesian coordinates (Sec. 1);
2. mass-centered coordinates — translational coordinate = **total** CoM, internal coordinates referred to the **nuclear** CoM (Sec. 2);
3. body-fixed (Eckart) Cartesian coordinates — Euler angles introduced, nuclear shape still Cartesian (Sec. 3.1–3.3);
4. normal-mode coordinates — the body-fixed frame with the shape parametrized by $Q_b$ (Sec. 3).

No Routhian is used anywhere in this file: at every stage the Lagrangian carries **all** degrees of freedom (electrons included) in velocity form, and the Hamiltonian is its full Legendre transform. Where this differs from the paper's mixed (Routhian) route in Sec. 3, the difference is flagged; the final Hamiltonians are identical, because conjugate momenta do not depend on the route taken.

Throughout: $N_{\mathrm{nucl}}$ nuclei ($j$, masses $M_j$), $N_{\mathrm{elec}}$ electrons ($i$, mass $m_e$),

$$
M_N=\sum_j M_j,\qquad M_{\mathrm{tot}}=M_N+N_{\mathrm{elec}}\,m_e ,
$$

and $V$ is the Coulomb potential of eq. (1.2), a function of interparticle separations only.

---

## 1. Lab-frame Cartesian coordinates

**Coordinates.** $\{\mathbf R_j,\mathbf r_i\}$, all in the space-fixed frame. No constraints.

**Lagrangian (velocity form).**

$$
L \;=\; \sum_j \tfrac12 M_j\,|\dot{\mathbf R}_j|^2 \;+\; \sum_i \tfrac12 m_e\,|\dot{\mathbf r}_i|^2 \;-\; V(\mathbf R_j,\mathbf r_i).
$$

**Momenta.** $\;\mathbf P_j=\partial L/\partial\dot{\mathbf R}_j=M_j\dot{\mathbf R}_j$, $\;\mathbf p_i=m_e\dot{\mathbf r}_i$.

**Hamiltonian (momentum form)** — eq. (1.1):

$$
H \;=\; \sum_j \frac{\mathbf P_j^{\,2}}{2M_j} \;+\; \sum_i \frac{\mathbf p_i^{\,2}}{2m_e} \;+\; V(\mathbf R_j,\mathbf r_i).
$$

*Notes.* Kinetic energy diagonal in both forms; the momentum–velocity map is trivial. This is the only stage where that is true.

---

## 2. Mass-centered coordinates (translation = total CoM, internal reference = nuclear CoM)

**Coordinates** — eqs. (2.1)–(2.3):

$$
\mathbf R_{\mathrm{cm}}=\frac{1}{M_{\mathrm{tot}}}\Bigl(\sum_j M_j\mathbf R_j+m_e\sum_i\mathbf r_i\Bigr),\qquad
\mathbf R'_j=\mathbf R_j-\mathbf R_{\mathrm{ncm}},\qquad
\mathbf r'_i=\mathbf r_i-\mathbf R_{\mathrm{ncm}},
$$

with $\mathbf R_{\mathrm{ncm}}=\tfrac1{M_N}\sum_j M_j\mathbf R_j$ the nuclear CoM.

**Constraints.** Coordinate: $\sum_j M_j\mathbf R'_j=0$ (2.4), hence also $\sum_j M_j\dot{\mathbf R}'_j=0$. Momentum (dual gauge): $\sum_j\mathbf P'_j=0$. Electrons unconstrained.

**Lagrangian (velocity form).**

$$
L \;=\; \tfrac12 M_{\mathrm{tot}}\,|\dot{\mathbf R}_{\mathrm{cm}}|^2
\;+\; \sum_j \tfrac12 M_j\,|\dot{\mathbf R}'_j|^2
\;+\; \sum_i \tfrac12 m_e\,|\dot{\mathbf r}'_i|^2
\;-\; \frac{m_e^2}{2M_{\mathrm{tot}}}\Bigl|\sum_i\dot{\mathbf r}'_i\Bigr|^2
\;-\; V(\mathbf R'_j,\mathbf r'_i).
$$

**Momenta.** $\;\mathbf P_{\mathrm{cm}}=M_{\mathrm{tot}}\dot{\mathbf R}_{\mathrm{cm}}$, $\;\mathbf P'_j=M_j\dot{\mathbf R}'_j$ (automatically satisfying the dual gauge), and — eq. (3.5) —

$$
\mathbf p'_i=m_e\dot{\mathbf r}'_i-\frac{m_e^2}{M_{\mathrm{tot}}}\sum_l\dot{\mathbf r}'_l
\qquad\Longleftrightarrow\qquad
\dot{\mathbf r}'_i=\frac{\mathbf p'_i}{m_e}+\frac{1}{M_N}\sum_l\mathbf p'_l .
$$

**Hamiltonian (momentum form)** — eq. (2.12):

$$
H \;=\; \frac{\mathbf P_{\mathrm{cm}}^{\,2}}{2M_{\mathrm{tot}}}
\;+\; \sum_j\frac{\mathbf P_j^{\prime\,2}}{2M_j}
\;+\; \sum_i\frac{\mathbf p_i^{\prime\,2}}{2m_e}
\;+\; \frac{1}{2M_N}\Bigl(\sum_i\mathbf p'_i\Bigr)^{\!2}
\;+\; V(\mathbf R'_j,\mathbf r'_i).
$$

*Notes.*

- **The mass polarization flips sign and mass between the two forms:** in momentum form it is $+\tfrac{1}{2M_N}(\sum_i\mathbf p'_i)^2$; in velocity form it is $-\tfrac{m_e^2}{2M_{\mathrm{tot}}}(\sum_i\dot{\mathbf r}'_i)^2$. Both statements encode the same physics — the non-trivial momentum–velocity map above converts one into the other under the Legendre transform.
- **Reduced-mass check** ($N_{\mathrm{elec}}=1$): the velocity form gives $\tfrac12 m_e(1-\tfrac{m_e}{M_{\mathrm{tot}}})|\dot{\mathbf r}'|^2=\tfrac12\mu_e|\dot{\mathbf r}'|^2$ and the momentum form $\mathbf p'^2\!/2\mu_e$, with $\mu_e=m_eM_N/M_{\mathrm{tot}}$ — the familiar electron reduced mass.
- The CoM decouples in both forms, with the correct total mass $M_{\mathrm{tot}}$, for any $\mathbf P_{\mathrm{cm}}$.

---

## 3. Body-fixed (Eckart) Cartesian coordinates

**Coordinates.** $\mathbf R_{\mathrm{cm}}$; Euler angles $\Omega=(\phi,\theta,\chi)$ defining the frame rotation $\mathbf S(\Omega)$; body-fixed Cartesian positions for nuclei **and** electrons:

$$
\mathbf R'_j=\mathbf S(\Omega)\,\bar{\mathbf R}_j,\qquad
\mathbf r'_i=\mathbf S(\Omega)\,\bar{\mathbf r}_i,\qquad
\boldsymbol\omega=\mathsf E(\Omega)\,\dot{\boldsymbol\Omega}.
$$

The nuclear shape is *not* yet parametrized — the $\bar{\mathbf R}_j$ remain Cartesian vectors.

**Constraints.** Bookkeeping: $3$ Euler angles $+\,3N_{\mathrm{nucl}}$ components of $\bar{\mathbf R}_j$ must equal the $3N_{\mathrm{nucl}}-3$ physical internal nuclear degrees of freedom, so six conditions are needed — precisely the Eckart conditions (Sec. 3.1):

$$
\underbrace{\sum_j M_j\bar{\mathbf R}_j=0}_{\text{translational (BF image of stage-2 constraint)}},\qquad
\underbrace{\sum_j M_j\,\bar{\mathbf R}_j^{\mathrm{eq}}\times\bar{\mathbf R}_j=0}_{\text{rotational (gauge-fixes the three Euler angles)}} ,
$$

with velocity versions $\sum_j M_j\dot{\bar{\mathbf R}}_j=0$, $\sum_j M_j\bar{\mathbf R}_j^{\mathrm{eq}}\times\dot{\bar{\mathbf R}}_j=0$, and momentum dual $\sum_j\bar{\mathbf P}_j=0$. Electrons unconstrained.

**Lagrangian (velocity form).** Rotations preserve norms, so substituting $\dot{\mathbf R}'_j=\mathbf S(\boldsymbol\omega\times\bar{\mathbf R}_j+\dot{\bar{\mathbf R}}_j)$, $\dot{\mathbf r}'_i=\mathbf S(\boldsymbol\omega\times\bar{\mathbf r}_i+\dot{\bar{\mathbf r}}_i)$ into the stage-2 Lagrangian:

$$
L \;=\; \tfrac12 M_{\mathrm{tot}}|\dot{\mathbf R}_{\mathrm{cm}}|^2
+\tfrac12\sum_j M_j\bigl|\boldsymbol\omega\times\bar{\mathbf R}_j+\dot{\bar{\mathbf R}}_j\bigr|^2
+\tfrac{m_e}{2}\sum_i\bigl|\boldsymbol\omega\times\bar{\mathbf r}_i+\dot{\bar{\mathbf r}}_i\bigr|^2
-\frac{m_e^2}{2M_{\mathrm{tot}}}\bigl|\boldsymbol\omega\times\bar{\mathbf r}+\dot{\bar{\mathbf r}}\bigr|^2
-V(\bar{\mathbf R}_j,\bar{\mathbf r}_i),
$$

with $\bar{\mathbf r}\equiv\sum_i\bar{\mathbf r}_i$; expanded via (3.9c),

$$
L \;=\; \tfrac12 M_{\mathrm{tot}}|\dot{\mathbf R}_{\mathrm{cm}}|^2
+\tfrac12\,\boldsymbol\omega^{T}\bigl(\mathbf I_N+\mathbf I_e-\mathbf I_X\bigr)\boldsymbol\omega
+\boldsymbol\omega\cdot\bigl(\mathbf L_N^{\mathrm{vib}}+\mathbf L_e-\mathbf L_X\bigr)
+\tfrac12\sum_j M_j|\dot{\bar{\mathbf R}}_j|^2
+\tfrac{m_e}{2}\sum_i|\dot{\bar{\mathbf r}}_i|^2
-\frac{m_e^2}{2M_{\mathrm{tot}}}|\dot{\bar{\mathbf r}}|^2
-V ,
$$

where $\mathbf I_N,\mathbf L_N^{\mathrm{vib}}$ are built from the $\bar{\mathbf R}_j$ exactly as in the definitions below (stage 4), only without the $Q$-parametrization.

**Momenta.** All three families are the body-fixed images of the stage-2 momenta:

$$
\bar{\mathbf P}_j=\frac{\partial L}{\partial\dot{\bar{\mathbf R}}_j}
=M_j\bigl(\boldsymbol\omega\times\bar{\mathbf R}_j+\dot{\bar{\mathbf R}}_j\bigr)=\mathbf S^{T}\mathbf P'_j ,
\qquad
\bar{\mathbf p}_i
=m_e\bigl(\boldsymbol\omega\times\bar{\mathbf r}_i+\dot{\bar{\mathbf r}}_i\bigr)
-\frac{m_e^2}{M_{\mathrm{tot}}}\bigl(\boldsymbol\omega\times\bar{\mathbf r}+\dot{\bar{\mathbf r}}\bigr)=\mathbf S^{T}\mathbf p'_i ,
$$

and the angular momentum $\mathbf J=\partial L/\partial\boldsymbol\omega=(\mathbf I_N+\mathbf I_e-\mathbf I_X)\boldsymbol\omega+\mathbf L_N^{\mathrm{vib}}+\mathbf L_e-\mathbf L_X$ is **not independent**: it satisfies the identity

$$
\mathbf J=\sum_j\bar{\mathbf R}_j\times\bar{\mathbf P}_j+\underbrace{\sum_i\bar{\mathbf r}_i\times\bar{\mathbf p}_i}_{\mathbf L_{\mathrm{elec}}} ,
$$

the momentum-side counterpart of the rotational Eckart gauge (three redundant coordinates $\Rightarrow$ three dependent momenta).

**Hamiltonian (momentum form).** Because the kinetic energy is a rotational scalar, $H$ is *form-invariant* under body-fixing — it is the stage-2 Hamiltonian with barred variables:

$$
H \;=\; \frac{\mathbf P_{\mathrm{cm}}^{\,2}}{2M_{\mathrm{tot}}}
\;+\; \sum_j\frac{\bar{\mathbf P}_j^{\,2}}{2M_j}
\;+\; \sum_i\frac{\bar{\mathbf p}_i^{\,2}}{2m_e}
\;+\; \frac{1}{2M_N}\Bigl(\sum_i\bar{\mathbf p}_i\Bigr)^{\!2}
\;+\; V(\bar{\mathbf R}_j,\bar{\mathbf r}_i).
$$

*Notes.*

- **Body-fixing alone changes nothing in $H$:** no $\boldsymbol\mu$, no explicit $\mathbf J$, no rotational kernel. The rotational kinetic energy is *hidden inside* the $\bar{\mathbf P}_j$ (each contains $M_j\,\boldsymbol\omega\times\bar{\mathbf R}_j$). Extracting an explicit rotation–vibration split requires parametrizing the shape — that is stage 4's job.
- What body-fixing **does** buy at this stage: the potential becomes $V(\bar{\mathbf R}_j,\bar{\mathbf r}_i)$, manifestly independent of $\Omega$ [eq. (3.8)] — rotational symmetry is explicit, and the Euler angles are cyclic.
- The identity $\mathbf J=\sum_j\bar{\mathbf R}_j\times\bar{\mathbf P}_j+\mathbf L_{\mathrm{elec}}$ is the seed of the key identity (3.19): once $\dot{\bar{\mathbf R}}_j$ is parametrized by $\dot Q_b$, its nuclear part becomes $\mathbf I_N\boldsymbol\omega+\mathbf L_N^{\mathrm{vib}}$ and inversion produces the Watson $\boldsymbol\mu$.

---

## 4. Normal-mode coordinates

**Coordinates.** $\mathbf R_{\mathrm{cm}}$; Euler angles $\Omega=(\phi,\theta,\chi)$ fixing $\mathbf S(\Omega)$ by the Eckart conditions; normal coordinates $Q_b$ ($b=1,\dots,3N_{\mathrm{nucl}}-6$); body-fixed electron coordinates $\bar{\mathbf r}_i$:

$$
\mathbf R'_j=\mathbf S(\Omega)\,\bar{\mathbf R}_j(Q),\qquad
\bar{\mathbf R}_j(Q)=\bar{\mathbf R}_j^{\mathrm{eq}}+\frac{1}{\sqrt{M_j}}\sum_b\mathbf l_{jb}Q_b,\qquad
\mathbf r'_i=\mathbf S(\Omega)\,\bar{\mathbf r}_i .
$$

The body-fixed angular velocity is $\boldsymbol\omega=\mathsf E(\Omega)\dot{\boldsymbol\Omega}$ [eqs. (3.2), (4.10)]. The Eckart/Sayvetz conditions are consumed by the chart: there are no residual constraints among $(\Omega,Q_b,\bar{\mathbf r}_i)$.

**Ingredient definitions** [eqs. (3.9)–(3.12), (3.14), (3.17a), (3.20)], with $\bar{\mathbf r}\equiv\sum_i\bar{\mathbf r}_i$:

$$
\mathbf I_N=\sum_j M_j(\bar R_j^2\mathbb 1-\bar{\mathbf R}_j\bar{\mathbf R}_j^T),\qquad
\mathbf I_e=m_e\sum_i(\bar r_i^2\mathbb 1-\bar{\mathbf r}_i\bar{\mathbf r}_i^T),\qquad
\mathbf I_X=\frac{m_e^2}{M_{\mathrm{tot}}}(\bar r^2\mathbb 1-\bar{\mathbf r}\bar{\mathbf r}^T),
$$

$$
\mathbf L_N^{\mathrm{vib}}=\sum_j M_j\bar{\mathbf R}_j\times\dot{\bar{\mathbf R}}_j,\qquad
\mathbf L_e=m_e\sum_i\bar{\mathbf r}_i\times\dot{\bar{\mathbf r}}_i,\qquad
\mathbf L_X=\frac{m_e^2}{M_{\mathrm{tot}}}\,\bar{\mathbf r}\times\dot{\bar{\mathbf r}},
$$

$$
\boldsymbol\sigma_b=\sum_a\boldsymbol\zeta_{ab}Q_a,\quad
\boldsymbol\zeta_{ab}=\sum_j\mathbf l_{ja}\times\mathbf l_{jb},\qquad
\boldsymbol\pi=\sum_b\boldsymbol\sigma_bP_b,\qquad
\boldsymbol\mu=(\mathbf I_N-\boldsymbol\sigma\boldsymbol\sigma^{\mathsf T})^{-1},\qquad
\mathbf L_{\mathrm{elec}}=\sum_i\bar{\mathbf r}_i\times\bar{\mathbf p}_i .
$$

**Lagrangian (velocity form, all degrees of freedom — no Routhian).** This is the stage-3 Lagrangian with the shape parametrized, $\bar{\mathbf R}_j\to\bar{\mathbf R}_j(Q)$ and $\dot{\bar{\mathbf R}}_j=\sum_b\mathbf l_{jb}\dot Q_b/\sqrt{M_j}$ (so $\tfrac12\sum_jM_j|\dot{\bar{\mathbf R}}_j|^2=\tfrac12\sum_b\dot Q_b^2$ by orthonormality of the $\mathbf l_{jb}$):

$$
L \;=\; \tfrac12 M_{\mathrm{tot}}|\dot{\mathbf R}_{\mathrm{cm}}|^2
+\tfrac12\sum_j M_j\bigl|\boldsymbol\omega\times\bar{\mathbf R}_j+\dot{\bar{\mathbf R}}_j\bigr|^2
+\tfrac{m_e}{2}\sum_i\bigl|\boldsymbol\omega\times\bar{\mathbf r}_i+\dot{\bar{\mathbf r}}_i\bigr|^2
-\frac{m_e^2}{2M_{\mathrm{tot}}}\bigl|\boldsymbol\omega\times\bar{\mathbf r}+\dot{\bar{\mathbf r}}\bigr|^2
-V(Q,\bar{\mathbf r}_i),
$$

or, expanded via the squared-term identity (3.9c) and eqs. (3.9d)–(3.9e):

$$
L \;=\; \tfrac12 M_{\mathrm{tot}}|\dot{\mathbf R}_{\mathrm{cm}}|^2
+\tfrac12\,\boldsymbol\omega^{T}\bigl(\mathbf I_N+\mathbf I_e-\mathbf I_X\bigr)\boldsymbol\omega
+\boldsymbol\omega\cdot\bigl(\mathbf L_N^{\mathrm{vib}}+\mathbf L_e-\mathbf L_X\bigr)
+\tfrac12\sum_b\dot Q_b^2
+\tfrac{m_e}{2}\sum_i|\dot{\bar{\mathbf r}}_i|^2
-\frac{m_e^2}{2M_{\mathrm{tot}}}|\dot{\bar{\mathbf r}}|^2
-V .
$$

**Momenta.** With $p_{\Omega_s}=\partial L/\partial\dot\Omega_s=(\mathsf E^{T}\mathbf J)_s$ defining the body-frame angular momentum $\mathbf J$ (conjugate to $\boldsymbol\omega$ in the quasi-velocity sense):

$$
\mathbf J=\frac{\partial L}{\partial\boldsymbol\omega}
=\bigl(\mathbf I_N+\mathbf I_e-\mathbf I_X\bigr)\boldsymbol\omega+\mathbf L_N^{\mathrm{vib}}+\mathbf L_e-\mathbf L_X ,
\qquad
P_b=\frac{\partial L}{\partial\dot Q_b}=\dot Q_b+\boldsymbol\omega\cdot\boldsymbol\sigma_b ,
$$

$$
\bar{\mathbf p}_i=\frac{\partial L}{\partial\dot{\bar{\mathbf r}}_i}
=m_e\bigl(\boldsymbol\omega\times\bar{\mathbf r}_i+\dot{\bar{\mathbf r}}_i\bigr)
-\frac{m_e^2}{M_{\mathrm{tot}}}\bigl(\boldsymbol\omega\times\bar{\mathbf r}+\dot{\bar{\mathbf r}}\bigr)
\qquad\text{[BF image of (3.5)]} .
$$

Substituting $\bar{\mathbf p}_i$ into $\mathbf L_{\mathrm{elec}}$ gives eq. (3.18), $\mathbf L_{\mathrm{elec}}=(\mathbf I_e-\mathbf I_X)\boldsymbol\omega+(\mathbf L_e-\mathbf L_X)$, whence the key identity (3.19):

$$
\mathbf J-\mathbf L_{\mathrm{elec}}=\mathbf I_N\,\boldsymbol\omega+\mathbf L_N^{\mathrm{vib}} .
$$

**Hamiltonian (momentum form)** — $H=\mathbf P_{\mathrm{cm}}^2/2M_{\mathrm{tot}}+H_{\mathrm{int}}$ with $H_{\mathrm{int}}$ the Watson-form eq. (3.22):

$$
H_{\mathrm{int}} \;=\;
\tfrac12\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}-\boldsymbol\pi\bigr)^{T}\boldsymbol\mu(Q)\,\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}-\boldsymbol\pi\bigr)
\;+\;\tfrac12\sum_b P_b^2
\;+\;\sum_i\frac{\bar{\mathbf p}_i^{\,2}}{2m_e}
\;+\;\frac{1}{2M_N}\Bigl(\sum_i\bar{\mathbf p}_i\Bigr)^{\!2}
\;+\;V_{NN}(Q)+V_{Ne}(Q,\bar{\mathbf r}_i)+V_{ee}(\bar{\mathbf r}_i).
$$

*Notes.*

- **The electron inertias appear in $L$ but cancel out of $H$.** The velocity-form rotational kernel is $\mathbf I_N+\mathbf I_e-\mathbf I_X$, yet the momentum-form kernel is the purely nuclear $\boldsymbol\mu=(\mathbf I_N-\boldsymbol\sigma\boldsymbol\sigma^{\mathsf T})^{-1}$: by (3.18)–(3.19) the entire electronic and mass-polarization content of $\mathbf J$ is absorbed into the kinematic shift $\mathbf J\to\mathbf J-\mathbf L_{\mathrm{elec}}$ — the "algebraic miracle" of §3.6. The paper's Routhian route never lets $\mathbf I_e,\mathbf I_X$ arise in the first place; the full-Lagrangian route used here makes them appear and then cancel. Either way, the same (3.22).
- **The mass polarization changes appearance again:** velocity form $-\tfrac{m_e^2}{2M_{\mathrm{tot}}}|\boldsymbol\omega\times\bar{\mathbf r}+\dot{\bar{\mathbf r}}|^2$ (negative, coupled to rotation), momentum form $+\tfrac{1}{2M_N}(\sum_i\bar{\mathbf p}_i)^2$ (positive, a rotational scalar).
- **No Wilson $G$-matrix:** the vibrational metric of rectilinear normal coordinates is the identity, so the vibrational block is the bare $\tfrac12\sum_bP_b^2$; all rotation–vibration coupling lives in $\boldsymbol\pi$ and $\boldsymbol\mu$.
- The CoM term is a decoupled spectator with the correct mass $M_{\mathrm{tot}}$, exactly as at stage 2.

---

## Summary table

| Stage | Kinetic energy, velocity form (in $L$) | Kinetic energy, momentum form (in $H$) | Mass polarization sits in | Constraints |
|---|---|---|---|---|
| 1. Cartesian | $\sum_j\tfrac12 M_j\dot{\mathbf R}_j^2+\sum_i\tfrac12 m_e\dot{\mathbf r}_i^2$ | $\sum_j\tfrac{\mathbf P_j^2}{2M_j}+\sum_i\tfrac{\mathbf p_i^2}{2m_e}$ | — (absent) | none |
| 2. Mass-centered | $\tfrac12 M_{\mathrm{tot}}\dot{\mathbf R}_{\mathrm{cm}}^2+\sum_j\tfrac12 M_j\dot{\mathbf R}_j^{\prime 2}+\sum_i\tfrac12 m_e\dot{\mathbf r}_i^{\prime 2}-\tfrac{m_e^2}{2M_{\mathrm{tot}}}(\sum_i\dot{\mathbf r}'_i)^2$ | $\tfrac{\mathbf P_{\mathrm{cm}}^2}{2M_{\mathrm{tot}}}+\sum_j\tfrac{\mathbf P_j^{\prime 2}}{2M_j}+\sum_i\tfrac{\mathbf p_i^{\prime 2}}{2m_e}+\tfrac{(\sum_i\mathbf p'_i)^2}{2M_N}$ | $-\tfrac{m_e^2}{2M_{\mathrm{tot}}}(\Sigma\dot{\mathbf r}')^2$ vs $+\tfrac{(\Sigma\mathbf p')^2}{2M_N}$ | $\sum_jM_j\mathbf R'_j=0$, $\sum_j\mathbf P'_j=0$ |
| 3. Body-fixed Cartesian | $\tfrac12 M_{\mathrm{tot}}\dot{\mathbf R}_{\mathrm{cm}}^2+\tfrac12\boldsymbol\omega^T(\mathbf I_N{+}\mathbf I_e{-}\mathbf I_X)\boldsymbol\omega+\boldsymbol\omega\!\cdot\!(\mathbf L_N^{\mathrm{vib}}{+}\mathbf L_e{-}\mathbf L_X)+\tfrac12\sum_jM_j\dot{\bar{\mathbf R}}_j^2+\tfrac{m_e}{2}\sum_i\dot{\bar{\mathbf r}}_i^2-\tfrac{m_e^2}{2M_{\mathrm{tot}}}\dot{\bar{\mathbf r}}^2$ | $\tfrac{\mathbf P_{\mathrm{cm}}^2}{2M_{\mathrm{tot}}}+\sum_j\tfrac{\bar{\mathbf P}_j^2}{2M_j}+\sum_i\tfrac{\bar{\mathbf p}_i^2}{2m_e}+\tfrac{(\sum_i\bar{\mathbf p}_i)^2}{2M_N}$ | in $\mathbf I_X,\mathbf L_X,T_X^{\mathrm{mp}}$ vs $+\tfrac{(\Sigma\bar{\mathbf p})^2}{2M_N}$ | Eckart pair $\sum_jM_j\bar{\mathbf R}_j=0$, $\sum_jM_j\bar{\mathbf R}_j^{\mathrm{eq}}\!\times\!\bar{\mathbf R}_j=0$; $\sum_j\bar{\mathbf P}_j=0$; $\mathbf J=\Sigma\bar{\mathbf R}\times\bar{\mathbf P}+\mathbf L_{\mathrm{elec}}$ |
| 4. Normal-mode | $\tfrac12 M_{\mathrm{tot}}\dot{\mathbf R}_{\mathrm{cm}}^2+\tfrac12\boldsymbol\omega^T(\mathbf I_N{+}\mathbf I_e{-}\mathbf I_X)\boldsymbol\omega+\boldsymbol\omega\!\cdot\!(\mathbf L_N^{\mathrm{vib}}{+}\mathbf L_e{-}\mathbf L_X)+\tfrac12\sum_b\dot Q_b^2+\tfrac{m_e}{2}\sum_i\dot{\bar{\mathbf r}}_i^2-\tfrac{m_e^2}{2M_{\mathrm{tot}}}\dot{\bar{\mathbf r}}^2$ | $\tfrac{\mathbf P_{\mathrm{cm}}^2}{2M_{\mathrm{tot}}}+\tfrac12(\mathbf J{-}\mathbf L_{\mathrm{elec}}{-}\boldsymbol\pi)^T\boldsymbol\mu(\mathbf J{-}\mathbf L_{\mathrm{elec}}{-}\boldsymbol\pi)+\tfrac12\sum_bP_b^2+\sum_i\tfrac{\bar{\mathbf p}_i^2}{2m_e}+\tfrac{(\sum_i\bar{\mathbf p}_i)^2}{2M_N}$ | in $\mathbf I_X,\mathbf L_X,T_X^{\mathrm{mp}}$ vs $+\tfrac{(\Sigma\bar{\mathbf p})^2}{2M_N}$ | none (Eckart conditions consumed by the chart) |

The through-line of the table: the mass polarization is *absent* in Cartesian coordinates, appears as a *negative velocity-form correction* $\propto m_e^2/M_{\mathrm{tot}}$ or a *positive momentum-form term* $\propto 1/M_N$ at stage 2, and from stage 3 onward additionally entangles with rotation in velocity form (through $\mathbf I_X,\mathbf L_X$) while remaining a clean rotational scalar in momentum form — which is precisely why the paper carries the electrons in momentum form from stage 2 onward. Note also how the Hamiltonian's *form* is inert through stages 2 → 3 (body-fixing is a rotation, and the kinetic energy is a rotational scalar): the familiar Watson structure — rotational kernel $\boldsymbol\mu$, Coriolis $\boldsymbol\pi$, bare vibrational $\tfrac12\Sigma P_b^2$ — only materializes at stage 4, when the shape is parametrized and the rotation–vibration split of $\dot{\bar{\mathbf R}}_j$ becomes available.

---

## Momenta: canonical, mechanical, and quasi

"Non-canonical" means two independent things across the four stages, and the
two tables below keep them apart:

- **canonical vs *mechanical*** (mass × velocity). A momentum is *canonical* if
  it is the variable conjugate to a coordinate, $\{q,p\}=1$ (quantized as
  $-i\hbar\,\partial_q$); *mechanical* if it equals mass × velocity. They
  coincide at stage 1, then **split**: $\mathbf p'_i\ne m_e\dot{\mathbf r}'_i$
  (mass-polarization gap, from referring electrons to the nuclear CoM), and
  $P_b\ne\dot Q_b$ (Coriolis gap, from the rotating frame).
- **canonical vs *quasi* (frame-projected / non-holonomic)**. The true
  canonical momenta conjugate to the Euler angles are $p_\phi,p_\theta,p_\chi$.
  The physically useful body-frame $\mathbf J$, vibrational $\boldsymbol\pi$, and
  angular *velocity* $\boldsymbol\omega$ are coordinate-dependent combinations of
  the canonical momenta (or, for $\boldsymbol\omega$, a velocity) — not
  canonical. This is the origin of the anomalous
  $\{J_\alpha,J_\beta\}=-\epsilon_{\alpha\beta\gamma}J_\gamma$.

### Table A — translational / internal / vibrational / electronic momenta (canonical vs mechanical)

| Stage | coord. $q$ | canonical momentum $p=\partial L/\partial\dot q$ | $\{q,p\}$ | $=$ mass × velocity? | quantum operator | Hermitian on |
|---|---|---|---|---|---|---|
| 1 | $\mathbf R_j$ | $\mathbf P_j=M_j\dot{\mathbf R}_j$ | $\delta_{jk}\delta_{\alpha\beta}$ | **yes** | $-i\hbar\nabla_{\mathbf R_j}$ | flat |
| 1 | $\mathbf r_i$ | $\mathbf p_i=m_e\dot{\mathbf r}_i$ | $\delta_{ik}\delta_{\alpha\beta}$ | **yes** | $-i\hbar\nabla_{\mathbf r_i}$ | flat |
| 2 | $\mathbf R_{\mathrm{cm}}$ | $\mathbf P_{\mathrm{cm}}=M_{\mathrm{tot}}\dot{\mathbf R}_{\mathrm{cm}}$ | $\delta_{\alpha\beta}$ | **yes** (= total mom., conserved) | $-i\hbar\nabla_{\mathbf R_{\mathrm{cm}}}$ | flat |
| 2 | $\mathbf R'_j$ | $\mathbf P'_j=M_j\dot{\mathbf R}'_j$ | $\delta_{jk}\delta_{\alpha\beta}$† | **yes** | $-i\hbar\nabla_{\mathbf R'_j}$ | flat |
| 2 | $\mathbf r'_i$ | $\mathbf p'_i=m_e\dot{\mathbf r}'_i-\tfrac{m_e^2}{M_{\mathrm{tot}}}\sum_l\dot{\mathbf r}'_l$ | $\delta_{ik}\delta_{\alpha\beta}$ | **NO** — mass-polarization gap | $-i\hbar\nabla_{\mathbf r'_i}$ | flat |
| 3 | $\bar{\mathbf R}_j$ | $\bar{\mathbf P}_j=M_j(\boldsymbol\omega\times\bar{\mathbf R}_j+\dot{\bar{\mathbf R}}_j)=\mathbf S^{\mathsf T}\mathbf P'_j$ | $\delta_{jk}\delta_{\alpha\beta}$† | **yes** (transported vel.) | $-i\hbar\nabla_{\bar{\mathbf R}_j}$ | curved‡ |
| 3, 4 | $\Omega_s$ | $p_{\Omega_s}=(\mathsf E^{\mathsf T}\mathbf J)_s$ | $\delta_{st}$ | — (rotational; Table B) | $-i\hbar\partial_{\Omega_s}$ | **not** alone§ |
| 4 | $Q_b$ | $P_b=\dot Q_b+\boldsymbol\omega\!\cdot\!\boldsymbol\sigma_b$ | $\delta_{bc}$ | **NO** — Coriolis gap ($\ne\dot Q_b$) | $-i\hbar\partial_{Q_b}$ | flat $dQ$ (after reweighting) |
| 3, 4 | $\bar{\mathbf r}_i$ | $\bar{\mathbf p}_i=\mathbf S^{\mathsf T}\mathbf p'_i$ | $\delta_{ik}\delta_{\alpha\beta}$ | **NO** — inherits mass-pol. gap | $-i\hbar\nabla_{\bar{\mathbf r}_i}$ | flat $\prod d^3\bar r$ |

† on the constraint surface, with the dual gauge $\sum_j\mathbf P'_j=0$
(resp. $\sum_j\bar{\mathbf P}_j=0$). ‡ curved measure
$\sin\theta\,d\Omega\,\prod_j d^3\bar R_j$ before the shape is parametrized.
§ $\hat p_{\Omega_s}=-i\hbar\partial_{\Omega_s}$ is **not** Hermitian on
$\sin\theta\,d\Omega$; only the combination $\hat J_\alpha$ is (Eq. 4.27).

### Table B — rotational objects (angular velocity & angular momenta, body frame)

| Object | definition | canonical? | algebra $\{\,,\}$ (classical) / $[\,,\,]$ (quantum) | quantum operator | role |
|---|---|---|---|---|---|
| $\boldsymbol\omega$ | $\mathsf E(\Omega)\dot{\boldsymbol\Omega}$ | **no** — quasi-*velocity* | — | (velocity, not quantized directly) | conjugate to $\mathbf J_N$ in $L$ |
| $p_{\Omega_s}$ | $\partial L/\partial\dot\Omega_s=(\mathsf E^{\mathsf T}\mathbf J)_s$ | **yes** (true conjugate) | $\{\Omega_s,p_{\Omega_t}\}=\delta_{st}$ | $-i\hbar\partial_{\Omega_s}$ | rarely used directly |
| $\mathbf J$ | $J_\alpha=(\mathsf E^{-1})_{s\alpha}p_{\Omega_s}$ | **no** — frame-projected | $\{J_\alpha,J_\beta\}=-\epsilon_{\alpha\beta\gamma}J_\gamma$; $[\hat J_\alpha,\hat J_\beta]=-i\hbar\epsilon_{\alpha\beta\gamma}\hat J_\gamma$ (**anomalous**) | $-i\hbar(\mathsf E^{-1})_{s\alpha}\partial_{\Omega_s}$ (Eq. 4.24) | total ang. mom.; $\hat J_\alpha$ Hermitian on $\sin\theta\,d\Omega$ |
| $\mathbf J_N$ | $\mathbf I_N\boldsymbol\omega+\mathbf L_N^{\mathrm{vib}}=\mathbf J-\mathbf L_{\mathrm{elec}}$ | **no** — conjugate to quasi-vel. $\boldsymbol\omega$ | same anomalous family | $\hat J_\alpha-\hat L_{\mathrm{elec},\alpha}$ | "nuclear" ang. mom. |
| $\boldsymbol\pi$ | $\sum_b\boldsymbol\sigma_b(Q)\,P_b$ | **no** — coord-dependent combo of $P_b$ | $[\hat P_b,\hat\pi_\alpha]=-i\hbar\,\zeta^{\alpha}_{bc}\hat P_c\ne0$ | $\sum_b\sigma_{\alpha b}\hat P_b$ | vibrational (Coriolis) ang. mom. |
| $\mathbf L_{\mathrm{elec}}$ | $\sum_i\bar{\mathbf r}_i\times\bar{\mathbf p}_i$ | derived from a **canonical** pair | $\{L_\alpha,L_\beta\}=+\epsilon_{\alpha\beta\gamma}L_\gamma$; $[\hat L_\alpha,\hat L_\beta]=+i\hbar\epsilon_{\alpha\beta\gamma}\hat L_\gamma$ (**normal**) | $-i\hbar\,\epsilon_{\alpha\beta\gamma}\sum_i\bar r_{i\beta}\partial_{\bar r_{i\gamma}}$ | the well-behaved one |
| $\mathbf L_N^{\mathrm{vib}}$ | $\sum_jM_j\bar{\mathbf R}_j\times\dot{\bar{\mathbf R}}_j=\boldsymbol\pi-\boldsymbol\sigma\boldsymbol\sigma^{\mathsf T}\boldsymbol\omega$ | **no** — velocity-form | — | — | intermediate; eliminated in (3.20) |

Cross-sector (anti)commutators — each operator acts on a different coordinate
set unless noted:
$$
[\hat J_\alpha,\hat L_{\mathrm{elec},\beta}]=0,\qquad
[\hat J_\alpha,\hat\pi_\beta]=0,\qquad
[\hat J_\alpha,\hat P_b]=0,\qquad
[\hat P_b,\hat\pi_\alpha]=-i\hbar\,\zeta^{\alpha}_{bc}\hat P_c .
$$

*Reading the two tables.*

- **Every entry marked canonical is quantized by the naive $-i\hbar\,\partial_q$;**
  Hermiticity on the correct measure is the separate job of the Podolsky
  reweighting (see `tutorial_watson_pseudopotential.md`). The quasi objects of
  Table B are *derived* from the canonical momenta — never quantized on their own.
- **Two gaps, two origins.** The *mechanical* gaps
  ($\mathbf p'_i\ne m_e\dot{\mathbf r}'_i$, $P_b\ne\dot Q_b$) come from the
  coordinate *choice*; the *quasi* gap ($\mathbf J,\boldsymbol\pi,\boldsymbol\omega$)
  comes from using a rotating frame.
- **Only $\mathbf J$ is anomalous; $\mathbf L_{\mathrm{elec}}$ is not.** Both live
  in the body frame, but $\mathbf J$ is the space-fixed angular momentum
  *projected through* the angle-dependent $\mathsf E(\Omega)$ (⇒ anomalous sign),
  whereas $\mathbf L_{\mathrm{elec}}$ is assembled *directly* from the canonical
  pair $(\bar{\mathbf r}_i,\bar{\mathbf p}_i)$ (⇒ normal sign). Sign derivation:
  `tutorial_quantization.md` §5.
- **The anomaly never bites $H$:** in $H_{\mathrm{int}}$ the only place
  $\{\hat J,\hat J\}$ could enter is contracted with the symmetric
  $\mu_{\alpha\beta}$, where it vanishes regardless of sign.
