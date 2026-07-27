# The rigid-rotor approximation: freezing the shape in the full molecular Hamiltonian

This note takes the two central results of the paper — the exact **classical**
internal Hamiltonian (3.22) and its **quantized Watson form** (4.41) — and asks
what survives under the **rigid-rotor approximation**: treat the molecule as a
rigid body that *rotates but does not vibrate*. We apply the approximation to
each Hamiltonian in turn and find

1. the classical rotor collapses to Euler's rigid body,
   $\tfrac12\boldsymbol\omega^{T}\mathsf I^{\mathrm{eq}}\boldsymbol\omega$;
2. the quantum rotor collapses to the textbook asymmetric top
   $A\hat J_a^2+B\hat J_b^2+C\hat J_c^2$;
3. the **Watson pseudopotential** — the one genuinely quantum term in (4.41),
   with no classical analogue — quietly *disappears*, for a reason that is a
   direct payoff of the coordinate-native derivation.

**Assumed background.** You have met the classical (3.22) and quantum (4.41)
Hamiltonians (both are reproduced below), and you know what the Watson
reciprocal-inertia tensor $\mathrm{\mu}$, the Coriolis vectors
$\boldsymbol\sigma_b$, and the vibrational angular momentum $\boldsymbol\pi$ are.
No new machinery is introduced — we only *switch things off*.

---

## 1. The two starting Hamiltonians

Everything below flows from the two boxed results of the paper. Both are written
in the body-fixed Eckart frame, in normal coordinates $Q_b$, with the total
centre-of-mass motion already removed.

**Classical internal Hamiltonian (3.22):**

$$
H_{\mathrm{int}}=\underbrace{\tfrac12\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}-\boldsymbol\pi\bigr)^{T}\mathrm{\mu}(Q)\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}-\boldsymbol\pi\bigr)}_{\text{rotation + Coriolis}}
\;+\;\underbrace{\tfrac12\sum_b P_b^2}_{\text{vibration}}
\;+\;\underbrace{\sum_i\frac{\bar{\mathbf p}_i^{\,2}}{2m_e}+\frac{1}{2M_N}\Bigl(\sum_i\bar{\mathbf p}_i\Bigr)^2}_{H_e}
\;+\;V .
$$

**Quantum Watson Hamiltonian (4.41),** acting on
$L^2(\mathcal C,\,d\mu_W)$ with the flat Watson measure
$d\mu_W=\sin\theta\,d\phi\,d\theta\,d\chi\,\prod_b dQ_b\,\prod_i d^3\bar{\mathbf r}_i$:

$$
\hat H_{\mathrm{int}}=\tfrac12\sum_{\alpha\beta}\bigl(\hat J_\alpha-\hat L_{\mathrm{elec},\alpha}-\hat\pi_\alpha\bigr)\mu_{\alpha\beta}(Q)\bigl(\hat J_\beta-\hat L_{\mathrm{elec},\beta}-\hat\pi_\beta\bigr)
+\tfrac12\sum_b\hat P_b^2
+\hat H_e
\;\underbrace{-\;\frac{\hbar^2}{8}\sum_\alpha\mu_{\alpha\alpha}(Q)}_{\text{Watson pseudopotential }\hat U}
+\,V .
$$

The two are *term-for-term identical* except for one thing: the quantum
Hamiltonian carries the extra piece

$$
\hat U(Q)=-\frac{\hbar^2}{8}\sum_\alpha\mu_{\alpha\alpha}(Q)
\qquad\text{(4.38),}
$$

a pure $\hbar^2$ $c$-number with **no classical counterpart**. It was not put in
by hand: it fell out of quantizing the *vibrational* kinetic energy on the curved
shape-space measure $\sqrt\gamma\,dQ$ (with $\gamma=\det\mathsf I'(Q)$), when the
$Q$-derivatives were symmetrized. Keep an eye on it — it is the term whose fate
is most interesting under the approximation.

The ingredients that carry the shape dependence are

$$
\boldsymbol\sigma_b(Q)=\sum_a\boldsymbol\zeta_{ab}\,Q_a,\qquad
\boldsymbol\pi=\sum_b\boldsymbol\sigma_b\,P_b\ \ (\hat\pi_\alpha=\textstyle\sum_b\sigma_{\alpha b}\hat P_b),\qquad
\mathrm{\mu}(Q)=\bigl(\mathsf I_N(Q)-\mathrm{\sigma}\mathrm{\sigma}^{T}\bigr)^{-1},
$$

with $\boldsymbol\sigma_b$ the **Coriolis vectors** (3.17a), linear in the
displacements $Q_a$, and $\mathrm{\sigma}$ the $3\times f$ matrix whose columns
are the $\boldsymbol\sigma_b$.

---

## 2. What "rigid rotor" means

Physically: pretend the nuclei are frozen into their equilibrium geometry and
glued together, so the molecule can only **reorient** in space — it cannot bend,
stretch, or twist. All vibrational amplitudes are set to zero and stay there.

Precisely, the approximation is the pair of statements

$$
\boxed{\;Q_b\equiv 0\quad\text{(clamp the shape at equilibrium)}\qquad\text{and}\qquad P_b\equiv 0\ \ (\hat P_b\equiv 0)\quad\text{(no vibrational motion).}\;}
$$

Everything else follows mechanically. Because the Coriolis vectors are **linear**
in $Q$ (that is the content of (3.17a)),

$$
\boldsymbol\sigma_b(0)=\sum_a\boldsymbol\zeta_{ab}\,\underbrace{Q_a}_{0}=0
\;\;\Longrightarrow\;\;
\mathrm{\sigma}\mathrm{\sigma}^{T}\big|_{Q=0}=0 ,
$$

so the four shape-dependent objects reduce to constants:

| object | full | rigid ($Q\to0,\ P\to0$) |
|---|---|---|
| vibrational KE | $\tfrac12\sum_b P_b^2$ | $0$ |
| vibrational ang. mom. | $\boldsymbol\pi=\sum_b\boldsymbol\sigma_b P_b$ | $0$ (both $\boldsymbol\sigma_b\to0$ **and** $P_b\to0$) |
| reciprocal inertia | $\mathrm{\mu}(Q)=(\mathsf I_N-\mathrm{\sigma}\mathrm{\sigma}^{T})^{-1}$ | $\mathrm{\mu}^{\mathrm{eq}}\equiv(\mathsf I^{\mathrm{eq}})^{-1}$, constant |
| potentials | $V_{NN}(Q),\,V_{Ne}(Q,\bar{\mathbf r})$ | $V_{NN}^{\mathrm{eq}},\,V_{Ne}(0,\bar{\mathbf r})$ |

Here $\mathsf I^{\mathrm{eq}}\equiv\mathsf I_N(0)$ is the **equilibrium nuclear
inertia tensor** — an ordinary constant $3\times3$ matrix — and
$\mathrm{\mu}^{\mathrm{eq}}=(\mathsf I^{\mathrm{eq}})^{-1}$ is its inverse. Note
that with $\mathrm{\sigma}\mathrm{\sigma}^{T}\to0$ the *modified* inertia
$\mathsf I'=\mathsf I_N-\mathrm{\sigma}\mathrm{\sigma}^{T}$ and the *bare* inertia
$\mathsf I_N$ coincide at equilibrium: at the rigid level there is no distinction
between them.

---

## 3. The classical rigid rotor

Drop the vibrational KE, set $\boldsymbol\pi\to0$, and freeze
$\mathrm{\mu}\to\mathrm{\mu}^{\mathrm{eq}}$ in (3.22):

$$
\boxed{\;H_{\mathrm{rr}}=\tfrac12\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}\bigr)^{T}\mathrm{\mu}^{\mathrm{eq}}\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}\bigr)+H_e^{\mathrm{eq}}+V^{\mathrm{eq}} .\;}
$$

The vibrational sector is gone entirely; what remains is a rotation of a rigid
frame (still dragging the electrons along through the shift
$\mathbf J\to\mathbf J-\mathbf L_{\mathrm{elec}}$) plus the electronic Hamiltonian
and the potential evaluated at the frozen geometry.

**Recognising Euler's rigid body.** Recall the key identity (3.19),
$\mathbf J-\mathbf L_{\mathrm{elec}}=\mathsf I_N(Q)\,\boldsymbol\omega+\mathbf L_N$,
where $\mathbf L_N=\boldsymbol\pi-\mathrm{\sigma}\mathrm{\sigma}^{T}\boldsymbol\omega$
is the nuclear vibrational angular momentum. At the rigid level
$\boldsymbol\pi\to0$ and $\mathrm{\sigma}\mathrm{\sigma}^{T}\to0$, so
$\mathbf L_N\to0$ and

$$
\mathbf J-\mathbf L_{\mathrm{elec}}=\mathsf I^{\mathrm{eq}}\,\boldsymbol\omega .
$$

Substituting back, the rotational kernel becomes the familiar rigid-body kinetic
energy,

$$
\tfrac12\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}\bigr)^{T}\mathrm{\mu}^{\mathrm{eq}}\bigl(\mathbf J-\mathbf L_{\mathrm{elec}}\bigr)
=\tfrac12\,\boldsymbol\omega^{T}\mathsf I^{\mathrm{eq}}\boldsymbol\omega ,
$$

i.e. Euler's $\tfrac12\sum_\alpha I_\alpha^{\mathrm{eq}}\omega_\alpha^2$ in the
principal-axis frame. If the electrons are in a closed-shell state so that
$\mathbf L_{\mathrm{elec}}$ averages to zero (see §6), $H_e^{\mathrm{eq}}+V^{\mathrm{eq}}$
becomes a constant electronic energy and we are left with nothing but a rigid top
spinning in free space. This is the classical picture the name "rigid rotor"
refers to.

---

## 4. The quantum rigid rotor

Now do the *same* three deletions in the quantum Hamiltonian (4.41):
$\tfrac12\sum_b\hat P_b^2\to0$, $\hat\pi_\alpha\to0$, and
$\mu_{\alpha\beta}(Q)\to\mu^{\mathrm{eq}}_{\alpha\beta}$ (constant). The result is

$$
\boxed{\;\hat H_{\mathrm{rr}}=\tfrac12\sum_{\alpha\beta}\bigl(\hat J_\alpha-\hat L_{\mathrm{elec},\alpha}\bigr)\mu^{\mathrm{eq}}_{\alpha\beta}\bigl(\hat J_\beta-\hat L_{\mathrm{elec},\beta}\bigr)+\hat H_e^{\mathrm{eq}}+V^{\mathrm{eq}}\;}
$$

acting now on the much smaller space
$L^2\bigl(\sin\theta\,d\phi\,d\theta\,d\chi\bigr)\otimes L^2\bigl(\prod_i d^3\bar{\mathbf r}_i\bigr)$:
the $\prod_b dQ_b$ factor of the measure has been removed along with the
vibrational coordinates.

**Where did the pseudopotential go?** Freezing $Q\to0$ turns
$\hat U(Q)=-\tfrac{\hbar^2}{8}\sum_\alpha\mu_{\alpha\alpha}(Q)$ into

$$
\hat U(0)=-\frac{\hbar^2}{8}\sum_\alpha\mu^{\mathrm{eq}}_{\alpha\alpha}
=-\frac{\hbar^2}{8}\operatorname{tr}\mathrm{\mu}^{\mathrm{eq}},
$$

a single **constant number**. A constant added to a Hamiltonian shifts every
eigenvalue by the same amount; it has no effect on spectra, transitions, or
dynamics, and is absorbed into the zero of energy. So in the rigid rotor the
Watson pseudopotential **contributes nothing physical** — it has effectively
disappeared. Section 5 explains why this is not a coincidence.

That leaves the naked rotational kernel
$\tfrac12\mu^{\mathrm{eq}}_{\alpha\beta}(\hat J_\alpha-\hat L_{\mathrm{elec},\alpha})(\hat J_\beta-\hat L_{\mathrm{elec},\beta})$
plus the electronic Hamiltonian at the equilibrium geometry — no vibrational
derivatives, no shape-dependent operators, no $\hbar^2$ residue.

---

## 5. Why the pseudopotential disappears — and why freezing commutes with quantizing

The disappearance of $\hat U$ is the cleanest possible illustration of *where the
pseudopotential came from in the first place*.

Recall the division of labour established in the paper:

- The **rotational block** (4.29) is Hermitian on the rotational measure
  $\sin\theta\,d\phi\,d\theta\,d\chi$ **with no pseudopotential**. The reason is
  geometric: the rotation generators $\hat J_\alpha$ are divergence-free on the
  Euler-angle sphere (the $\sin\theta$ weight is exactly the one they preserve),
  so symmetrizing $\tfrac12\mu_{\alpha\beta}\hat J_\alpha\hat J_\beta$ produces no
  leftover $c$-number.
- The **entire** pseudopotential $\hat U$ came from the *vibrational* block —
  specifically from the $Q$-derivatives of the shape Jacobian $\sqrt\gamma(Q)$
  when the bare $\hat P_b$ were made Hermitian on $\sqrt\gamma\,dQ$.

The rigid-rotor approximation switches off exactly the machinery that *generates*
$\hat U$: with no vibrational coordinates there is no $\gamma(Q)$, no
$Q$-derivatives to symmetrize, and hence no pseudopotential — only the frozen
number $\hat U(0)$ is left behind, and even that is inert.

This has a sharp consequence for the **order of operations**. There are two ways
to build a rigid-rotor operator:

- **Freeze, then quantize.** Impose $Q=0$ *classically* first, reducing the
  configuration space to the Euler angles (times the electrons), then quantize on
  $\sin\theta\,d\Omega$. Because the rotational block needs no pseudopotential,
  this gives exactly
  $\tfrac12\mu^{\mathrm{eq}}_{\alpha\beta}\hat J_\alpha\hat J_\beta$.
- **Quantize, then freeze.** Quantize the full problem to (4.41) — pseudopotential
  and all — then set $Q=0$. This gives the same rotational operator, plus the
  irrelevant constant $\hat U(0)$.

$$
\boxed{\;\text{freeze}\circ\text{quantize}\;=\;\text{quantize}\circ\text{freeze}\quad(\text{up to a physically inert constant}).\;}
$$

The two routes **commute**. This is *not* automatic: quantizing a system with a
frozen (constrained) curvilinear coordinate generically produces a
constraint-induced pseudopotential — the well-known ordering ambiguity of
constrained quantization — and freezing would then fail to commute with
quantizing. Here it does commute, precisely because the surviving degree of
freedom is *rotation*, whose block was pseudopotential-free from the start. The
rigid rotor is "clean" for the same reason the rotational block was clean:
the unimodular geometry of the rotation group.

---

## 6. Born–Oppenheimer and the fate of $\mathbf L_{\mathrm{elec}}$

There is a genuine asymmetry in §§2–4 worth pausing on. The rigid-rotor
approximation killed the vibrational angular momentum $\boldsymbol\pi$ outright,
yet the **electronic** angular momentum $\mathbf L_{\mathrm{elec}}$ is still
sitting in the boxed rotor of §4. Most textbooks, by contrast, write the rigid
rotor as the bare $\mathbf J^2/2I$ with no electronic term at all. The resolution
is that $\mathbf L_{\mathrm{elec}}$ is removed by a *different* approximation —
Born–Oppenheimer — that we have not yet made.

### 6.1 Why $\boldsymbol\pi$ dies but $\mathbf L_{\mathrm{elec}}$ survives

"Rigid rotor" is a statement about the **nuclear** sector alone: freeze the shape,
$Q_b=0$, $\hat P_b=0$. The vibrational angular momentum is a nuclear-shape object,
$\boldsymbol\pi=\sum_b\boldsymbol\sigma_b\hat P_b$, and it vanishes *twice over* —
once because $\hat P_b\to0$, once because $\boldsymbol\sigma_b(0)\to0$. But
$\hat{\mathbf L}_{\mathrm{elec}}=\sum_i\bar{\mathbf r}_i\times\hat{\bar{\mathbf p}}_i$
is built from the **electrons**, which keep moving no matter how rigidly the
nuclei are clamped. Freezing the geometry does nothing to it. Removing
$\mathbf L_{\mathrm{elec}}$ therefore requires a separate approximation acting on
the *electronic* degrees of freedom — and that is exactly what Born–Oppenheimer
supplies. The rotor of §4 is the honest **pre-BO** object: electrons are still
fully dynamical and rotationally coupled through the shift
$\hat{\mathbf J}\to\hat{\mathbf J}-\hat{\mathbf L}_{\mathrm{elec}}$.

### 6.2 What BO does: the $(\hat{\mathbf J}-\hat{\mathbf L}_{\mathrm{elec}})^2$ decomposition

BO selects a single electronic state $|\phi_n\rangle$, the solution of the
clamped-nuclei electronic problem at the frozen geometry,
$\hat H_{\mathrm{el}}^{\mathrm{eq}}|\phi_n\rangle=E_n^{\mathrm{eq}}|\phi_n\rangle$,
and evaluates the nuclear operator in that state. Because $\hat{\mathbf J}$
(differential in the Euler angles) and $\hat{\mathbf L}_{\mathrm{elec}}$
(differential in the body-frame electron coordinates) act on **disjoint
variables**, they commute and the rotational kernel expands with no ordering
subtlety:

$$
\tfrac12(\hat{\mathbf J}-\hat{\mathbf L}_{\mathrm{elec}})^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}(\hat{\mathbf J}-\hat{\mathbf L}_{\mathrm{elec}})
=\underbrace{\tfrac12\hat{\mathbf J}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf J}}_{\text{pure rotation}}
\;-\;\underbrace{\hat{\mathbf J}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf L}_{\mathrm{elec}}}_{\text{rotation–electronic}}
\;+\;\underbrace{\tfrac12\hat{\mathbf L}_{\mathrm{elec}}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf L}_{\mathrm{elec}}}_{\text{electronic}} .
$$

Taking the electronic expectation $\langle\phi_n|\cdots|\phi_n\rangle$ (the
BO/adiabatic projection) sends the three terms to three very different places:

| term | after $\langle\phi_n\vert\cdot\vert\phi_n\rangle$ | fate |
|---|---|---|
| $\tfrac12\hat{\mathbf J}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf J}$ | untouched ($\hat{\mathbf J}$ acts on Euler angles; $\langle\phi_n\vert\phi_n\rangle=1$) | **the rotor**, $\mathbf J^2/2I$ |
| $-\hat{\mathbf J}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\langle\phi_n\vert\hat{\mathbf L}_{\mathrm{elec}}\vert\phi_n\rangle$ | $\propto\langle\hat{\mathbf L}_{\mathrm{elec}}\rangle_n$ | **vanishes** (non-degenerate state, §6.3) |
| $\tfrac12\langle\phi_n\vert\hat{\mathbf L}_{\mathrm{elec}}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf L}_{\mathrm{elec}}\vert\phi_n\rangle$ | a $c$-number | constant, absorbed into $E_n^{\mathrm{eq}}$ |

Only the first term keeps any dependence on $\hat{\mathbf J}$. So after BO the
rotational Hamiltonian is $\tfrac12\hat{\mathbf J}^{\mathsf T}\mathrm{\mu}^{\mathrm{eq}}\hat{\mathbf J}$
plus constants — the bare $\mathbf J^2/2I$ of the textbooks. **This is the missing
step:** the books quietly assume BO (so the projection has been done) *and* a
non-degenerate electronic ground state (so the middle term is zero).

### 6.3 Quenching: why the cross term vanishes

The rotation–electronic term dies by a symmetry, not an estimate. For a
**non-degenerate** spatial electronic state, time-reversal invariance of the
electronic Hamiltonian lets us choose $\phi_n$ **real**. But the angular-momentum
operator $\hat{\mathbf L}_{\mathrm{elec}}=-i\hbar\sum_i\bar{\mathbf r}_i\times\nabla_i$
is purely *imaginary*. Hence

$$
\langle\phi_n|\hat{\mathbf L}_{\mathrm{elec}}|\phi_n\rangle
=(\text{real }\phi_n)\times(-i\hbar)\times(\text{real integral})
=\ \text{purely imaginary},
$$

while the expectation value of a Hermitian operator must be **real**. A number
that is simultaneously real and imaginary is $0$:

$$
\boxed{\;\langle\phi_n|\hat{\mathbf L}_{\mathrm{elec}}|\phi_n\rangle=0\qquad(\text{non-degenerate state}).\;}
$$

This is the **quenching of orbital angular momentum** — exact at this level, and
the precise reason a closed-shell $^1\Sigma$ (or $^1A$) molecule has no
first-order electronic contribution to its rotational Hamiltonian.

### 6.4 When $\mathbf L_{\mathrm{elec}}$ does *not* disappear

The clean $\mathbf J^2/2I$ is special to the non-degenerate, adiabatic case. It
fails in two ways:

- **Open-shell or orbitally degenerate states** ($^2\Pi$, $^3\Sigma$, radicals,
  Renner–Teller and Jahn–Teller systems): $\phi_n$ *cannot* be chosen real, the
  quenching argument breaks, $\langle\hat{\mathbf L}_{\mathrm{elec}}\rangle\neq0$,
  and the cross term survives as genuine **rotation–electronic coupling**
  ($\Lambda$-doubling, spin-rotation, the various Hund's coupling cases). Here the
  bare rigid rotor is simply wrong.
- **Second order, even for a closed shell**: BO discards the *off-diagonal*
  couplings $\langle\phi_n|\hat{\mathbf L}_{\mathrm{elec}}|\phi_m\rangle$
  ($m\neq n$) that mix electronic states. Restoring them at second order gives
  small corrections — the rotational $g$-factor (the molecule's magnetic response
  to rotation) and a tiny renormalization of $A,B,C$. So the statement that
  $\mathbf L_{\mathrm{elec}}$ disappears is exact only to leading adiabatic order.

In short: **rigid-rotor kills $\boldsymbol\pi$; Born–Oppenheimer (plus a
non-degenerate state) kills $\mathbf L_{\mathrm{elec}}$.** They are independent
approximations, and only after *both* does the Hamiltonian reduce to the bare
$\mathbf J^2/2I$ found in most textbooks.

---

## 7. Landing it: the asymmetric top

Choose the body frame to be the **principal-axis frame** of the equilibrium
geometry, in which $\mathsf I^{\mathrm{eq}}=\operatorname{diag}(I_a,I_b,I_c)$ and
therefore $\mathrm{\mu}^{\mathrm{eq}}=\operatorname{diag}(1/I_a,1/I_b,1/I_c)$.
Assume the standard case of §6 — a non-degenerate (closed-shell) electronic
ground state under Born–Oppenheimer, so $\mathbf L_{\mathrm{elec}}$ is quenched
and the electronic energy $E_0^{\mathrm{eq}}$ is a constant we drop. What remains
is the pure nuclear rotor.

**The asymmetric-top Hamiltonian.** With $\mathrm{\mu}^{\mathrm{eq}}$ diagonal and
$\hat{\mathbf L}_{\mathrm{elec}}$ gone,

$$
\boxed{\;\hat H_{\mathrm{rr}}=\frac{\hat J_a^2}{2I_a}+\frac{\hat J_b^2}{2I_b}+\frac{\hat J_c^2}{2I_c}\;=\;A\,\hat{j}_a^2+B\,\hat{j}_b^2+C\,\hat{j}_c^2,\qquad A=\frac{\hbar^2}{2I_a},\ \ B=\frac{\hbar^2}{2I_b},\ \ C=\frac{\hbar^2}{2I_c}.\;}
$$

This is the standard rigid **asymmetric top**. The first form uses the
paper's operators $\hat J_\alpha$ (which carry a factor $\hbar$); factoring that
out, $\hat J_\alpha=\hbar\,\hat{j}_\alpha$ with the *dimensionless* body-fixed
components $\hat{j}_\alpha$ (spectrum $\hat{\mathbf j}^2=J(J+1)$), gives the
spectroscopists' form with **rotational constants** $A\ge B\ge C$ (largest for
the smallest moment, $I_a\le I_b\le I_c$). The $\hat{j}_a,\hat{j}_b,\hat{j}_c$ are
the *body-fixed* components (the Euler-angle operators of (4.24)); they obey the
anomalous commutation relations $[\hat{j}_a,\hat{j}_b]=-i\,\hat{j}_c$ and act on
the Wigner basis $|J,K,M\rangle$. Two familiar special cases fall out
immediately:

- **Symmetric top** ($I_a=I_b\ne I_c$, so $A=B$):
  $\hat H=B\,\hat{\mathbf j}^2+(C-B)\hat{j}_c^2$, with eigenvalues
  $B\,J(J+1)+(C-B)K^2$.
- **Spherical top** ($A=B=C$): $B\,J(J+1)$, $(2J+1)^2$-fold degenerate.
- **Linear molecule**: the moment about the molecular axis vanishes,
  $I_c\to0$ so $C\to\infty$; that forces $K=0$ (no rotation about the axis), one
  Euler angle ($\chi$) becomes redundant, and the spectrum collapses to the
  diatomic $B\,J(J+1)$.

**What the rigid rotor threw away.** Reading (3.22)/(4.41) backwards, the terms we
deleted are exactly the leading corrections to the rigid top:

- the $Q$-dependence of $\mathrm{\mu}(Q)$ — expanding
  $\mathrm{\mu}(Q)=\mathrm{\mu}^{\mathrm{eq}}+(\partial\mathrm{\mu}/\partial Q_b)Q_b+\dots$
  restores **centrifugal distortion** (the rotor stretches as it spins);
- the vibrational angular momentum $\boldsymbol\pi$ — restores **Coriolis
  coupling** between rotation and vibration;
- the vibrational kinetic energy $\tfrac12\sum_b\hat P_b^2$ (with the potential's
  $Q$-dependence) — restores the **vibrational levels** themselves;
- and the Watson pseudopotential $\hat U(Q)$ — restores the small $\hbar^2$
  shift that accompanies large-amplitude shape change.

The rigid rotor is thus the $Q=0$, $\hat P=0$ vertex of the full theory; the rest
of (4.41) is the systematic vibration–rotation expansion around it.
