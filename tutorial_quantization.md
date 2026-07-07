# A guide to the two confusing steps in Section 4

This note explains, at "first-year quantum mechanics" level, the two things in
Section 4 / Appendix B (version `353fa80`) that look like black magic:

1. **The decomposition** $\ \hat J_\alpha \mu_{\alpha\beta} \hat J_\beta
   = \tfrac12 \mu_{\alpha\beta}\{\hat J_\alpha,\hat J_\beta\}
   + \tfrac12 \mu_{\alpha\beta}[\hat J_\alpha,\hat J_\beta]$ — why is this
   *allowed*, and what is it *for*?
2. **The anomalous commutator** $\ [\hat J_\alpha^{\rm BF},\hat J_\beta^{\rm BF}]
   = -i\hbar\,\epsilon_{\alpha\beta\gamma}\hat J_\gamma^{\rm BF}$ — where does
   the wrong sign come from, and why doesn't it break anything?

Along the way it answers the background question: *why does quantizing the
rovibrational Hamiltonian require any machinery at all* (Podolsky, symmetrized
products, …) instead of just writing $p \to -i\hbar\,\partial$?

**Assumed background:** you know what a commutator is, what a Hermitian
operator is, and how to integrate by parts. No differential geometry needed.

**Notation** (matching the paper):

| Symbol | Meaning | Acts on |
|---|---|---|
| $\hat J_\alpha$ | body-fixed (BF) total angular momentum | Euler angles $\phi,\theta,\chi$ |
| $\hat L_\alpha$ | BF electronic angular momentum | electron coordinates $\bar{\mathbf r}_i$ |
| $\hat P_b = -i\hbar\,\partial_{Q_b}$ | vibrational momentum | normal coordinates $Q_b$ |
| $\hat\pi_\alpha = \sum_b \sigma_{b\alpha}(Q)\hat P_b$ | vibrational angular momentum | $Q_b$ |
| $\mu_{\alpha\beta}(Q)$ | Watson reciprocal inertia, **symmetric**, real, depends on $Q$ only | (multiplication) |

Greek indices $\alpha,\beta,\gamma \in \{1,2,3\}$ are BF Cartesian directions;
repeated indices are summed.

---

## 1. Why quantization needs care at all

### 1.1 Classically, order doesn't matter. Quantum mechanically, it does.

The classical rovibrational kinetic energy is (Eq. 4.10a)

$$
T_{\rm rovib}
= \tfrac12 \sum_{\alpha\beta} A_\alpha\, \mu_{\alpha\beta}(Q)\, A_\beta
  + \tfrac12 \sum_b P_b^2,
\qquad
A_\alpha \equiv J_\alpha - L_{{\rm elec},\alpha} - \pi_\alpha .
$$

Classically, $J_\alpha$, $\mu_{\alpha\beta}(Q)$, $P_b$ are just **numbers**, so

$$
J_\alpha\,\mu_{\alpha\beta}\,\pi_\beta
= \pi_\beta\,\mu_{\alpha\beta}\,J_\alpha
= \mu_{\alpha\beta}\,J_\alpha\,\pi_\beta
\quad\text{(all the same thing).}
$$

Quantum mechanically these are **operators**, and operators need not commute.
The single classical expression corresponds to *several inequivalent* quantum
operators, differing by terms of order $\hbar$. "Quantizing" means picking one
of them — and the choice must satisfy one non-negotiable requirement:

> **The Hamiltonian must be Hermitian**, $\hat H^\dagger = \hat H$.
> Otherwise its eigenvalues need not be real and probability is not conserved.
> Hermiticity is the constraint that removes (most of) the ordering freedom.

Quick reminder: the adjoint $\hat O^\dagger$ is defined by
$\langle \phi | \hat O \psi\rangle = \langle \hat O^\dagger \phi | \psi \rangle$
for all states, where $\langle\phi|\psi\rangle = \int \phi^*\psi \,(\text{measure})$.
Two facts we use constantly:

- **(Rule 1)** $(\hat X \hat Y)^\dagger = \hat Y^\dagger \hat X^\dagger$ — the order *reverses*.
- **(Rule 2)** $\hat X\hat Y = \hat Y\hat X + [\hat X,\hat Y]$ — sliding one
  operator past another costs a commutator.

### 1.2 Hermiticity depends on the *measure* — a toy example

Here is the part rookies usually haven't seen. Whether $-i\hbar\,\partial$ is
Hermitian **depends on which integral defines your inner product.**

Take a free particle in 2D and switch to polar coordinates $(r,\theta)$. The
area element is $dx\,dy = r\,dr\,d\theta$ — the factor $r$ is the Jacobian.
So the inner product is

$$
\langle\phi|\psi\rangle = \int_0^\infty\!\!\int_0^{2\pi} \phi^*\,\psi\; r\,dr\,d\theta .
$$

Is $\hat p_r \equiv -i\hbar\,\partial_r$ Hermitian on this space? Integrate by
parts in $r$ (boundary terms vanish for normalizable states):

$$
\int_0^\infty \phi^* \left(-i\hbar\,\partial_r\psi\right) r\,dr
= \int_0^\infty \left(-i\hbar\,\partial_r\phi\right)^* \psi\; r\,dr
\;+\; \underbrace{\int_0^\infty (-i\hbar)\,\phi^*\psi\; dr}_{\text{from } \partial_r r = 1} .
$$

The extra term means $\hat p_r^\dagger = \hat p_r - \tfrac{i\hbar}{r} \neq \hat p_r$.
**The naive momentum is not Hermitian**, purely because the measure carries an
$r$. The fix is to shove half of the Jacobian's logarithmic derivative into
the operator:

$$
\hat p_r^{\rm herm} = -i\hbar\left(\partial_r + \frac{1}{2r}\right)
= -i\hbar\, r^{-1/2}\,\partial_r\, r^{1/2},
$$

which you can check *is* Hermitian on the measure $r\,dr$.

**This is the whole content of the Podolsky prescription.** In the molecular
problem, the coordinates (Euler angles, normal coordinates) are curvilinear
like $(r,\theta)$, and the Jacobian of the coordinate change is
$\sin\theta\,\sqrt{\det \mathbf I'(Q)}$ (Eq. 4.4) instead of $r$. Podolsky's
formula (Eq. 4.3) is nothing more than "build the kinetic operator so that it
is Hermitian on the *correct* measure," done once and for all in any
coordinates. When you then rescale the wavefunction to get rid of the
$Q$-dependent Jacobian (Eq. 4.5, exactly like the $r^{1/2}$ above), the
leftover $c$-number is the **Watson pseudopotential**
$\hat U = -\tfrac{\hbar^2}{8}\sum_\alpha \mu_{\alpha\alpha}$.

### 1.3 Division of labour

So there are **two separate jobs**, and it pays never to confuse them:

| Job | Tool | Output |
|---|---|---|
| Make the operator Hermitian on the curved measure | Podolsky / rescaling $\psi_W = (\det\mathbf I')^{1/4}\psi$ | Watson pseudopotential $\hat U(Q)$ |
| Order the non-commuting factors inside $\hat A\,\mu\,\hat A$ | commutator algebra (Appendix B) | Hermitian operators, **no** extra $c$-number |

The decomposition you asked about lives entirely in the **second** job. Let's
do it.

---

## 2. The identity $\hat X\hat Y = \tfrac12\{\hat X,\hat Y\} + \tfrac12[\hat X,\hat Y]$

### 2.1 Why you are *allowed* to write it

Because it is an **identity** — you add and subtract the same thing:

$$
\hat X\hat Y
= \tfrac12\hat X\hat Y + \tfrac12\hat Y\hat X
\;+\; \tfrac12\hat X\hat Y - \tfrac12\hat Y\hat X
= \tfrac12\{\hat X,\hat Y\} + \tfrac12[\hat X,\hat Y].
$$

That's the entire proof. Nothing was assumed about $\hat X$, $\hat Y$. It is
not an approximation, not a quantization postulate, not specific to angular
momentum. It is the operator version of splitting a function into even + odd
parts, or a matrix into symmetric + antisymmetric parts:

$$
M_{ij} = \underbrace{\tfrac12(M_{ij}+M_{ji})}_{\text{symmetric}}
       + \underbrace{\tfrac12(M_{ij}-M_{ji})}_{\text{antisymmetric}} .
$$

### 2.2 Why you would *want* to write it

Because the two pieces have opposite behaviour under the adjoint. Let
$\hat X$, $\hat Y$ be Hermitian. Using Rule 1:

$$
\{\hat X,\hat Y\}^\dagger = (\hat X\hat Y)^\dagger + (\hat Y\hat X)^\dagger
= \hat Y\hat X + \hat X\hat Y = +\{\hat X,\hat Y\}
\qquad\text{(Hermitian)},
$$

$$
[\hat X,\hat Y]^\dagger = \hat Y\hat X - \hat X\hat Y = -[\hat X,\hat Y]
\qquad\text{(anti-Hermitian)}.
$$

So the split is the operator analogue of writing a complex number as
$z = \mathrm{Re}\,z + i\,\mathrm{Im}\,z$:

- the **anticommutator** is the "real part" — manifestly a legitimate,
  Hermitian piece of a Hamiltonian;
- the **commutator** is the "imaginary part" — the *only* place where
  non-Hermiticity, ordering ambiguity, or $\hbar$-dependent surprises can
  hide. (Recall $[\hat x,\hat p]=i\hbar$: commutators of Hermitian operators
  are always $i\hbar$ times something Hermitian, hence "small," of quantum
  origin.)

**The split is a diagnostic tool.** You isolate the dangerous part, look at
it, and — in the lucky cases — watch it vanish. That's exactly what happens
next.

---

## 3. Applying it to $\hat J_\alpha \mu_{\alpha\beta} \hat J_\beta$ (term T1)

### Step 0: pull $\mu$ out front

$\mu_{\alpha\beta}(Q)$ is a function of the normal coordinates $Q$ only, while
$\hat J_\alpha$ differentiates Euler angles only. Operators living in
different coordinate sectors commute, so $[\hat J_\alpha, \mu_{\alpha\beta}]=0$
and

$$
\hat J_\alpha\, \mu_{\alpha\beta}\, \hat J_\beta
= \mu_{\alpha\beta}\, \hat J_\alpha \hat J_\beta .
$$

(This innocent step is *why* T1 is easy and the cross terms with $\hat\pi$
are not: $\hat\pi$ contains $\hat P_b$, which does **not** commute with
$\mu(Q)$.)

### Step 1: split

$$
\mu_{\alpha\beta}\,\hat J_\alpha \hat J_\beta
= \tfrac12\,\mu_{\alpha\beta}\{\hat J_\alpha,\hat J_\beta\}
+ \tfrac12\,\mu_{\alpha\beta}[\hat J_\alpha,\hat J_\beta] .
$$

Still an identity — nothing assumed yet.

### Step 2: the commutator piece dies by symmetry

Insert the algebra (anomalous sign, from Section 5 below — but watch: the
argument won't even care about the sign):

$$
\tfrac12\,\mu_{\alpha\beta}[\hat J_\alpha,\hat J_\beta]
= \tfrac12\,\mu_{\alpha\beta}\,(-i\hbar\,\epsilon_{\alpha\beta\gamma})\,\hat J_\gamma
= -\tfrac{i\hbar}{2}\,\Bigl(\textstyle\sum_{\alpha\beta}\mu_{\alpha\beta}\,\epsilon_{\alpha\beta\gamma}\Bigr)\hat J_\gamma .
$$

Now the punchline: **a symmetric object contracted with an antisymmetric
object over the same index pair is zero.** Watch it happen explicitly for
$\gamma = 3$ (only $\alpha\beta \in \{12, 21\}$ contribute, since
$\epsilon_{\alpha\beta 3}$ vanishes otherwise):

$$
\sum_{\alpha\beta}\mu_{\alpha\beta}\,\epsilon_{\alpha\beta 3}
= \mu_{12}\,\underbrace{\epsilon_{123}}_{+1} + \mu_{21}\,\underbrace{\epsilon_{213}}_{-1}
= \mu_{12} - \mu_{21} = 0,
$$

because $\mu_{12} = \mu_{21}$ (the inertia tensor is symmetric). Same for
$\gamma = 1, 2$. In general:

$$
\sum_{\alpha\beta}\mu_{\alpha\beta}\epsilon_{\alpha\beta\gamma}
\;\overset{\alpha\leftrightarrow\beta}{=}\;
\sum_{\alpha\beta}\mu_{\beta\alpha}\epsilon_{\beta\alpha\gamma}
= -\sum_{\alpha\beta}\mu_{\alpha\beta}\epsilon_{\alpha\beta\gamma}
\quad\Longrightarrow\quad
\text{it equals minus itself, hence } 0 .
$$

### Conclusion for T1

$$
\boxed{\;\hat J_\alpha\,\mu_{\alpha\beta}\,\hat J_\beta
= \tfrac12\,\mu_{\alpha\beta}\,\{\hat J_\alpha,\hat J_\beta\}\;}
$$

- It is **Hermitian exactly as written** (you can also see this directly:
  $(\hat J_\alpha\mu_{\alpha\beta}\hat J_\beta)^\dagger
  = \hat J_\beta\mu_{\alpha\beta}\hat J_\alpha$, relabel
  $\alpha\leftrightarrow\beta$, use $\mu$ symmetric — you get back the
  original).
- It produces **no extra $c$-number** — so the Watson pseudopotential does
  *not* come from here.
- The argument used only that $\mu$ is symmetric and that
  $[\hat J,\hat J] \propto \epsilon\,\hat J$. The **sign** of the algebra
  (normal $+i\hbar$ or anomalous $-i\hbar$) is irrelevant: the same result
  holds verbatim for the electronic term $\hat L\mu\hat L$ (T2), which has
  the normal sign.

**So, to your question "why can you decompose":** you can always decompose —
it's free. The reason the *paper* does it is that the decomposition cleanly
separates "the classical-looking, safely Hermitian part" (anticommutator)
from "the potentially dangerous quantum part" (commutator), and then symmetry
kills the dangerous part. It is a proof of *innocence*, not a computation of
something new.

---

## 4. Why the cross terms *must* be symmetrized (T5, T6)

For $J\mu J$ the ordering question was almost trivial. Here is where it gets
real. Consider the classical rotation–vibration coupling $-2\,J_\alpha
\mu_{\alpha\beta}\pi_\beta$. Quantum mechanically, is
$\hat J_\alpha \mu_{\alpha\beta}\hat\pi_\beta$ Hermitian? Take the adjoint
(Rule 1, everything individually Hermitian — $\hat\pi^\dagger = \hat\pi$ is
proved in Eq. B.6):

$$
(\hat J_\alpha \mu_{\alpha\beta}\hat\pi_\beta)^\dagger
= \hat\pi_\beta\, \mu_{\alpha\beta}\, \hat J_\alpha .
$$

Can we shuffle this back to the original order? $\hat J$ commutes with
everything here, fine — but $\hat\pi_\beta$ does **not** commute with
$\mu_{\alpha\beta}(Q)$:

$$
[\hat\pi_\beta, \mu_{\alpha\beta}]
= \sum_b \sigma_{b\beta}\,[\hat P_b, \mu_{\alpha\beta}]
= -i\hbar \sum_b \sigma_{b\beta}\,\partial_b \mu_{\alpha\beta} \;\neq\; 0 .
$$

So $(\hat J\mu\hat\pi)^\dagger \neq \hat J\mu\hat\pi$: **a single ordering is
not Hermitian**, and writing it alone would be a wrong (non-Hermitian)
Hamiltonian. But notice what the adjoint *is*: after relabeling
$\alpha\leftrightarrow\beta$ and using $\mu$ symmetric,

$$
(\hat J_\alpha \mu_{\alpha\beta}\hat\pi_\beta)^\dagger
= \hat\pi_\alpha\, \mu_{\alpha\beta}\, \hat J_\beta
\quad\text{— exactly the other ordering.}
$$

And for *any* operator $\hat X$, the combination $\hat X + \hat X^\dagger$ is
automatically Hermitian. So the correct quantization of the classical
$2\,J\mu\pi$ is the **symmetric sum of both orderings**:

$$
\hat J_\alpha \mu_{\alpha\beta}\hat\pi_\beta
+ \hat\pi_\alpha \mu_{\alpha\beta}\hat J_\beta
= \hat J_\alpha\,\{\mu_{\alpha\beta},\hat\pi_\beta\}
\qquad\text{(Hermitian by construction).}
$$

This is the same "real part" logic as Section 2 — symmetrize to get the
Hermitian combination — but now it's *forced*, not optional. Conveniently,
the classical expression already comes as the sum of the two (equal) classical
orderings $J\mu\pi + \pi\mu J$, so the symmetric quantum sum is its natural
image.

**Does this introduce arbitrariness?** In principle, different Hermitian
orderings of the same classical expression can differ by $O(\hbar^2)$
potential-like terms — quantization is genuinely not unique. The consistency
check of Appendix B is that the ordering adopted (the one inherited from the
Podolsky operator) produces **no leftover pure $c$-number from the kernel**:
every term in the expanded operators (Eqs. B.10, B.13) still carries at least
one momentum. Hence the *entire* $c$-number content of the theory is the
Watson pseudopotential, and it comes from the measure (job 1), not from
ordering (job 2). Clean division of labour, as promised.

---

## 5. The anomalous commutator, from scratch

### 5.1 The setup: two frames, one operator

- **Space-fixed (SF) frame:** lab axes, glued to the wall.
- **Body-fixed (BF) frame:** axes glued to the molecule; their orientation
  relative to the lab is given by the Euler angles
  $\Omega = (\phi,\theta,\chi)$, packaged in the rotation matrix
  $\mathbf S(\Omega)$ ($\det \mathbf S = 1$, $\mathbf S^{-1} = \mathbf S^{\mathsf T}$).

The total angular momentum is one physical observable; the two sets of
components are related by projection (Eq. 4.6):

$$
\hat J_\alpha^{\rm BF} = S_{\beta\alpha}(\Omega)\; \hat J_\beta^{\rm SF} .
$$

Here is the **crucial observation**, and the entire origin of the anomaly:

> $\hat J^{\rm SF}$ is a differential operator in the Euler angles, and
> $S_{\beta\alpha}(\Omega)$ is a *function* of those same Euler angles.
> So $\hat J^{\rm BF}$ is a product of two objects that **do not commute
> with each other.** When you commute two $\hat J^{\rm BF}$'s, you must
> differentiate the $\mathbf S$'s too — and those extra terms overwhelm and
> flip the normal algebra.

### 5.2 The two ingredients

**Ingredient 1 — SF components are normal.** In SF coordinates nothing exotic
is happening; the usual angular momentum algebra holds:

$$
[\hat J_\gamma^{\rm SF}, \hat J_\delta^{\rm SF}]
= +i\hbar\,\epsilon_{\gamma\delta\lambda}\,\hat J_\lambda^{\rm SF} .
$$

**Ingredient 2 — how $\hat J^{\rm SF}$ acts on $\mathbf S$.** Angular momentum
generates rotations: for any *vector* observable $\hat V_\delta$ (SF
components), $[\hat J_\gamma^{\rm SF}, \hat V_\delta] = i\hbar\,
\epsilon_{\gamma\delta\lambda}\hat V_\lambda$. Now fix the BF index $\beta$ of
$S_{\delta\beta}$ and read the column as a vector: it is the $\beta$-th body
axis written in SF components — an SF vector. Hence

$$
[\hat J_\gamma^{\rm SF}, S_{\delta\beta}]
= i\hbar\,\epsilon_{\gamma\delta\lambda}\, S_{\lambda\beta} .
$$

The rotation acts on the **SF index** ($\delta$) of $\mathbf S$; the BF index
$\beta$ is a spectator.

We'll also need two small technical facts:

- **(a)** Because $\det\mathbf S = 1$, three direction cosines contract to a
  single epsilon:
  $\epsilon_{\gamma\delta\lambda}S_{\gamma\alpha}S_{\delta\beta}S_{\lambda\nu}
  = \epsilon_{\alpha\beta\nu}$. (This is just the statement "the determinant
  is the totally antisymmetric contraction," i.e. proper rotations preserve
  cross products.)
- **(b)** Orthogonality inverts the projection:
  $\hat J_\lambda^{\rm SF} = S_{\lambda\nu} \hat J_\nu^{\rm BF}$.

### 5.3 The computation, every step

Write $\hat J \equiv \hat J^{\rm SF}$ for brevity. Expand using the product
rule for commutators, $[\hat A\hat B, \hat C\hat D]$ with the $\mathbf S$'s
being commuting functions:

$$
[\hat J_\alpha^{\rm BF}, \hat J_\beta^{\rm BF}]
= [S_{\gamma\alpha}\hat J_\gamma,\; S_{\delta\beta}\hat J_\delta]
= \underbrace{S_{\gamma\alpha}S_{\delta\beta}\,[\hat J_\gamma,\hat J_\delta]}_{\text{(A) bare SF algebra}}
+ \underbrace{S_{\gamma\alpha}\,[\hat J_\gamma, S_{\delta\beta}]\,\hat J_\delta}_{\text{(B) }\hat J\text{ hits the right }\mathbf S}
- \underbrace{S_{\delta\beta}\,[\hat J_\delta, S_{\gamma\alpha}]\,\hat J_\gamma}_{\text{(C) }\hat J\text{ hits the left }\mathbf S} .
$$

(If terms (B) and (C) were absent — if $\mathbf S$ were a constant matrix —
you'd get the normal algebra back. They are not absent.)

Insert the two ingredients. All three terms acquire a common factor
$i\hbar\,\epsilon_{\gamma\delta\lambda}$ (for (C), use
$-\epsilon_{\delta\gamma\lambda} = +\epsilon_{\gamma\delta\lambda}$):

$$
[\hat J_\alpha^{\rm BF}, \hat J_\beta^{\rm BF}]
= i\hbar\,\epsilon_{\gamma\delta\lambda}\Bigl(
  S_{\gamma\alpha}S_{\delta\beta}\,\hat J_\lambda
+ S_{\gamma\alpha}S_{\lambda\beta}\,\hat J_\delta
+ S_{\delta\beta}S_{\lambda\alpha}\,\hat J_\gamma \Bigr).
$$

Now convert every remaining $\hat J^{\rm SF}$ back to BF using fact (b),
$\hat J_\lambda = S_{\lambda\nu}\hat J_\nu^{\rm BF}$, so that each term
contains *three* $S$'s, and collapse them with fact (a):

$$
\epsilon_{\gamma\delta\lambda}S_{\gamma\alpha}S_{\delta\beta}S_{\lambda\nu}\,\hat J_\nu^{\rm BF}
= \epsilon_{\alpha\beta\nu}\,\hat J_\nu^{\rm BF},
\qquad
\epsilon_{\gamma\delta\lambda}S_{\gamma\alpha}S_{\lambda\beta}S_{\delta\nu}\,\hat J_\nu^{\rm BF}
= \epsilon_{\alpha\nu\beta}\,\hat J_\nu^{\rm BF},
$$
$$
\epsilon_{\gamma\delta\lambda}S_{\delta\beta}S_{\lambda\alpha}S_{\gamma\nu}\,\hat J_\nu^{\rm BF}
= \epsilon_{\nu\beta\alpha}\,\hat J_\nu^{\rm BF} .
$$

(In each case, just match which slot of the epsilon each BF index lands in.)
Add the three epsilons:

$$
\epsilon_{\alpha\beta\nu} + \epsilon_{\alpha\nu\beta} + \epsilon_{\nu\beta\alpha}
= \epsilon_{\alpha\beta\nu} - \epsilon_{\alpha\beta\nu} - \epsilon_{\alpha\beta\nu}
= -\,\epsilon_{\alpha\beta\nu},
$$

since each of the last two is one index-swap away from the first. Result:

$$
\boxed{\;[\hat J_\alpha^{\rm BF}, \hat J_\beta^{\rm BF}]
= -\,i\hbar\,\epsilon_{\alpha\beta\nu}\,\hat J_\nu^{\rm BF}.\;}
$$

**Bookkeeping of the sign:** the bare SF algebra (A) contributed $+1$ unit of
$\epsilon$; the two derivative-of-$\mathbf S$ terms (B) and (C) contributed
$-1$ each. Total $+1-1-1 = -1$. The anomaly is literally "the frame turns
under you twice as hard as the operator algebra turns."

### 5.4 Intuition for the minus sign

- **Passive vs. active.** $\hat J$ generates rotations of the *system*.
  But a BF component is measured along an axis *glued to the system*.
  Rotating the system while reading components off axes that co-rotate is
  equivalent to rotating the reference frame the *opposite* way — an active
  rotation seen passively. Opposite rotation ⇒ opposite sign in the algebra.
- **Left vs. right actions (for the mathematically curious).** SF components
  generate left multiplication on the rotation group, BF components generate
  right multiplication. The two actions have Lie algebras differing by a
  sign, and they mutually commute — which is why, additionally,
  $[\hat J^{\rm SF}_\alpha, \hat J^{\rm BF}_\beta] = 0$ (a fact that looks
  bizarre until you hear it this way).

### 5.5 Why the anomaly is harmless

1. **Same Casimir:** $\sum_\alpha (\hat J^{\rm BF}_\alpha)^2 =
   \sum_\alpha (\hat J^{\rm SF}_\alpha)^2 = \hat{\mathbf J}^2$ (use
   orthogonality of $\mathbf S$). Same total angular momentum, same spectrum
   $\hbar^2 J(J+1)$.
2. **Ladder operators swap roles:** with the anomalous sign,
   $\hat J^{\rm BF}_1 \pm i\hat J^{\rm BF}_2$ *lower/raise* the BF projection
   quantum number $k$ instead of raising/lowering it. A pure relabeling; the
   set of states $|J, k, m\rangle$ is untouched.
3. **In the Hamiltonian it never bites:** as shown in Section 3, the only
   place the $\hat J$ algebra could enter the kinetic operator is via
   $\mu_{\alpha\beta}[\hat J_\alpha, \hat J_\beta]$, and that dies against the
   symmetric $\mu$ *regardless of the sign*.

### 5.6 Why $\hat L_{\rm elec}$ is *not* anomalous

The electronic angular momentum is built **directly in the BF frame** from BF
electron coordinates and momenta,
$\hat L_\alpha = \sum_i (\bar{\mathbf r}_i \times \hat{\bar{\mathbf p}}_i)_\alpha$.
No $\mathbf S(\Omega)$ factor appears — nothing to differentiate — so the
ordinary $\hat x\hat p$ algebra gives the normal sign
$[\hat L_\alpha,\hat L_\beta] = +i\hbar\,\epsilon_{\alpha\beta\gamma}\hat L_\gamma$.
The anomaly is a property of *how you obtained the components* (projection
through an angle-dependent frame), not of the BF frame itself. And since
$\hat J^{\rm BF}$ touches only Euler angles while $\hat L$ touches only
electron coordinates, $[\hat J^{\rm BF}_\alpha, \hat L_\beta] = 0$.

---

## 6. FAQ

**Q: Is the $\tfrac12\{,\} + \tfrac12[,]$ decomposition an approximation?**
No. It's the identity "add and subtract $\tfrac12 \hat Y\hat X$." Always exact,
always allowed. Its value is diagnostic: anticommutator = manifestly Hermitian
part, commutator = the only place trouble can hide.

**Q: So where does the ordering choice actually happen?**
Not in T1/T2/T4 (all factors mutually commute or the commutator dies by
symmetry — no choice to make). It happens in the terms containing $\hat\pi$
(T5, T6), where Hermiticity *forces* the symmetric sum of the two orderings;
and globally in adopting the Podolsky/Laplace–Beltrami operator, which fixes
the measure and produces the Watson pseudopotential.

**Q: Could a different Hermitian ordering give different physics?**
Different orderings can differ by $O(\hbar^2)$ potential-like terms — the
classical limit alone does not pick a unique quantum theory. The
Laplace–Beltrami choice is the standard one: it is coordinate-independent
(built only from the metric), reduces correctly as $\hbar\to 0$, and Appendix B
verifies that with this choice the kernel ordering contributes *no* stray
$c$-numbers — the only $O(\hbar^2)$ potential is Watson's
$-\tfrac{\hbar^2}{8}\sum_\alpha\mu_{\alpha\alpha}$, from the measure.

**Q: Does the anomalous sign change any observable?**
No. Same $\hat{\mathbf J}^2$, same spectrum, ladder operators relabeled, and
the sign never survives contraction with the symmetric $\mu_{\alpha\beta}$ in
the Hamiltonian.

**Q: One-sentence summary of each confusion?**
- *Decomposition:* splitting $\hat J_\alpha\hat J_\beta$ into anticommutator +
  commutator is free algebra; it's done to show the commutator part is killed
  by the symmetry of $\mu$, so $\hat J\mu\hat J$ is Hermitian as written with
  no hidden extra terms.
- *Anomalous algebra:* BF components are SF components multiplied by
  Euler-angle-dependent direction cosines; commuting them, the derivatives of
  the direction cosines contribute $-2$ units of $\epsilon$ against the bare
  $+1$, flipping the sign.

---

*Note:* on the current `main` branch, Section 4 was rewritten in a
"coordinate-native" way that avoids invoking the anomalous BF algebra
altogether (the ordering results are obtained directly in terms of
$\partial_{Q_b}$ and Euler-angle operators). This tutorial explains the
classic route as it stands at commit `353fa80`, which is also the standard
textbook presentation (Watson 1968; Bunker & Jensen).
