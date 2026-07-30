# Why $U_{\rm vib}+U_{\rm Cor}$ collapse into Watson's $-\frac{\hbar^2}{8}\sum_\alpha\mu_{\alpha\alpha}$

This note explains, at "first-year quantum mechanics" level, the most magical
single step in Section 4 / Appendix B: the reweighting of the vibrational
kinetic operator leaves behind a messy pile of derivatives of a determinant
[Eqs. (4.34), (B.18), together (B.19)],

$$
U \;=\; \frac{\hbar^{2}}{8}\Bigl[
\partial_b\bigl(g^{QQ}_{bc}\,\partial_c\ln\gamma\bigr)
+\tfrac14\,g^{QQ}_{bc}\,(\partial_b\ln\gamma)(\partial_c\ln\gamma)\Bigr],
\qquad
\gamma=\det\mathsf I',\quad
g^{QQ}=\mathbb 1+\mathrm{\sigma}^{T}\mathrm{\mu}\mathrm{\sigma},
$$

and this pile — second derivatives of $\ln\det$, squares of first derivatives,
Coriolis dressings — collapses into **one clean trace**:

$$
\boxed{\;U=-\frac{\hbar^{2}}{8}\sum_{\alpha=1}^{3}\mu_{\alpha\alpha}(Q)\;}
\qquad\text{— Watson's pseudopotential, Eq. (4.38).}
$$

Why on earth should that happen? The answer has two parts, and both deserve
the word *miracle* until you see the proof:

- **Miracle 1 (all derivatives become constants).** The Eckart conditions
  force the modified inertia into Watson's exact closed form
  $\mathsf I'=A\,(\mathsf I^{\rm eq})^{-1}A$ with $A$ *linear* in $Q$ — so every
  derivative of $\gamma$ and $\mathrm{\mu}$ reduces to constant matrices.
- **Miracle 2 (the assembly telescopes).** When the five groups of $U$ are
  evaluated with the Eckart *sum rules*, every complicated term cancels
  pairwise, and the sole survivor is $-\operatorname{tr}\mathrm{\mu}$.

**Nothing in this note is quoted.** Every relation used before the final
assembly is derived here from the definitions: the inertia derivative, the
completeness relation, all five sum rules, the closed form, and the five
groups. Section 13 checks all of it numerically.

The plan:

| § | Content |
|---|---|
| 0 | Notation and conventions |
| 1 | Warm-up: where the residue comes from (a 1D computation you can do on a napkin) |
| 2 | The kinematic input, from scratch: $\mathsf I_N$, $\mathrm{\sigma}$, $\mathrm{\zeta}$, $\mathsf I'$ |
| 3 | The residue in full and its five-group split |
| 4 | **The inertia derivative $\mathsf I'_{N,b}$, derived** — and the Eckart-symmetry lemma |
| 5 | The Eckart completeness relation, derived |
| 6 | The two workhorse $\epsilon$-contractions |
| 7 | The five sum rules, each proven in full |
| 8 | Miracle 1: Watson's closed form |
| 9 | Jacobi's formula → all derivatives are constants |
| 10 | The five groups, evaluated |
| 11 | The two remaining $V$-traces |
| 12 | Miracle 2: the cancellation table |
| 13 | A runnable Python check (no packages needed) |
| 14 | FAQ |

**Assumed background:** product rule, integration by parts, matrix
multiplication, determinants and inverses of small matrices, the Levi-Civita
symbol. No differential geometry, no group theory.

---

## 0. Notation and conventions

Repeated indices are summed unless stated otherwise.

| Index | Ranges over |
|---|---|
| $j,k$ | nuclei, $1\dots N_{\rm nucl}$ |
| $b,c,d$ | vibrational modes, $1\dots f=3N_{\rm nucl}-6$ |
| $\alpha,\beta,\gamma,\delta,\rho,\tau,\lambda,\nu$ and $p,q,r,w,x,y,z,m,n$ | Cartesian, $1,2,3$ |

$\epsilon_{\alpha\beta\gamma}$ is the Levi-Civita symbol; $\mathbb 1$ is the
$3\times3$ identity. All $3\times3$ matrices below are written without indices
where possible; $M^{T}$ is the transpose.

| Symbol | Meaning |
|---|---|
| $\mathbf e_j=\sqrt{M_j}\,\bar{\mathbf R}^{\rm eq}_j$ | mass-weighted equilibrium positions |
| $\mathbf u_j(Q)=\mathbf e_j+\sum_b\mathbf l_{jb}Q_b$ | mass-weighted positions at shape $Q$ |
| $\mathbf l_{jb}$ | Eckart mode vectors: orthonormal, $\perp$ translations and rotations |
| $\mathsf I_N(Q)$ | nuclear inertia tensor at shape $Q$ |
| $\mathrm{\sigma}(Q)$, $\zeta^{\alpha}_{bc}$ | Coriolis vectors and the constant Coriolis coefficients |
| $\mathsf I'=\mathsf I_N-\mathrm{\sigma}\mathrm{\sigma}^{T}$ | **Coriolis-modified** inertia (prime, *no* subscript) |
| $\mathsf I^{\rm eq}=\mathsf I_N(0)=\mathsf I'(0)$ | equilibrium inertia |
| $\mathrm{\mu}=(\mathsf I')^{-1}$, $\mathrm{\mu}^{\rm eq}=(\mathsf I^{\rm eq})^{-1}$ | reciprocal inertias |
| $\gamma=\det\mathsf I'$ | shape Jacobian |
| $\displaystyle \mathsf I'_{N,b}\equiv\frac{\partial\mathsf I_N}{\partial Q_b}\bigg\|_{Q=0}$ | **constant** inertia derivatives (§4) |
| $S=\sum_j\mathbf e_j\mathbf e_j^{T}$ | equilibrium second-moment matrix, $\mathsf I^{\rm eq}=(\operatorname{tr}S)\mathbb 1-S$ |
| $\mathsf C_b=\sum_j\mathbf e_j\mathbf l_{jb}^{T}$ | equilibrium–displacement overlap (§4.2) |
| $s_b=\sum_j\mathbf e_j\!\cdot\!\mathbf l_{jb}=\operatorname{tr}\mathsf C_b$ | its trace |
| $A(Q)=\mathsf I^{\rm eq}+\tfrac12\sum_b\mathsf I'_{N,b}\,Q_b$ | Watson's linear matrix |
| $\mathsf B_b=A^{-1}\mathsf I'_{N,b}$ | the combination every derivative produces |
| $\Lambda_b=\partial_b\ln\gamma$ | log-derivative of the Jacobian |
| $Y=\sum_j\mathbf u_j\mathbf e_j^{T}$, $\ \boldsymbol\kappa$, $\ V[M]^{\gamma}$ | auxiliaries defined in §4.4, §7.2, §10.1 |

> **⚠ Read the prime carefully.** Two different objects wear a prime, and the
> subscript is what tells them apart:
>
> - $\mathsf I'$ — prime alone, **no** subscript — is the *Coriolis-modified*
>   inertia $\mathsf I_N-\mathrm{\sigma}\mathrm{\sigma}^{T}$. The prime is part
>   of the name; it is **not** a derivative.
> - $\mathsf I'_{N,b}$ — prime *plus* the subscript $N,b$ — is the
>   $\partial/\partial Q_b$ **derivative of $\mathsf I_N$**, evaluated at
>   $Q=0$. Here the prime *is* the derivative mark and the subscript says
>   "of $\mathsf I_N$, with respect to mode $b$".
>
> They are never the same object and never appear in the same slot:
> $\mathsf I'$ is a function of $Q$; $\mathsf I'_{N,b}$ is a constant matrix.
> A bare $\mathsf I'$ in this note *always* means the modified inertia.

**Why not $a^b$?** The paper and Watson's 1968 original write these constants
$a^{b}_{\alpha\beta}$ (Watson: $a^{k}_{\alpha\beta}$). That is standard in the
vibration–rotation literature, but it reads badly: in $a^{b}$ the letter $a$
is a *symbol name*, not an index, while $b$ *is* an index — so the expression
looks like a rank-2 tensor component when it is really "the $b$-th matrix in a
list". $\mathsf I'_{N,b}$ says what the object is. **Dictionary to the paper:**

$$
\mathsf I'_{N,b}\;=\;a^{b}\ \text{[Eq. (B.21)]},\qquad
\mathsf C_b=E^{b}\ \text{[§B.6, Step 1]},\qquad
\mathsf B_b=B^{b},\qquad
\mathsf I^{\rm eq}=\mathsf I^{e},\qquad
\mathrm{\mu}^{\rm eq}=\mathrm{\mu}^{e}.
$$

The mode label is a **subscript on every object** in this note
($\mathsf I'_{N,b}$, $\mathsf C_b$, $\mathsf B_b$, $s_b$, $\Lambda_b$,
$\sigma_{\alpha b}$, $Q_b$); Cartesian indices are Greek and always in the
subscript slot *after* it, e.g. $(\mathsf I'_{N,b})_{\alpha\beta}$. Superscripts
are reserved for Cartesian labels of families of matrices ($V[M]^{\gamma}$) and
for $\zeta^{\alpha}_{bc}$, which the literature fixes.

---

## 1. Warm-up: where does the residue come from? (1D)

Take a single coordinate $Q$ with a curved measure $w(Q)\,dQ$, and the
kinetic operator that is Hermitian on it (this is Eq. (4.5) in one dimension
with $g^{QQ}=1$, $\sqrt g=w$):

$$
\hat T=-\frac{\hbar^{2}}{2}\,\frac1w\,\partial\bigl(w\,\partial\,\cdot\bigr).
$$

Move to the flat measure: $\psi_W=w^{1/2}\psi$, so operators get conjugated,
$\hat T_W=w^{1/2}\,\hat T\,w^{-1/2}$. Compute it with nothing but the product
rule, one derivative at a time ($'\equiv\partial_Q$):

$$
\partial\bigl(w^{-1/2}\psi_W\bigr)
= w^{-1/2}\Bigl[\psi_W'-\tfrac12\bigl(\ln w\bigr)'\psi_W\Bigr],
$$

$$
\frac1w\,\partial\Bigl(w\cdot w^{-1/2}\bigl[\psi_W'-\tfrac12(\ln w)'\psi_W\bigr]\Bigr)
= w^{-1/2}\Bigl[\psi_W''
\;-\;\tfrac14\bigl((\ln w)'\bigr)^{2}\psi_W
\;-\;\tfrac12(\ln w)''\,\psi_W\Bigr]
$$

(the two cross terms $\pm\tfrac12(\ln w)'\psi_W'$ cancel — check it!). So

$$
\hat T_W=-\frac{\hbar^{2}}{2}\,\partial^{2}
\;+\;\underbrace{\frac{\hbar^{2}}{8}\Bigl[\,2(\ln w)''+\bigl((\ln w)'\bigr)^{2}\Bigr]}_{\text{the residue}} .
$$

**This is the entire mechanism.** The flat-measure operator is the naive one
*plus a multiplicative leftover* built from derivatives of $\ln w$. In the
molecule, $w=\sqrt\gamma$ per mode and there are $f$ modes plus the Coriolis
dressing $\mathrm{\sigma}^T\mathrm{\mu}\mathrm{\sigma}$; substituting
$w=\gamma^{1/2}$ above gives exactly the paper's per-mode residue
$\frac{\hbar^2}{8}[\partial_b^2\ln\gamma+\tfrac14(\partial_b\ln\gamma)^2]$
of Eq. (4.34). Nothing about *this* step is special to molecules. The magic
is what the sum over modes does next.

---

## 2. The kinematic input, from scratch

Everything in this note is built from four definitions and three conditions.
Here they are, with the elementary consequences we will use constantly.

### 2.1 Shape coordinates and the Eckart conditions

The nuclei sit at body-frame positions $\bar{\mathbf R}_j$; mass-weight them,
$\mathbf u_j\equiv\sqrt{M_j}\,\bar{\mathbf R}_j$, and expand about equilibrium
in *rectilinear* normal coordinates:

$$
\boxed{\;\mathbf u_j(Q)=\mathbf e_j+\sum_b\mathbf l_{jb}\,Q_b\;}
\qquad
\mathbf e_j=\sqrt{M_j}\,\bar{\mathbf R}^{\rm eq}_j ,
\qquad
\frac{\partial\mathbf u_j}{\partial Q_b}=\mathbf l_{jb}\ \text{(constant).}
\tag{2.1}
$$

The last statement is the whole reason rectilinear coordinates work: the
displacement vectors do not themselves depend on $Q$. The $3N_{\rm nucl}$
numbers $\{l_{jb,\alpha}\}$ obey **three** conditions:

$$
\textbf{(E1) translational:}\quad \sum_j\sqrt{M_j}\,\mathbf l_{jb}=0,
\qquad\text{together with}\quad \sum_j\sqrt{M_j}\,\mathbf e_j=0
\ \ \text{(centre of mass at the origin);}
\tag{2.2}
$$

$$
\textbf{(E2) rotational:}\quad \sum_j\mathbf e_j\times\mathbf l_{jb}=0;
\qquad\qquad
\textbf{(E3) orthonormality:}\quad \sum_j\mathbf l_{jb}\!\cdot\!\mathbf l_{jc}=\delta_{bc}.
\tag{2.3}
$$

(E1) and (E2) are the Eckart conditions — they *define* the body frame; (E3)
is a normalisation choice. Note that (E1) and (2.1) give immediately

$$
\sum_j\sqrt{M_j}\,\mathbf u_j(Q)=0\qquad\text{for every }Q,
\tag{2.4}
$$

i.e. the body frame origin stays at the nuclear centre of mass at every shape.
That innocuous fact is what kills the translational term of the completeness
relation in *every* sum rule of §7.

### 2.2 Inertia, Coriolis vectors, modified inertia

The nuclear inertia tensor is, in mass-weighted variables,

$$
\mathsf I_N(Q)=\sum_j\Bigl[\,|\mathbf u_j|^{2}\,\mathbb 1-\mathbf u_j\mathbf u_j^{T}\Bigr],
\qquad\text{i.e.}\quad
(\mathsf I_N)_{\alpha\beta}=\sum_j\bigl[u_j^2\,\delta_{\alpha\beta}-u_{j\alpha}u_{j\beta}\bigr].
\tag{2.5}
$$

(The familiar $\sum_jM_j(R_j^2\delta_{\alpha\beta}-R_{j\alpha}R_{j\beta})$ with
$\mathbf u_j=\sqrt{M_j}\bar{\mathbf R}_j$.) The Coriolis vectors are the
vibrational angular momenta per unit mode velocity,

$$
\boldsymbol\sigma_b(Q)\equiv\sum_j\mathbf u_j\times\mathbf l_{jb},
\qquad
\sigma_{\alpha b}=\epsilon_{\alpha\gamma\delta}\sum_ju_{j\gamma}\,l_{jb,\delta},
\tag{2.6}
$$

and the **modified inertia** that appears everywhere in the Watson Hamiltonian
is

$$
\mathsf I'(Q)\equiv\mathsf I_N-\mathrm{\sigma}\mathrm{\sigma}^{T},
\qquad
\mathrm{\mu}\equiv(\mathsf I')^{-1},\qquad\gamma\equiv\det\mathsf I' .
\tag{2.7}
$$

Here $\mathrm{\sigma}$ is the $3\times f$ matrix with entries $\sigma_{\alpha b}$,
so $(\mathrm{\sigma}\mathrm{\sigma}^{T})_{\alpha\beta}=\sum_b\sigma_{\alpha b}\sigma_{\beta b}$.

### 2.3 $\mathrm{\sigma}$ is *linear* in $Q$; the constants $\zeta$

Substitute (2.1) into (2.6):

$$
\sigma_{\alpha b}
=\underbrace{\Bigl(\sum_j\mathbf e_j\times\mathbf l_{jb}\Bigr)_\alpha}_{=0\ \text{by (E2)}}
+\sum_cQ_c\Bigl(\sum_j\mathbf l_{jc}\times\mathbf l_{jb}\Bigr)_\alpha
\qquad\Longrightarrow\qquad
\boxed{\;\sigma_{\alpha b}=\sum_c\zeta^{\alpha}_{cb}\,Q_c\;}
\tag{2.8}
$$

with the **constant** Coriolis coefficients

$$
\zeta^{\alpha}_{bc}\equiv\Bigl(\sum_j\mathbf l_{jb}\times\mathbf l_{jc}\Bigr)_\alpha
=\epsilon_{\alpha\gamma\delta}\sum_jl_{jb,\gamma}\,l_{jc,\delta},
\qquad
\zeta^{\alpha}_{bc}=-\zeta^{\alpha}_{cb},\qquad \zeta^{\alpha}_{bb}=0 .
\tag{2.9}
$$

Two consequences we will need:

$$
\partial_b\,\sigma_{\alpha c}=\zeta^{\alpha}_{bc},
\qquad\qquad
\sum_b\partial_b\,\sigma_{\alpha b}=\sum_b\zeta^{\alpha}_{bb}=0 .
\tag{2.10}
$$

The second one — the Coriolis field is *divergence-free in the elementary
sense* — is what makes $\hat\pi_\alpha$ Hermitian without symmetrisation
[(B.5)] and, in §3, lets us drop one whole term from the product rule.

Finally, since $\mathrm{\sigma}$ is linear and $\mathsf I_N$ quadratic in $Q$,
$\mathsf I'$ is **exactly quadratic** in $Q$ — no truncation anywhere in this
note. And at $Q=0$, $\mathrm{\sigma}=0$, so

$$
\mathsf I'(0)=\mathsf I_N(0)=\mathsf I^{\rm eq}=(\operatorname{tr}S)\mathbb 1-S,
\qquad S=\sum_j\mathbf e_j\mathbf e_j^{T},
\tag{2.11}
$$

where the last equality is (2.5) at $Q=0$ with $\sum_j|\mathbf e_j|^2=\operatorname{tr}S$.

---

## 3. The residue in full, split into five groups

Adding the identity part (4.34) and the dressed part (B.18) gives (B.19),
quoted at the top. Now expand it. Write $G_{bc}\equiv(\mathrm{\sigma}^{T}\mathrm{\mu}\mathrm{\sigma})_{bc}
=\sigma_{\alpha b}\mu_{\alpha\beta}\sigma_{\beta c}$, so
$g^{QQ}_{bc}=\delta_{bc}+G_{bc}$, and introduce the 3-vector

$$
v_\alpha\equiv\sum_c\sigma_{\alpha c}\Lambda_c
\qquad\Longrightarrow\qquad
G_{bc}\Lambda_c=\sigma_{\alpha b}\,(\mathrm{\mu}\mathbf v)_\alpha .
\tag{3.1}
$$

The identity part of (B.19) gives $N_0$ and $N_4$'s $\delta_{bc}$ share
directly. For the dressed part, apply the product rule to
$\partial_b\bigl[\sigma_{\alpha b}(\mathrm{\mu}\mathbf v)_\alpha\bigr]$:

$$
\partial_b\bigl[\sigma_{\alpha b}(\mathrm{\mu}\mathbf v)_\alpha\bigr]
=\underbrace{\Bigl(\sum_b\partial_b\sigma_{\alpha b}\Bigr)}_{=0\ \text{by (2.10)}}(\mathrm{\mu}\mathbf v)_\alpha
+\sigma_{\alpha b}\,\partial_b(\mathrm{\mu}\mathbf v)_\alpha ,
$$

and expand the surviving term with
$\partial_bv_\beta=\zeta^{\beta}_{bc}\Lambda_c+\sigma_{\beta c}\,\partial_b\Lambda_c$
[from (3.1) and (2.10)]:

$$
\sigma_{\alpha b}\,\partial_b(\mathrm{\mu}\mathbf v)_\alpha
=\underbrace{\sigma_{\alpha b}(\partial_b\mathrm{\mu})_{\alpha\beta}v_\beta}_{N_1}
+\underbrace{\sigma_{\alpha b}\mu_{\alpha\beta}\zeta^{\beta}_{bc}\Lambda_c}_{N_2}
+\underbrace{\sigma_{\alpha b}\mu_{\alpha\beta}\sigma_{\beta c}\,\partial_b\Lambda_c}_{N_3\,=\,G_{bc}\partial_b\Lambda_c}.
$$

Collecting, and writing the $\tfrac14 g^{QQ}\Lambda\Lambda$ piece as
$\tfrac14\sum_b\Lambda_b^2+\tfrac14\mathbf v^{T}\mathrm{\mu}\mathbf v$, we get
Eq. (B.20):

$$
\frac{8}{\hbar^{2}}\,U
=\underbrace{\sum_b\Bigl[\partial_b\Lambda_b+\tfrac14\Lambda_b^2\Bigr]}_{N_{0}\ \text{(identity)}}
+\underbrace{\sigma_{\alpha b}(\partial_b\mathrm{\mu})_{\alpha\beta}v_\beta}_{N_{1}}
+\underbrace{\sigma_{\alpha b}\mu_{\alpha\beta}\zeta^{\beta}_{bc}\Lambda_c}_{N_{2}}
+\underbrace{G_{bc}\,\partial_b\Lambda_c}_{N_{3}}
+\underbrace{\tfrac14\,\mathbf v^{T}\mathrm{\mu}\mathbf v}_{N_{4}} .
\tag{3.2}
$$

Our job: evaluate all five and show the sum is $-\operatorname{tr}\mathrm{\mu}$.
Everything hinges on being able to *compute* $\Lambda_b=\partial_b\ln\det\mathsf I'$
— which brings us to the first miracle.

---

## 4. The inertia derivative $\mathsf I'_{N,b}$

This is the object the whole construction is built on, so we derive it in
full — nothing here is quoted.

### 4.1 Differentiate $\mathsf I_N$

Start from the definition (2.5) and differentiate with respect to $Q_b$,
using $\partial\mathbf u_j/\partial Q_b=\mathbf l_{jb}$ from (2.1). Do it
index by index. For the first term,

$$
\frac{\partial}{\partial Q_b}\Bigl[\sum_j u_j^{2}\Bigr]\delta_{\alpha\beta}
=\frac{\partial}{\partial Q_b}\Bigl[\sum_j u_{j\lambda}u_{j\lambda}\Bigr]\delta_{\alpha\beta}
=2\Bigl(\sum_j u_{j\lambda}\,l_{jb,\lambda}\Bigr)\delta_{\alpha\beta}
=2\Bigl(\sum_j\mathbf u_j\!\cdot\!\mathbf l_{jb}\Bigr)\delta_{\alpha\beta},
$$

and for the second, by the product rule on $u_{j\alpha}u_{j\beta}$,

$$
\frac{\partial}{\partial Q_b}\Bigl[\sum_ju_{j\alpha}u_{j\beta}\Bigr]
=\sum_j\bigl[l_{jb,\alpha}u_{j\beta}+u_{j\alpha}l_{jb,\beta}\bigr].
$$

Hence, at **any** shape $Q$,

$$
\frac{\partial\mathsf I_N}{\partial Q_b}
=2\Bigl(\sum_j\mathbf u_j\!\cdot\!\mathbf l_{jb}\Bigr)\mathbb 1
-\sum_j\bigl[\mathbf l_{jb}\mathbf u_j^{T}+\mathbf u_j\mathbf l_{jb}^{T}\bigr].
\tag{4.1}
$$

Now set $Q=0$, where $\mathbf u_j=\mathbf e_j$. Define the **overlap matrix**

$$
\boxed{\;\mathsf C_b\equiv\sum_j\mathbf e_j\mathbf l_{jb}^{T},
\qquad (\mathsf C_b)_{\alpha\beta}=\sum_je_{j\alpha}\,l_{jb,\beta},
\qquad s_b\equiv\sum_j\mathbf e_j\!\cdot\!\mathbf l_{jb}=\operatorname{tr}\mathsf C_b\;}
\tag{4.2}
$$

(the trace identity is immediate: $\operatorname{tr}\mathsf C_b=\sum_j e_{j\alpha}l_{jb,\alpha}$).
Since $\sum_j\mathbf l_{jb}\mathbf e_j^{T}=\mathsf C_b^{T}$, Eq. (4.1) at $Q=0$ reads

$$
\boxed{\;
\mathsf I'_{N,b}\;\equiv\;\frac{\partial\mathsf I_N}{\partial Q_b}\bigg\|_{Q=0}
\;=\;2s_b\,\mathbb 1-\mathsf C_b-\mathsf C_b^{T}\;}
\tag{4.3}
$$

**— this is the relation the paper states as (B.21).** In components,
$(\mathsf I'_{N,b})_{\alpha\beta}=\sum_j\bigl[2(\mathbf e_j\!\cdot\!\mathbf l_{jb})\delta_{\alpha\beta}
-l_{jb,\alpha}e_{j\beta}-e_{j\alpha}l_{jb,\beta}\bigr]$. Note that (4.3) is
manifestly **symmetric** and **constant** (it is built entirely from
$\mathbf e_j$ and $\mathbf l_{jb}$, neither of which depends on $Q$), whereas
(4.1) at general $Q$ is *not* constant — this is precisely why the definition
fixes $Q=0$.

### 4.2 The Eckart-symmetry lemma: $\mathsf C_b$ is symmetric

Take the rotational Eckart condition (E2) and write it in components. Using
$(\mathbf a\times\mathbf b)_\gamma=\epsilon_{\gamma\alpha\beta}a_\alpha b_\beta$,

$$
0=\Bigl(\sum_j\mathbf e_j\times\mathbf l_{jb}\Bigr)_\gamma
=\epsilon_{\gamma\alpha\beta}\sum_je_{j\alpha}l_{jb,\beta}
=\epsilon_{\gamma\alpha\beta}\,(\mathsf C_b)_{\alpha\beta}.
\tag{4.4}
$$

Contracting an $\epsilon$ with a matrix extracts exactly its antisymmetric
part (and nothing else): $\epsilon_{\gamma\alpha\beta}M_{\alpha\beta}
=\epsilon_{\gamma\alpha\beta}\tfrac12(M_{\alpha\beta}-M_{\beta\alpha})$, and this
vanishes for all $\gamma$ iff $M=M^{T}$. Hence

$$
\boxed{\;\text{(E2)}\iff \mathsf C_b=\mathsf C_b^{T}\;}
\tag{4.5}
$$

**the rotational Eckart condition is *nothing but* the statement that the
equilibrium–displacement overlap is a symmetric matrix.** That is the seed of
everything downstream.

### 4.3 The lemma in the form we use

Put (4.5) into (4.3):

$$
\mathsf I'_{N,b}=2s_b\mathbb 1-2\,\mathsf C_b
\qquad\Longleftrightarrow\qquad
\boxed{\;\mathsf C_b=s_b\,\mathbb 1-\tfrac12\,\mathsf I'_{N,b}\;}
\tag{4.6}
$$

and take a trace of the left form, using $\operatorname{tr}\mathsf C_b=s_b$ and
$\operatorname{tr}\mathbb 1=3$:

$$
\operatorname{tr}\mathsf I'_{N,b}=6s_b-2s_b=4s_b
\qquad\Longrightarrow\qquad
\boxed{\;s_b=\tfrac14\operatorname{tr}\mathsf I'_{N,b}\;}
\tag{4.7}
$$

Innocent-looking — but this is the seed of everything: *the entire overlap
matrix between equilibrium and displacement is fixed by the single constant
matrix $\mathsf I'_{N,b}$.* Every appearance of $\mathbf e_j$ and
$\mathbf l_{jb}$ from now on will be traded for $\mathsf I'_{N,b}$ through (4.6).

### 4.4 Two corollaries: $A(Q)$ and $Y(Q)$

Define Watson's linear matrix

$$
\boxed{\;A(Q)\equiv\mathsf I^{\rm eq}+\tfrac12\sum_b\mathsf I'_{N,b}\,Q_b\;}
\qquad\Longrightarrow\qquad
\partial_bA=\tfrac12\,\mathsf I'_{N,b}\ \ \text{(constant)},\qquad A=A^{T}.
\tag{4.8}
$$

*First corollary.* Directly from (4.8),

$$
\sum_bQ_b\,\mathsf I'_{N,b}=2\bigl(A-\mathsf I^{\rm eq}\bigr).
\tag{4.9}
$$

Taking a trace and using (4.7) plus $\operatorname{tr}\mathsf I^{\rm eq}=3\operatorname{tr}S-\operatorname{tr}S=2\operatorname{tr}S$:

$$
\sum_bs_bQ_b=\tfrac14\operatorname{tr}\!\sum_bQ_b\mathsf I'_{N,b}
=\tfrac12\operatorname{tr}\bigl(A-\mathsf I^{\rm eq}\bigr)
=\tfrac12\operatorname{tr}A-\operatorname{tr}S .
\tag{4.10}
$$

*Second corollary.* Define $Y\equiv\sum_j\mathbf u_j\mathbf e_j^{T}$, i.e.
$Y_{\alpha\beta}=\sum_ju_{j\alpha}e_{j\beta}$. Expanding with (2.1) and using
the symmetry (4.5),

$$
Y=S+\sum_bQ_b\sum_j\mathbf l_{jb}\mathbf e_j^{T}
=S+\sum_bQ_b\,\mathsf C_b^{T}
=S+\sum_bQ_b\,\mathsf C_b .
$$

Now insert (4.6), then (4.9) and (4.10), then $\mathsf I^{\rm eq}=(\operatorname{tr}S)\mathbb 1-S$:

$$
Y=S+\Bigl(\tfrac12\operatorname{tr}A-\operatorname{tr}S\Bigr)\mathbb 1-\bigl(A-\mathsf I^{\rm eq}\bigr)
=S+\tfrac12(\operatorname{tr}A)\mathbb 1-(\operatorname{tr}S)\mathbb 1-A+(\operatorname{tr}S)\mathbb 1-S ,
$$

$$
\boxed{\;Y=\sum_j\mathbf u_j\mathbf e_j^{T}=\tfrac12(\operatorname{tr}A)\,\mathbb 1-A\;}
\qquad\Longrightarrow\qquad
\operatorname{tr}Y=\tfrac32\operatorname{tr}A-\operatorname{tr}A=\tfrac12\operatorname{tr}A .
\tag{4.11}
$$

$Y$ is symmetric — a fact we use repeatedly to kill $\epsilon$-contractions.

---

## 5. The Eckart completeness relation

Every sum rule in §7 is a contraction over the mode index $b$. The one tool
that evaluates such contractions is completeness. Here it is, derived.

### 5.1 The basis

The mass-weighted displacement space is $\mathbb R^{3N_{\rm nucl}}$: a vector
$\mathbf w$ has components $w_{j\alpha}$, with inner product
$\langle\mathbf w,\mathbf w'\rangle=\sum_{j\alpha}w_{j\alpha}w'_{j\alpha}$.
It is spanned by three families:

| family | count | components |
|---|---|---|
| vibrations | $f=3N_{\rm nucl}-6$ | $l_{jb,\alpha}$ |
| translations | 3 | $T^{(\nu)}_{j\alpha}=\sqrt{M_j/M_N}\;\delta_{\nu\alpha}$, $\ M_N=\sum_jM_j$ |
| rotations | 3 | $t^{(\gamma)}_{j\alpha}=(\hat{\mathbf n}_\gamma\times\mathbf e_j)_\alpha=\epsilon_{\gamma\rho\alpha}e_{j\rho}$ |

*Mutual orthogonality.* Vibrations $\perp$ translations is (E1); vibrations
$\perp$ rotations is (E2) [since $\sum_{j\alpha}t^{(\gamma)}_{j\alpha}l_{jb,\alpha}
=\epsilon_{\gamma\rho\alpha}\sum_je_{j\rho}l_{jb,\alpha}=0$ by (4.4)];
translations $\perp$ rotations because
$\sum_j\sqrt{M_j}\,(\hat{\mathbf n}_\gamma\times\mathbf e_j)
=\hat{\mathbf n}_\gamma\times\sum_j\sqrt{M_j}\,\mathbf e_j=0$ by the
centre-of-mass condition (2.2).

*Translations are already orthonormal:*
$\sum_{j\alpha}T^{(\nu)}_{j\alpha}T^{(\nu')}_{j\alpha}
=\delta_{\nu\nu'}\sum_jM_j/M_N=\delta_{\nu\nu'}$.

### 5.2 The rotational Gram matrix is $\mathsf I^{\rm eq}$

Rotations are **not** normalised. Their Gram matrix is, using
$\epsilon_{\gamma\rho\alpha}\epsilon_{\delta\tau\alpha}
=\delta_{\gamma\delta}\delta_{\rho\tau}-\delta_{\gamma\tau}\delta_{\rho\delta}$
(the BAC–CAB identity in index form):

$$
\sum_{j\alpha}t^{(\gamma)}_{j\alpha}t^{(\delta)}_{j\alpha}
=\epsilon_{\gamma\rho\alpha}\epsilon_{\delta\tau\alpha}\sum_je_{j\rho}e_{j\tau}
=\bigl(\delta_{\gamma\delta}\delta_{\rho\tau}-\delta_{\gamma\tau}\delta_{\rho\delta}\bigr)S_{\rho\tau}
=\delta_{\gamma\delta}\operatorname{tr}S-S_{\delta\gamma}
=\mathsf I^{\rm eq}_{\gamma\delta}.
\tag{5.1}
$$

**The equilibrium inertia is literally the Gram matrix of the rotational
directions.** That is why $\mathrm{\mu}^{\rm eq}=(\mathsf I^{\rm eq})^{-1}$ — and
not some other matrix — appears in the completeness relation.

### 5.3 The projector onto a non-orthonormal family

For any linearly independent set $\{\mathbf w_A\}$ with Gram matrix
$g_{AB}=\langle\mathbf w_A,\mathbf w_B\rangle$, the orthogonal projector onto
their span is $P=\sum_{AB}\mathbf w_A\,(g^{-1})_{AB}\,\mathbf w_B^{T}$. (Check:
$P\mathbf w_C=\sum_{AB}\mathbf w_A(g^{-1})_{AB}g_{BC}=\mathbf w_C$, and $P$ is
symmetric, so it is *the* projector.) With $g=\mathsf I^{\rm eq}$ from (5.1),
the rotational projector has components

$$
P^{\rm rot}_{j\alpha,k\beta}
=t^{(\gamma)}_{j\alpha}\,\mu^{\rm eq}_{\gamma\delta}\,t^{(\delta)}_{k\beta}
=\epsilon_{\gamma\rho\alpha}e_{j\rho}\;\mu^{\rm eq}_{\gamma\delta}\;\epsilon_{\delta\tau\beta}e_{k\tau},
$$

and the translational one, being orthonormal, is simply
$P^{\rm tr}_{j\alpha,k\beta}=\sqrt{M_jM_k}/M_N\,\delta_{\alpha\beta}$.

### 5.4 The relation

The three families together span all of $\mathbb R^{3N_{\rm nucl}}$ and are
mutually orthogonal, so their projectors sum to the identity
$\delta_{jk}\delta_{\alpha\beta}$. Solving for the vibrational block:

$$
\boxed{\;
\sum_b l_{jb,\alpha}\,l_{kb,\beta}
=\delta_{jk}\delta_{\alpha\beta}
-\underbrace{\frac{\sqrt{M_jM_k}}{M_N}\delta_{\alpha\beta}}_{\text{translations}}
-\underbrace{\epsilon_{\gamma\rho\alpha}e_{j\rho}\;\mu^{\rm eq}_{\gamma\delta}\;\epsilon_{\delta\tau\beta}e_{k\tau}}_{\text{rotations}}\;}
\tag{5.2}
$$

— Eq. (B.24). Two working facts follow immediately and are used *every single
time* below:

- **Fact (i) — translations never contribute.** Every coefficient we ever
  contract (5.2) against is built from $\mathbf e_j$, $\mathbf u_j$, or
  $\mathbf l_{jc}$, and all three obey $\sum_j\sqrt{M_j}(\cdot)_j=0$ — by
  (2.2), (2.4), and (E1) respectively. The translational term dies on contact.
- **Fact (ii) — the rotational term always produces $\tfrac12\mathsf I'_{N,c}$
  or $A$.** Each side of the rotational term carries one $\epsilon\,e_j$
  factor, which meets the $\epsilon$ of whatever object is being contracted.
  The resulting $\epsilon\epsilon$ pair always evaluates to one of exactly two
  constant matrices. That is §6.

---

## 6. The two workhorse $\epsilon$-contractions

Both use the single identity

$$
\epsilon_{\alpha\gamma\delta}\,\epsilon_{q\rho\delta}
=\delta_{\alpha q}\delta_{\gamma\rho}-\delta_{\alpha\rho}\delta_{\gamma q}.
\tag{6.1}
$$

**W1 (against an $\mathbf l$).** Contract (6.1) with $\sum_jl_{jc,\gamma}e_{j\rho}
=(\mathsf C_c^{T})_{\gamma\rho}=(\mathsf C_c)_{\rho\gamma}$:

$$
\epsilon_{\alpha\gamma\delta}\,\epsilon_{q\rho\delta}\sum_jl_{jc,\gamma}\,e_{j\rho}
=\delta_{\alpha q}\underbrace{\sum_jl_{jc,\rho}e_{j\rho}}_{=\,s_c}
-\underbrace{\sum_jl_{jc,q}e_{j\alpha}}_{=\,(\mathsf C_c)_{\alpha q}}
=s_c\delta_{\alpha q}-(\mathsf C_c)_{\alpha q},
$$

and inserting the lemma (4.6), $\mathsf C_c=s_c\mathbb 1-\tfrac12\mathsf I'_{N,c}$:

$$
\boxed{\;\epsilon_{\alpha\gamma\delta}\,\epsilon_{q\rho\delta}\sum_jl_{jc,\gamma}\,e_{j\rho}
=\tfrac12\,(\mathsf I'_{N,c})_{\alpha q}\;}
\tag{6.2}
$$

**W2 (against a $\mathbf u$).** Identically, with $\sum_ju_{j\gamma}e_{j\rho}=Y_{\gamma\rho}$:

$$
\epsilon_{\alpha\gamma\delta}\,\epsilon_{q\rho\delta}\sum_ju_{j\gamma}\,e_{j\rho}
=\delta_{\alpha q}\operatorname{tr}Y-Y_{q\alpha}
$$

and inserting (4.11), $Y=\tfrac12(\operatorname{tr}A)\mathbb 1-A$ with
$\operatorname{tr}Y=\tfrac12\operatorname{tr}A$:

$$
\boxed{\;\epsilon_{\alpha\gamma\delta}\,\epsilon_{q\rho\delta}\sum_ju_{j\gamma}\,e_{j\rho}
=\tfrac12(\operatorname{tr}A)\delta_{\alpha q}-\Bigl[\tfrac12(\operatorname{tr}A)\delta_{\alpha q}-A_{q\alpha}\Bigr]
=A_{\alpha q}\;}
\tag{6.3}
$$

That is the whole toolbox. **Everything below is (5.2) + (6.2) + (6.3) +
bookkeeping.**

---

## 7. The five sum rules, proven in full

The recipe is fixed: *insert completeness (5.2); drop translations by Fact
(i); evaluate the rotational term with W1 or W2.*

### 7.1 The $\zeta\zeta$ rule

Write both $\zeta$'s with (2.9) and contract over the mode index $d$:

$$
\sum_d\zeta^{\alpha}_{bd}\zeta^{\beta}_{cd}
=\epsilon_{\alpha\gamma\delta}\,\epsilon_{\beta\lambda\nu}
\sum_{jk}l_{jb,\gamma}\,l_{kc,\lambda}
\Bigl[\underbrace{\textstyle\sum_dl_{jd,\delta}\,l_{kd,\nu}}_{\text{insert (5.2)}}\Bigr].
$$

- **$\delta_{jk}\delta_{\delta\nu}$ part.** Contract the two epsilons with
  (6.1): $\epsilon_{\alpha\gamma\delta}\epsilon_{\beta\lambda\delta}
  =\delta_{\alpha\beta}\delta_{\gamma\lambda}-\delta_{\alpha\lambda}\delta_{\gamma\beta}$.
  Then use orthonormality (E3), $\sum_j\mathbf l_{jb}\!\cdot\!\mathbf l_{jc}=\delta_{bc}$:
  $$\delta_{\alpha\beta}\delta_{bc}-\sum_jl_{jb,\beta}\,l_{jc,\alpha}.$$
- **Translation part:** zero, by Fact (i).
- **Rotation part.** It factorises into a $j$-side and a $k$-side, each of
  which is exactly W1:
  $$
  -\Bigl[\epsilon_{\alpha\gamma\delta}\epsilon_{\gamma'\rho\delta}\sum_jl_{jb,\gamma}e_{j\rho}\Bigr]
  \mu^{\rm eq}_{\gamma'\delta'}
  \Bigl[\epsilon_{\beta\lambda\nu}\epsilon_{\delta'\tau\nu}\sum_kl_{kc,\lambda}e_{k\tau}\Bigr]
  =-\tfrac14\bigl(\mathsf I'_{N,b}\,\mathrm{\mu}^{\rm eq}\,\mathsf I'_{N,c}\bigr)_{\alpha\beta},
  $$
  where the last step used $\mathsf I'_{N,c}=(\mathsf I'_{N,c})^{T}$.

$$
\boxed{\;\sum_d\zeta^{\alpha}_{bd}\,\zeta^{\beta}_{cd}
=\delta_{\alpha\beta}\delta_{bc}-\sum_jl_{jb,\beta}\,l_{jc,\alpha}
-\tfrac14\bigl(\mathsf I'_{N,b}\,\mathrm{\mu}^{\rm eq}\,\mathsf I'_{N,c}\bigr)_{\alpha\beta}\;}
\tag{7.1}
$$

— Eq. (B.26).

### 7.2 The $\mathsf I'\mathsf I'$ rule

Use $\mathsf I'_{N,b}=2s_b\mathbb 1-2\mathsf C_b$ from (4.6):

$$
\sum_b(\mathsf I'_{N,b})_{yz}(\mathsf I'_{N,b})_{wx}
=4\Bigl(\sum_bs_b^2\Bigr)\delta_{yz}\delta_{wx}
-4\delta_{yz}\sum_bs_b(\mathsf C_b)_{wx}
-4\delta_{wx}\sum_bs_b(\mathsf C_b)_{yz}
+4\sum_b(\mathsf C_b)_{yz}(\mathsf C_b)_{wx}.
$$

Three contractions are needed; all three are the same recipe.

**(a) $\sum_bs_b\mathsf C_b=S$.** Write $s_b=\sum_k\mathbf e_k\!\cdot\!\mathbf l_{kb}$ and insert (5.2):

$$
\sum_bs_b(\mathsf C_b)_{yz}=\sum_{jk}e_{jy}\,e_{kq}\underbrace{\sum_bl_{kb,q}\,l_{jb,z}}_{(5.2)} .
$$

The $\delta_{jk}\delta_{qz}$ part gives $\sum_je_{jy}e_{jz}=S_{yz}$;
translations die; the rotation part contains
$\sum_ke_{kq}e_{k\rho}\epsilon_{\gamma\rho q}=\epsilon_{\gamma\rho q}S_{q\rho}=0$
because $S$ is symmetric. Hence

$$
\sum_bs_b\,\mathsf C_b=S .
\tag{7.2}
$$

**(b) $\sum_bs_b^2=\operatorname{tr}S$.** Take the trace of (7.2) and use
$\operatorname{tr}\mathsf C_b=s_b$.

**(c) $\sum_b\mathsf C_b\otimes\mathsf C_b$.** Same insertion:

$$
\sum_b(\mathsf C_b)_{yz}(\mathsf C_b)_{wx}=\sum_{jk}e_{jy}e_{kw}\underbrace{\sum_bl_{jb,z}l_{kb,x}}_{(5.2)}
=\delta_{zx}S_{yw}-\bigl[\epsilon_{\gamma\rho z}S_{\rho y}\bigr]\mu^{\rm eq}_{\gamma\delta}\bigl[\epsilon_{\delta\tau x}S_{\tau w}\bigr].
$$

Introduce the shorthand (any symmetric $M$; note $V[\mathbb 1]=0$)

$$
V[M]^{\gamma}_{yz}\equiv\epsilon_{\gamma\rho y}M_{\rho z}+\epsilon_{\gamma\rho z}M_{\rho y}.
\tag{7.3}
$$

The left-hand side of (c) is symmetric under $y\!\leftrightarrow\!z$ and under
$w\!\leftrightarrow\!x$ (because $\mathsf C_b$ is symmetric), so we may
symmetrise the right-hand side — which turns $\delta_{zx}S_{yw}$ into the
four-term average and $4\,\epsilon S\,\mu^{\rm eq}\,\epsilon S$ into
$V[S]\mu^{\rm eq}V[S]$. Assembling:

$$
\boxed{\;
\begin{aligned}
\sum_b(\mathsf I'_{N,b})_{yz}(\mathsf I'_{N,b})_{wx}
=\;&4(\operatorname{tr}S)\delta_{yz}\delta_{wx}-4S_{wx}\delta_{yz}-4S_{yz}\delta_{wx}\\
&+\delta_{yw}S_{zx}+\delta_{yx}S_{zw}+\delta_{zw}S_{yx}+\delta_{zx}S_{yw}
-V[S]^{\gamma}_{yz}\,\mu^{\rm eq}_{\gamma\delta}\,V[S]^{\delta}_{wx}
\end{aligned}\;}
\tag{7.4}
$$

— Eq. (B.28).

### 7.3 The $\mathrm{\sigma}\mathsf I'$ rule

Write $\sigma_{\alpha b}=\epsilon_{\alpha\lambda\nu}\sum_ju_{j\lambda}l_{jb,\nu}$
and $(\mathsf I'_{N,b})_{\gamma\delta}=\sum_k[2(\mathbf e_k\!\cdot\!\mathbf l_{kb})\delta_{\gamma\delta}
-l_{kb,\gamma}e_{k\delta}-e_{k\gamma}l_{kb,\delta}]$. Three terms.

**Term $\propto\delta_{\gamma\delta}$.** Its $\delta_{jk}$ part is
$2\delta_{\gamma\delta}\epsilon_{\alpha\lambda\nu}\sum_ju_{j\lambda}e_{j\nu}
=2\delta_{\gamma\delta}\epsilon_{\alpha\lambda\nu}Y_{\lambda\nu}=0$ since $Y$
is symmetric (4.11); its rotation part contains
$\sum_ke_{k\lambda}e_{k\tau}\epsilon_{\delta'\tau\lambda}=0$, again by symmetry
of $S$. **The whole term vanishes.**

**Terms $-l_{kb,\gamma}e_{k\delta}$ and $-e_{k\gamma}l_{kb,\delta}$.** Take the
first; the second is its $\gamma\!\leftrightarrow\!\delta$ mirror. Insert (5.2):

- $\delta_{jk}$ part: $-\epsilon_{\alpha\lambda\gamma}\sum_ju_{j\lambda}e_{j\delta}
  =-\epsilon_{\alpha\lambda\gamma}Y_{\lambda\delta}$; with (4.11) the
  $\tfrac12(\operatorname{tr}A)\mathbb 1$ piece drops ($\epsilon$ of $\delta$) and we
  are left with $+\epsilon_{\alpha\lambda\gamma}A_{\lambda\delta}$.
- rotation part: the $j$-side is exactly **W2** (6.3), giving $A_{\alpha\gamma'}$;
  the $k$-side gives $\epsilon_{\delta'\tau\gamma}S_{\tau\delta}$. Together:
  $(A\mathrm{\mu}^{\rm eq})_{\alpha\delta'}\epsilon_{\delta'\tau\gamma}S_{\tau\delta}$.

Adding the mirror term ($\gamma\!\leftrightarrow\!\delta$) recombines both
halves into the $V$-brackets of (7.3):

$$
\boxed{\;\sum_b\sigma_{\alpha b}\,\mathsf I'_{N,b}
=V[A]^{\alpha}+(A\mathrm{\mu}^{\rm eq})_{\alpha p}\,V[S]^{p}\;}
\tag{7.5}
$$

— Eq. (B.29).

### 7.4 The $\mathrm{\sigma}\zeta$ rule

$$
\sum_b\sigma_{\alpha b}\zeta^{\beta}_{bc}
=\epsilon_{\alpha\lambda\nu}\epsilon_{\beta m n}\sum_{jk}u_{j\lambda}\,l_{kc,n}
\Bigl[\underbrace{\textstyle\sum_bl_{jb,\nu}l_{kb,m}}_{(5.2)}\Bigr].
$$

- **$\delta_{jk}\delta_{\nu m}$ part.** $\epsilon_{\alpha\lambda\nu}\epsilon_{\beta\nu n}
  =-\epsilon_{\alpha\lambda\nu}\epsilon_{\beta n\nu}
  =-(\delta_{\alpha\beta}\delta_{\lambda n}-\delta_{\alpha n}\delta_{\lambda\beta})$, so it gives
  $$-\delta_{\alpha\beta}\sum_j\mathbf u_j\!\cdot\!\mathbf l_{jc}+\sum_ju_{j\beta}l_{jc,\alpha}.$$
  The first sum is elementary from (2.1), (4.2) and (E3):
  $$\sum_j\mathbf u_j\!\cdot\!\mathbf l_{jc}=\sum_j\mathbf e_j\!\cdot\!\mathbf l_{jc}
  +\sum_bQ_b\sum_j\mathbf l_{jb}\!\cdot\!\mathbf l_{jc}=s_c+Q_c .$$
- **Translation part:** zero.
- **Rotation part.** The $j$-side is W2 $\to A_{\alpha\gamma'}$; the $k$-side is
  W1 after one sign flip, $\epsilon_{\beta m n}\epsilon_{\delta'\tau m}
  =-\epsilon_{\beta n m}\epsilon_{\delta'\tau m}$, giving
  $-\tfrac12(\mathsf I'_{N,c})_{\beta\delta'}$. With the overall minus of the
  rotation term the two signs cancel:
  $+\tfrac12(A\mathrm{\mu}^{\rm eq}\mathsf I'_{N,c})_{\alpha\beta}$.

$$
\boxed{\;\sum_b\sigma_{\alpha b}\,\zeta^{\beta}_{bc}
=\sum_ju_{j\beta}\,l_{jc,\alpha}-\delta_{\alpha\beta}\bigl(s_c+Q_c\bigr)
+\tfrac12\bigl(A\mathrm{\mu}^{\rm eq}\mathsf I'_{N,c}\bigr)_{\alpha\beta}\;}
\tag{7.6}
$$

— Eq. (B.30).

### 7.5 The auxiliaries

**(A1) $\sum_cs_c\,\mathsf I'_{N,c}=2\,\mathsf I^{\rm eq}$.** Immediate from
(4.6), (7.2) and (b): $\sum_cs_c\mathsf I'_{N,c}=2(\sum_cs_c^2)\mathbb 1-2\sum_cs_c\mathsf C_c
=2(\operatorname{tr}S)\mathbb 1-2S=2\mathsf I^{\rm eq}$.

**(A2) $\sum_c(s_c+Q_c)\Lambda_c=6$.** Anticipating $\Lambda_c=\operatorname{tr}(A^{-1}\mathsf I'_{N,c})$
from §9, combine (A1) with (4.9):

$$
\sum_c(s_c+Q_c)\Lambda_c
=\operatorname{tr}\Bigl(A^{-1}\bigl[\underbrace{2\mathsf I^{\rm eq}}_{\text{(A1)}}
+\underbrace{2(A-\mathsf I^{\rm eq})}_{(4.9)}\bigr]\Bigr)
=2\operatorname{tr}\bigl(A^{-1}A\bigr)=6 .
\tag{7.7}
$$

*That* is where the "mysterious clean number 6" of §10.5 comes from: it is
$2\times$ the dimension of space, and nothing else.

**(A3) the vector $\boldsymbol\ell_j$.** Define $\ell_{j\alpha}\equiv\sum_c\Lambda_c\,l_{jc,\alpha}$.
Using $\Lambda_c=\operatorname{tr}(A^{-1}\mathsf I'_{N,c})=2s_c\operatorname{tr}A^{-1}
-2\sum_k(A^{-1}\mathbf e_k)\!\cdot\!\mathbf l_{kc}$ [from (4.6)] and inserting
(5.2) in both resulting contractions:

- $\sum_cs_c\,l_{jc,\alpha}=e_{j\alpha}$ — the $\delta_{jk}$ part; the rotation
  part carries $\epsilon_{\gamma\rho y}S_{y\rho}=0$.
- $\sum_c\sum_k(A^{-1}\mathbf e_k)_y\,l_{kc,y}\,l_{jc,\alpha}
  =(A^{-1}\mathbf e_j)_\alpha-\epsilon_{\delta\tau\alpha}e_{j\tau}(\mathrm{\mu}^{\rm eq}\boldsymbol\kappa)_\delta$,
  where the rotation part produced $\sum_k(A^{-1}\mathbf e_k)_ye_{k\rho}=(A^{-1}S)_{y\rho}$
  and hence, by the definition below, $\epsilon_{\gamma\rho y}(A^{-1}S)_{y\rho}=\kappa_\gamma$.

With the **antisymmetry vector**

$$
\kappa_\gamma\equiv\epsilon_{\gamma\rho\tau}\,(SA^{-1})_{\rho\tau}
=\tfrac12\operatorname{tr}\bigl(A^{-1}V[S]^{\gamma}\bigr),
\tag{7.8}
$$

$$
\boxed{\;\boldsymbol\ell_j
=2\bigl[(\operatorname{tr}A^{-1})\,\mathbf e_j-A^{-1}\mathbf e_j\bigr]
+2\,\epsilon_{\delta\tau\alpha}e_{j\tau}(\mathrm{\mu}^{\rm eq}\boldsymbol\kappa)_\delta\;}
\tag{7.9}
$$

— Eq. (B.40), first half.

**Note the pattern:** the rotational corrections only ever produce
$\tfrac12\mathsf I'_{N,c}$, $A$, $S$, and $\mathrm{\mu}^{\rm eq}$ — the same four
constant objects, over and over. That is why the assembly below has any chance
of closing.

---

## 8. Miracle 1: Watson's closed form $\mathsf I'=A\,\mathrm{\mu}^{\rm eq}A$

Expand $\mathsf I_N$ and $\mathrm{\sigma}\mathrm{\sigma}^{T}$ exactly — both are
just quadratic polynomials in $Q$. For $\mathsf I_N$, substitute (2.1) into (2.5):

$$
\sum_j|\mathbf u_j|^2=\operatorname{tr}S+2\sum_bs_bQ_b+\sum_{bc}Q_bQ_c\,\delta_{bc},
\qquad
\sum_j\mathbf u_j\mathbf u_j^{T}=S+\sum_bQ_b(\mathsf C_b+\mathsf C_b^{T})
+\sum_{bc}Q_bQ_c\sum_j\mathbf l_{jb}\mathbf l_{jc}^{T},
$$

so that, using $2s_b\mathbb 1-\mathsf C_b-\mathsf C_b^{T}=\mathsf I'_{N,b}$ (4.3),

$$
\mathsf I_{N}
=\mathsf I^{\rm eq}+\sum_b\mathsf I'_{N,b}\,Q_b
+\sum_{bc}Q_bQ_c\Bigl[\delta_{bc}\mathbb 1-\textstyle\sum_j\mathbf l_{jb}\mathbf l_{jc}^{T}\Bigr],
\qquad
(\mathrm{\sigma}\mathrm{\sigma}^{T})_{\alpha\beta}
=\sum_{bc}Q_bQ_c\sum_d\zeta^{\alpha}_{bd}\zeta^{\beta}_{cd},
\tag{8.1}
$$

the second from (2.8). Subtract, and insert the $\zeta\zeta$ rule (7.1) into
the quadratic term:

$$
\Bigl[\delta_{bc}\delta_{\alpha\beta}-\sum_jl_{jb,\alpha}l_{jc,\beta}\Bigr]
-\Bigl[\delta_{\alpha\beta}\delta_{bc}-\sum_jl_{jb,\beta}l_{jc,\alpha}
-\tfrac14(\mathsf I'_{N,b}\mathrm{\mu}^{\rm eq}\mathsf I'_{N,c})_{\alpha\beta}\Bigr].
$$

The $\delta_{bc}\delta_{\alpha\beta}$ terms cancel outright. The two
$\sum_j l\,l$ terms cancel *under the $b\!\leftrightarrow\!c$ symmetrisation
that the factor $Q_bQ_c$ enforces* — relabelling $b\!\leftrightarrow\!c$ in the
second turns $\sum_jl_{jc,\beta}l_{jb,\alpha}$ into the first. The **only**
survivor at quadratic order is $\tfrac14\mathsf I'_{N,b}\mathrm{\mu}^{\rm eq}\mathsf I'_{N,c}$.
So

$$
\mathsf I'=\mathsf I^{\rm eq}+\sum_b\mathsf I'_{N,b}Q_b
+\tfrac14\sum_{bc}Q_bQ_c\,\mathsf I'_{N,b}\mathrm{\mu}^{\rm eq}\mathsf I'_{N,c},
$$

which is *precisely* the expansion of $A\mathrm{\mu}^{\rm eq}A$ with
$A=\mathsf I^{\rm eq}+\tfrac12\sum_b\mathsf I'_{N,b}Q_b$ — multiply out and use
$\mathsf I^{\rm eq}\mathrm{\mu}^{\rm eq}=\mathbb 1$ on the constant and linear
orders. Hence, **exactly, to all orders**:

$$
\boxed{\;\mathsf I'(Q)=A(Q)\,\mathrm{\mu}^{\rm eq}\,A(Q),
\qquad A(Q)=\mathsf I^{\rm eq}+\tfrac12\sum_b\mathsf I'_{N,b}\,Q_b\;}
\qquad\text{— Eq. (4.37).}
\tag{8.2}
$$

Corollaries you get for free (check each by multiplying out; all used later):

$$
\mathrm{\mu}=A^{-1}\mathsf I^{\rm eq}A^{-1},\qquad
\mathrm{\mu} A=A^{-1}\mathsf I^{\rm eq},\qquad
A\mathrm{\mu} A=\mathsf I^{\rm eq},\qquad
\mathrm{\mu} A\mathrm{\mu}^{\rm eq}=A^{-1},\qquad
\gamma=\frac{(\det A)^{2}}{\det\mathsf I^{\rm eq}} .
\tag{8.3}
$$

Two further identities we will need in §10 and §12, both one-liners from (8.3)
and $\mathsf I^{\rm eq}=(\operatorname{tr}S)\mathbb 1-S$:

$$
\operatorname{tr}(\mathrm{\mu} A)=\operatorname{tr}(A^{-1}\mathsf I^{\rm eq})
=\operatorname{tr}S\operatorname{tr}A^{-1}-\operatorname{tr}(A^{-1}S),
\qquad
\operatorname{tr}\mathrm{\mu}=\operatorname{tr}(A^{-2}\mathsf I^{\rm eq})
=\operatorname{tr}S\operatorname{tr}A^{-2}-\operatorname{tr}(A^{-2}S).
\tag{8.4}
$$

---

## 9. Jacobi's formula, and why all derivatives are now constants

**Jacobi's formula** (freshman proof): for a symmetric positive-definite
matrix $M(Q)$, diagonalize $M=O\,{\rm diag}(m_i)\,O^{T}$; then
$\ln\det M=\sum_i\ln m_i=\operatorname{tr}\ln M$, and differentiating,

$$
\partial_b\ln\det M=\operatorname{tr}\bigl(M^{-1}\,\partial_bM\bigr).
$$

Apply it to $\gamma=(\det A)^2/\det\mathsf I^{\rm eq}$, with
$\partial_bA=\tfrac12\mathsf I'_{N,b}$ *constant*:

$$
\boxed{\;\Lambda_b=\partial_b\ln\gamma=2\operatorname{tr}\bigl(A^{-1}\partial_bA\bigr)
=\operatorname{tr}\bigl(A^{-1}\mathsf I'_{N,b}\bigr)=\operatorname{tr}\mathsf B_b,
\qquad
\partial_c\Lambda_b=-\tfrac12\operatorname{tr}\bigl(\mathsf B_b\mathsf B_c\bigr)\;}
\tag{9.1}
$$

using $\partial_cA^{-1}=-A^{-1}(\partial_cA)A^{-1}$ (differentiate
$A^{-1}A=\mathbb 1$), and likewise, from $\mathrm{\mu}=A^{-1}\mathsf I^{\rm eq}A^{-1}$,

$$
\partial_b\mathrm{\mu}=-\tfrac12\bigl(\mathsf B_b\,\mathrm{\mu}+\mathrm{\mu}\,\mathsf B_b^{T}\bigr),
\qquad \mathsf B_b\equiv A^{-1}\mathsf I'_{N,b},\quad \mathsf B_b^{T}=\mathsf I'_{N,b}A^{-1}.
\tag{9.2}
$$

**This is Miracle 1 in action:** $\Lambda_b$, its derivatives, and
$\partial_b\mathrm{\mu}$ — every derivative in the five groups is now an
algebraic expression in the constant matrices $\mathsf I'_{N,b}$ and the linear
matrix $A(Q)$. **No derivative remains to be taken.**

---

## 10. The five groups, evaluated

Shorthand for the rest of the note: $\bar A\equiv A^{-1}$,
$t_1\equiv\operatorname{tr}\bar A$, $t_2\equiv\operatorname{tr}\bar A^{2}$.

### 10.1 $N_0$ — two contractions of the $\mathsf I'\mathsf I'$ rule

By (9.1), $N_0=\sum_b\bigl[-\tfrac12\operatorname{tr}(\mathsf B_b\mathsf B_b)+\tfrac14(\operatorname{tr}\mathsf B_b)^2\bigr]$.
Both sums are (7.4) contracted with two factors of $\bar A$:

$$
\sum_b\operatorname{tr}(\mathsf B_b\mathsf B_b)=\bar A_{yz}\bar A_{wx}\sum_b(\mathsf I'_{N,b})_{zw}(\mathsf I'_{N,b})_{xy},
\qquad
\sum_b(\operatorname{tr}\mathsf B_b)^2=\bar A_{yz}\bar A_{wx}\sum_b(\mathsf I'_{N,b})_{zy}(\mathsf I'_{N,b})_{xw}.
$$

Carry them out term by term — each Kronecker delta of (7.4) closes a pair of
$\bar A$'s into a trace, and the last term of (7.4) becomes either a
$V$-trace or, when both $\bar A$'s close on their *own* $V$, the square
$\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa$ via (7.8):

$$
\sum_b\operatorname{tr}(\mathsf B_b\mathsf B_b)
=4\operatorname{tr}S\,t_2-6\operatorname{tr}(\bar A^{2}S)+2t_1\operatorname{tr}(\bar AS)
-\mu^{\rm eq}_{\gamma\delta}\operatorname{tr}\bigl(\bar AV[S]^{\gamma}\bar AV[S]^{\delta}\bigr),
$$

$$
\sum_b(\operatorname{tr}\mathsf B_b)^2
=4\operatorname{tr}S\,t_1^{2}-8t_1\operatorname{tr}(\bar AS)+4\operatorname{tr}(\bar A^{2}S)
-4\,\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa .
\tag{10.1}
$$

Combining with the coefficients $-\tfrac12$ and $+\tfrac14$:

$$
\boxed{\;
N_0=\operatorname{tr}S\bigl[t_1^{2}-2t_2\bigr]
+4\operatorname{tr}(\bar A^{2}S)-3\,t_1\operatorname{tr}(\bar AS)
+\tfrac12\mu^{\rm eq}_{\gamma\delta}\operatorname{tr}\bigl(\bar AV[S]^{\gamma}\bar AV[S]^{\delta}\bigr)
-\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa\;}
\tag{10.2}
$$

— Eq. (B.35).

### 10.2 $N_1=0$ — a divergence lemma

By (9.2), $N_1$ needs $\sum_b\sigma_{\alpha b}(\mathsf B_b\mathrm{\mu})_{\alpha\beta}$
and $\sum_b\sigma_{\alpha b}(\mathrm{\mu}\mathsf B_b^{T})_{\alpha\beta}$, both fixed
by the $\mathrm{\sigma}\mathsf I'$ rule (7.5). Two elementary facts do all the work:

- *An $\epsilon$ contracted with a symmetric matrix vanishes.*
- $S\mathrm{\mu}^{\rm eq}=(\operatorname{tr}S)\mathrm{\mu}^{\rm eq}-\mathbb 1$ is
  **symmetric** [from $S=(\operatorname{tr}S)\mathbb 1-\mathsf I^{\rm eq}$], and
  $\mathrm{\mu}A=\bar A\mathsf I^{\rm eq}=(\operatorname{tr}S)\bar A-\bar AS$ has
  antisymmetric part $\boldsymbol\kappa$ by (7.8).

**First contraction.** $\bar A_{\alpha\gamma}V[A]^{\alpha}_{\gamma\delta}
=\epsilon_{\alpha\rho\gamma}\bar A_{\alpha\gamma}A_{\rho\delta}
+\epsilon_{\alpha\rho\delta}(\bar AA)_{\alpha\rho}=0+0$, and
$\bar A_{\alpha\gamma}(A\mathrm{\mu}^{\rm eq})_{\alpha p}=\mu^{\rm eq}_{\gamma p}$,
so the second piece of (7.5) gives $\mu^{\rm eq}_{\gamma p}V[S]^{p}_{\gamma\delta}
=\epsilon_{p\rho\delta}(S\mathrm{\mu}^{\rm eq})_{\rho p}=0$. Hence
$\sum_b\sigma_{\alpha b}(\mathsf B_b\mathrm{\mu})_{\alpha\beta}=0$.

**Second contraction.** $\mu_{\alpha\gamma}V[A]^{\alpha}_{\gamma\delta}
=\epsilon_{\alpha\rho\delta}(\mathrm{\mu}A)_{\alpha\rho}=+\kappa_\delta$, while
$\mu_{\alpha\gamma}(A\mathrm{\mu}^{\rm eq})_{\alpha p}=(\mathrm{\mu}A\mathrm{\mu}^{\rm eq})_{\gamma p}=\bar A_{\gamma p}$
by (8.3), giving $\bar A_{\gamma p}V[S]^{p}_{\gamma\delta}
=\epsilon_{p\rho\delta}(S\bar A)_{\rho p}=-\kappa_\delta$. **The pair cancels.**

$$
\boxed{\;\sum_b\sigma_{\alpha b}(\partial_b\mathrm{\mu})_{\alpha\beta}=0
\iff \sum_b\partial_b(\mathrm{\mu}\mathrm{\sigma})_{\alpha b}=0
\quad\Longrightarrow\quad N_1=0\;}
\tag{10.3}
$$

The three vector fields $(\mathrm{\mu}\mathrm{\sigma})_{\alpha b}$ are
**divergence-free** on mode space — a pretty fact in its own right.

### 10.3 $N_4$

Contract (7.5) with $\bar A$ to get $v_\alpha=\sum_b\sigma_{\alpha b}\Lambda_b$.
The $V[A]$ part dies, $\operatorname{tr}(\bar AV[A]^{\alpha})=2\epsilon_{\alpha\rho\delta}(A\bar A)_{\rho\delta}=0$;
the $V[S]$ part gives $\operatorname{tr}(\bar AV[S]^{p})=2\kappa_p$ by (7.8). So

$$
v_\alpha=2(A\mathrm{\mu}^{\rm eq}\boldsymbol\kappa)_\alpha
\qquad\Longrightarrow\qquad
\boxed{\;N_4=\tfrac14\mathbf v^{T}\mathrm{\mu}\mathbf v
=\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\underbrace{(A\mathrm{\mu}A)}_{=\,\mathsf I^{\rm eq}}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa
=\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa\;}
\tag{10.4}
$$

which **cancels the last term of (10.2) on the spot.**

### 10.4 $N_3$

By (9.1), $N_3=-\tfrac12G_{bc}\operatorname{tr}(\mathsf B_b\mathsf B_c)
=-\tfrac12\mu_{\alpha\beta}\operatorname{tr}\bigl(\bar AW_\alpha\bar AW_\beta\bigr)$
with $W_\alpha\equiv\sum_b\sigma_{\alpha b}\mathsf I'_{N,b}=V[A]^{\alpha}+(A\mathrm{\mu}^{\rm eq})_{\alpha p}V[S]^{p}$
from (7.5). Expand the square; the two mixed terms are equal under
$\alpha\!\leftrightarrow\!\beta$ ($\mathrm{\mu}$ symmetric), and the coefficients
collapse by (8.3): $\mu_{\alpha\beta}(A\mathrm{\mu}^{\rm eq})_{\beta p}=\bar A_{\alpha p}$
and $(A\mathrm{\mu}^{\rm eq})^{T}\mathrm{\mu}(A\mathrm{\mu}^{\rm eq})
=\mathrm{\mu}^{\rm eq}(A\mathrm{\mu}A)\mathrm{\mu}^{\rm eq}=\mathrm{\mu}^{\rm eq}$. Hence

$$
\boxed{\;
N_3=-\tfrac12\mu_{\alpha\beta}\operatorname{tr}\bigl(\bar AV[A]^{\alpha}\bar AV[A]^{\beta}\bigr)
-\bar A_{\alpha p}\operatorname{tr}\bigl(\bar AV[A]^{\alpha}\bar AV[S]^{p}\bigr)
-\tfrac12\mu^{\rm eq}_{pq}\operatorname{tr}\bigl(\bar AV[S]^{p}\bar AV[S]^{q}\bigr)\;}
\tag{10.5}
$$

— Eq. (B.38). Its **last term cancels the $\mathrm{\mu}^{\rm eq}$-quadratic term
of (10.2)**. All $\mathrm{\mu}^{\rm eq}$-dependence is now gone from the running
total.

### 10.5 $N_2$

Contract the $\mathrm{\sigma}\zeta$ rule (7.6) with $\mu_{\alpha\beta}\Lambda_c$.
Three pieces.

**Piece 1** $=\mu_{\alpha\beta}\sum_ju_{j\beta}\ell_{j\alpha}$ with
$\boldsymbol\ell_j$ from (7.9). Using $Y=\tfrac12(\operatorname{tr}A)\mathbb 1-A$
(4.11) and $\operatorname{tr}(\bar A\mathrm{\mu}A)=\operatorname{tr}\mathrm{\mu}$:

- $2t_1\operatorname{tr}(\mathrm{\mu}Y)=t_1\operatorname{tr}A\operatorname{tr}\mathrm{\mu}-2t_1\operatorname{tr}(\mathrm{\mu}A)$;
- $-2\operatorname{tr}(\bar A\mathrm{\mu}Y)=-\operatorname{tr}A\operatorname{tr}(\bar A\mathrm{\mu})+2\operatorname{tr}\mathrm{\mu}$;
- the $\epsilon$-piece: $\epsilon_{\delta\tau\alpha}(\mathrm{\mu}Y)_{\alpha\tau}=+\kappa_\delta$
  (the $\operatorname{tr}A$ part is symmetric and drops; the rest is $-\epsilon(\mathrm{\mu}A)=+\kappa$
  as in §10.2), giving $+2\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa$.

**Piece 2** $=-\operatorname{tr}\mathrm{\mu}\sum_c(s_c+Q_c)\Lambda_c=-6\operatorname{tr}\mathrm{\mu}$
by the auxiliary (7.7).

**Piece 3** $=\tfrac12\operatorname{tr}(\mathrm{\mu}A\mathrm{\mu}^{\rm eq}P)
=\tfrac12\operatorname{tr}(\bar AP)$ with $P\equiv\sum_c\Lambda_c\mathsf I'_{N,c}$, using
$\mathrm{\mu}A\mathrm{\mu}^{\rm eq}=\bar A$. But
$\operatorname{tr}(\bar AP)=\sum_c\Lambda_c\operatorname{tr}(\bar A\mathsf I'_{N,c})=\sum_c\Lambda_c^{2}$,
already computed in (10.1). So Piece 3 $=2\operatorname{tr}S\,t_1^{2}-4t_1\operatorname{tr}(\bar AS)
+2\operatorname{tr}(\bar A^{2}S)-2\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa$.

The two $\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa$ terms cancel
between Pieces 1 and 3, and
$-2t_1\operatorname{tr}(\mathrm{\mu}A)=-2\operatorname{tr}S\,t_1^{2}+2t_1\operatorname{tr}(\bar AS)$
by (8.4) removes the remaining $\operatorname{tr}S\,t_1^2$:

$$
\boxed{\;
N_2=\operatorname{tr}A\bigl[t_1\operatorname{tr}\mathrm{\mu}-\operatorname{tr}(\bar A\mathrm{\mu})\bigr]
-2\,t_1\operatorname{tr}(\bar AS)-4\operatorname{tr}\mathrm{\mu}+2\operatorname{tr}(\bar A^{2}S)\;}
\tag{10.6}
$$

— Eq. (B.41).

---

## 11. The two remaining $V$-traces

Both are evaluated with the $\epsilon\epsilon$ **pair identity** (any symmetric $M$),
obtained by contracting the $3\times3$ determinant expansion of
$\epsilon_{iab}\epsilon_{jcd}$ with $M_{ij}$:

$$
M_{ij}\,\epsilon_{iab}\,\epsilon_{jcd}
=\operatorname{tr}M\,(\delta_{ac}\delta_{bd}-\delta_{ad}\delta_{bc})
-M_{ac}\delta_{bd}+M_{bc}\delta_{ad}+M_{ad}\delta_{bc}-M_{bd}\delta_{ac},
\tag{11.1}
$$

together with the one simplification that makes everything collapse:

$$
\bigl(\bar AV[A]^{\alpha}\bigr)_{\gamma\delta}
=\epsilon_{\alpha\rho\lambda}\bar A_{\gamma\lambda}A_{\rho\delta}
+\epsilon_{\alpha\rho\delta}\underbrace{(\bar AA)_{\gamma\rho}}_{=\,\delta_{\gamma\rho}}
=\epsilon_{\alpha\rho\lambda}\bar A_{\gamma\lambda}A_{\rho\delta}+\epsilon_{\alpha\gamma\delta}.
\tag{11.2}
$$

**(11a).** Squaring (11.2) and contracting with $\mu_{\alpha\beta}$ gives four
terms. The pure-$\epsilon$ term is $\mu_{\alpha\beta}\epsilon_{\alpha\gamma\delta}\epsilon_{\beta\delta\gamma}=-2\operatorname{tr}\mathrm{\mu}$;
the two cross terms are equal and each evaluates by (11.1) to
$\operatorname{tr}\mathrm{\mu}[\operatorname{tr}A\,t_1-1]-t_1\operatorname{tr}(\mathrm{\mu}A)-\operatorname{tr}A\operatorname{tr}(\bar A\mathrm{\mu})$;
the doubly-$K$ term simplifies because $\bar A A=\mathbb 1$ collapses its tensor
factor to $\delta\delta$, leaving $-2\operatorname{tr}\mathrm{\mu}$. Adding, and
using (8.4) for $\operatorname{tr}(\mathrm{\mu}A)$:

$$
\mu_{\alpha\beta}\operatorname{tr}\bigl(\bar AV[A]^{\alpha}\bar AV[A]^{\beta}\bigr)
=2\operatorname{tr}A\bigl[t_1\operatorname{tr}\mathrm{\mu}-\operatorname{tr}(\bar A\mathrm{\mu})\bigr]
-6\operatorname{tr}\mathrm{\mu}-2t_1\bigl[\operatorname{tr}S\,t_1-\operatorname{tr}(\bar AS)\bigr].
\tag{11.3}
$$

**(11b).** The same expansion against $\bar AV[S]^{p}$, contracted with
$\bar A_{\alpha p}$. Here two collapses do the work:
$\sum_qA_{\rho q}(\bar AS)_{q m}=S_{\rho m}$ and $\sum_qA_{\rho q}\bar A_{qn}=\delta_{\rho n}$.
Each of the four terms is one application of (11.1); they come out to
$-t_1\operatorname{tr}(\bar AS)+\operatorname{tr}(\bar A^{2}S)$ (twice) and
$\operatorname{tr}S[t_1^{2}-t_2]-2t_1\operatorname{tr}(\bar AS)+2\operatorname{tr}(\bar A^{2}S)$ (twice):

$$
\bar A_{\alpha p}\operatorname{tr}\bigl(\bar AV[A]^{\alpha}\bar AV[S]^{p}\bigr)
=2\operatorname{tr}S\bigl[t_1^{2}-t_2\bigr]
-6\,t_1\operatorname{tr}(\bar AS)+6\operatorname{tr}(\bar A^{2}S).
\tag{11.4}
$$

— Eqs. (B.43)–(B.44).

---

## 12. Miracle 2: the cancellation table

Add $N_0+N_1+N_2+N_3+N_4$. First, three **pairwise cancellations** of the
ugly objects:

| ugly object | appears in | cancels against |
|---|---|---|
| $\boldsymbol\kappa^{T}\mathrm{\mu}^{\rm eq}\boldsymbol\kappa$ | $N_0$ (with $-$) | $N_4$ (with $+$) |
| $\mu^{\rm eq}\text{-weighted }V[S]V[S]$ traces | $N_0$ (with $+\tfrac12$) | $N_3$ (with $-\tfrac12$) |
| $\operatorname{tr}A\,[t_1\operatorname{tr}\mathrm{\mu}-\operatorname{tr}(\bar A\mathrm{\mu})]$ | $N_2$ (with $+$) | $-\tfrac12\times$(11.3) inside $N_3$ (with $-$) |

Everything that could still remember the fine structure of the molecule —
the equilibrium reciprocal inertia $\mathrm{\mu}^{\rm eq}$, the antisymmetry
vector $\boldsymbol\kappa$, the trace of $A$ itself — is now gone. What is
left is a handful of ordinary traces; collect their coefficients:

| trace | from $N_0$ | from $N_2$ | from $-\tfrac12$(11.3) | from $-$(11.4) | total |
|---|---|---|---|---|---|
| $\operatorname{tr}S\;t_1^{2}$ | $+1$ | $0$ | $+1$ | $-2$ | $0$ |
| $\operatorname{tr}S\;t_2$ | $-2$ | $0$ | $0$ | $+2$ | $0$ |
| $\operatorname{tr}(\bar A^{2}S)$ | $+4$ | $+2$ | $0$ | $-6$ | $0$ |
| $t_1\operatorname{tr}(\bar AS)$ | $-3$ | $-2$ | $-1$ | $+6$ | $0$ |
| $\operatorname{tr}\mathrm{\mu}$ | $0$ | $-4$ | $+3$ | $0$ | $\mathbf{-1}$ |

Four columns of coefficients, four zeros, one survivor:

$$
\frac{8}{\hbar^{2}}\,U=-\operatorname{tr}\mathrm{\mu}
\qquad\Longrightarrow\qquad
\boxed{\;U=-\frac{\hbar^{2}}{8}\sum_{\alpha}\mu_{\alpha\alpha}(Q).\;}
$$

That is Watson's pseudopotential. Notice what *kind* of result this is: not
an approximation, not a leading order — an exact algebraic identity, true at
every shape $Q$, for every molecule, provided the coordinates are rectilinear
Eckart normal coordinates.

**Physical footnotes.**

- $U$ is $O(\hbar^{2})$ and purely multiplicative — a genuine (tiny) potential
  energy on the shape space, with no classical counterpart.
- It came entirely from the *measure* ($\sqrt\gamma$ in the volume element).
  The ordering of the kernel operators contributes no $c$-number at all
  (Appendix B.1–B.5); clean division of labour.
- The collapse *needs* the Eckart conditions and rectilinear normal
  coordinates: for general curvilinear internal coordinates the closed form
  (8.2) fails and extra potential-like terms survive.

---

## 13. Check it yourself (pure Python, no packages)

The script builds a bent triatomic (H$_2$O-like masses), constructs
Eckart-compliant mode vectors, picks a random shape $Q$, and then verifies
**every boxed relation in this note** — the inertia derivative (4.3), the
lemma (4.6), completeness (5.2), the contractions (6.2)–(6.3), all five sum
rules of §7, the closed form (8.2), Jacobi (9.1), and each of the five groups
$N_0\dots N_4$ — before confirming that their sum is $-\operatorname{tr}\mathrm{\mu}$.
The last block evaluates the residue (B.19) itself by finite differences.
Expected output: every algebraic check $\lesssim10^{-14}$, and the two
finite-difference numbers agreeing to $\sim10^{-5}$.

```python
"""Every relation used in the collapse (B.19) -> U = -(hbar^2/8) tr(mu). No packages."""
import random, math
random.seed(1)
R3 = range(3)

# ---- tiny 3x3 linear algebra ----
def mm(A,B):  return [[sum(A[i][k]*B[k][j] for k in R3) for j in R3] for i in R3]
def T(A):     return [[A[j][i] for j in R3] for i in R3]
def tr(A):    return sum(A[i][i] for i in R3)
def det3(A):  return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                      + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
def inv3(A):
    d = det3(A)
    C = [[(A[(i+1)%3][(j+1)%3]*A[(i+2)%3][(j+2)%3]-A[(i+1)%3][(j+2)%3]*A[(i+2)%3][(j+1)%3]) for j in R3] for i in R3]
    return [[C[j][i]/d for j in R3] for i in R3]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a,b):   return sum(x*y for x,y in zip(a,b))
def eps(i,j,k): return (i-j)*(j-k)*(k-i)//2          # Levi-Civita

# ---- equilibrium geometry (masses O,H,H), shifted to the nuclear centre of mass ----
M   = [16.0, 1.0, 1.0]
Req = [[0.0,0.0,0.065],[0.0,0.757,-0.520],[0.0,-0.757,-0.520]]
MN  = sum(M)
com = [sum(M[j]*Req[j][a] for j in R3)/MN for a in R3]
e   = [[math.sqrt(M[j])*(Req[j][a]-com[a]) for a in R3] for j in R3]      # e_j

# ---- 3 mode vectors l_jb: orthonormal, orthogonal to translations & rotations ----
ext = []
for a in R3:                                                 # translations
    v=[0.0]*9
    for j in R3: v[3*j+a]=math.sqrt(M[j]/MN)
    ext.append(v)
for g in R3:                                                 # rotations  n_g x e_j
    n=[0.0]*3; n[g]=1.0
    v=[]
    for j in R3: v+=cross(n,e[j])
    ext.append(v)
def gram_schmidt(vs):
    out=[]
    for v in vs:
        w=list(v)
        for u in out:
            c=dot(w,u); w=[wi-c*ui for wi,ui in zip(w,u)]
        n=math.sqrt(dot(w,w))
        if n>1e-12: out.append([wi/n for wi in w])
    return out
full = gram_schmidt(ext + [[random.uniform(-1,1) for _ in range(9)] for _ in range(6)])
L    = [[[full[6+b][3*j+a] for a in R3] for j in R3] for b in R3]         # l[b][j][alpha]

# ---- the constants of Sec. 4:  S, I^eq, s_b, C_b, I'_{N,b}, zeta ----
S    = [[sum(e[j][a]*e[j][b] for j in R3) for b in R3] for a in R3]
Ieq  = [[tr(S)*(a==b)-S[a][b] for b in R3] for a in R3]
mueq = inv3(Ieq)
s    = [sum(dot(e[j],L[b][j]) for j in R3) for b in R3]
C    = [[[sum(e[j][x]*L[b][j][y] for j in R3) for y in R3] for x in R3] for b in R3]
dI   = [[[2*s[b]*(x==y)-C[b][x][y]-C[b][y][x] for y in R3] for x in R3] for b in R3]   # I'_{N,b}
zeta = [[[sum(cross(L[b][j],L[c][j])[al] for j in R3) for c in R3] for b in R3] for al in R3]

def rep(tag,val,ref=0.0): print(f"  {tag:46s} {abs(val-ref):.2e}")
def V(Mx,g,y,z): return sum(eps(g,r,y)*Mx[r][z]+eps(g,r,z)*Mx[r][y] for r in R3)

print("Sec. 4  (inertia derivative)")
rep("C_b symmetric  [rotational Eckart]", max(abs(C[b][x][y]-C[b][y][x]) for b in R3 for x in R3 for y in R3))
rep("C_b = s_b 1 - (1/2) I'_{N,b}",
    max(abs(C[b][x][y]-(s[b]*(x==y)-0.5*dI[b][x][y])) for b in R3 for x in R3 for y in R3))
rep("s_b = (1/4) tr I'_{N,b}", max(abs(s[b]-0.25*tr(dI[b])) for b in R3))

print("Sec. 5-6  (completeness and the two contractions)")
err=0.0
for j in R3:
  for k in R3:
    for al in R3:
      for be in R3:
        rot=sum(eps(g,r,al)*e[j][r]*mueq[g][d]*eps(d,t,be)*e[k][t] for g in R3 for r in R3 for d in R3 for t in R3)
        err=max(err,abs(sum(L[b][j][al]*L[b][k][be] for b in R3)
                        -((j==k)*(al==be)-math.sqrt(M[j]*M[k])/MN*(al==be)-rot)))
rep("completeness  sum_b l_jb,a l_kb,B", err)
rep("W1  eps eps sum_j l_jc e_j = (1/2) I'_{N,c}",
    max(abs(sum(eps(al,g,d)*eps(q,r,d)*L[c][j][g]*e[j][r] for g in R3 for d in R3 for r in R3 for j in R3)
            -0.5*dI[c][al][q]) for c in R3 for al in R3 for q in R3))

# ---- pick a shape Q and build everything there ----
Q0  = [random.uniform(-0.25,0.25) for _ in R3]
u   = [[e[j][a]+sum(L[b][j][a]*Q0[b] for b in R3) for a in R3] for j in R3]
A   = [[Ieq[x][y]+0.5*sum(dI[b][x][y]*Q0[b] for b in R3) for y in R3] for x in R3]
Ai  = inv3(A)
IN  = [[sum(dot(u[j],u[j]) for j in R3)*(a==b)-sum(u[j][a]*u[j][b] for j in R3) for b in R3] for a in R3]
sig = [[sum(cross(u[j],L[b][j])[al] for j in R3) for b in R3] for al in R3]
Ip  = [[IN[a][b]-sum(sig[a][c]*sig[b][c] for c in R3) for b in R3] for a in R3]
mu  = inv3(Ip)
Y   = [[sum(u[j][a]*e[j][b] for j in R3) for b in R3] for a in R3]
SA  = mm(S,Ai)
kap = [sum(eps(g,p,q)*SA[p][q] for p in R3 for q in R3) for g in R3]

rep("W2  eps eps sum_j u_j e_j = A",
    max(abs(sum(eps(al,g,d)*eps(q,r,d)*u[j][g]*e[j][r] for g in R3 for d in R3 for r in R3 for j in R3)
            -A[al][q]) for al in R3 for q in R3))
rep("Y = (1/2)(tr A) 1 - A", max(abs(Y[x][y]-(0.5*tr(A)*(x==y)-A[x][y])) for x in R3 for y in R3))

print("Sec. 7  (the five sum rules)")
rep("zeta-zeta", max(abs(sum(zeta[al][b][d]*zeta[be][c][d] for d in R3)
      -((al==be)*(b==c)-sum(L[b][j][be]*L[c][j][al] for j in R3)
        -0.25*mm(dI[b],mm(mueq,dI[c]))[al][be]))
      for al in R3 for be in R3 for b in R3 for c in R3))
rep("I'I'", max(abs(sum(dI[b][y][z]*dI[b][w][x] for b in R3)
      -(4*tr(S)*(y==z)*(w==x)-4*S[w][x]*(y==z)-4*S[y][z]*(w==x)
        +(y==w)*S[z][x]+(y==x)*S[z][w]+(z==w)*S[y][x]+(z==x)*S[y][w]
        -sum(V(S,g,y,z)*mueq[g][d]*V(S,d,w,x) for g in R3 for d in R3)))
      for y in R3 for z in R3 for w in R3 for x in R3))
AM=mm(A,mueq)
rep("sigma-I'", max(abs(sum(sig[al][b]*dI[b][g][d] for b in R3)
      -(V(A,al,g,d)+sum(AM[al][p]*V(S,p,g,d) for p in R3)))
      for al in R3 for g in R3 for d in R3))
rep("sigma-zeta", max(abs(sum(sig[al][b]*zeta[be][b][c] for b in R3)
      -(sum(u[j][be]*L[c][j][al] for j in R3)-(al==be)*(s[c]+Q0[c])
        +0.5*mm(AM,dI[c])[al][be]))
      for al in R3 for be in R3 for c in R3))
rep("sum_c s_c I'_{N,c} = 2 I^eq",
    max(abs(sum(s[c]*dI[c][x][y] for c in R3)-2*Ieq[x][y]) for x in R3 for y in R3))

print("Sec. 8-9  (closed form, Jacobi)")
rep("I' = A mu^eq A", max(abs(Ip[x][y]-mm(A,mm(mueq,A))[x][y]) for x in R3 for y in R3))
rep("gamma = (det A)^2 / det I^eq", (det3(Ip)-det3(A)**2/det3(Ieq))/det3(Ip))
B   = [mm(Ai,dI[b]) for b in R3]
Lam = [tr(B[b]) for b in R3]
rep("sum_c (s_c + Q_c) Lambda_c = 6", sum((s[c]+Q0[c])*Lam[c] for c in R3), 6.0)
muk = [sum(mueq[g][d]*kap[d] for d in R3) for g in R3]
rep("ell_j identity", max(abs(sum(L[c][j][al]*Lam[c] for c in R3)
      -(2*(tr(Ai)*e[j][al]-sum(Ai[al][g]*e[j][g] for g in R3))
        +2*sum(eps(d,t,al)*e[j][t]*muk[d] for d in R3 for t in R3)))
      for j in R3 for al in R3))

print("Sec. 10-12  (the five groups and the collapse)")
t1=tr(Ai); t2=tr(mm(Ai,Ai)); AS=tr(mm(Ai,S)); A2S=tr(mm(mm(Ai,Ai),S))
kmk=sum(kap[g]*mueq[g][d]*kap[d] for g in R3 for d in R3)
VV =sum(mueq[g][d]*sum(Ai[p][q]*V(S,g,q,r)*Ai[r][w]*V(S,d,w,p)
        for p in R3 for q in R3 for r in R3 for w in R3) for g in R3 for d in R3)
T43=sum(mu[a1][a2]*sum(Ai[p][q]*V(A,a1,q,r)*Ai[r][w]*V(A,a2,w,p)
        for p in R3 for q in R3 for r in R3 for w in R3) for a1 in R3 for a2 in R3)
T44=sum(Ai[a1][a2]*sum(Ai[p][q]*V(A,a1,q,r)*Ai[r][w]*V(S,a2,w,p)
        for p in R3 for q in R3 for r in R3 for w in R3) for a1 in R3 for a2 in R3)
G   = [[sum(sig[al][b]*mu[al][be]*sig[be][c] for al in R3 for be in R3) for c in R3] for b in R3]
v   = [sum(sig[al][c]*Lam[c] for c in R3) for al in R3]
dmu = [[[-0.5*(mm(B[b],mu)[x][y]+mm(mu,T(B[b]))[x][y]) for y in R3] for x in R3] for b in R3]
N0 = sum(-0.5*tr(mm(B[b],B[b]))+0.25*Lam[b]**2 for b in R3)
N1 = sum(sig[al][b]*dmu[b][al][be]*v[be] for al in R3 for be in R3 for b in R3)
N2 = sum(sig[al][b]*mu[al][be]*zeta[be][b][c]*Lam[c] for al in R3 for be in R3 for b in R3 for c in R3)
N3 = sum(G[b][c]*(-0.5*tr(mm(B[b],B[c]))) for b in R3 for c in R3)
N4 = 0.25*sum(v[a1]*mu[a1][a2]*v[a2] for a1 in R3 for a2 in R3)
rep("N0 = closed form", N0-(tr(S)*(t1*t1-2*t2)+4*A2S-3*t1*AS+0.5*VV-kmk))
rep("N1 = 0", N1)
rep("N2 = closed form", N2-(tr(A)*(t1*tr(mu)-tr(mm(Ai,mu)))-2*t1*AS-4*tr(mu)+2*A2S))
rep("N3 = closed form", N3-(-0.5*T43-T44-0.5*VV))
rep("N4 = kappa mu^eq kappa", N4-kmk)
rep("V-trace (11.3)", T43-(2*tr(A)*(t1*tr(mu)-tr(mm(Ai,mu)))-6*tr(mu)-2*t1*(tr(S)*t1-AS)))
rep("V-trace (11.4)", T44-(2*tr(S)*(t1*t1-t2)-6*t1*AS+6*A2S))
rep("N0+N1+N2+N3+N4 = -tr mu", N0+N1+N2+N3+N4+tr(mu))

# ---- and the residue (B.19) itself, by finite differences ----
def geom(Q):
    uu = [[e[j][a]+sum(L[b][j][a]*Q[b] for b in R3) for a in R3] for j in R3]
    I0 = [[sum(dot(uu[j],uu[j]) for j in R3)*(a==b)-sum(uu[j][a]*uu[j][b] for j in R3) for b in R3] for a in R3]
    sg = [[sum(cross(uu[j],L[b][j])[al] for j in R3) for b in R3] for al in R3]
    Iq = [[I0[a][b]-sum(sg[a][c]*sg[b][c] for c in R3) for b in R3] for a in R3]
    m  = inv3(Iq)
    g  = [[(b==c)+sum(sg[al][b]*m[al][be]*sg[be][c] for al in R3 for be in R3) for c in R3] for b in R3]
    return Iq, m, det3(Iq), g
h = 1e-5
def d(fun,i,Q=Q0):
    Qp=list(Q); Qp[i]+=h; Qm=list(Q); Qm[i]-=h
    fp,fm = fun(Qp),fun(Qm)
    if isinstance(fp,float): return (fp-fm)/(2*h)
    return [[(fp[a][b]-fm[a][b])/(2*h) for b in range(len(fp[0]))] for a in range(len(fp))]
lng, gQQ = (lambda Q: math.log(geom(Q)[2])), (lambda Q: geom(Q)[3])
_,_,_,g = geom(Q0)
Lfd  = [d(lng,b) for b in R3]
dLfd = [[d(lambda q,c=c: d(lng,c,q), b) for c in R3] for b in R3]
dgfd = [d(gQQ,b) for b in R3]
Ufd  = sum(dgfd[b][b][c]*Lfd[c] + g[b][c]*dLfd[b][c] + 0.25*g[b][c]*Lfd[b]*Lfd[c] for b in R3 for c in R3)
print("Sec. 3  (the residue itself, finite differences)")
print(f"  (8/hbar^2) U from (B.19) : {Ufd:.10f}")
print(f"  -tr(mu)                  : {-tr(mu):.10f}")
```

Actual output:

```
Sec. 4  (inertia derivative)
  C_b symmetric  [rotational Eckart]             2.21e-16
  C_b = s_b 1 - (1/2) I'_{N,b}                   1.11e-16
  s_b = (1/4) tr I'_{N,b}                        1.11e-16
Sec. 5-6  (completeness and the two contractions)
  completeness  sum_b l_jb,a l_kb,B              3.33e-16
  W1  eps eps sum_j l_jc e_j = (1/2) I'_{N,c}    1.67e-16
  W2  eps eps sum_j u_j e_j = A                  1.11e-16
  Y = (1/2)(tr A) 1 - A                          2.22e-16
Sec. 7  (the five sum rules)
  zeta-zeta                                      4.44e-16
  I'I'                                           8.88e-16
  sigma-I'                                       5.69e-16
  sigma-zeta                                     2.22e-16
  sum_c s_c I'_{N,c} = 2 I^eq                    4.44e-16
Sec. 8-9  (closed form, Jacobi)
  I' = A mu^eq A                                 1.11e-16
  gamma = (det A)^2 / det I^eq                   1.17e-16
  sum_c (s_c + Q_c) Lambda_c = 6                 1.78e-15
  ell_j identity                                 1.33e-15
Sec. 10-12  (the five groups and the collapse)
  N0 = closed form                               8.88e-16
  N1 = 0                                         4.99e-18
  N2 = closed form                               3.18e-15
  N3 = closed form                               7.70e-16
  N4 = kappa mu^eq kappa                         8.67e-18
  V-trace (11.3)                                 1.07e-14
  V-trace (11.4)                                 1.78e-15
  N0+N1+N2+N3+N4 = -tr mu                        0.00e+00
Sec. 3  (the residue itself, finite differences)
  (8/hbar^2) U from (B.19) : -4.2032179515
  -tr(mu)                  : -4.2032133067
```

Try changing the geometry, the masses, or the random seed for the mode
vectors: the agreement persists, because the identity is exact for *any*
Eckart-compliant set — the mode vectors don't even need to diagonalize a
force field.

---

## 14. FAQ

**Q: Why $\mathsf I'_{N,b}$ instead of the literature's $a^{b}$?**
Because in $a^{b}$ the letter $a$ is a *name* and $b$ is an *index*, so the
symbol looks like a tensor component when it is a list of matrices — the very
confusion the Cartesian subscripts then compound, $a^{b}_{\alpha\beta}$.
$\mathsf I'_{N,b}$ names the object: the $Q_b$-derivative of $\mathsf I_N$.
See §0 for the dictionary to the paper's notation and Watson's.

**Q: Isn't $\mathsf I'_{N,b}$ confusable with $\mathsf I'$?**
Only if you ignore the subscript, which is exactly what distinguishes them.
$\mathsf I'$ alone is the Coriolis-modified inertia — a $Q$-dependent matrix
whose prime is part of its name. $\mathsf I'_{N,b}$ carries the subscript
$N,b$, is constant, and is a derivative of $\mathsf I_N$ (not of $\mathsf I'$).
The two never occupy the same slot in any equation in this note. See the
warning box in §0.

**Q: Where exactly does each Eckart condition enter?**
The translational conditions ($\sum_j\sqrt{M_j}\,\mathbf l_{jb}=0$) silently
kill the translation term of completeness (5.2) in every sum rule — Fact (i).
The rotational condition ($\sum_j\mathbf e_j\times\mathbf l_{jb}=0$) *is* the
symmetry of $\mathsf C_b$ (§4.2) — the seed of the closed form. Orthonormality
of the $\mathbf l_b$'s gives the $\delta_{bc}$'s and $\sum_j\mathbf u_j\cdot\mathbf l_{jc}=s_c+Q_c$.
All three are used; none is optional.

**Q: Is the closed form $\mathsf I'=A\mathrm{\mu}^{\rm eq}A$ an expansion?**
No — $\mathsf I'$ is exactly quadratic in $Q$ (inertia of points moving
linearly, minus a Coriolis term that is quadratic by (2.8)), and
$A\mathrm{\mu}^{\rm eq}A$ is exactly quadratic too; §8 matches all three orders
exactly. It holds at arbitrarily large $Q$ (as long as $A$ stays invertible).

**Q: Why does $U$ come out negative (an attractive well)?**
$\operatorname{tr}\mathrm{\mu}=\sum_\alpha 1/I'_\alpha$-like is positive,
so $U<0$ everywhere; it is larger where the moments of inertia are small
(light, compact molecules). For H$_2$O it is of order $-10\,{\rm cm}^{-1}$ —
small but spectroscopically visible.

**Q: The five groups contain $\operatorname{tr}A$, $\boldsymbol\kappa$,
$\mathrm{\mu}^{\rm eq}$… why must they cancel?**
Because the answer they add up to, $-\operatorname{tr}\mathrm{\mu}$, doesn't
contain them — but that is hindsight. The honest statement: the cancellations
are consequences of the *same* completeness relation (5.2) that generated the
terms in the first place. The Eckart frame doesn't just simplify the classical
Hamiltonian; it fine-tunes the quantum residue too.

**Q: What if I use curvilinear internal coordinates (bond lengths, angles)?**
Then $\partial\mathbf u_j/\partial Q_b$ is no longer constant — Eq. (2.1)
fails at the first step — so $\mathsf I_N$ is not quadratic, $\mathsf I'$ is no
longer $A\mathrm{\mu}^{\rm eq}A$ with linear $A$, Miracle 1 fails, and the
residue does *not* collapse to a trace. Extra pseudopotential terms survive.
Watson's clean result is a special property of rectilinear normal coordinates.

**Q: One-sentence summary?**
The reweighting residue looks like an arbitrary function of the shape, but the
Eckart conditions force $\det\mathsf I'=(\det A)^2/\det\mathsf I^{\rm eq}$ with
$A$ linear in $Q$, so every derivative becomes a constant matrix, and the
completeness sum rules make all the complicated terms cancel in pairs —
leaving exactly $-\frac{\hbar^2}{8}\operatorname{tr}\mathrm{\mu}$.

---

*Version note:* equation numbers of the form (B.xx) and (4.xx) refer to the
current paper (Podolsky sandwich (4.32′), closed form (4.37), Appendix §B.6
with Eqs. (B.20)–(B.46)); numbers of the form (n.m) are internal to this note.
The paper still writes the inertia derivatives $a^{b}$ and the equilibrium
inertia $\mathsf I^{e}$ — see the dictionary in §0. The companion tutorial
`tutorial_coriolis_ordering.md` explains why the *Coriolis* block leaves no
residue at all; `tutorial_quantization.md` covers the older
ordering/anomalous-algebra questions.
