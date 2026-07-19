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
\gamma=\det\mathbf I',\quad
g^{QQ}=\mathbb 1+\boldsymbol\sigma^{T}\boldsymbol\mu\boldsymbol\sigma,
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
  $\mathbf I'=A\,(\mathbf I^{e})^{-1}A$ with $A$ *linear* in $Q$ — so every
  derivative of $\gamma$ and $\boldsymbol\mu$ reduces to constant matrices.
- **Miracle 2 (the assembly telescopes).** When the five groups of $U$ are
  evaluated with the Eckart *sum rules*, every complicated term cancels
  pairwise, and the sole survivor is $-\operatorname{tr}\boldsymbol\mu$.

The plan:

1. Warm-up: where the residue comes from (a 1D computation you can do on a napkin).
2. The residue in full and its five-group split.
3. Miracle 1: the Eckart-symmetry lemma, one sum rule proven in full, and
   Watson's closed form.
4. Jacobi's formula from scratch → all derivatives are constants.
5. The toolbox: the remaining sum rules (quoted, with proof idea).
6. The five groups evaluated.
7. Miracle 2: the cancellation table.
8. A runnable Python check (no packages needed).
9. FAQ.

**Assumed background:** product rule, integration by parts, matrix
multiplication, determinants and inverses of small matrices. No differential
geometry, no group theory.

**Notation** (matching the paper; repeated indices summed):

| Symbol | Meaning |
|---|---|
| $\mathbf e_j=\sqrt{M_j}\,\bar{\mathbf R}^{\rm eq}_j$ | mass-weighted equilibrium positions |
| $\mathbf u_j(Q)=\mathbf e_j+\sum_b\mathbf l_{jb}Q_b$ | mass-weighted positions at shape $Q$ |
| $\mathbf l_{jb}$ | Eckart mode vectors: orthonormal, $\perp$ translations and rotations |
| $\mathbf I_N(Q)$, $\boldsymbol\sigma(Q)$, $\mathbf I'=\mathbf I_N-\boldsymbol\sigma\boldsymbol\sigma^{T}$ | inertia, Coriolis vectors, modified inertia |
| $\boldsymbol\mu=(\mathbf I')^{-1}$, $\gamma=\det\mathbf I'$ | reciprocal inertia; shape Jacobian |
| $S=\sum_j\mathbf e_j\mathbf e_j^{T}$, $\mathbf I^{e}=(\operatorname{tr}S)\mathbb 1-S$, $\boldsymbol\mu^{e}=(\mathbf I^{e})^{-1}$ | equilibrium constants |
| $a^{b}=\partial\mathbf I_N/\partial Q_b\big|_{Q=0}$ | **constant** inertia derivatives |
| $A(Q)=\mathbf I^{e}+\tfrac12\sum_b a^{b}Q_b$ | Watson's linear matrix |
| $\Lambda_b=\partial_b\ln\gamma$ | log-derivative of the Jacobian |

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
dressing $\boldsymbol\sigma^T\boldsymbol\mu\boldsymbol\sigma$; substituting
$w=\gamma^{1/2}$ above gives exactly the paper's per-mode residue
$\frac{\hbar^2}{8}[\partial_b^2\ln\gamma+\tfrac14(\partial_b\ln\gamma)^2]$
of Eq. (4.34). Nothing about *this* step is special to molecules. The magic
is what the sum over modes does next.

---

## 2. The residue in full, split into five groups

Adding the identity part (4.34) and the dressed part (B.18) gives (B.19),
quoted at the top. Expand it with the product rule, using
$\zeta^{\alpha}_{bb}=0$ (antisymmetry of the Coriolis constants) to drop the
derivative of the explicit $\sigma_{\alpha b}$. With
$v_\alpha\equiv\sum_c\sigma_{\alpha c}\Lambda_c$, the result is Eq. (B.20):

$$
\frac{8}{\hbar^{2}}\,U
=\underbrace{\sum_b\Bigl[\partial_b\Lambda_b+\tfrac14\Lambda_b^2\Bigr]}_{N_{0}\ \text{(identity)}}
+\underbrace{\sigma_{\alpha b}(\partial_b\boldsymbol\mu)_{\alpha\beta}v_\beta}_{N_{1}}
+\underbrace{\sigma_{\alpha b}\mu_{\alpha\beta}\zeta^{\beta}_{bc}\Lambda_c}_{N_{2}}
+\underbrace{G_{bc}\,\partial_b\Lambda_c}_{N_{3}}
+\underbrace{\tfrac14\,\mathbf v^{T}\boldsymbol\mu\mathbf v}_{N_{4}} .
$$

Our job: evaluate all five and show the sum is $-\operatorname{tr}\boldsymbol\mu$.
Everything hinges on being able to *compute* $\Lambda_b=\partial_b\ln\det\mathbf I'$
— which brings us to the first miracle.

---

## 3. Miracle 1: Watson's closed form $\mathbf I'=A\,\boldsymbol\mu^{e}A$

### 3.1 The Eckart-symmetry lemma (two lines)

Define $E^{b}_{\alpha\beta}\equiv\sum_je_{j\alpha}l_{jb,\beta}$. The
*rotational Eckart condition* says $\sum_j\mathbf e_j\times\mathbf l_{jb}=0$;
component-wise that is $\epsilon_{\gamma\alpha\beta}E^{b}_{\alpha\beta}=0$,
i.e. **the antisymmetric part of $E^b$ vanishes: $E^b$ is symmetric.**
Comparing with the explicit inertia derivative
$a^{b}=2(\sum_j\mathbf e_j\!\cdot\!\mathbf l_{jb})\mathbb 1-E^{b}-E^{bT}$ and
taking a trace [Eq. (B.22)]:

$$
E^{b}=s_b\,\mathbb 1-\tfrac12\,a^{b},
\qquad
s_b\equiv\sum_j\mathbf e_j\!\cdot\!\mathbf l_{jb}=\tfrac14\operatorname{tr}a^{b}.
$$

Innocent-looking — but this is the seed of everything: *the overlap matrix
between equilibrium and displacement is entirely fixed by the inertia
derivative $a^b$.*

### 3.2 One sum rule, proven in full (the template)

The mode vectors $\mathbf l_b$ are $f=3N-6$ orthonormal vectors in the
$3N$-dimensional mass-weighted displacement space. The other six directions
are the three translations and the three rotations
$t^{\gamma}_{j\beta}=(\hat{\mathbf n}_\gamma\times\mathbf e_j)_\beta$, whose
Gram matrix is the equilibrium inertia:
$\sum_j\mathbf t^{\gamma}_j\!\cdot\!\mathbf t^{\delta}_j
=\hat{\mathbf n}_\gamma^{T}\bigl[\textstyle\sum_j(e_j^2\mathbb 1-\mathbf e_j\mathbf e_j^{T})\bigr]\hat{\mathbf n}_\delta
=I^{e}_{\gamma\delta}$ (one line, using the BAC–CAB identity). **Completeness**
of the whole set then reads [Eq. (B.24)]:

$$
\sum_b l_{jb,\alpha}\,l_{kb,\beta}
=\delta_{jk}\delta_{\alpha\beta}
-\underbrace{\frac{\sqrt{M_jM_k}}{M_N}\delta_{\alpha\beta}}_{\text{translations}}
-\underbrace{\epsilon_{\gamma\rho\alpha}e_{j\rho}\,\mu^{e}_{\gamma\delta}\,\epsilon_{\delta\sigma\beta}e_{k\sigma}}_{\text{rotations}} .
$$

Two working facts: *(i)* every coefficient we ever contract with is built from
$\mathbf e$, $\mathbf u$, or $\mathbf l$, all orthogonal to translations —
so the translation term **never contributes**; *(ii)* the rotation term always
meets an $\epsilon$-contraction that evaluates, via
$\epsilon_{\alpha\gamma\delta}\epsilon_{\gamma''\rho\delta}
=\delta_{\alpha\gamma''}\delta_{\gamma\rho}-\delta_{\alpha\rho}\delta_{\gamma\gamma''}$
and the lemma of §3.1, to one of two constant matrices [Eq. (B.25)]:
$\tfrac12a^{c}$ (when contracted against an $\mathbf l$) or $A$ (against a
$\mathbf u$).

Now the **$\zeta\zeta$ sum rule**. Write
$\zeta^{\alpha}_{bd}=\epsilon_{\alpha\gamma\delta}\sum_jl_{jb,\gamma}l_{jd,\delta}$
and contract two of them over the mode index $d$, inserting completeness:

$$
\sum_d\zeta^{\alpha}_{bd}\zeta^{\beta}_{cd}
=\epsilon_{\alpha\gamma\delta}\,\epsilon_{\beta\varepsilon\varphi}
\sum_{jk}l_{jb,\gamma}\,l_{kc,\varepsilon}
\Bigl[\underbrace{\textstyle\sum_dl_{jd,\delta}l_{kd,\varphi}}_{\text{completeness}}\Bigr].
$$

- the $\delta_{jk}\delta_{\delta\varphi}$ part: contract the two epsilons
  ($\epsilon_{\alpha\gamma\delta}\epsilon_{\beta\varepsilon\delta}
  =\delta_{\alpha\beta}\delta_{\gamma\varepsilon}-\delta_{\alpha\varepsilon}\delta_{\gamma\beta}$),
  then use orthonormality $\sum_j\mathbf l_{jb}\!\cdot\!\mathbf l_{jc}=\delta_{bc}$:
  you get $\delta_{\alpha\beta}\delta_{bc}-\sum_jl_{jb,\beta}l_{jc,\alpha}$;
- the translation part: dies (fact *(i)*);
- the rotation part: each side is an "$\epsilon\,\mathbf l\,\mathbf e$"
  contraction $=\tfrac12a$ (fact *(ii)*), leaving
  $-\tfrac14(a^{b}\boldsymbol\mu^{e}a^{c})_{\alpha\beta}$.

$$
\boxed{\;\sum_d\zeta^{\alpha}_{bd}\zeta^{\beta}_{cd}
=\delta_{\alpha\beta}\delta_{bc}-\sum_jl_{jb,\beta}l_{jc,\alpha}
-\tfrac14\bigl(a^{b}\boldsymbol\mu^{e}a^{c}\bigr)_{\alpha\beta}\;}
\qquad\text{— Eq. (B.26).}
$$

Every other sum rule in the toolbox (§5) is proven by *exactly this recipe*:
insert completeness, drop translations, evaluate the rotational correction
with fact *(ii)*.

### 3.3 The closed form

Expand $\mathbf I_N$ and $\boldsymbol\sigma\boldsymbol\sigma^{T}$ exactly —
both are just quadratic polynomials in $Q$ [Eq. (B.31)]:

$$
\mathbf I_{N}
=\mathbf I^{e}+\sum_ba^{b}Q_b
+\sum_{bc}Q_bQ_c\Bigl[\delta_{bc}\mathbb 1-\textstyle\sum_j\mathbf l_{jb}\mathbf l_{jc}^{T}\Bigr],
\qquad
\boldsymbol\sigma\boldsymbol\sigma^{T}
=\sum_{bc}Q_bQ_c\,\textstyle\sum_d\boldsymbol\zeta_{bd}\boldsymbol\zeta_{cd}^{T} .
$$

Subtract, and insert the boxed $\zeta\zeta$ rule: the $\delta_{bc}\mathbb 1$
terms cancel, the two $\sum_j\mathbf l\,\mathbf l^{T}$ terms cancel under the
$b\leftrightarrow c$ symmetrization that $Q_bQ_c$ enforces, and the *only*
survivor at quadratic order is $\tfrac14a^{b}\boldsymbol\mu^{e}a^{c}$ — which
is precisely the quadratic term of $A\boldsymbol\mu^{e}A$. Matching the
constant and linear orders is immediate. Hence, **exactly, to all orders**:

$$
\boxed{\;\mathbf I'(Q)=A(Q)\,\boldsymbol\mu^{e}\,A(Q),
\qquad A(Q)=\mathbf I^{e}+\tfrac12\sum_ba^{b}Q_b\;}
\qquad\text{— Eq. (4.37).}
$$

Corollaries you get for free (check each by multiplying out; all used later):

$$
\boldsymbol\mu=A^{-1}\mathbf I^{e}A^{-1},\qquad
A\boldsymbol\mu A=\mathbf I^{e},\qquad
\boldsymbol\mu A\boldsymbol\mu^{e}=A^{-1},\qquad
\gamma=\frac{(\det A)^{2}}{\det\mathbf I^{e}} .
$$

---

## 4. Jacobi's formula, and why all derivatives are now constants

**Jacobi's formula** (freshman proof): for a symmetric positive-definite
matrix $M(Q)$, diagonalize $M=O\,{\rm diag}(m_i)\,O^{T}$; then
$\ln\det M=\sum_i\ln m_i=\operatorname{tr}\ln M$, and differentiating,

$$
\partial_b\ln\det M=\operatorname{tr}\bigl(M^{-1}\,\partial_bM\bigr).
$$

Apply it to $\gamma=(\det A)^2/\det\mathbf I^{e}$, with
$\partial_bA=\tfrac12a^{b}$ *constant*:

$$
\Lambda_b=\partial_b\ln\gamma=2\operatorname{tr}\bigl(A^{-1}\partial_bA\bigr)
=\operatorname{tr}\bigl(A^{-1}a^{b}\bigr)\equiv\operatorname{tr}B^{b},
\qquad
\partial_c\Lambda_b
=-\tfrac12\operatorname{tr}\bigl(B^{b}B^{c}\bigr),
$$

using $\partial_c A^{-1}=-A^{-1}(\partial_cA)A^{-1}$ (differentiate
$A^{-1}A=\mathbb 1$). **This is Miracle 1 in action:** $\Lambda_b$, its
derivatives, $\partial_b\boldsymbol\mu=-\tfrac12(B^{b}\boldsymbol\mu
+\boldsymbol\mu B^{bT})$ — every derivative in the five groups is now an
algebraic expression in the constant matrices $a^{b}$ and the linear matrix
$A(Q)$. No derivative remains to be taken.

---

## 5. The toolbox: the remaining sum rules (quoted)

Proven by the §3.2 recipe — insert completeness, drop translations, evaluate
the rotational correction — in Appendix §B.6 (Eqs. B.28–B.30, B.39–B.40).
With the shorthand
$V[M]^{\alpha}_{\gamma\delta}\equiv\epsilon_{\alpha\rho\gamma}M_{\rho\delta}
+\epsilon_{\alpha\rho\delta}M_{\rho\gamma}$ and
$\kappa_\gamma\equiv\epsilon_{\gamma\sigma\tau}(SA^{-1})_{\sigma\tau}$
(the antisymmetric part of $SA^{-1}$):

| Rule | Statement | Where proven |
|---|---|---|
| $aa$ | $\sum_ba^{b}\otimes a^{b}=[\text{terms in }S,\mathbb 1]-V[S]\boldsymbol\mu^{e}V[S]$ | (B.28) |
| $\sigma a$ | $\sum_b\sigma_{\alpha b}a^{b}=V[A]^{\alpha}+(A\boldsymbol\mu^{e})_{\alpha\varepsilon}V[S]^{\varepsilon}$ | (B.29) |
| $\sigma\zeta$ | $\sum_b\sigma_{\alpha b}\zeta^{\beta}_{bc}=\sum_ju_{j\beta}l_{jc,\alpha}-\delta_{\alpha\beta}(s_c{+}Q_c)+\tfrac12(A\boldsymbol\mu^{e}a^{c})_{\alpha\beta}$ | (B.30) |
| auxiliaries | $\sum_cs_ca^{c}=2\mathbf I^{e}$, $\ \sum_c(s_c{+}Q_c)\Lambda_c=6$ | (B.39) |

Note the pattern: the rotational corrections only ever produce $\tfrac12a$,
$A$, $S$, and $\boldsymbol\mu^{e}$ — the same four constant objects, over and
over. That is why the assembly below has any chance of closing.

---

## 6. The five groups, evaluated

Using §4 and the toolbox (full algebra: Appendix §B.6, Steps 5–10). Write
$t_1=\operatorname{tr}A^{-1}$, and remember
$\operatorname{tr}\boldsymbol\mu
=\operatorname{tr}(A^{-2}\mathbf I^{e})
=\operatorname{tr}S\operatorname{tr}A^{-2}-\operatorname{tr}(A^{-2}S)$.

- **$N_0$** [(B.35)]: contract the $aa$ rule twice with $A^{-1}$:
  $$
  N_0=\operatorname{tr}S\bigl[t_1^{2}-2\operatorname{tr}A^{-2}\bigr]
  +4\operatorname{tr}(A^{-2}S)-3\,t_1\operatorname{tr}(A^{-1}S)
  +\tfrac12\mu^{e}_{\gamma\delta}\operatorname{tr}\bigl(A^{-1}V[S]^{\gamma}A^{-1}V[S]^{\delta}\bigr)
  -\boldsymbol\kappa^{T}\boldsymbol\mu^{e}\boldsymbol\kappa .
  $$
- **$N_1=0$** [(B.36)]: a pretty lemma — the vector fields
  $(\boldsymbol\mu\boldsymbol\sigma)_{\alpha b}$ are **divergence-free** on
  mode space, $\sum_b\partial_b(\boldsymbol\mu\boldsymbol\sigma)_{\alpha b}=0$.
  (Proof: apply the $\sigma a$ rule; every piece is an $\epsilon$-contraction
  of a *symmetric* matrix — e.g. $S\boldsymbol\mu^{e}$ is symmetric because
  $S$ and $\boldsymbol\mu^{e}$ are both functions of $S$ — plus a pair
  $\pm\kappa$ that cancels.)
- **$N_4$** [(B.37)]: the $\sigma a$ rule gives
  $v_\alpha=2(A\boldsymbol\mu^{e}\boldsymbol\kappa)_\alpha$, and with
  $A\boldsymbol\mu A=\mathbf I^{e}$,
  $$
  N_4=\tfrac14\mathbf v^{T}\boldsymbol\mu\mathbf v
  =\boldsymbol\kappa^{T}\boldsymbol\mu^{e}\boldsymbol\kappa .
  $$
- **$N_3$** [(B.38)]: with $W_\alpha\equiv\sum_b\sigma_{\alpha b}a^{b}$
  ($\sigma a$ rule) and $\boldsymbol\mu A\boldsymbol\mu^{e}=A^{-1}$,
  $$
  N_3=-\tfrac12\mu_{\alpha\beta}\operatorname{tr}\bigl(A^{-1}V[A]^{\alpha}A^{-1}V[A]^{\beta}\bigr)
  -A^{-1}_{\alpha\varepsilon}\operatorname{tr}\bigl(A^{-1}V[A]^{\alpha}A^{-1}V[S]^{\varepsilon}\bigr)
  -\tfrac12\mu^{e}_{\varepsilon\varphi}\operatorname{tr}\bigl(A^{-1}V[S]^{\varepsilon}A^{-1}V[S]^{\varphi}\bigr).
  $$
- **$N_2$** [(B.41)]: the $\sigma\zeta$ rule plus the auxiliaries (this is
  where the mysterious clean number $\sum_c(s_c+Q_c)\Lambda_c=6$ earns its keep):
  $$
  N_2=\operatorname{tr}A\bigl[t_1\operatorname{tr}\boldsymbol\mu
  -\operatorname{tr}(A^{-1}\boldsymbol\mu)\bigr]
  -2\,t_1\operatorname{tr}(A^{-1}S)-4\operatorname{tr}\boldsymbol\mu
  +2\operatorname{tr}(A^{-2}S).
  $$

The two remaining $V$-traces in $N_3$ evaluate, via the standard
$\epsilon\epsilon$-pair identity [(B.42)] and $A^{-1}A=\mathbb 1$, to
[(B.43)–(B.44)]:

$$
\mu_{\alpha\beta}\operatorname{tr}(A^{-1}V[A]^{\alpha}A^{-1}V[A]^{\beta})
=2\operatorname{tr}A\bigl[t_1\operatorname{tr}\boldsymbol\mu-\operatorname{tr}(A^{-1}\boldsymbol\mu)\bigr]
-6\operatorname{tr}\boldsymbol\mu-2t_1\bigl[\operatorname{tr}S\,t_1-\operatorname{tr}(A^{-1}S)\bigr],
$$
$$
A^{-1}_{\alpha\varepsilon}\operatorname{tr}(A^{-1}V[A]^{\alpha}A^{-1}V[S]^{\varepsilon})
=2\operatorname{tr}S\bigl[t_1^{2}-\operatorname{tr}A^{-2}\bigr]
-6\,t_1\operatorname{tr}(A^{-1}S)+6\operatorname{tr}(A^{-2}S).
$$

---

## 7. Miracle 2: the cancellation table

Add $N_0+N_1+N_2+N_3+N_4$. First, three **pairwise cancellations** of the
ugly objects:

| ugly object | appears in | cancels against |
|---|---|---|
| $\boldsymbol\kappa^{T}\boldsymbol\mu^{e}\boldsymbol\kappa$ | $N_0$ (with $-$) | $N_4$ (with $+$) |
| $\mu^{e}\text{-weighted }V[S]V[S]$ traces | $N_0$ (with $+\tfrac12$) | $N_3$ (with $-\tfrac12$) |
| $\operatorname{tr}A\,[t_1\operatorname{tr}\boldsymbol\mu-\operatorname{tr}(A^{-1}\boldsymbol\mu)]$ | $N_2$ (with $+$) | $-\tfrac12\times$(B.43) inside $N_3$ (with $-$) |

Everything that could still remember the fine structure of the molecule —
the equilibrium reciprocal inertia $\boldsymbol\mu^{e}$, the antisymmetry
vector $\boldsymbol\kappa$, the trace of $A$ itself — is now gone. What is
left is a handful of ordinary traces; collect their coefficients:

| trace | from $N_0$ | from $N_2$ | from $-\tfrac12$(B.43) | from $-$(B.44) | total |
|---|---|---|---|---|---|
| $\operatorname{tr}S\;t_1^{2}$ | $+1$ | $0$ | $+1$ | $-2$ | $0$ |
| $\operatorname{tr}S\operatorname{tr}A^{-2}$ | $-2$ | $0$ | $0$ | $+2$ | $0$ |
| $\operatorname{tr}(A^{-2}S)$ | $+4$ | $+2$ | $0$ | $-6$ | $0$ |
| $t_1\operatorname{tr}(A^{-1}S)$ | $-3$ | $-2$ | $-1$ | $+6$ | $0$ |
| $\operatorname{tr}\boldsymbol\mu$ | $0$ | $-4$ | $+3$ | $0$ | $\mathbf{-1}$ |

Four columns of coefficients, four zeros, one survivor:

$$
\frac{8}{\hbar^{2}}\,U=-\operatorname{tr}\boldsymbol\mu
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
  (4.37) fails and extra potential-like terms survive.

---

## 8. Check it yourself (pure Python, no packages)

The script below builds a bent triatomic (H$_2$O-like masses), constructs
Eckart-compliant mode vectors, picks a random shape $Q$, evaluates the residue
(B.19) by finite differences, and compares with $-\operatorname{tr}\boldsymbol\mu$.
It also verifies Watson's closed form to machine precision. Expected output:
the first two numbers agree to $\sim10^{-5}$ (finite-difference noise), the
closed-form error is $\sim10^{-16}$.

```python
"""U_vib + U_Cor = -(hbar^2/8) tr(mu):  numerical check  [(B.19) -> (4.38)]"""
import random, math
random.seed(1)

# ---- tiny 3x3 linear algebra ----
def mm(A,B):  return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A):    return sum(A[i][i] for i in range(len(A)))
def det3(A):  return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                      + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
def inv3(A):
    d = det3(A)
    C = [[(A[(i+1)%3][(j+1)%3]*A[(i+2)%3][(j+2)%3]-A[(i+1)%3][(j+2)%3]*A[(i+2)%3][(j+1)%3]) for j in range(3)] for i in range(3)]
    return [[C[j][i]/d for j in range(3)] for i in range(3)]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a,b):   return sum(x*y for x,y in zip(a,b))

# ---- equilibrium geometry (masses O,H,H), shifted to the nuclear centre of mass ----
M   = [16.0, 1.0, 1.0]
Req = [[0.0,0.0,0.065],[0.0,0.757,-0.520],[0.0,-0.757,-0.520]]
com = [sum(M[j]*Req[j][a] for j in range(3))/sum(M) for a in range(3)]
e   = [[math.sqrt(M[j])*(Req[j][a]-com[a]) for a in range(3)] for j in range(3)]   # e_j

# ---- 3 vibration vectors l_jb: orthonormal, orthogonal to translations & rotations ----
ext = []
for a in range(3):                                   # translations
    v=[0.0]*9
    for j in range(3): v[3*j+a]=math.sqrt(M[j]/sum(M))
    ext.append(v)
for g in range(3):                                   # rotations  n_g x e_j
    n=[0.0]*3; n[g]=1.0
    v=[]
    for j in range(3): v+=cross(n,e[j])
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
L    = [[[full[6+b][3*j+a] for a in range(3)] for j in range(3)] for b in range(3)]  # l[b][j][alpha]

# ---- geometry at a shape Q:  I', mu, gamma, g^QQ ----
def geom(Q):
    u   = [[e[j][a]+sum(L[b][j][a]*Q[b] for b in range(3)) for a in range(3)] for j in range(3)]
    IN  = [[sum(dot(u[j],u[j]) for j in range(3))*(a==b) - sum(u[j][a]*u[j][b] for j in range(3)) for b in range(3)] for a in range(3)]
    sig = [[sum(cross(u[j],L[b][j])[al] for j in range(3)) for b in range(3)] for al in range(3)]
    Ip  = [[IN[a][b]-sum(sig[a][c]*sig[b][c] for c in range(3)) for b in range(3)] for a in range(3)]
    mu  = inv3(Ip)
    g   = [[(b==c)+sum(sig[al][b]*mu[al][be]*sig[be][c] for al in range(3) for be in range(3)) for c in range(3)] for b in range(3)]
    return Ip, mu, det3(Ip), g

Q0, h = [random.uniform(-0.25,0.25) for _ in range(3)], 1e-5
def d(fun,i,Q=Q0):                                   # finite-difference d/dQ_i
    Qp=list(Q); Qp[i]+=h; Qm=list(Q); Qm[i]-=h
    fp,fm = fun(Qp),fun(Qm)
    if isinstance(fp,float): return (fp-fm)/(2*h)
    return [[(fp[a][b]-fm[a][b])/(2*h) for b in range(len(fp[0]))] for a in range(len(fp))]

lng  = lambda Q: math.log(geom(Q)[2])
gQQ  = lambda Q: geom(Q)[3]
Ip,mu,gam,g = geom(Q0)
Lam  = [d(lng,b) for b in range(3)]                          #  Lambda_b
dLam = [[d(lambda q,c=c: d(lng,c,q), b) for c in range(3)] for b in range(3)]
dg   = [d(gQQ,b) for b in range(3)]

# ---- the residue (B.19), times 8/hbar^2, versus  -tr(mu) ----
U = sum(dg[b][b][c]*Lam[c] + g[b][c]*dLam[b][c] + 0.25*g[b][c]*Lam[b]*Lam[c]
        for b in range(3) for c in range(3))
print("residue (8/hb^2) U :", U)
print("-tr(mu)            :", -tr(mu))

# ---- bonus: Watson's closed form (4.37),  I' = A (I^e)^-1 A ----
Ie = [[sum(dot(e[j],e[j]) for j in range(3))*(a==b) - sum(e[j][a]*e[j][b] for j in range(3)) for b in range(3)] for a in range(3)]
aC = [[[2*sum(dot(e[j],L[b][j]) for j in range(3))*(x==y)
        - sum(L[b][j][x]*e[j][y]+e[j][x]*L[b][j][y] for j in range(3)) for y in range(3)] for x in range(3)] for b in range(3)]
A  = [[Ie[x][y]+0.5*sum(aC[b][x][y]*Q0[b] for b in range(3)) for y in range(3)] for x in range(3)]
AmA = mm(A,mm(inv3(Ie),A))
print("max |I' - A Ie^-1 A|:", max(abs(Ip[x][y]-AmA[x][y]) for x in range(3) for y in range(3)))
print("Lambda_b vs tr(A^-1 a^b):", max(abs(Lam[b]-tr(mm(inv3(A),aC[b]))) for b in range(3)))
```

Try changing the geometry, the masses, or the random seed for the mode
vectors: the agreement persists, because the identity is exact for *any*
Eckart-compliant set — the mode vectors don't even need to diagonalize a
force field.

---

## 9. FAQ

**Q: Where exactly does each Eckart condition enter?**
The translational conditions ($\sum_j\sqrt{M_j}\,\mathbf l_{jb}=0$) silently
kill the translation term of completeness in every sum rule. The rotational
condition ($\sum_j\mathbf e_j\times\mathbf l_{jb}=0$) is the symmetry of
$E^{b}$ (§3.1) — the seed of the closed form. Orthonormality of the
$\mathbf l_b$'s gives the $\delta_{bc}$'s. All three are used; none is
optional.

**Q: Is the closed form $\mathbf I'=A\boldsymbol\mu^{e}A$ an expansion?**
No — $\mathbf I'$ is exactly quadratic in $Q$ (inertia of points moving
linearly), and $A\boldsymbol\mu^{e}A$ is exactly quadratic too; §3.3 matches
all three orders exactly. It holds at arbitrarily large $Q$ (as long as $A$
stays invertible).

**Q: Why does $U$ come out negative (an attractive well)?**
$\operatorname{tr}\boldsymbol\mu=\sum_\alpha 1/I'_\alpha$-like is positive,
so $U<0$ everywhere; it is larger where the moments of inertia are small
(light, compact molecules). For H$_2$O it is of order $-10\,{\rm cm}^{-1}$ —
small but spectroscopically visible.

**Q: The five groups contain $\operatorname{tr}A$, $\boldsymbol\kappa$,
$\boldsymbol\mu^{e}$… why must they cancel?**
Because the answer they add up to, $-\operatorname{tr}\boldsymbol\mu$, doesn't
contain them — but that is hindsight. The honest statement: the cancellations
are consequences of the *same* completeness relation that generated the terms
in the first place. The Eckart frame doesn't just simplify the classical
Hamiltonian; it fine-tunes the quantum residue too.

**Q: What if I use curvilinear internal coordinates (bond lengths, angles)?**
Then $\mathbf I'$ is no longer $A\boldsymbol\mu^{e}A$ with linear $A$,
Miracle 1 fails, and the residue does *not* collapse to a trace — extra
pseudopotential terms survive. Watson's clean result is a special property of
rectilinear normal coordinates.

**Q: One-sentence summary?**
The reweighting residue looks like an arbitrary function of the shape, but the
Eckart conditions force $\det\mathbf I'=(\det A)^2/\det\mathbf I^{e}$ with $A$
linear in $Q$, so every derivative becomes a constant matrix, and the
completeness sum rules make all the complicated terms cancel in pairs —
leaving exactly $-\frac{\hbar^2}{8}\operatorname{tr}\boldsymbol\mu$.

---

*Version note:* equation numbers refer to the current paper (Podolsky sandwich
(4.32′), closed form (4.37), Appendix §B.6 with Eqs. (B.20)–(B.46)). The
companion tutorial `tutorial_coriolis_ordering.md` explains why the *Coriolis*
block leaves no residue at all; `tutorial_quantization.md` covers the older
ordering/anomalous-algebra questions.
