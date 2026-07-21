# From Lab-Frame Kinetic Energy to Nuclear-CoM Coordinates — A Step-by-Step Guide

**Goal of this tutorial.** We start with the plain kinetic-energy Lagrangian of a molecule written in lab-frame Cartesian coordinates,

$$
L \;=\; \frac{1}{2}\sum_j M_j\dot{\mathbf R}_j^{\,2} \;+\; \frac{1}{2}\,m_e\sum_i\dot{\mathbf r}_i^{\,2},
$$

where $\mathbf R_j$ are the positions of the nuclei (mass $M_j$, index $j$) and $\mathbf r_i$ are the positions of the electrons (mass $m_e$, index $i$), all measured from one fixed point in the lab.

We will rewrite $L$ in new coordinates chosen so that:

1. **One coordinate is the *true* center of mass of the whole molecule** (nuclei *and* electrons) — this is the one that moves in a straight line at constant velocity when no outside force acts, exactly like the center of mass you met in introductory mechanics.
2. **All the other coordinates describe the internal shape** of the molecule (how the nuclei and electrons sit relative to one another), referred to the **center of mass of the nuclei only**.

No physics is added or removed — this is just a relabeling of the same $3(N_{\text{nucl}}+N_{\text{elec}})$ numbers. But the relabeling is extremely useful: it will let the total translational motion of the molecule drop out of the internal dynamics completely, so that everything a chemist cares about (bond lengths, vibrations, electron motion) can be studied in a frame that isn't flying across the lab.

Every algebra step below is written out in full — nothing is skipped.

---

## Step 1 — Meet the players

| symbol | meaning |
|---|---|
| $\mathbf R_j$ | lab-frame position of nucleus $j$ |
| $M_j$ | mass of nucleus $j$ |
| $\mathbf r_i$ | lab-frame position of electron $i$ |
| $m_e$ | electron mass (same for every electron) |
| $M_N \equiv \sum_j M_j$ | total **nuclear** mass |
| $M_{\text{tot}} \equiv M_N + N_{\text{elec}}\,m_e$ | total mass of the whole molecule |

$N_{\text{elec}}$ is just the number of electrons, so $N_{\text{elec}}\,m_e = \sum_i m_e$.

---

## Step 2 — Two different "centers of mass"

There are two natural averaging points here, and it matters which one we use for what.

**(a) The center of mass of the nuclei only** — call it $\mathbf R_{\text{ncm}}$:

$$
\mathbf R_{\text{ncm}} \;=\; \frac{1}{M_N}\sum_j M_j\mathbf R_j .
$$

This ignores the electrons entirely. Chemists like this point because it is a good stand-in for "where the molecular skeleton is," independent of the (much lighter, much faster) electrons.

**(b) The center of mass of the *entire* molecule** — call it $\mathbf R_{\text{cm}}$:

$$
\mathbf R_{\text{cm}} \;=\; \frac{1}{M_{\text{tot}}}\Bigl(\sum_j M_j\mathbf R_j + m_e\sum_i\mathbf r_i\Bigr).
$$

This is the point that Newton's laws guarantee moves in a straight line at constant velocity for an isolated molecule (total momentum is conserved). It is *not* the same point as $\mathbf R_{\text{ncm}}$, because the electrons pull the true center of mass slightly away from the nuclei-only center of mass.

**The key design choice of this whole derivation:** we will use $\mathbf R_{\text{cm}}$ (total CoM) as our translational coordinate — because that is the one with the clean, force-free motion — but we will measure everyone's *internal* position relative to $\mathbf R_{\text{ncm}}$ (nuclear CoM) — because that is the point chemists actually want as their internal reference. Using two different reference points sounds odd at first, but it is exactly what makes the final answer clean, as we'll see.

---

## Step 3 — Define the new coordinates

Keep $\mathbf R_{\text{cm}}$ as one coordinate. For everybody else — nuclei *and* electrons — define a position measured relative to the nuclear CoM:

$$
\mathbf R'_j \;=\; \mathbf R_j - \mathbf R_{\text{ncm}}, \qquad
\mathbf r'_i \;=\; \mathbf r_i - \mathbf R_{\text{ncm}}.
$$

So our full new coordinate set is $\bigl(\mathbf R_{\text{cm}},\ \{\mathbf R'_j\},\ \{\mathbf r'_i\}\bigr)$ — exactly as many numbers as we started with ($\{\mathbf R_j\},\{\mathbf r_i\}$), just relabeled.

---

## Step 4 — A free geometric fact we'll need

Straight from the definition of $\mathbf R_{\text{ncm}}$ in Step 2(a):

$$
\sum_j M_j \mathbf R'_j \;=\; \sum_j M_j\bigl(\mathbf R_j - \mathbf R_{\text{ncm}}\bigr)
\;=\; \sum_j M_j\mathbf R_j \;-\; M_N\mathbf R_{\text{ncm}}
\;=\; M_N\mathbf R_{\text{ncm}} - M_N\mathbf R_{\text{ncm}} \;=\; 0.
$$

This isn't an extra assumption — it's automatically true for *any* configuration, just because of how $\mathbf R_{\text{ncm}}$ was defined (a "mass-weighted average sits at zero relative to itself"). Differentiating in time,

$$
\sum_j M_j \dot{\mathbf R}'_j = 0 . \tag{4.1}
$$

Keep this identity in your pocket — it is what makes a cross-term vanish in Step 6.

---

## Step 5 — The "train" trick: writing old velocities in terms of new ones

To rewrite $L$ we need $\dot{\mathbf R}_j$ and $\dot{\mathbf r}_i$ in terms of the new variables. First, invert Step 3:

$$
\mathbf R_j = \mathbf R'_j + \mathbf R_{\text{ncm}}, \qquad \mathbf r_i = \mathbf r'_i + \mathbf R_{\text{ncm}}.
$$

Both equations still contain $\mathbf R_{\text{ncm}}$, which is *not* one of our new coordinates ($\mathbf R_{\text{cm}}$ is). So we need to re-express $\mathbf R_{\text{ncm}}$ itself using the new coordinates. Start from the definition of $\mathbf R_{\text{cm}}$ in Step 2(b) and substitute $\mathbf R_j = \mathbf R'_j+\mathbf R_{\text{ncm}}$, $\mathbf r_i = \mathbf r'_i+\mathbf R_{\text{ncm}}$:

$$
\mathbf R_{\text{cm}} = \frac{1}{M_{\text{tot}}}\Bigl(\sum_j M_j[\mathbf R'_j+\mathbf R_{\text{ncm}}] + m_e\sum_i[\mathbf r'_i+\mathbf R_{\text{ncm}}]\Bigr).
$$

Split the sum:

$$
\mathbf R_{\text{cm}} = \frac{1}{M_{\text{tot}}}\Bigl(\underbrace{\sum_j M_j\mathbf R'_j}_{=\,0\ \text{by (4.1)}} + M_N\mathbf R_{\text{ncm}} + m_e\sum_i\mathbf r'_i + N_{\text{elec}}m_e\mathbf R_{\text{ncm}}\Bigr).
$$

The first term is zero by Step 4. The two $\mathbf R_{\text{ncm}}$ pieces combine to $(M_N+N_{\text{elec}}m_e)\mathbf R_{\text{ncm}} = M_{\text{tot}}\mathbf R_{\text{ncm}}$. So:

$$
\mathbf R_{\text{cm}} = \mathbf R_{\text{ncm}} + \frac{m_e}{M_{\text{tot}}}\sum_i\mathbf r'_i
\quad\Longrightarrow\quad
\mathbf R_{\text{ncm}} = \mathbf R_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\sum_i\mathbf r'_i. \tag{5.1}
$$

This says something physically simple: the true center of mass sits slightly off from the nuclear center of mass, displaced by the electrons' own (mass-weighted) spread $\sum_i\mathbf r'_i$.

Now substitute (5.1) back into $\mathbf R_j = \mathbf R'_j+\mathbf R_{\text{ncm}}$ and $\mathbf r_i=\mathbf r'_i+\mathbf R_{\text{ncm}}$:

$$
\mathbf R_j = \mathbf R'_j + \mathbf R_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\sum_l\mathbf r'_l, \qquad
\mathbf r_i = \mathbf r'_i + \mathbf R_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\sum_l\mathbf r'_l.
$$

Differentiate with respect to time:

$$
\dot{\mathbf R}_j = \dot{\mathbf R}'_j + \mathbf v, \qquad \dot{\mathbf r}_i = \dot{\mathbf r}'_i + \mathbf v,
\qquad\text{where}\qquad
\mathbf v \;\equiv\; \dot{\mathbf R}_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\sum_l\dot{\mathbf r}'_l . \tag{5.2}
$$

**This is the whole trick.** Every single particle's lab velocity — nucleus or electron — is its own *internal* velocity ($\dot{\mathbf R}'_j$ or $\dot{\mathbf r}'_i$) plus one and the same extra velocity $\mathbf v$, shared by everybody. Think of a molecule as passengers on a moving train: each passenger walks around inside the train car at their own internal velocity, but every passenger also shares the train's common drift velocity $\mathbf v$. That $\mathbf v$ is *not* simply $\dot{\mathbf R}_{\text{cm}}$ — it has a small correction from the electrons, which is exactly the "mass polarization" effect we're about to uncover.

---

## Step 6 — Substitute into $L$ and expand

Plug (5.2) into the original Lagrangian:

$$
L = \frac{1}{2}\sum_j M_j\bigl(\dot{\mathbf R}'_j+\mathbf v\bigr)^2 + \frac{1}{2}m_e\sum_i\bigl(\dot{\mathbf r}'_i+\mathbf v\bigr)^2 .
$$

Expand each square, $(\mathbf a+\mathbf b)^2 = \mathbf a^2+2\mathbf a\cdot\mathbf b+\mathbf b^2$:

$$
L = \frac{1}{2}\sum_j M_j\Bigl[\dot{\mathbf R}_j^{\prime\,2} + 2\dot{\mathbf R}'_j\cdot\mathbf v + \mathbf v^2\Bigr]
+ \frac{1}{2}m_e\sum_i\Bigl[\dot{\mathbf r}_i^{\prime\,2} + 2\dot{\mathbf r}'_i\cdot\mathbf v + \mathbf v^2\Bigr].
$$

Group term by term:

$$
L = \underbrace{\frac{1}{2}\sum_j M_j\dot{\mathbf R}_j^{\prime\,2}}_{\text{nuclear internal KE}}
+ \underbrace{\frac{1}{2}m_e\sum_i\dot{\mathbf r}_i^{\prime\,2}}_{\text{electron internal KE}}
+ \underbrace{\Bigl(\sum_j M_j\dot{\mathbf R}'_j\Bigr)\cdot\mathbf v}_{\text{(A)}}
+ \underbrace{m_e\Bigl(\sum_i\dot{\mathbf r}'_i\Bigr)\cdot\mathbf v}_{\text{(B)}}
+ \underbrace{\frac{1}{2}\bigl(M_N+N_{\text{elec}}m_e\bigr)\mathbf v^2}_{\text{(C)}\ =\ \frac12 M_{\text{tot}}\mathbf v^2}.
$$

**Term (A) vanishes immediately**, by the identity (4.1) we saved in Step 4: $\sum_j M_j\dot{\mathbf R}'_j=0$.

We're left with

$$
L = \frac{1}{2}\sum_j M_j\dot{\mathbf R}_j^{\prime\,2} + \frac{1}{2}m_e\sum_i\dot{\mathbf r}_i^{\prime\,2}
+ \frac{1}{2}M_{\text{tot}}\mathbf v^2 + m_e\Bigl(\sum_i\dot{\mathbf r}'_i\Bigr)\cdot\mathbf v . \tag{6.1}
$$

---

## Step 7 — Unpack $\mathbf v$ and watch the last cross-term cancel

Everything in (6.1) is already in the coordinates we want, *except* $\mathbf v$, which by (5.2) is $\mathbf v = \dot{\mathbf R}_{\text{cm}} - \tfrac{m_e}{M_{\text{tot}}}\mathbf S$, where for brevity we write $\mathbf S\equiv\sum_i\dot{\mathbf r}'_i$ (the total internal electron velocity).

**Square $\mathbf v$:**

$$
\mathbf v^2 = \Bigl(\dot{\mathbf R}_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\mathbf S\Bigr)^{\!2}
= \dot{\mathbf R}_{\text{cm}}^{\,2} - \frac{2m_e}{M_{\text{tot}}}\,\dot{\mathbf R}_{\text{cm}}\cdot\mathbf S + \frac{m_e^2}{M_{\text{tot}}^2}\mathbf S^2 .
$$

Multiply by $\tfrac12 M_{\text{tot}}$:

$$
\frac{1}{2}M_{\text{tot}}\mathbf v^2 = \frac{1}{2}M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2}
\;-\; m_e\,\dot{\mathbf R}_{\text{cm}}\cdot\mathbf S
\;+\; \frac{m_e^2}{2M_{\text{tot}}}\mathbf S^2 . \tag{7.1}
$$

**Now the last term of (6.1),** $m_e\,\mathbf S\cdot\mathbf v$:

$$
m_e\,\mathbf S\cdot\mathbf v = m_e\,\mathbf S\cdot\Bigl(\dot{\mathbf R}_{\text{cm}} - \frac{m_e}{M_{\text{tot}}}\mathbf S\Bigr)
= m_e\,\mathbf S\cdot\dot{\mathbf R}_{\text{cm}} - \frac{m_e^2}{M_{\text{tot}}}\mathbf S^2 . \tag{7.2}
$$

**Add (7.1) and (7.2).** The two $\pm m_e\dot{\mathbf R}_{\text{cm}}\cdot\mathbf S$ pieces are equal and opposite — they cancel exactly:

$$
\frac{1}{2}M_{\text{tot}}\mathbf v^2 + m_e\,\mathbf S\cdot\mathbf v
= \frac{1}{2}M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2}
\;\underbrace{-\,m_e\dot{\mathbf R}_{\text{cm}}\cdot\mathbf S + m_e\dot{\mathbf R}_{\text{cm}}\cdot\mathbf S}_{=\,0}
\;+\;\frac{m_e^2}{2M_{\text{tot}}}\mathbf S^2 - \frac{m_e^2}{M_{\text{tot}}}\mathbf S^2
$$
$$
= \frac{1}{2}M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2} - \frac{m_e^2}{2M_{\text{tot}}}\mathbf S^2 . \tag{7.3}
$$

---

## Step 8 — The final Lagrangian

Insert (7.3) back into (6.1), and restore $\mathbf S=\sum_i\dot{\mathbf r}'_i$:

$$
\boxed{\;
L \;=\; \frac{1}{2}M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2}
\;+\; \frac{1}{2}\sum_j M_j\dot{\mathbf R}_j^{\prime\,2}
\;+\; \frac{1}{2}m_e\sum_i\dot{\mathbf r}_i^{\prime\,2}
\;-\; \frac{m_e^2}{2M_{\text{tot}}}\Bigl(\sum_i\dot{\mathbf r}'_i\Bigr)^{\!2}
\;}
$$

That's it — four terms, and every one of them has a clear physical meaning:

1. $\dfrac{1}{2}M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2}$ — **free translational motion of the whole molecule**, carrying the *correct total mass* $M_{\text{tot}}$ (nuclei + electrons together). This term involves *only* $\mathbf R_{\text{cm}}$ — it does not mix with anything internal.
2. $\dfrac{1}{2}\sum_j M_j\dot{\mathbf R}_j^{\prime\,2}$ — internal kinetic energy of the nuclei, relative to their own center of mass.
3. $\dfrac{1}{2}m_e\sum_i\dot{\mathbf r}_i^{\prime\,2}$ — internal kinetic energy of the electrons, relative to the nuclear center of mass.
4. $-\dfrac{m_e^2}{2M_{\text{tot}}}\Bigl(\sum_i\dot{\mathbf r}'_i\Bigr)^{2}$ — the **electron mass-polarization term**. It's a small correction (suppressed by $m_e/M_{\text{tot}}$, a tiny ratio) that couples the electrons to each other through their common recoil against the nuclei. It appeared because we referred the electrons to the *nuclear* CoM rather than the *total* CoM — see the box below.

If you also want a potential energy $V$ in the full molecular Lagrangian, it can simply be added back on: since $V$ only ever depends on the *distances* between particles (Coulomb interactions), and every pairwise distance is left untouched by this whole change of coordinates ($\mathbf R'_j-\mathbf R'_k=\mathbf R_j-\mathbf R_k$, and likewise for any other pair), $V$ just carries over unchanged as $-V(\mathbf R'_j,\mathbf r'_i)$. None of the algebra above needed to touch it.

---

## Step 9 — Why also go to momentum form?

The boxed $L$ from Step 8 is a perfectly good, complete description — but it's written in terms of *velocities* $\dot{\mathbf R}_{\text{cm}},\dot{\mathbf R}'_j,\dot{\mathbf r}'_i$. In quantum mechanics we don't quantize velocities; we quantize **momenta** (a velocity becomes an operator only *through* the momentum, $\hat{\mathbf p}=-i\hbar\nabla$). So before this can be turned into a Hamiltonian/Schrödinger-equation problem, $L$ has to be rewritten as a function of momenta instead of velocities. That rewriting is a standard, mechanical recipe called the **Legendre transform**, and it always proceeds in the same three moves:

1. **Define** each new momentum as $p_q \equiv \partial L/\partial\dot q$ for every coordinate $q$.
2. **Invert** those relations to get every velocity as a function of the momenta.
3. **Assemble** $H \equiv \sum_q p_q\dot q - L$ and substitute the inverted velocities, so that only momenta (and coordinates) remain.

We now do exactly this to the boxed $L$ of Step 8.

---

## Step 10 — Define the conjugate momenta

Differentiate $L$ with respect to each velocity in turn, holding everything else fixed.

**Translational momentum.** Only the first term of $L$ involves $\dot{\mathbf R}_{\text{cm}}$:

$$
\mathbf P_{\text{cm}} \equiv \frac{\partial L}{\partial\dot{\mathbf R}_{\text{cm}}} = M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}.
$$

**Nuclear-internal momentum.** Only the second term involves $\dot{\mathbf R}'_j$, and only through nucleus $j$ itself:

$$
\mathbf P'_j \equiv \frac{\partial L}{\partial\dot{\mathbf R}'_j} = M_j\dot{\mathbf R}'_j.
$$

Both of these are the ordinary "mass times velocity" you already expect.

**Electron momentum.** This one takes a little more care, because *two* terms of $L$ involve $\dot{\mathbf r}'_i$: the electron kinetic-energy term (through index $i$ alone), and the mass-polarization term $-\tfrac{m_e^2}{2M_{\text{tot}}}\mathbf S^2$, where $\mathbf S=\sum_l\dot{\mathbf r}'_l$ — and $\mathbf S$ contains $\dot{\mathbf r}'_i$ too, since it's a sum over *every* electron including this one. Since $\mathbf S$ is a plain sum, $\partial\mathbf S/\partial\dot{\mathbf r}'_i=1$ (one unit vector's worth, for each Cartesian component), so by the chain rule $\partial(\mathbf S^2)/\partial\dot{\mathbf r}'_i = 2\mathbf S$. Putting both pieces together:

$$
\mathbf p'_i \equiv \frac{\partial L}{\partial\dot{\mathbf r}'_i} = m_e\dot{\mathbf r}'_i - \frac{m_e^2}{M_{\text{tot}}}\mathbf S. \tag{10.1}
$$

Notice this is **not** simply $m_e\dot{\mathbf r}'_i$: the electron's canonical momentum is offset by a term built from the *total* internal electron velocity $\mathbf S$. Physically, this is the recoil bookkeeping again — when *any* electron moves, it nudges the nuclear CoM (through the total-CoM constraint hidden in $\mathbf v$ of Step 5), and *every* electron's coordinate is measured against that same nuclear CoM, so each electron's momentum has to carry a small correction shared by all of them.

---

## Step 11 — Invert: velocities in terms of momenta

$\mathbf P_{\text{cm}}$ and $\mathbf P'_j$ are already solved:

$$
\dot{\mathbf R}_{\text{cm}} = \frac{\mathbf P_{\text{cm}}}{M_{\text{tot}}}, \qquad \dot{\mathbf R}'_j = \frac{\mathbf P'_j}{M_j}.
$$

The electron relation (10.1) still has the unknown $\mathbf S$ tangled up in it, so we do the same move as in Step 5: sum (10.1) over all electrons to isolate $\mathbf S$.

$$
\sum_i\mathbf p'_i = m_e\mathbf S - N_{\text{elec}}\frac{m_e^2}{M_{\text{tot}}}\mathbf S = m_e\mathbf S\Bigl(1-\frac{N_{\text{elec}}m_e}{M_{\text{tot}}}\Bigr).
$$

Since $M_{\text{tot}}=M_N+N_{\text{elec}}m_e$, the bracket is $\bigl(M_{\text{tot}}-N_{\text{elec}}m_e\bigr)/M_{\text{tot}}=M_N/M_{\text{tot}}$, so

$$
\sum_i\mathbf p'_i = \frac{m_eM_N}{M_{\text{tot}}}\mathbf S
\quad\Longrightarrow\quad
\mathbf S = \frac{M_{\text{tot}}}{m_eM_N}\sum_l\mathbf p'_l. \tag{11.1}
$$

Now solve (10.1) for $\dot{\mathbf r}'_i$ and substitute (11.1):

$$
\dot{\mathbf r}'_i = \frac{\mathbf p'_i}{m_e} + \frac{m_e}{M_{\text{tot}}}\mathbf S
= \frac{\mathbf p'_i}{m_e} + \frac{m_e}{M_{\text{tot}}}\cdot\frac{M_{\text{tot}}}{m_eM_N}\sum_l\mathbf p'_l
\;\Longrightarrow\;
\boxed{\;\dot{\mathbf r}'_i = \frac{\mathbf p'_i}{m_e} + \frac{1}{M_N}\sum_l\mathbf p'_l\;}. \tag{11.2}
$$

A clean result: each electron's internal velocity is its own momentum divided by $m_e$, plus one shared "recoil" term $\tfrac1{M_N}\sum_l\mathbf p'_l$ — the momentum-space mirror image of the shared drift velocity $\mathbf v$ we met back in Step 5.

---

## Step 12 — Assemble $\sum p\dot q$

The Legendre-transform recipe says $H=\sum_qp_q\dot q - L$. Compute the sum term by term, using Step 11:

$$
\mathbf P_{\text{cm}}\cdot\dot{\mathbf R}_{\text{cm}} = \frac{\mathbf P_{\text{cm}}^{\,2}}{M_{\text{tot}}}, \qquad
\sum_j\mathbf P'_j\cdot\dot{\mathbf R}'_j = \sum_j\frac{\mathbf P_j^{\prime\,2}}{M_j}.
$$

For the electrons, dot (11.2) into $\mathbf p'_i$ and sum over $i$:

$$
\sum_i\mathbf p'_i\cdot\dot{\mathbf r}'_i
= \sum_i\mathbf p'_i\cdot\Bigl(\frac{\mathbf p'_i}{m_e}+\frac{1}{M_N}\sum_l\mathbf p'_l\Bigr)
= \frac{1}{m_e}\sum_i\mathbf p_i^{\prime\,2} + \frac{1}{M_N}\Bigl(\sum_i\mathbf p'_i\Bigr)^{\!2},
$$

where in the last step $\sum_i\mathbf p'_i\cdot\sum_l\mathbf p'_l=\bigl(\sum_i\mathbf p'_i\bigr)^2$ because it's the same sum dotted into itself. Collecting everything:

$$
\sum_qp_q\dot q = \frac{\mathbf P_{\text{cm}}^{\,2}}{M_{\text{tot}}} + \sum_j\frac{\mathbf P_j^{\prime\,2}}{M_j} + \frac{1}{m_e}\sum_i\mathbf p_i^{\prime\,2} + \frac{1}{M_N}\Bigl(\sum_i\mathbf p'_i\Bigr)^{\!2}. \tag{12.1}
$$

---

## Step 13 — Subtract $L$ and simplify: the Hamiltonian

We still need $L$ itself written in momentum variables, so substitute (11.2) into the boxed Lagrangian of Step 8. The translational and nuclear terms are immediate:

$$
\frac12M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}^{\,2} = \frac{\mathbf P_{\text{cm}}^{\,2}}{2M_{\text{tot}}}, \qquad
\frac12\sum_jM_j\dot{\mathbf R}_j^{\prime\,2} = \sum_j\frac{\mathbf P_j^{\prime\,2}}{2M_j}.
$$

For the electron sector, write $\mathbf c\equiv\tfrac1{M_N}\sum_l\mathbf p'_l$ for short, so (11.2) reads $\dot{\mathbf r}'_i=\mathbf p'_i/m_e+\mathbf c$ and $\mathbf S=\sum_i\dot{\mathbf r}'_i = \tfrac1{m_e}\sum_i\mathbf p'_i+N_{\text{elec}}\mathbf c = \mathbf c\bigl(\tfrac{M_N}{m_e}+N_{\text{elec}}\bigr)=\dfrac{M_{\text{tot}}}{m_e}\,\mathbf c$ (using $\sum_i\mathbf p'_i=M_N\mathbf c$ and $M_{\text{tot}}=M_N+N_{\text{elec}}m_e$). Then, term by term:

$$
\frac12m_e\sum_i\dot{\mathbf r}_i^{\prime\,2} = \frac12m_e\sum_i\Bigl(\frac{\mathbf p'_i}{m_e}+\mathbf c\Bigr)^{\!2}
= \frac{1}{2m_e}\sum_i\mathbf p_i^{\prime\,2} + \mathbf c\cdot\sum_i\mathbf p'_i + \frac12N_{\text{elec}}m_e\mathbf c^2
= \frac{1}{2m_e}\sum_i\mathbf p_i^{\prime\,2} + M_N\mathbf c^2 + \frac12N_{\text{elec}}m_e\mathbf c^2,
$$

$$
-\frac{m_e^2}{2M_{\text{tot}}}\mathbf S^2 = -\frac{m_e^2}{2M_{\text{tot}}}\cdot\frac{M_{\text{tot}}^2}{m_e^2}\mathbf c^2 = -\frac{M_{\text{tot}}}{2}\mathbf c^2 .
$$

Adding these two, the $\mathbf c^2$ coefficients are $M_N+\tfrac12N_{\text{elec}}m_e-\tfrac12M_{\text{tot}} = M_N+\tfrac12N_{\text{elec}}m_e-\tfrac12(M_N+N_{\text{elec}}m_e)=\tfrac12M_N$, so every $\mathbf c^2$ term collapses to $\tfrac12M_N\mathbf c^2=\tfrac1{2M_N}\bigl(\sum_i\mathbf p'_i\bigr)^2$ (using $M_N\mathbf c=\sum_i\mathbf p'_i$). Hence

$$
L\bigl(\dot q(p)\bigr) = \frac{\mathbf P_{\text{cm}}^{\,2}}{2M_{\text{tot}}} + \sum_j\frac{\mathbf P_j^{\prime\,2}}{2M_j} + \frac{1}{2m_e}\sum_i\mathbf p_i^{\prime\,2} + \frac{1}{2M_N}\Bigl(\sum_i\mathbf p'_i\Bigr)^{\!2}. \tag{13.1}
$$

Comparing (13.1) to (12.1), the electron and nuclear/translational pieces of (12.1) are each exactly **twice** the corresponding piece of (13.1) — as they must be, since kinetic energy is quadratic in the velocities. Now finish the recipe, $H=\sum p\dot q - L$, i.e. (12.1) $-$ (13.1):

$$
\boxed{\;
H \;=\; \frac{\mathbf P_{\text{cm}}^{\,2}}{2M_{\text{tot}}}
\;+\; \sum_j\frac{\mathbf P_j^{\prime\,2}}{2M_j}
\;+\; \sum_i\frac{\mathbf p_i^{\prime\,2}}{2m_e}
\;+\; \frac{1}{2M_N}\Bigl(\sum_i\mathbf p'_i\Bigr)^{\!2}
\;}
$$

and, if a potential was carried along, simply add $+V(\mathbf R'_j,\mathbf r'_i)$ (a Legendre transform never touches terms with no velocity dependence, so $V$ passes straight through with a flipped sign relative to its appearance in $L$).

Every term keeps the meaning it had in $L$: free translation of the whole molecule with mass $M_{\text{tot}}$; ordinary kinetic energy for each nucleus and each electron; and the mass-polarization term, which has *changed sign and form* — in $L$ it was a *negative* cross term built from velocities, $-\tfrac{m_e^2}{2M_{\text{tot}}}(\sum\dot{\mathbf r}')^2$; in $H$ it is a *positive* term built from momenta, $+\tfrac1{2M_N}(\sum\mathbf p')^2$. Both describe the same physics — the Legendre transform simply presents it from the other side.

---

## Step 14 — Sanity checks

A good derivation should reduce to something you already trust in simple limits.

- **No electrons at all** ($N_{\text{elec}}=0$): then $M_{\text{tot}}=M_N$, term 3 and term 4 vanish, and $L=\tfrac12 M_N\dot{\mathbf R}_{\text{cm}}^2+\tfrac12\sum_jM_j\dot{\mathbf R}_j'^2$ — exactly the ordinary "CoM + relative motion" split you'd write for a bare cluster of nuclei.
- **One electron** ($N_{\text{elec}}=1$, so there's only one $\mathbf r'$ and the sum in term 4 has a single entry): terms 3 and 4 combine to
$$
\frac12 m_e\dot{\mathbf r}'^{\,2} - \frac{m_e^2}{2M_{\text{tot}}}\dot{\mathbf r}'^{\,2}
= \frac12 m_e\Bigl(1-\frac{m_e}{M_{\text{tot}}}\Bigr)\dot{\mathbf r}'^{\,2}
= \frac12\mu_e\dot{\mathbf r}'^{\,2},
\qquad \mu_e \equiv \frac{m_eM_N}{M_{\text{tot}}},
$$
which is exactly the familiar **reduced mass** $\mu_e$ from the two-body problem you meet in introductory quantum chemistry (e.g. the hydrogen atom, where $\mu_e = m_eM_{\text{nucleus}}/(m_e+M_{\text{nucleus}})$).
- **Total momentum is conserved automatically:** because $L$ has no $\mathbf R_{\text{cm}}$ dependence at all (only $\dot{\mathbf R}_{\text{cm}}$ appears), $\mathbf R_{\text{cm}}$ is a *cyclic coordinate* — its conjugate momentum $\mathbf P_{\text{cm}}=M_{\text{tot}}\dot{\mathbf R}_{\text{cm}}$ is conserved for any motion, and the translational and internal parts of the dynamics **never talk to each other**. That was the whole point of the exercise.
- **The Hamiltonian's one-electron limit reproduces the same reduced mass.** For $N_{\text{elec}}=1$ the last two terms of the boxed $H$ collapse to a single momentum, $\dfrac{\mathbf p'^{\,2}}{2m_e}+\dfrac{\mathbf p'^{\,2}}{2M_N}=\dfrac{\mathbf p'^{\,2}}{2}\Bigl(\dfrac1{m_e}+\dfrac1{M_N}\Bigr)=\dfrac{\mathbf p'^{\,2}}{2\mu_e}$, with $\mu_e=m_eM_N/(m_e+M_N)=m_eM_N/M_{\text{tot}}$ — the *same* $\mu_e$ found for the velocity-form $L$ above. Velocity form and momentum form must always agree on such limits; here they do.
- **The mass-polarization term flips character between $L$ and $H$**, but never disappears: $-\tfrac{m_e^2}{2M_{\text{tot}}}(\sum\dot{\mathbf r}')^2$ in velocity form becomes $+\tfrac1{2M_N}(\sum\mathbf p')^2$ in momentum form. Neither form is "more correct" — they're the same physics, related by the Legendre transform of Steps 9–13.

---

## Recap in one sentence

By describing the molecule with *one* coordinate for the true, total center of mass, and *all the rest* relative to the center of mass of the nuclei alone, the kinetic energy splits cleanly into "the whole molecule drifting through space" (with the right total mass) plus "everything happening inside it" (nuclear vibration/rotation and electron motion, joined by one small mass-polarization correction) — with **no cross-terms and no approximation** — and the same clean split survives, term for term, when the Legendre transform rewrites the whole thing in momentum form.

---

*This matches the Hamiltonian derived by a different (canonical, generating-function) route in `sections/02_mass_centered_coordinates.tex` and the "Stage 2" entry of `markdown/lagrangians_hamiltonians.md` — reassuring, since conjugate momenta don't depend on which route you take to find them.*
