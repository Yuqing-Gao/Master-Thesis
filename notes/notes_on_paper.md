# Notes on the paper by Mandel, Brezina & Vaněk

## Contents

- [1 Basic concepts, backgrounds](#1-basic-concepts-backgrounds)
  - [1.1 Background and motivation](#11-background-and-motivation)
  - [1.2 The main problem: fixed pattern energy minimization - the scalar cases](#12-the-main-problem-fixed-pattern-energy-minimization---the-scalar-cases)
    - [1.2.1 The allowed update directions](#121-the-allowed-update-directions)
    - [1.2.2 The main algorithm](#122-the-main-algorithm)
    - [1.2.3 1D scalar case example](#123-1d-scalar-case-example)
    - [1.2.4 2D scalar case example](#124-2d-scalar-case-example)
  - [1.3 The block cases](#13-the-block-cases)
- [2 The Problem I Need to Solve](#2-the-problem-i-need-to-solve)
  - [2.1 How to decide $`\{N^{(k)}\}`$](#21-how-to-decide-nk)
    - [2.1.1 Coupling strength](#211-coupling-strength)
    - [2.1.2 Connection score](#212-connection-score)
    - [2.1.3 Energy score](#213-energy-score)
  - [2.2 How to decide $`N^{(k)} \mapsto P^{(k)}`$](#22-how-to-decide-nk-mapsto-pk)

## 1 Basic concepts, backgrounds

### 1.1 Background and motivation

Given a linear system $`Ax=b`$, where $`A\in\mathbb{R}^{n\times n}`$ denotes the fine level matrix. Here $`A`$ is assumed to be SPD. Write the prolongation as $`P=[p_1,\ldots,p_m]\in\mathbb{R}^{n\times m}`$. For an arbitrary coarse level vector $`x_c`$, $`Px_c`$ is a linear combination of $`p_1,\ldots,p_m`$, i.e. $`\operatorname{range}(P)=\operatorname{span}\{p_1,\ldots,p_m\}`$. The coarse correction in AMG is $`x\leftarrow x+Pe`$, where $`e\in\mathbb{R}^m`$. The error that the smoother cannot eliminate should be well approximated by a linear combination of the column vectors of $`P`$. It does not have to be exactly represented in every case. Thus the choice of $`P`$ is important.

The authors state that a good $`P`$ should satisfy 3 mutual restraint conditions:

- minimization of energy: the energy $`p_j^TAp_j`$ of $`p_j`$ should be as little as possible. Later the authors minimize the sum of these energies.
- preservation of null space: On the interior nodes, $`\operatorname{span}(P)`$ should contain the zero-energy modes. In the scalar case, for instance the diffusion problem, the zero-energy modes are constanst vectors, say multiples of $`\mathbf{1}`$. In linear elasticity, the zero-energy modes are rigid-body translations and rotations, whose discrete representations form the matrix $`B`$.
- limited overlap: the supports of every $`p_j`$ should not be too large, and different supports should not overlap too much. Otherwise $`P`$ and $`A_c=P^TAP`$ would be too dense.

We check first if the $`P`$ in unsmoothed and smoothed aggregation multigrid satisfy those conditions.

**unsmoothed aggregation** Assume the fine level nodes are divided into $`m`$ non-overlapping, nonempty aggregates $`\mathcal{A}_1,\ldots,\mathcal{A}_m`$ which cover all fine nodes. Define the tentative prolongator

```math
\widehat{P}_{ij}=\begin{cases}
1,&i\in\mathcal{A}_j,\\
0,&\text{otherwise.}
\end{cases}
```

The tentative prolongator satisfies $`\widehat{P}\mathbf{1}_m=\mathbf{1}_n`$ (there is exactly one 1 in every row of $`\widehat{P}`$). Thus $`\mathbf{1}_n\in\operatorname{range}(\widehat{P})`$. For the actual SPD Dirichlet matrix, this vector is not in $`\ker A`$; it comes from the constant mode of the problem without Dirichlet constraints.

But the energy of $`\hat{p}_j`$ can be large, because its values change fast at aggregate boundaries. We want to lower the energy by smoothing.

**smoothed aggregation** Let $`P=S\widehat{P}`$, where $`S=q(D^{-1}A)`$, $`q(0)=1`$, $`q`$ is a polynomial, and $`D`$ is a suitable SPD matrix. A simple choice is $`D=\operatorname{diag}(A)`$. If $`A\mathbf{1}_n=0`$, then

```math
\begin{aligned}
P\mathbf{1}_m&=S\widehat{P}\mathbf{1}_m\\
&=S\mathbf{1}_n\\
&=q(D^{-1}A)\mathbf{1}_n\\
&=q(0)\mathbf{1}_n&& (A\mathbf{1}_n=0)\\
&=\mathbf{1}_n&& (q(0)=1).
\end{aligned}
```

*Here $`A`$ is assumed to be SPD, which contradicts $`A\mathbf{1}_n=0`$. So this is not a global proof for the SPD Dirichlet problem.*

Under the assumptions of Remark 4.5, the first iteration of Algorithm 1 starting from the unsmoothed prolongator $`\widehat{P}`$ reduces to

```math
P=(I-\omega A)\widehat{P}\implies P\mathbf{1}_m=\mathbf{1}_n-\omega A\mathbf{1}_n.
```

Thus a row with $`(A\mathbf{1}_n)_i=0`$ still has row sum 1 after this step. The authors use $`(f_i)`$ to choose the rows where they require this.

But smoothing can enlarge the supports. Repeated smoothing may enlarge them again and again.

### 1.2 The main problem: fixed pattern energy minimization - the scalar cases

The authors’ idea is to give an allowed pattern first, and then lower the energy without going outside this pattern. One useful choice is larger than the pattern of $`\widehat{P}`$: it allows the aggregates and their immediate neighbours.

What is $`N`$? It is a given 0–1 matrix with the same size as $`P`$, and $`n_{ij}=0\implies p_{ij}=0`$. So $`N`$ is an upper bound on the support. $`n_{ij}=1`$ only means this position is allowed; $`p_{ij}`$ can still be zero. The authors define a linear space

```math
\mathcal{N}=\{P\in\mathbb{R}^{n\times m}:p_{ij}=0\text{ if }n_{ij}=0\}.
```

The problem they want to solve is

```math
\begin{aligned}
\min_P\ J(P)&=\frac12\sum_{j=1}^m p_j^TAp_j=\frac12(P,AP),\\
\text{s.t.}\quad&P\in\mathcal{N},\\
&\sum_j p_{ij}=1\quad\text{if }f_i=1.
\end{aligned}
```

Here $`(X,Y)=\sum_{i,j}x_{ij}y_{ij}=\operatorname{tr}(X^TY)`$ is the Frobenius inner product, which is a number. In contrast, $`X*Y=(x_{ij}y_{ij})`$ is the term-by-term product, which is still a matrix.

The vector $`F=(f_i)`$ is used for the Dirichlet boundary:

```math
f_i=\begin{cases}
0,&\text{node }i\text{ belongs to an element with a Dirichlet constraint,}\\
1,&\text{otherwise.}
\end{cases}
```

So it does not only ask whether node $`i`$ itself is a Dirichlet node. A free node next to the Dirichlet boundary can also have $`f_i=0`$. If $`f_i=1`$, its row sum must be 1. If $`f_i=0`$, there is no row-sum condition on this row.

Since $`A`$ is SPD, $`J(P)`$ is strictly convex. Assuming that there is a feasible $`P_0`$, the problem has a unique minimizer $`P^*`$.

#### 1.2.1 The allowed update directions
The authors want to use projected gradient descent. The feasible set is

```math
\mathcal{F}=\left\{P\in\mathcal{N}:\sum_j p_{ij}=1\text{ if }f_i=1\right\}.
```

It is not a linear space in general. So when updating $`P`$, we need to avoid breaking the row-sum-1 condition. The authors introduce the subspace

```math
\mathcal{Z}=\left\{Q\in\mathcal{N}:\sum_j q_{ij}=0\text{ if }f_i=1\right\}\subseteq\mathcal{N}.
```

This is the space of allowed changes.

The gradient of $`J`$ in the full matrix space is $`AP`$. The authors prove that its gradient in $`\mathcal{N}`$ is $`N*(AP)`$. For any $`Q\in\mathcal{N}`$,

```math
\left.\frac{d}{dt}J(P+tQ)\right|_{t=0}=(Q,AP)=(Q,N*(AP)).
```

The term-by-term product makes the gradient zero at forbidden entries, so $`N*(AP)\in\mathcal{N}`$. But this is not enough to keep the row sums. We still need to project it onto $`\mathcal{Z}`$.

The authors write this projection as $`Z:S\mapsto Q`$, where $`S\in\mathcal{N}`$ is just an input matrix and $`Q\in\mathcal{Z}`$ is the output. This $`S`$ is not the smoother from the previous paragraph. The row formula is

```math
q_{i*}=s_{i*}-f_i\frac{(s_{i*},n_{i*})}{(n_{i*},n_{i*})}n_{i*}.
```

$`n_{i*}`$ is the $`i`$-th row of $`N`$, which has 1 at the allowed entries. So $`(s_{i*},n_{i*})`$ is the sum of the allowed entries in row $`i`$, and

```math
(n_{i*},n_{i*})=\sum_j n_{ij}^2=\#\{j:n_{ij}=1\}.
```

The fraction is therefore the average of the allowed entries of this row. If $`f_i=1`$, subtract this average from every allowed entry, and the row sum becomes zero. If $`f_i=0`$, do not change the row. A constrained row must have at least one allowed entry. For an unconstrained empty row, the projection returns zero without evaluating this fraction.

#### 1.2.2 The main algorithm 
Given a feasible $`P_0`$, the authors rewrite the problem as its first-order condition: find $`P`$ such that

```math
P-P_0\in\mathcal{Z},\qquad Z(N*(AP))=0.
```

Since $`J`$ is strictly convex on the affine feasible set, this condition is sufficient too. The unique point satisfying it is $`P^*`$.

Algorithm 1 is then:

1. Find a feasible $`P_0`$ and choose $`0<\omega<2/\lambda_{\max}(A)`$.
2. Update $`P_{t+1}=P_t-\omega Z(N*(AP_t))`$, for $`t=0,1,\ldots`$.

I use $`t`$ for the iteration number here, since $`n`$ already means the number of fine nodes.

In Remark 4.5, starting from $`\widehat{P}`$, with the one-neighbour pattern and the stated boundary conditions, the first iteration gives

```math
P_1=(I-\omega A)\widehat{P}.
```

#### 1.2.3 1D scalar case example

choose thirteen positive edge weights

```math
(w_0, \ldots, w_{12})
=
(1, 1, 0.25, 1, 2, 0.5, 0.8, 1, 3, 1, 0.4, 1, 1).
```

They describe the weighted chain

```math
\partial_L
\overset{1}{-} 1
\overset{1}{-} 2
\overset{0.25}{-} 3
\overset{1}{-} 4
\overset{2}{-} 5
\overset{0.5}{-} 6
\overset{0.8}{-} 7
\overset{1}{-} 8
\overset{3}{-} 9
\overset{1}{-} 10
\overset{0.4}{-} 11
\overset{1}{-} 12
\overset{1}{-} \partial_R.
```

Those weights are to describe the coupling strength. The matrix $`A \in \mathbb{R}^{12 \times 12}`$ is therefore defined by

```math
a_{ii} = w_{i-1} + w_i,
\qquad
a_{i,i+1} = a_{i+1,i} = -w_i.
```

And the full matrix is

```math
A =
\left(
\begin{array}{rrrrrrrrrrrr}
2  & -1    & 0     & 0  & 0    & 0    & 0    & 0  & 0  & 0    & 0    & 0 \\
-1 & 1.25  & -0.25 & 0  & 0    & 0    & 0    & 0  & 0  & 0    & 0    & 0 \\
0  & -0.25 & 1.25  & -1 & 0    & 0    & 0    & 0  & 0  & 0    & 0    & 0 \\
0  & 0     & -1    & 3  & -2   & 0    & 0    & 0  & 0  & 0    & 0    & 0 \\
0  & 0     & 0     & -2 & 2.5  & -0.5 & 0    & 0  & 0  & 0    & 0    & 0 \\
0  & 0     & 0     & 0  & -0.5 & 1.3  & -0.8 & 0  & 0  & 0    & 0    & 0 \\
0  & 0     & 0     & 0  & 0    & -0.8 & 1.8  & -1 & 0  & 0    & 0    & 0 \\
0  & 0     & 0     & 0  & 0    & 0    & -1   & 4  & -3 & 0    & 0    & 0 \\
0  & 0     & 0     & 0  & 0    & 0    & 0    & -3 & 4  & -1   & 0    & 0 \\
0  & 0     & 0     & 0  & 0    & 0    & 0    & 0  & -1 & 1.4  & -0.4 & 0 \\
0  & 0     & 0     & 0  & 0    & 0    & 0    & 0  & 0  & -0.4 & 1.4  & -1 \\
0  & 0     & 0     & 0  & 0    & 0    & 0    & 0  & 0  & 0    & -1   & 2
\end{array}
\right).
```

The edge weights are only used to generate this test matrix. In a real AMG problem, $`A`$ would already be given.

**fixed aggregates**

The twelve nodes are artificially divided into six aggregates:

```math
\begin{aligned}
\mathcal{A}_1 &= \{1,2\}, &
\mathcal{A}_2 &= \{3,4\}, &
\mathcal{A}_3 &= \{5,6\}, \\
\mathcal{A}_4 &= \{7,8\}, &
\mathcal{A}_5 &= \{9,10\}, &
\mathcal{A}_6 &= \{11,12\}.
\end{aligned}
```

Let $`a(i)`$ denote the index of the aggregate containing node $`i`$. Then

```math
(a(1), \ldots, a(12))
=
(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6).
```

The tentative unsmoothed prolongator is

```math
\widehat{P} =
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}.
```

For instance, row 5 is $`(0,0,1,0,0,0)`$, which means that fine node 5 initially uses only the coarse variable corresponding to aggregate $`\mathcal{A}_3`$.

**boundary indicator**

The scalar boundary indicator is

```math
f = (0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)^T.
```

According to the definition of $`(f_i)`$ in the paper, rows 2 to 11 must satisfy the row-sum constraint

```math
\sum_j p_{ij} = 1.
```

Rows 1 and 12 touch a Dirichlet boundary element, so the row-sum constraint is not imposed there.

**initial pattern**

The initial pattern is

```math
N^{(0)} = \operatorname{pattern}(\widehat{P}).
```

Thus

```math
N^{(0)} =
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}.
```

**construction of all candidate positions**

For every nonzero entry $`a_{ik}`$, node $`i`$ is directly coupled to node $`k`$. If node $`k`$ belongs to aggregate $`\mathcal{A}_j`$, then node $`i`$ is allowed to use coarse variable $`j`$ in the full one-neighbour pattern.

Therefore,

```math
N^{\mathrm{full}}_{ij} = 1
```

if

```math
j = a(i)
```

or if there exists $`k \in \mathcal{A}_j`$ such that

```math
a_{ik} \ne 0.
```

For instance, row 5 of $`A`$ has nonzero positions

```math
k = 4,\ 5,\ 6.
```

These nodes belong to aggregates

```math
\mathcal{A}_2,\ \mathcal{A}_3,\ \mathcal{A}_3.
```

Therefore row 5 of the full pattern allows columns 2 and 3:

```math
N^{\mathrm{full}}_{5*} = (0,1,1,0,0,0).
```

The original pattern allowed only column 3, so

```math
(5,2)
```

is a new candidate position.

The full one-neighbour pattern is

```math
N^{\mathrm{full}} =
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
1 & \textcolor{blue}{1} & 0 & 0 & 0 & 0 \\
\textcolor{blue}{1} & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & \textcolor{blue}{1} & 0 & 0 & 0 \\
0 & \textcolor{blue}{1} & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & \textcolor{blue}{1} & 0 & 0 \\
0 & 0 & \textcolor{blue}{1} & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & \textcolor{blue}{1} & 0 \\
0 & 0 & 0 & \textcolor{blue}{1} & 1 & 0 \\
0 & 0 & 0 & 0 & 1 & \textcolor{blue}{1} \\
0 & 0 & 0 & 0 & \textcolor{blue}{1} & 1 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}.
```

The complete candidate set is

```math
\begin{aligned}
\mathcal{C} = \{&
(2,2), (3,1), (4,3), (5,2), (6,4), (7,3), \\
&
(8,5), (9,4), (10,6), (11,5)\},
\end{aligned}
```

which depends only on the nonzero structure of $`A`$ and the aggregate assignment.

#### 1.2.4 2D scalar case example

I used a public case（the rail371 example with 371 degrees of freedom） from [Zenodo](https://doi.org/10.5281/zenodo.5113560). I read the matrices $`S`$ and $`M_\Gamma`$ and form

```math
A=S+M_\Gamma.
```

Here $`S`$ is the diffusion stiffness matrix and $`M_\Gamma`$ contains the boundary heat-transfer contribution. 

I build the aggregates from the nonzero connections of $`A`$, then visit the
nodes in increasing index order. If a node has not been assigned,
start a new aggregate with it and add up to three unassigned neighbors
with the largest absolute entries in its row. This gives 122 aggregates, each containing at most 4 nodes.

As in the 1D example, the tentative prolongator has one 1 in each row:

```math
\widehat P_{ij}=
\begin{cases}
1,&j=a(i),\\
0,&\text{otherwise},
\end{cases}
\qquad
N^{(0)}=\operatorname{pattern}(\widehat P).
```

The model has no Dirichlet boundary conditions. I set $`f_i=1`$ for every row, so all rows of $`P`$ must preserve the row sum 1. 

### 1.3 The block cases

t.b.c.

## 2 The Problem I Need to Solve

I need to find a continous (here “continuous” means a sequence of nested discrete patterns) path from unsmoothed to smoothed aggregation by defining a priori interpolation patterns

```math
N^{(0)}\subseteq N^{(1)}\subseteq\cdots\subseteq N^{(K)}.
```

The starting pattern is the pattern of $`\widehat{P}`$. The final pattern should allow smoothed aggregation. In between, we allow some new entries earlier and some later.

This is how I currently separate the tasks:

- Design a rule to choose where nonzeros are allowed and in which order the candidates enter, i.e. how to decide $`\{N^{(k)}\}`$.
- Given a pattern $`N^{(k)}`$, get $`P^{(k)}`$.

### 2.1 How to decide $`\{N^{(k)}\}`$

#### 2.1.1 Coupling strength
Let $`S_i`$ be the set of strong couplings. Then

```math
j\in S_i \Longleftrightarrow -a_{ij} \ge \theta\max_{k\neq i}(-a_{ik}),\quad j\neq i, a_{ij}< 0.
```

Let

```math
c_{ij} = \frac{-a_{ij}}{\max_{k\neq i}(-a_{ik})},\quad j\neq i.
```

Then $`c_{ij}\ge\theta`$ means the connection between fine grid points $`i`$ and $`j`$ is strong. 

Now consider the connection between fine grid points $`i`$ and aggregate $`\mathcal{A}_j`$. I think a natural choice of criterion could be to check if $`t_{ij}:=\max_{k\in\mathcal{A}_j}c_{ik}\ge\theta`$. That means if $`i`$ is strongly connected to at least one node of this aggregates, then the corresponding candidate is allowed to enter. Then lower the threshold gradually to 0, we can get a nested patterns.

In practice we don't need to lower the threshold. We sort all candidate positions by decreasing $`t_{ij}`$ and add them one at a time. Note that only negative off-diagonal entries contribute to the strength. All other entries are assigned zero strength. If a row has no negative off-diagonal entries, all its strengths are set to zero.

#### 2.1.2 Connection score

The coupling strength only looks at single entries of $`A`$. I want a score which looks at the energy $`J`$ directly and everything should be computed from the initial prolongator $`\widehat P`$ .

Recall

```math
J(\widehat P)=\frac12\operatorname{tr}(\widehat P^TA\widehat P).
```

Initially, row $`i`$ has value 1 in column $`c`$ and 0 elsewhere. Once position $`(i,j)`$ is allowed, an amount $`\delta`$ will be transferred to it while keeping the row sum = 1:

```math
P(i,j)=\delta,\qquad P(i,c)=1-\delta.
```

Here $`c`$ = `agg(i)`, the original corase column of $`i`$. Take one column $`p`$ of $`\widehat P`$. Adding $`\delta`$ to its $`i`$-th coefficient gives $`p+\delta e_i`$, where $`e_i`$ has a single 1 in its $`i`$-th entry. Expanding the entire expression directly gives

```math
\frac12(p+\delta e_i)^TA(p+\delta e_i)
= \frac12p^TAp+\frac{\delta}{2}p^TAe_i+\frac{\delta}{2}e_i^TAp+\frac{\delta^2}{2}e_i^TAe_i.
```

Since $`A`$ is symmetric, the two middle terms are equal. Also, $`e_i^TAe_i=A_{ii}`$, so the change is

```math
\delta(Ap)_i+\frac12A_{ii}\delta^2.
```

Denote the gradient of $`J`$ at $`\widehat P`$ as $`G=A\widehat P`$. Now calculate the changes in both columns together.

- Adding $`\delta`$ to column $`j`$ contributes

  ```math
  \delta G(i,j)+\frac12A_{ii}\delta^2.
  ```

- Subtracting $`\delta`$ from column $`c`$ contributes

  ```math
  -\delta G(i,c)+\frac12A_{ii}\delta^2.
  ```

Adding these contributions gives

```math
\boxed{
\Delta J(\delta)
=\delta[G(i,j)-G(i,c)]+A_{ii}\delta^2.
}
```

Taking derivative w.r.t. $`\delta`$ and set it to 0, the minimizer is used as the connection score:

```math
\boxed{s^{\mathrm{conn}}_{ij}:=\delta^*_{ij}
=\frac{G(i,c)-G(i,j)}{2A_{ii}}.}
```

In the 1D scalar example, every candidate row has exactly one edge leaving its original aggregate. Thus

```math
G(i,j)=-w, \qquad G(i,c)=w,
```

Which leads to

```math
s^{\mathrm{conn}}_{ij}
=\frac{w-(-w)}{2A_{ii}}
=\frac{w}{A_{ii}}
=\frac{|G(i,j)|}{A_{ii}}.
```

#### 2.1.3 Energy score

Recall

```math
\Delta J(\delta)
=\delta[G(i,j)-G(i,c)]+A_{ii}\delta^2.
```

Define

```math
\boxed{
s^{\mathrm{en}}_{ij}
=-\Delta J(\delta^*)
=\frac{[G(i,c)-G(i,j)]^2}{4A_{ii}},
}
```

which describes the largest energy decrease obtained by adjusting the coefficients in columns $`c`$ and $`j`$ of row $`i`$.

### 2.2 How to decide $`N^{(k)} \mapsto P^{(k)}`$

t.b.c.