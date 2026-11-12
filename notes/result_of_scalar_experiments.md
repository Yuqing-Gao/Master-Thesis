# Result of the 1D and 2D experiments

## 1 The 1D case

### 1.1 The scores of the 10 candidates

| candidate | edge | $w$ | $a_{ii}$ | coupling strength | connection score | energy score | true decrease if added alone|
|---|---|---|---|---|---|---|---|
| $(2,2)$ | 2–3 | 0.25 | 1.25 | 0.25 | 0.2000 | 0.0500 | 0.2500 |
| $(3,1)$ | 2–3 | 0.25 | 1.25 | 0.25 | 0.2000 | 0.0500 | 0.0500 |
| $(4,3)$ | 4–5 | 2 | 3 | 1 | 0.6667 | 1.3333 | 1.3333 |
| $(5,2)$ | 4–5 | 2 | 2.5 | 1 | 0.8000 | 1.6000 | 1.6000 |
| $(6,4)$ | 6–7 | 0.8 | 1.3 | 1 | 0.6154 | 0.4923 | 0.4923 |
| $(7,3)$ | 6–7 | 0.8 | 1.8 | 0.8 | 0.4444 | 0.3556 | 0.3556 |
| $(8,5)$ | 8–9 | 3 | 4 | 1 | 0.7500 | 2.2500 | 2.2500 |
| $(9,4)$ | 8–9 | 3 | 4 | 1 | 0.7500 | 2.2500 | 2.2500 |
| $(10,6)$ | 10–11 | 0.4 | 1.4 | 0.4 | 0.2857 | 0.1143 | 0.1143 |
| $(11,5)$ | 10–11 | 0.4 | 1.4 | 0.4 | 0.2857 | 0.1143 | 0.3674 |

"True decrease" means: minimal energy on $N^{(0)}$ ($6.95$) minus minimal energy on $N^{(0)}$ plus this one position.

The 3 orders that follow the scores:

- coupling: $(4,3),(5,2),(6,4),(8,5),(9,4),(7,3),(10,6),(11,5),(2,2),(3,1)$
- connection: $(5,2),(8,5),(9,4),(4,3),(6,4),(7,3),(10,6),(11,5),(2,2),(3,1)$
- energy: $(8,5),(9,4),(5,2),(4,3),(6,4),(7,3),(10,6),(11,5),(2,2),(3,1)$

For each budget $b$, best_J is the minimum energy over all subsets containing exactly $b$ candidate positions. For each ranking, ranked_J is the minimum energy obtained on the pattern formed by its first $b$ candidates.

| Budget $b$ | Exhaustive `best_J` | Coupling `ranked_J` | Connection `ranked_J` | Energy `ranked_J` |
|---:|---:|---:|---:|---:|
| 0 | 6.950000 | 6.950000 | 6.950000 | 6.950000 |
| 1 | 4.700000 | 5.616667 | 5.350000 | 4.700000 |
| 2 | 3.100000 | 5.235714 | 3.100000 | 4.378571 |
| 3 | 2.305396 | 4.522493 | 2.778571 | 2.778571 |
| 4 | 1.866621 | 2.272493 | 2.664286 | 2.664286 |
| 5 | 1.616621 | 1.951064 | 1.951064 | 1.951064 |
| 6 | 1.376866 | 1.804707 | 1.804707 | 1.804707 |
| 7 | 1.242075 | 1.538709 | 1.538709 | 1.538709 |
| 8 | 1.135305 | 1.326835 | 1.326835 | 1.326835 |
| 9 | 1.076835 | 1.076835 | 1.076835 | 1.076835 |
| 10 | 1.024610 | 1.024610 | 1.024610 | 1.024610 |

### 1.2 Interactions between candidates

In this case
$$
a_{ii} = w_{i-1} + w_i,
\qquad
a_{i,i+1} = a_{i+1,i} = -w_i.
$$

so

$$
\begin{aligned}
(P^T A P)_{jj}
&= \sum_{i=1}^{12} a_{ii}P_{ij}^2
 + \sum_{i=1}^{11} a_{i,i+1}P_{ij}P_{i+1,j}
 + \sum_{i=1}^{11} a_{i+1,i}P_{i+1,j}P_{ij}\\
&= \sum_{i=1}^{12}(w_{i-1}+w_i)P_{ij}^2
 - 2\sum_{i=1}^{11}w_i P_{ij}P_{i+1,j}\\
&= w_0P_{1j}^2
 + \sum_{i=1}^{11}w_i
   (P_{ij}^2+P_{i+1,j}^2-2P_{ij}P_{i+1,j})
 + w_{12}P_{12,j}^2\\
&= w_0P_{1j}^2
 + \sum_{i=1}^{11}w_i(P_{i+1,j}-P_{ij})^2
 + w_{12}P_{12,j}^2.
\end{aligned}
$$

Then
$$
\begin{aligned}
J(P)
&= \frac12\operatorname{tr}(P^TAP)\\
&=\frac12\sum_{j=1}^{6}(P^TAP)_{jj}\\
&=\frac12\sum_{j=1}^{6}
\left[
w_0P_{1j}^2
+\sum_{i=1}^{11}w_i(P_{i+1,j}-P_{ij})^2
+w_{12}P_{12,j}^2
\right].
\end{aligned}
$$
Thus the energy contribution of an edge $i$ to $i+1$ can be written as
$$
E_{i,i+1}(P)=\frac{w_i}{2}
\sum_{j=1}^{6}(P_{i+1,j}-P_{ij})^2,
\qquad i=1,\ldots,11.
\tag{1}
$$

Now look at the candidates $(8,5)$ and $(9,4)$. Let $t$ and $u$ be the amount of transfer after these 2 candidates are open.

| $i$ | $P_{i4}$ |$P_{i5}$ |
| :---: | :---: | :---: |
| 7  | $1$   | $0$   |
| 8  | $1-t$ | $t$   |
| 9  | $u$   | $1-u$ |
| 10 | $0$   | $1$   |

The edges that are influnced are 7-8, 8-9, 9-10. By $(1)$, their energy contributions are $t^2$, $3(1-t-u)^2$, $u^2$. Before those 2 candidates open, the contributions are $0,3,0$. Thus $\Delta J = -6t-6u+4t^2+4u^2+6tu$. Take the partial derivatives and set them to 0, there are
$$
\min\Delta J = - \frac{18}{7} \approx -2.57143.
$$
But only one of $(8,5)$ and $(9,4)$ are open, the energy decreses by 2.2500. So the joint reduction is smaller than the sum of the individual reductions. I think it's because there are a $+6tu >0$. $\Delta J$ is negative, $6tu$ would make the value of reduction smaller.  

Then look at the candidates $(5,2)$ and $(6,4)$. Similarly

| $i$ |$P_{i2}$ |$P_{i3}$ |$P_{i4}$ |
| :---: | :---: | :---: | :---: |
| 4 | $1$ | $0$ | $0$ |
| 5 | $t$ | $1-t$ | $0$ |
| 6 | $0$ | $1-u$ | $u$ |
| 7 | $0$ | $0$ | $1$ |

By $(1)$ we have $E_{4,5}=2(1-t)^2$, $E_{6,7}=0.8(1-u)^2$, $E_{5,6}=\frac{0.5}{2}\left[t^2+(t-u)^2+u^2\right]=0.5(t^2+u^2-tu).$ Thus $\Delta J = -4t-1.6u+2.5t^2+1.3u^2-0.5tu$. Take the partial derivatives and set them to 0, there are
$$
\min\Delta J = - \frac{608}{255} \approx -2.38431.
$$

When they are open individually, the sum of reduction is $1.6000 + 0.4923 \approx 2.0923$. Because there are $-0.5tu$, they would help reduce $J$ together.

### 1.3 Explanation of the tables in 1.1

When budget = 1, energy score chooses $(8,5)$ and gets the best reduction of energy. Because for this interior candidate, it is exactly how much the energy decreases if we add only 1 candidate. So at this budget energy score is the best.

When budget = 2, energy score chooses $(9,4)$. By their interaction it will lower the sum of reduction of energy, which we have discussed in 1.2. Connection score chooses $(5,2)$ and $(8,5)$, which do not interact, so it happens to choose the best combination. The choice of coupling $(4,3)$, $(5,2)$ are also at the 2 ends of 1 edge, so there are also overlap in the reduction.

When budget = 4, Coupling finds a better combination. The positions that Connection and Energy choose are the same under this budget. They are $\{(4,3),(5,2),(8,5),(9,4)\},$ which are at the ends of the 2 edges that have the biggest weights. Both of those pairs make the value of reduction smaller. The set that Coupling chooses is $\{(4,3),(5,2),(8,5),(6,4)\}.$ As we have discussed in 1.2, $(6,4)$ and $(5,2)$ make the value of reduction bigger. Thus it is more effective.

### 1.4 A short conclusion

The interactions between candidates affect the energy reduction. Joint gains can be smaller or larger than the sum of individual gains. Those scores that used in the example evaluate candidates independently and do not account for these interactions.

## 2 The 2D case (rail, $n=371$)

In this example there are too many candidates thus cannot find `best_J` exhaustively. I recorded the result and asked chatgpt to draw a comparison figure, see `experiments/result/rail371_energy_comparison.svg`.

The connection and energy scores perform better at small budgets, but coupling strength attains lower energies from budget 162. I think a possible explanation may be: both connection and energy scores use AP^hat (the initial gradient of energy) to compute the scores and rank the candidates. After more candidates are open, this gradient may be "out of date"? That is, as the pattern grows, candidate interactions change the actual marginal gains, so the initial scores may become less representative of these gains.