#import "47851.typ": *
#show: lecture.with([Assignment 2], 7)
#set heading(numbering: "Q1.")

= Eigenvalues and SDPs (15 points)

*Part A.* Let $a_i tilde cal(N)(mu_i, Sigma_i)$ be multivariate normal random vectors. Recast the following problem as an SDP assuming $1\/2 <= eta < 1$:
$ max_x c^T x quad "subject to" quad Pr(a_i^T x <= b_i) >= eta quad (i = 1, dots, m). $

*Part B.* Recast the following problem as an SDP:
$ min_X (lambda_(max)(X)) / (lambda_(min)(X)) quad "subject to" quad ⟨A_i, X⟩ = b_i quad (i = 1, dots, m), quad X succ.eq 0. $

*Part C.* Recast the following problem as an SDP:
$ min_(X in SS^n) sum_(j=1)^k lambda_j (X) quad "subject to" quad ⟨A_i, X⟩ = b_i quad (i = 1, dots, m), $
where $lambda_1 (X) >= dots >= lambda_n (X)$ are the eigenvalues.

= Lovász's bound on Shannon capacity (15 points)

Let $alpha(G)$ be the maximum size of an independent set. Recall
$ Theta(G) = lim_(k -> oo) root(k, alpha(G^k)) = sup_(k >= 1) root(k, alpha(G^k)).
$

*Part A.* Let $J$ be the all-ones matrix and define
$ nu (G) = max_X {tr(J X) : tr(X) = 1, space X succ.eq 0, space X_(i j) = 0 "if" {i, j} in E}. $
Show that $nu (G) >= alpha(G)$.

*Part B.* Show that $alpha(G_A times.square G_B) >= alpha(G_A) alpha(G_B)$. The strong product has vertex set $V_A times V_B$. Distinct vertices $(i, k), (j, ell)$ are adjacent exactly when
$
  ({i, j} in E_A "or" i = j) quad "and" quad ({k, ell} in E_B "or" k = ell).
$

*Part C.* Show that $X, Y succ.eq 0$ implies $X times.o Y succ.eq 0$. For $A in RR^(m times n)$ and $B in RR^(p times q)$, the Kronecker product is
$
  A times.o B = mat(
    delim: "[",
    A_(11) B, A_(12) B, dots.c, A_(1n) B;
    dots.v, dots.v, dots.down, dots.v;
    A_(m 1) B, A_(m 2) B, dots.c, A_(m n) B;
  ).
$

*Part D.* Show that $nu (G_A times.square G_B) <= J(G_A) J(G_B)$.

_Hint:_ Use the SDP dual.

*Part E.* Show that $nu (G) >= Theta(G)$ and $sqrt(5) >= nu(C_5) >= Theta(C_5) >= sqrt(5)$.
