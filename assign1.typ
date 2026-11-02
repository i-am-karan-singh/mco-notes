#import "47851.typ": *
#show: lecture.with([Assignment 1], 6)
#counter(heading).update(0)
#set heading(numbering: "Q1.")

= Inner Products and Dual Norms (10 points)

*Part A.* Prove that a norm on a real vector space $V$ is induced by an inner product if and only if it satisfies the parallelogram identity
$ norm(x + y)^2 + norm(x - y)^2 = 2 norm(x)^2 + 2 norm(y)^2 quad forall x, y in V. $

*Part B.* For a norm $norm(dot)$ on $RR^n$, define its dual norm as
$ norm(x)_* = max_(norm(y) <= 1) y^top x. $
+ Calculate the dual norms of $norm(dot)_2$ and $norm(dot)_oo$, where $norm(x)_2 = sqrt(x^top x)$ and $norm(x)_oo = max_i abs(x_i)$.
+ Show that $ partial norm(dot)(x) = {y : norm(y)_* <= 1, space y^top x = norm(x)}. $

= Local Descent and Checking PSD (10 points)

*Part A.* Find a function $f: RR^n -> RR$ and a point $x in RR^n$ such that there is no (non-zero) descent direction at $x$, and yet $x$ is not a local minimizer.

*Part B.* Given $x in RR^n$, form a $(n+1) times (n+1)$-dimensional matrix
$
  X = mat(
    delim: "[",
    1, x_1, x_2, dots.c, x_n;
    x_1, x_1, 0, dots.c, 0;
    x_2, 0, x_2, dots.c, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x_n, 0, 0, dots.c, x_n;
  ).
$
Show that $X succ.eq 0$ if and only if $x_i >= 0, forall i in [n]$ and $sum_(i=1)^n x_i <= 1$.

= Quasiconvexity and Back (10 points)

*Part A.* Prove that $f(x) = log(sum_(i=1)^n a_i e^(x_i))$ where $a_i > 0$ is convex.

*Part B.* Call a function _kinda-convex_ if all its sublevel sets are convex.
+ Give a function that is kinda-convex but not convex.
+ Prove that, for any convex $f: RR -> RR$ and nondecreasing $g: RR -> RR$, $f(g(x))$ is kinda-convex.
+ Can such functions be made convex by stretching the $x$-axis? More precisely, given a kinda-convex $h: RR -> RR$, does there exist a differentiable, invertible $g: RR -> RR$ such that $h(g(x))$ is convex? Prove your claim.

*Challenge (Zero Credits):* Try to generalize (A) to $f(X_1,dots, X_n) = log tr (sum_(i=1)^n a_i e^(X_i))$. Matrix exponential simply exponentiates the eigenvalues in the spectral decomposition.
