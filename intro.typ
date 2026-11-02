#import "47851.typ": *
#show: lecture.with([Introduction and Complexity of Local Minimization], 1)

This is a fast-paced course that develops the basic toolkit of convex optimization. We will see what makes convexity an attractive lens to examine optimization problems, and develop efficient algorithms. A good part of the course will also be spent on (conceptual) applications in applied mathematics and engineering. In roughly five weeks, we will visit:

+ Convex sets and functions
+ Duality: Conic, Lagrangian, Fenchel
+ Structured conic optimization: theory, models, applications
+ Large-scale optimization: ellipsoid, first-order methods, acceleration, lower bounds

#outline()

= A Motivating Example

Consider a language or a script with $n$ alphabets.
Some letter pairs look alike.
If one were to communicate using alphabets alone, one must avoid the use of at least one alphabet in all such letter pairs.
How many letters are usable?

#align(center)[#image("assets/stable-set.svg")]

Consider a graph $G = (V, E)$ where $V$ is the set of letters and $(u,v) in E$ if and only if $u$ and $v$ look alike; this notion is symmetric. We are looking for the largest _stable set_.

$ alpha(G) = max_(S subset.eq V) { |S| : S "contains at most one endpoint of every edge" }. $

Other related graph-theoretic notions:

- Clique is a subset of vertices that are all pairwise connected.
- Vertex Cover is a subset of vertices such that for every edge, at least one endpoint is in it.

#exercise[
  Prove that $S$ is an independent set on $G$ iff $S$ is a clique on the edge-complement of $G$.
]

#exercise[
  Prove that $S$ is an independent set on $G$ iff $V - S$ is a vertex cover on $G$.
]

Another interesting fact is that the number of maximum (in fact, maximal too) independent sets is $3^(n/3) << 2^n$. Can you think of graphs for which this is attained?

A question of Claude Shannon (1956) concerns "worst-case information theory": How much more efficient are words than letters? If one were to communicate with $k$-letter words, $|V| = n^k$, and 2 words $(alpha_1 dots alpha_k)$, $(beta_1, dots beta_k)$ would be ambiguous if and only if they are not the same, yet for all positions $1 <= ell <= k$, $alpha_ell$ and $beta_ell$ are ambiguous (or the same). Drawing edges between such words, we get $G^k$. Graph theorists call it the "strong graph product," denoted by $times.square$, applied $k$ times.

#exercise[
  How do the adjacency matrices of $G$ and $G^k$ relate? Think of Kronecker products.
]

We expect $alpha(G^k)$ to scale exponentially in $k$; after all the vertex set does. So to normalize this effect, let us define
$ Theta(G) = lim_(k -> infinity) root(k, alpha(G^k)). $

#exercise[
  Prove that $alpha(G^k) >= (alpha(G))^k$.

  _Hint:_ Using unambiguous letters alone one can construct unambiguous words.
]

But in fact, $Theta(G)$ can be quite a bit greater than $alpha(G)$, indicating words are more communication-efficient than letters. Here was the concrete question Shannon posed.

$C_n$ is the cycle graph with $n$ vertices. Clearly $alpha(C_5) = 2$. How large is $Theta(C_5)$?

This was resolved by Laszlo Lovasz, in 1979, 23 years later in a result that is one of many origins of Semi-Definite Programming; this we will study, and you'll prove the result in one of the assignments!

#theorem(title: "Lovasz '79")[Using SDPs, we have
  $ Theta(C_5) = sqrt(5) > 2. $
]

An unresolved question even in the LLM age is: what is $Theta(C_7)$?

= Preliminaries and Conventions

The following are notions of _similarity_ and _length_ respectively.

#definition[
  $angle(dot, dot): V times V -> RR$, where $V$ is a vector space, is an inner product if
  + $angle(x,x) >= 0$ for all $x in V$ and $angle(x,x) = 0$ if and only if $x=0$.
  + $angle(x,y) = angle(y,x)$ for all $x,y in V$.
  + $angle(alpha x,y) = alpha angle(x,y)$ for all $x,y in V$ and $alpha in RR$.
  + $angle(x,y+z) = angle(x,y) + angle(x,z)$ for all $x,y,z in V$.
]

#definition[
  $norm(dot): V -> RR$, for vector space $V$, is a norm if
  + $norm(x) >= 0$ for all $x in V$ and $norm(x) = 0$ if and only if $x = 0$.
  + $norm(alpha x) = |alpha| norm(x)$ for all $x in V$ and $alpha in RR$.
  + $norm(x+y) <= norm(x) + norm(y)$ for all $x,y in V$.
]

If the second part of condition (1) fails in either definition, we call these semi inner products and semi norms respectively. Under such definitions, non-zero vectors can have zero _length_.

#lemma(title: "Cauchy-Schwarz")[
  For any (semi) inner product, for all $x,y in V$
  $ angle(x,y)^2 <= angle(x,x) angle(y,y). $
]

#proof[
  Consider that $norm(s x - t y)^2 >= 0$ for all $x,y in V$. We will fix $s,t$ later.
  $ s^2 angle(x,x) + t^2 angle(y,y) - 2 s t angle(x,y) >= 0 $
  $ => (s^2 angle(x,x) + t^2 angle(y,y)) / 2 >= s t angle(x,y). $
  Substitute $s=1/sqrt(angle(x,x))$, $t=1/sqrt(angle(y,y))$, then we are done unless $angle(x,x)$ or $angle(y,y)$ is zero.
]

#exercise[
  Complete the above proof for cases when $angle(x,x)$ or $angle(y,y)$ is zero.

  _Hint:_ Think about what happens if we choose arbitrarily large $s$ or $t$ in these cases.
]

#proposition[
  $f(x) = sqrt(angle(x,x))$ is a norm.
]

#proof[
  (1) & (2) are straightforward. For all $x,y in V$, we have
  $ f(x+y)^2 &= angle(x,x) + angle(y,y) + 2 angle(x,y) \
    &<= angle(x,x) + angle(y,y) + 2 f(x) f(y) quad ("C.S.") \
    &<= (f(x) + f(y))^2. $

]

== When do minimizers exist?

For practical optimization problems, a desirable property is that there is indeed a best decision, no matter how difficult or easy it is to find.

#definition[
  A set $S$ is closed if for any (convergent) sequence contained within it, the corresponding limit also lies in it.

  A set $S$ is bounded if $exists C in RR$ such that $forall x,y in S$, $norm(x-y) <= C$.

  Compact = bounded and closed.
]

A generic optimization problem looks like
$ cal(O)(f,K) equiv min_(x in K) f(x). $

#definition[
  $cal(O)(f,K)$ is feasible if $K$ is not empty.

  $cal(O)(f,K)$ is bounded below if $exists C in RR$, such that $forall x in K, f(x) >= C$.

  $cal(O)(f,K)$ is solvable if $exists x^* in K$ such that $forall x in K, f(x) >= f(x^*)$.
]

For linear programming:
$ "feasible & bounded below" => "solvable." $
But *false* for Convex programs: $min_(x in RR) e^(-x)$.

This is the cause of a lot of pain when stating any non-LP duality result.

#theorem[
  If $K$ is compact & non-empty and $f$ is continuous, then a minimizer exists.
]

This theorem is nice, but we would also like something that applies to unconstrained problems.

#theorem[
  If $K$ is closed & non-empty, $f$ is continuous with all its sub-level sets bounded, then a minimizer exists.
]

#definition[
  $alpha$-sublevel set of $f$ on $K$ for any $alpha in RR$ is $S_alpha (f) = {x in K : f(x) <= alpha}.$
]

= Convexity of Sets and Functions

Consider $V$ to be any vector space and $K$ to be any convex set so that adding elements makes sense.

#definition(title: "Convex Set")[
  $S subset.eq V$ is convex if $forall x,y in S$, $lambda in [0,1]$, $lambda x + (1-lambda)y in S$.
]

#definition(title: "Convex Function")[
  $f: K -> RR$ is convex if $forall x,y in K$, $lambda in [0,1]$,
  $ f(lambda x + (1-lambda)y) <= lambda f(x) + (1-lambda) f(y). $
]

The set ${lambda x + (1-lambda) y: lambda in [0,1]}$ is the line segment connecting $x$ and $y$. So convex sets are precisely those that contain all line segments for which corresponding endpoints are within the set. And convex functions are precisely those for which any linear interpolation overapproximates the function value. How are these related?

#definition(title: "Epigraph")[
  $op("epi")(f) = {(t,x): t >= f(x)}$.
]

#align(center)[#image("assets/epigraph.svg")]

#theorem[
  $f$ is a convex function if and only if $op("epi")(f)$ is a convex set.
]
#proof[
  $=>$ Take $(t,x), (s,y) in op("epi")(f)$. Then
  $ lambda (t,x) + (1-lambda)(s,y) = (lambda t + (1-lambda)s, lambda x + (1-lambda)y) in op("epi")(f) $
  since
  $ lambda t + (1-lambda)s &>= lambda f(x) + (1-lambda)f(y) quad "by membership" \
  &>= f(lambda x + (1-lambda)y). $

  $arrow.l.double$ Take $x,y in K$, $lambda in [0,1]$.
  Note that $(f(x),x), (f(y),y) in op("epi")(f)$.
  Hence $ (lambda f(x) + (1-lambda)f(y), lambda x + (1-lambda)y) in op("epi")(f). $
  Thus $lambda f(x) + (1-lambda)f(y) >= f(lambda x + (1-lambda)y)$.
]

One attraction of convexity is that it is a robust notion, that is, it is preserved under many (useful) operations. We will describe a few operations for convex sets below. You should be able to translate these to convex functions via the above link.

#exercise[
  Prove that Cartesian products, Minkowski sums, linear transformations and intersections of convex sets are convex.
]

If you have not seen this before, $A+B = {a+b: a in A, b in B}$ is a Minkowski sum.

#exercise[
  Construct a closed set $S subset.eq RR^n$, and a matrix $A in RR^(m times n)$ such that the linear transformation $A S = {A x: x in S}$ is not closed.

  _Hint:_ Linear transformations can be used to do projections onto subspaces.
]

#exercise[
  Prove $S_alpha (f)$ is a convex set $forall alpha in RR$ if $f$ is a convex function. Is the converse true?
]

= Local and Global Minimizers

#definition[
  $overline(x)$ is a local minimum for $f$ on $K$ if $exists epsilon > 0$ such that
  $ f(overline(x)) <= f(x) quad forall x in K inter B(overline(x),epsilon). $

  // Typo: “for x” in the next three definitions is “for f”.
  $overline(x)$ is a strict local minimum for $f$ on $K$ if $exists epsilon > 0$ such that
  $ f(overline(x)) < f(x) quad forall x in K inter B(overline(x),epsilon) - {overline(x)}. $

  $overline(x)$ is a global minimum for $f$ on $K$ if $forall x in K$
  $ f(overline(x)) <= f(x). $

  $overline(x)$ is a strict global minimum for $f$ in $K$ if
  $ forall x in K - {overline(x)}, quad f(overline(x)) < f(x). $
]

Here and elsewhere, $B(overline(x),epsilon) = {x: norm(x-overline(x)) <= epsilon}$.
Think of $norm(dot)$ as the Euclidean norm, but for finite-dimensional vector spaces, all norms are equivalent (sandwich each other up to constants), so it does not matter. These definitions apply to any (for example, non-convex) function. Also note that global minima are local minima.

#exercise[
  Give a non-constant function over $K = RR$ all of whose local minima are global _maxima_.
]

#proposition[
  For a convex function $f$, any local minimum is also a global minimum.
]
#proof[
  Take a local minimum $overline(x) in K$. If $overline(x)$ is not a global minimum, then $exists x^* in K$ such that $f(x^*) < f(overline(x))$.
  Note $x^*$ need not be a global minimum; a priori it is not clear if the latter exists.
  Now set $x = (1-lambda)overline(x) + lambda x^* = overline(x) + lambda (x^* - overline(x)) in K$;
  we will ensure $lambda in (0,1]$. We have
  $ f(x) <= (1-lambda)f(overline(x)) + lambda f(x^*) < f(overline(x)), $
  yet $x$ can be made arbitrarily close to $overline(x)$.
  For example, if we want the distance to be at most $epsilon$, we can set $lambda = min{1, epsilon/norm(overline(x)-x^*)}$.
]

#align(center)[#image("assets/local-global.svg")]

You might expect that finding a local minimum is somewhat easy: given that it imposes only local qualifications, maybe some greedy procedure produces these. If so, for convex functions, the above proposition translates them to global minima. This is a key attraction of convexity: doing well enough locally ensures global performance. But this intuition is misleading, as we show next.

= Computational Hardness of Local Minimization

#theorem[
  Given $overline(x) in RR^d$, checking if $overline(x)$ is a local minimum for a degree 4 polynomial is NP-hard.
]

#theorem[
  Given $overline(x) in Delta$, where $Delta = {x in RR^n: x >= 0, sum_(i=1)^n x_i <= 1}$, checking if $overline(x) in Delta$ is a local (constrained) minimum for a degree 2 polynomial is NP-hard.
]

We will prove a weaker result.

#theorem(title: "STRICT-4-HARD")[
  Checking if $x in RR^n$ is a local strict minimum for a degree 4 polynomial is NP-hard.
]

Note that these results concern certification of a given point. It can still be that other local minima are easy to find; this does happen for a few problem classes. But also finding something that is guaranteed to exist (via compactness, say) cannot be NP-hard, which concerns decision problems, not search. These problems are better translated to PPAD-hardness.

== Global problems can be difficult

We will build towards a proof of the last theorem. The easy bit is to believe that global minimization is hard; so let us start there.

#definition[
  *3SAT* is NP-complete and asks if $exists$ a boolean assignment of $x,y,z dots$ such that
  $ (x or overline(y) or z) and (overline(x) or y or overline(z)) and dots $
  is satisfied. Three literals in each clause are joined by ORs with negations allowed. Clauses themselves are joined by ANDs.

  *1-in-3SAT* is also NP-complete. It asks if there is a boolean assignment that makes exactly one of the ORs in each clause true.

  *$d$-POLYPOS* asks if given a degree-$d$ multivariate polynomial ($d$ is a fixed number, so no representational issues here), is it true that $forall x in RR^n, P(x) > 0$.
]

#claim[
  6-POLYPOS is NP-hard.
]
#proof[
  We will embed 3SAT in a degree 6-poly.
  $ P(x) &= sum_(i=1)^n (x_i (1-x_i))^2 \
    &+ underbrace(((x_1+(1-x_2)+x_3-1)(x_1+(1-x_2)+x_3-2)(x_1+(1-x_2)+x_3-3))^2, x or overline(y) or z) + dots $
  There is a degree $6$ term for each OR clause.

  $P(x) = 0$ if and only if $x_i in {0,1}$ and for each clause, at least one (signed) literal is true, e.g.,
  $ x_1 + (1-x_2) + x_3 in {1, 2, 3} $
  which is the same as 3SAT being satisfied.
]

#claim[
  4-POLYPOS is NP-hard.
]
#proof[
  Same proof, except we embed 3-in-1-SAT.
  This time we only need $(x_1+(1-x_2)+x_3-1)^2$.
]

We are butchering complexity-theoretic reductions a bit. We are using Cook/algorithmic reductions, whereas Karp/logical reductions are more standard. But Cook-NP-hard is plenty hard too.

== Embedding global hardness locally

For any degree-$d$ polynomial $P$, let $P_h (x,y) := y^d P(x/y)$ be its homogeneous counterpart. Note that $P_h (x,1) = P(x)$. Also note that this transformation is symbolic; so we can plug in $y=0$.

*Observation:* For any homogeneous polynomial $P_h$, $forall alpha >= 0$, $P_h (alpha x) = alpha^d P_h (x)$, where $d$ is the degree.

#claim[
  For a homogeneous poly $q$, $x=0$ is a strict local min $<=> q(x)>0$ for all $x != 0$.
]
#proof[
  $arrow.l.double$ $q(0)=0$ always. We are done.

  // Typo: exclude y=0 from the strict inequality.
  $=>$ $exists epsilon > 0$ such $forall y in B(0,epsilon) - {0}$, $q(y)>q(0)=0$.
  But then consider any $z in RR^n$, $z != 0$.
  $ q(z) = q(norm(z)/epsilon times z/norm(z) epsilon)
    = (norm(z)/epsilon)^d q(z/norm(z) epsilon) > 0 $
  because $z/norm(z) epsilon$ belongs to $B(0,epsilon)$.
]

The reason for considering homogenized polynomials is as follows: they display the same quantitative behavior in a neighborhood of the origin as they display elsewhere, up to suitable scaling. This characteristic can be used to hide global hardness in a small near-zero neighborhood; the latter is what local minima care about.

Finally, we know that testing positivity of degree $4$ polynomials is NP-hard. If this were also established for homogeneous equivalents, we would be done. Except this doesn't make sense directly because homogeneous polynomials evaluate to zero at zero.
To complete this translation, consider the following assertion for even degree polynomials:
$ P(x)>0 quad forall x quad <=> quad P_h (x,y)>0 quad forall (x,y) != 0. $

To be clear, we will fail to prove this as stated, but we will establish this for special cases which will nevertheless suffice for us.

Here, $arrow.l.double$ holds always as $P(x)=P_h (x,1)>0$.

$=>$ is tricky. Let's try the contrapositive.
$ exists (x,y) != 0 "such" P_h (x,y) <= 0. $
We wish to conclude $exists z != 0$, $P(z) <= 0$.
This is indeed fine if $y != 0$ incidentally.
$ P(z:= x/y) = (P_h (x,y))/y^d <= 0. $

But this doesn't work if $y=0$. In fact, consider the counterexample below.

#example[
  $P(x_1,x_2)>0$ everywhere, yet $P_h (1,0,0)=0$, where
  $ P(x_1,x_2) = x_1^2 + (1-x_1 x_2)^2, $
  $ P_h (x_1,x_2,y) = y^2 x_1^2 + (y^2-x_1 x_2)^2. $
]

The subtle observation is that the bad case never happens for the polynomial constructed for 1-in-3SAT.
Indeed
$ P_h (x,y) = sum_(i=1)^n (x_i (y-x_i))^2 + (y(x_1+(y-x_2)+x_3)-y^2)^2 + dots $
With $y=0$, the first term becomes $sum_(i=1)^n x_i^4$.
Hence $P_h (x,0) <= 0 => x_i = 0$ for all $i$.
But the only problematic case was when $exists (x,y) != 0$ and $y=0$ where $P_h (x,y) <= 0$.

So our proof of the hardness of certifying strict local minima is complete.
