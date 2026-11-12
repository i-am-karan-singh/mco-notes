#import "47851.typ": *
#show: lecture.with([Volumetric Algorithms], 4)

With this unit, we begin our study of convex optimization algorithms. The algorithms in this unit work the same way that you might look for apartment keys in a packed suitcase, by iterative discarding large chunks of space that you are sure do contain the keys. These algorithms are not entirely efficient in practice (and some not even in theory), but they are the most natural starting point.

#outline()

= The Center-of-Mass Algorithm

For any compact body $K subset.eq RR^n$, let $upright("vol")(K)=integral_K dif x$ be the volume, $overline(K)=integral_K 1/(upright("vol")(K)) x dif x$ be the center of mass. Another equivalent expression is $overline(K)= EE_(x tilde upright("Unif")(K))[x]$. For a compact, full-dimensional, convex set $K$, consider $min_(x in K) c^top x$; this can represent any convex program. (Why?)

*The Center-of-Mass Algorithm:*
+ Let $K_1=K$.
+ For $t=1 dots T$:
  - Compute $x_t = overline(K_t) = 1/(upright("vol")(K_t)) integral_(K_t) x dif x$.
  - Update $K_(t+1) = K_t inter {x: c^top x <= c^top x_t}$.
+ Output $overline(x)=arg min_(x in {x_1,dots,x_(T+1)}) f(x)$.

#theorem[
  If $max_(x,y in K) c^top (x-y) <=F$, then after $T$ steps, we have
  $ c^top overline(x) <= min_(x in K) c^top x + F (1-1/e)^(T/n). $
  In particular, if $T>=n log(F/epsilon)$, then $overline(x)$ is $epsilon$-optimal.
]

The following is the central result we will rely on. In words, any half-space through the center of mass of a convex body rejects at least $1/e$ fraction of the volume. With this interpretation, the COM algorithm has the same flavor as binary search.

#lemma(title: "Gr\u{fc}nbaum's Lemma")[
	For any convex, compact $K$, with COM $overline(x)$:
  $ forall c in RR^n, upright("vol")(K inter {x:c^top (x-overline(x))<=0}) >= 1/e upright("vol")(K). $
]

#align(center)[#image("assets/com-cuts.svg", width: 65%)]

#proof[
  Let $x^star=arg min_(x in K) f(x)$. Then, take $X_epsilon^star={(1-epsilon)x^star+epsilon x:x in K}=(1-epsilon)x^star+epsilon K$.
  Now, $upright("vol")(X_epsilon^star)=epsilon^n upright("vol")(K)$.
  Also, $X_epsilon^star subset.eq K=K_1$, by convexity of $K$.
  $ max_(x in X_epsilon^star)c^top x <= min_(x in K)c^top x+epsilon F, $
  by construction of $X_epsilon^star$. In words, $X_epsilon^star$ is a small set of points, all with good objective values. We will prove that although initially completely inside $K_1$, some of it must fall outside $K_t$ for large enough $t$. Whenever this first happens, $x_t$ must be better than some $x in X_epsilon^star$. To see this:
  $ upright("vol")(K_(t+1)) &<= (1-1/e) upright("vol")(K_t) & text("(Gr\u{fc}nbaum)") \
  &<= (1-1/e)^t upright("vol")(K_1) $
  Now, set $epsilon>(1-1/e)^(T/n)$. Then, $X_epsilon^star subset.eq K_1$, yet $upright("vol")(K_T)<upright("vol")(X_epsilon^star)$.
  Hence, $exists t in [T]$ and $x_epsilon^star in X_epsilon^star$ such that $x_epsilon^star in K_t$, $x_epsilon^star in.not K_(t+1)$.
  By construction, $c^top x_t<c^top x_epsilon^star<=min_(x in K)c^top x+epsilon F$.
]

== The Cone as an Extremal Example

What we did above is basically the standard way to convert volume-reduction arguments to function value suboptimality. Note that we lose a factor of dimension $n$ in this reduction. In the rest of this note, we will prove Gr\u{fc}nbaum's lemma.

#example[
	Let us say that the $(n-1)$-dimension volume of a sphere of radius $r$ is $C_(n-1) r^(n-1)$. Let $K$ be the right circular cone. Then:
	$
		upright("vol")(K)=integral_0^R C_(n-1) r^(n-1) dif r = C_(n-1)/n R^n\
		overline(K) = n/(C_(n-1) R^n) integral_0^R C_(n-1) r^n dif r = n/(n+1)R.
	$
	From this, we can calculate the volume retained in a cut parallel to the base.
	$
		upright("vol")(K inter {x:x_1<=overline(K)) = integral_0^(n/(n+1)R) C_(n-1)r^(n-1) dif r = C_(n-1)/n (n/(n+1))^n R^n\
		(upright("vol")(K^-))/(upright("vol")(K))=(n/(n+1))^n >= 1/e.
	$
	We will see that the cone (and this cut) is the worst case for Gr\u{fc}nbaum.
]

#align(center)[#image("assets/com-cone.svg", width: 65%)]

Although this is a specific example, we will use it as a proof strategy. We will reduce general convex bodies to right circular cones.

== The Brunn-Minkowski Inequality

#theorem(title: "Brunn-Minkowski Inequality")[
  For non-empty compact sets $A,B subset.eq RR^n$,
  $ upright("vol")(A+B)^(1/n)>=upright("vol")(A)^(1/n)+upright("vol")(B)^(1/n). $
]

#exercise[
	Prove the above statement when both $A$ and $B$ are axis aligned hyper-rectangles.
]

#proof[
  We will extend the above exercise to the case when $A$ and $B$ are unions of disjoint cuboids. By a limiting argument, this extends to compact bodies. Our induction hypothesis is that the stated inequality is true when $A$ and $B$ contain $n$ disjoint cuboids in total. Volume is translation invariant. Hence, shift the $x_1=0$ plane so that at least one cuboid in $A$ lies entirely above it.
  $ A_+ = A inter {x:x_1>=0}, quad A_-=A-A_+, $
  $ B_+ = B inter {x:x_1>=0}, quad B_-=B-B_+. $
  Translate $B$ along $x_1$ so that $(upright("vol")(A^+))/(upright("vol")(A)) = (upright("vol")(B^+))/(upright("vol")(B))$; such translation always exists due to the intermediate value theorem. Also both $A^+$ and $A^-$ have at least one fewer cuboid than $A$. (Why? This is not literally true, but true enough for good choices of $A$ and $x_1$.)

  #align(center)[#image("assets/brunn-minkowski-cuboids.svg", width: 100%)]

  Notice that $(A^++B^+) inter (A^-+B^-)=emptyset$ since $x_1=0$ separates the two sets, yet we also have that $(A^++B^+) union (A^-+B^-) subset.eq A+B$. Hence, using our induction hypothesis, we have
  $
  	upright("vol")(A+B) &>= upright("vol")(A^++B^+)+upright("vol")(A^-+B^-)\
   	&>= (upright("vol")(A^+)^(1/n)+upright("vol")(B^+)^(1/n))^n + (upright("vol")(A^-)^(1/n)+upright("vol")(B^-)^(1/n))^n \
    &= (upright("vol")(A^+)+upright("vol")(A^-)) (1 + ( (upright("vol")(B))/(upright("vol")(A)) )^(1/n))^n \
    &= (upright("vol")(A)^(1/n)+upright("vol")(B)^(1/n))^n.
  $
]

#corollary[
  $upright("vol")(K inter {x_1=alpha})^(1/(n-1))$ is concave in $alpha$ for any compact, convex set $K subset.eq RR^n$. Here, we mean $(n-1)$-dimensional volume, and this quantity is the effective radius.
]

#proof[
  Let $K_alpha=K inter {x_1=alpha}$. Note that $lambda K_alpha+(1-lambda)K_beta subset.eq K_(lambda alpha+(1-lambda)beta)$, $forall alpha,beta in RR$, $lambda in [0,1]$ since $K$ is convex. Then, we have that $ upright("vol")(K_(lambda alpha+(1-lambda)beta))^(1/(n-1)) >= lambda upright("vol")(K_alpha)^(1/(n-1))+(1-lambda)upright("vol")(K_beta)^(1/(n-1)). $
]

== Proof of Gr\u{fc}nbaum's lemma

Now, we are ready to complete Gr\u{fc}nbaum's lemma.

#proof(title: "Proof of Gr\u{fc}nbaum's lemma")[
  Without loss of generality, we can orient our convex body so that $x_1=0$ is the cutting hyperplane.

  Replace every slice of $K$ along $x_1$ axis with a $(n-1)$-dimensional sphere of equal $(n-1)$-dimensional volume. This step preserves volumes of both sections on either side of $x_1=0$; also $x_1$ coordinate of center-of-mass stays the same. So, it suffices to establish the claim for this new body.
  #align(center)[#image("assets/steiner-symmetrization.svg", width: 80%)]

  $ K_+=K inter {x:x_1>=0}, quad K_-=K inter {x:x_1<=0}. $
  First note, this new body is convex, since we didn't modify $upright("vol")(K inter {x_1=alpha})$ and $upright("vol")(K inter {x_1=alpha})^(1/(n-1))$ was concave in $alpha$ for the old, and hence, the new, body.

  Replace $K^+$ with a cone with the same spherical base as $K^+$, so that the cone and $K^+$ are equi-volume. Extend this cone in the negative $x_1$-region till this extension has volume equal to $upright("vol")(K^-)$, again always possible by the intermediate value theorem.

  These operations are volume preserving, but what happens to the center of mass? The claim is that it can only move rightwards. In other words, this transformation increases the $x_1$ coordinate of center of mass from $0$ to something non-negative. This is once again a consequence of concavity of $upright("vol")(K inter {x_1=alpha})^(1/(n-1))$ in $alpha$. Post this transformation, we have a perfect cone with non-negative $x_1$ coordinate of COM. Hence, taking $x_1^upright("COM")>=0$ to be $x_1$-coordinate of COM, we get
  $ (upright("vol")(K_+))/(upright("vol")(K)) = (upright("vol")(K inter {x_1>=0}))/(upright("vol")(K)) >= (upright("vol")(K inter {x_1>=x_1^upright("COM")}))/(upright("vol")(K)) >= 1/e. $
]

#align(center)[#image("assets/com-transport.svg", width: 65%)]

== The Trouble with the Center-of-Mass Algorithm

Computing the COM of a general convex body given a polyhedral description is \#P-hard. Nevertheless it is possible to get an $epsilon$-far point to COM in polynomial (in $1/epsilon,m,n$) time. This absolves the COM algorithm somewhat. But even the $epsilon$-approximation algorithm is almost as tedious as solving an LP in the first place. We will however see a better algorithm next, in fact the first poly-time algorithm for solving LPs, inspired by the COM algorithm.

= The Ellipsoid Method

The ellipsoid method can be thought of as a variant of the center-of-mass method, but one implementable in polynomial time. It was the first provably poly-time algorithm for LPs.

*Assumption:* $K$ is convex, compact and $r B_2 subset.eq K subset.eq R B_2$ where $B_2={x:norm(x)_2<=1}$.

#align(center)[#image("assets/ellipsoid-bounds.svg", width: 30%)]

*Ellipsoid Algorithm:*
+ Initialize $cal(E)_1=R B_2$.
+ For $t=1 dots T$:
  - Let $x_t$ be the center of $cal(E)_t$.
  - Is $x_t in K$? (Membership oracle)
  - If yes, construct an ellipse $cal(E)_(t+1)$ containing $cal(E)_t inter {x:c^top (x-x_t)<=0}$.
  - If no, ask for a halfspace $w_t$ such that $forall x in K$, $w_t^top (x-x_t)>=0$; thus all of $K$ is contained in it. Construct an ellipse $cal(E)_(t+1)$ containing $cal(E)_t inter {x:w_t^top (x-x_t)>=0}$. (Separation oracle)

#align(center)[#image("assets/ellipsoid-oracles.svg", width: 80%)]

If $K={x:A x<=b}$, then $x_t in K$ can be answered in linear time by checking all constraints. If $x_t in.not K$, then $exists i in [m]$ such $a_i^top x_t>b_i$. But $forall x in K,a_i^top x<=b_i$ implying $a_i^top (x-x_t)<=0, forall x in K$. This gives us the required separating hyperplane. This also works for functional constraints of the form $f(x)<=0$, where $f$ is a convex function. Because if $f(x_t)>0$, then for any feasible $x$, we have $ nabla f(x_t)^top (x_t -x)>=f(x_t)-f(x) > 0 $
which is the required separating hyperplane.

We've demonstrated that for LPs with polynomially many constraints, membership and separation oracles can be implemented efficiently. However, this is not the only case when this is possible. For certain structured LPs with exponentially/infinitely many constraints, ellipsoid is still a poly-time algorithm as long as separation/membership queries are efficiently answerable.

#lemma[
  For any ellipsoid $cal(E)_0$ with center $x_0 in RR^n$, and vector $w_0$, we can efficiently construct an ellipsoid $cal(E)_1$ containing $cal(E)_0 inter {x:w_0^top (x-x_0)>=0}$ satisfying
  $ upright("vol")(cal(E)_1)<=upright("vol")(cal(E)_0)e^(-1/(2(n+1))). $
]

#theorem[
  Let $overline(x)$ be a feasible point in ${x_1 dots x_T}$ achieving the smallest objective value. Then
  $ c^top overline(x)-min_(x in K)c^top x <= (F R)/r e^(-T/(2(n+1)n)), text("where") F=max_(x,y in K)c^top (x-y). $
  As long as $T>=2n(n+1)log((F R)/(r epsilon))$, $overline(x)$ is $epsilon$-optimal.
]

Notice that this is slower than COM by a factor of $n$, because volume reduction is $approx 1-1/n$ in each step, instead of a constant. This is the price ellipsoid method pays for efficient implementability.

#proof[
  We will follow the same recipe as that for COM. Let $x^star in arg min_(x in K)c^top x$, and again define $X_epsilon^star=(1-epsilon)x^star+epsilon K subset.eq K$.
  Now, we have $upright("vol")(X_epsilon^star)>=upright("vol")(epsilon r B_2)=C_n (epsilon r)^n$ for some $C_n$ such that $upright("vol")(B_2)=C_n$. Also, we have
  $ max_(x in X_epsilon^star)c^top x <= c^top x^star+max_(x in X_epsilon^star)c^top (x-x^star)<=c^top x^star+epsilon F. $
  By the volume reduction lemma, with repeated applications, we get
  $ upright("vol")(cal(E)_(T+1))<=upright("vol")(cal(E)_1)e^(-T/(2(n+1)))=C_n R^n e^(-T/(2(n+1))). $
  Choose $epsilon>R/r e^(-T/(2n(n+1)))$. Then since $upright("vol")(cal(E)_(T+1))<upright("vol")(X_epsilon^star)$, $exists t$ and $x_epsilon^star in X_epsilon^star$ such that $x_epsilon^star in cal(E)_t$, yet $x_epsilon^star in.not cal(E)_(t+1)$. Further, note that this can only happen on the YES branch, because all existing feasible points are retained on the NO branch. Also YES must happen at least once. (Why?) Hence, we get that $c^top x_t<c^top x_epsilon^star<=c^top x^star+epsilon F$.
]

== Volume Reduction with Ellipsoids

Now, we will finish up the proof of the volume reduction lemma. Note that this lemma was crucial for the optimality result. Also, such niceties don't work for some simpler shapes.

#align(center)[#image("assets/no-volume-reduction.svg", width: 80%)]

Smallest circle containing this hemicircle is the original circle itself. So no volume reduction for the circle method, and similarly, no volume reduction for the rectangle method.

#align(center)[#image("assets/ellipsoid-cut.svg", width: 40%)]

#proof[
  Consider a simpler case when we start with the unit ball $B_2={x:norm(x)_2<=1}$; $x_1>=0$ is halfspace. Any ellipse can be written as
  $ norm(x-x_0)_(H_0^(-1))^2=(x-x_0)^top H_0^(-1)(x-x_0)<=1 $
  where $x_0$ is the center with principal axes of lengths $sqrt(lambda_1),dots,sqrt(lambda_n)$ where $lambda_i$'s are the eigenvectors of $H_0$. (This is a generalization of $(x-x_0)^2/a^2+(y-y_0)^2/b^2<=1$.) Such an ellipse has volume $C_n sqrt(lambda_1 dots lambda_n)=C_n sqrt(det(H_0))$, because stretching a body along any axis by a factor two doubles its volume.

  Now, back to $cal(E)_1$. By symmetry, we take $x_1=t e_1$ as the center of $cal(E)_1$. This ellipse touches $e_1$ and $B_2 inter {x_1=0}$. So, our ansatz is $H=a e_1 e_1^top+b(I-e_1 e_1^top)$. This is an eigenvalue decomposition. Also, $H^(-1)=1/a e_1 e_1^top+1/b (I-e_1 e_1^top)$. So using these two extremal characterizations, we get
  $
  	(1-t)^2/a=1 => a=(1-t)^2, quad t^2/a+1/b=1 => b=1/(1-t^2/a)=(1-t)^2/(1-2t).\
   upright("volume")(cal(E)_1)=C_n sqrt(a b^(n-1))=C_n (1-t)^n/(1-2t)^((n-1)/2).
  $
  Minimizing this, for $t=1/(n+1)$, using $1+x<=e^x$, we get
  $ upright("volume")(cal(E)_1) &= C_n sqrt(a b^(n-1))=C_n (n/(n+1))(n^2/(n^2-1))^((n-1)/2) \
  &<= C_n e^(-1/(n+1)) e^((n-1)/(2(n^2-1)))=e^(-1/(2(n+1)))upright("volume")(cal(E)_0 equiv B_2). $

  Note that $(upright("volume")(cal(E)_1))/(upright("volume")(cal(E)_0))$ is invariant under rotations, because all volumes are, but also under stretching of any axes. Thus $(upright("volume")(cal(E)_1))/(upright("volume")(cal(E)_0))$ is invariant under any invertible coordinate transformation, since any invertible linear map $A=U Sigma V^top$ for orthogonal $U,V$ and diagonal $Sigma$ with positive entries by singular value decomposition. Thus our volume reduction result holds starting with any ellipse and half-space.
]

#exercise[Verify that $cal(E)_1 supset.eq cal(E)_0 inter {x_1>=0}$.]

#exercise[
	Another way to simplify this is by upper bounding both terms in
	$
		(1-t)^n/(1-2t)^((n-1)/2) = (n/(n+1))^n ((n+1)/(n-1))^((n-1)/2)
	$
	using $1+x<=e^x$. Explain why this fails despite also being a "first-order" expansion.
]

== Computational Implementation

Finally, although unnecessary for our proof, we have also constructed the smallest ellipse subject to the containment requirement; our upper bounds on its volume might have been a bit loose though. For a computationally explicit implementation, we extend this construction to the general case directly. So $ cal(E)_0={x:norm(x-x_0)_(H_0^(-1))^2<=1}$. We will construct $cal(E)_1 supset.eq cal(E)_0 inter {x:w^top (x-x_0)>=0}$.

Let $x=H_0^(1/2)y+x_0$. In $y$-space, $cal(E)_0={y:norm(y)_2^2<=1}$. Also
$ w^top (x-x_0)>=0 equiv (H_0^(1/2)w)^top y>=0 equiv ((H_0^(1/2)w)^top y)/norm(w)_(H_0)>=0. $

Now $cal(E)_1 = {y: norm(y-1/(n+1) (H_0^(1/2)w)/norm(w)_(H_0))_M^2 <= 1}$, where
$ M=((n+1)/n)^2 (H_0^(1/2)w w^top H_0^(1/2))/norm(w)_(H_0)^2 + (n^2-1)/n^2 (I-(H_0^(1/2)w w^top H_0^(1/2))/norm(w)_(H_0)^2). $

Going back the $x$-space, $cal(E)_1={x:norm(x-x_1)_(H_1^(-1))^2<=1}$, where
$
	x_1=x_0+1/(n+1) (H_0 w)/norm(w)_(H_0), \
	H_1=(((n+1)/n)^2 (w w^top)/norm(w)_(H_0)^2 + (n^2-1)/n^2 (H_0^(-1)-(w w^top)/norm(w)_(H_0)^2))^(-1).
$
