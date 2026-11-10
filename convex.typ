#import "47851.typ": *
#show: lecture.with([Convexity and Separation], 2)

So far, we have seen that an attractive property of convex functions is that local minima are also global minima. However, we have also seen that even local minima can be hard to find and certify generally. So we will begin by relaxing the definition of local minima itself, in a way that retains the local-to-global property while being tractable.

#outline()

= First-Order Conditions for Local Minimization

#definition[
	$d$ is a descent direction of $f$ at $x$ if there exists $overline(alpha) > 0$ such that $f(x+alpha d) < f(x)$ for all $alpha in (0,overline(alpha))$.
]

#claim[
	If $angle(nabla f(x),d) < 0$, then $d$ is a descent direction for $f$ at $x$.
]

#proof[
	By Taylor's theorem, we have
	$ f(x+epsilon d) = f(x) + angle(nabla f(x),d) epsilon + r(epsilon,x,d), $
	where $r(x,epsilon,d) = o(epsilon)$, that is, $lim_(epsilon -> 0^+) r(epsilon,x,d)slash epsilon = 0$.
	There exists $overline(epsilon)>0$ such that for all $0 < epsilon <= overline(epsilon)$, $|r(epsilon,x,d) slash epsilon| <= |angle(nabla f(x),d)| slash 2$. Then, we get
	$ (f(x+epsilon d)-f(x))/epsilon = angle(nabla f(x),d) + r(epsilon,x,d)/epsilon<= 1/2 angle(nabla f(x),d) < 0. $
]

#align(center)[#image("assets/descent-directions.svg")]

The converse is false. Consider $f(x)=-x^2$; then at $0$ any non-zero direction is a descent direction.

#proposition[
	$x$ is a local min for $f$ $=>$ $nabla f(x)=0$.]

#proof[
	Say $nabla f(x)!= 0$; then $angle(nabla f(x),-nabla f(x)) < 0$. Hence $-nabla f(x) slash norm(nabla f(x))$ is a descent direction.

	Hence there exists $overline(alpha)>0$ such that for all $0 < epsilon < overline(alpha) $, there exists $x' in B(x,epsilon)$, namely $ x' = x - epsilon (nabla f(x))/norm(nabla f(x)), $ such that $f(x') < f(x)$. Hence $x$ is not a local min.
]

To summarize all this, we have:
#align(center)[
	#diagram($
		x "is a local min" edge(=>) edge("dr", =>) & edge(=>) & nabla f(x)=0 \
		& exists.not "a descent direction" edge("ur", =>)
	$)
]

#exercise[
	Provide separate counterexamples for each of the three converses.
]

== The Constrained Case

Now, let's do the constrained case. We have dealt with non-convex functions so far. Coming to the constrained case, you might be tempted to also treat general (non-convex) sets. But, for example, if $K$ is the set of vertices on a hypercube, all feasible points are local minima, regardless of the objective function. To avoid such pathologies, it makes sense to only consider a convex set $K$.

#proposition[
	$overline(x)$ is a local min for $f$ on a convex set $K$ $=>$ $angle(nabla f(overline(x)),x-overline(x)) >= 0$ for all $x in K$.
] <necessary-first-order>

#proof[
	Say there exists $x in K$, $angle(nabla f(overline(x)),x-overline(x)) < 0$. Then consider $ y=(1-lambda)overline(x)+lambda x=overline(x)+lambda(x-overline(x)) in K. $
	Now, by Taylor's theorem,
	$f(y)=f(overline(x))+lambda angle(nabla f(overline(x)),x-overline(x))+r(lambda,x,overline(x))$ where $r(lambda,x,overline(x))=o(lambda)$.
	Then there exists $overline(lambda)>0$ such that for all $0<lambda<overline(lambda)$, $|r(lambda,x,overline(x))/lambda| < 1/2 |angle(nabla f(overline(x)),x-overline(x))|$.
	For all such $lambda$'s, we have $ f(y) <= f(overline(x))+lambda/2 angle(nabla f(overline(x)),x-overline(x)) < f(x). $
	Thus, there does not exist an $epsilon>0$ such that $f(y)>=f(overline(x)) forall y in B(overline(x),epsilon)$, because $lambda(x-overline(x))$ can be made arbitrarily small, yet $f(overline(x)+lambda(x-overline(x)))<f(overline(x))$.
]

The normal cone of $K$ at $x$ captures all the directions starting at $x$ along which one cannot progress while being restricted to $K$ in a very strong sense, that is, it is not even possible to move along directions which correlate with a member of the cone.

#definition[
	$N_K (overline(x)) = {d : angle(d,x-overline(x)) <= 0, forall x in K}$ is the normal cone of $K$ at $overline(x)$.
]

#align(center)[#image("assets/normal-cones.svg")]

With this definition, we get an alternative characterization:
$ angle(nabla f(overline(x)),x-overline(x)) >= 0, forall x in K equiv-nabla f(overline(x)) in N_K (overline(x)). $

In words, one would like to go in the direction of $-nabla f(overline(x))$ to decrease function value, but all directions that cast a positive shadow on $-nabla f(overline(x))$ are unapproachable at $overline(x)$ by moving in feasible directions.

= Positive Semidefinite Matrices
Now we have relaxed analogs of local minima which can be easily checked and found. Nevertheless, we have lost something in this relaxation. Can we at least come up with more stringent notions that help us certify local minima in special cases?

To do this, we will consider high-order derivatives, for which the definition below captures the central property of interest.

#definition[
	A symmetric matrix $S in RR^(n times n)$ is positive semi-definite if $x^top S x >= 0, forall x in RR^n$. In this case, we will write $S succ.eq 0$. Also $S_+^(n times n)$ is the space of such matrices.

	Similarly, we'll say $S$ is positive definite or $S succ 0$ if $forall x != 0, x^top S x > 0$.
]

#theorem(title: "Spectral Theorem")[
	Any symmetric matrix can be written as $ S = sum_(i=1)^n lambda_i u_i u_i^top $ where $lambda_i in RR$ are eigenvalues, and $u_i$'s are the corresponding eigenvectors (all orthonormal).
]

#theorem[
	The following statements are all equivalent:
	+ $S succ.eq 0$.
	+ All eigenvalues of $S$ are non-negative.
	+ $S = L L^top$ for some $L in RR^(n times k)$ for some $k$.
	+ There exist vectors $v_1,dots,v_n in RR^n$ such that $S_(i j)=v_i^top v_j, forall i,j in [n]$.
	+ Determinants of all principal submatrices are non-negative.
]

Recall that principal submatrices are those formed by only keeping (or equivalently, deleting) columns and rows corresponding to the same set of indices.

#proof[
	We will leave out (5), and form a chain from (1) through (4).

	(1) $=>$ (2) For all $x$, $x^top S x >= 0$. Take $x=u_i$. Then $0 <= u_i^top S u_i = lambda_i norm(u_i)^2$.

	For the next two steps, we reinterpret matrix multiplication two ways. Consider $A in RR^(n times k)$, $B in RR^(k times n)$; then the standard definition implies $ (A B)_(i j) = (i^("th") " row of " A)^top (j^("th") " column of " B). $
	But matrix multiplication can also be written via outer products as
	$ A B = sum_(i=1)^k (i^("th") " column of " A)(i^("th") " row of " B)^top. $

	(2) $=>$ (3) Hence, $L L^top = sum_(i=1)^n (i^("th") " column of " L)(i^("th") " column of " L)^T$. So to mirror the spectral decomposition, we can choose
	$ L = [sqrt(lambda_1) u_1, sqrt(lambda_2) u_2, dots, sqrt(lambda_n) u_n], $
	which is valid because eigenvalues are non-negative.

	(3) $=>$ (4) Here we appeal to the standard interpretation of matrix products. In particular, we have $ S_(i j) = (L L^top)_(i j) = (i^("th") " row of " L)^top (j^("th") " row of " L), $
	which means we can simply choose $v_i$ to be the $i^("th")$ row of $L$.

	(4) $=>$ (1) Note that $x^top S x = sum_(i,j) x_i x_j S_(i j) = sum_(i,j) x_i x_j v_i^top v_j = norm(sum_(i=1)^n x_i v_i)^2 >= 0$.
]

Let's try to make sense of (5) by proving its equivalence for 2x2 matrices. Take $ S = mat(a,b;b,c). $

(5) says $a,c >= 0$, $a c-b^2 >= 0$. Recall that $lambda$'s are roots of $ 0=|S-lambda I|=(a-lambda)(c-lambda)-b^2=lambda^2-(a+c)lambda+(a c-b^2). $
Note $lambda_1 lambda_2 = a c-b^2 >= 0$. Also $lambda_1+lambda_2=a+c >= 0$. Hence $lambda_1,lambda_2 >= 0 => S succ.eq 0$.

(1) says $S succ.eq 0 => e_1^top S e_1 = a >= 0, e_2^top S e_2 = c >= 0$. Also $lambda_1,lambda_2 >= 0 =>lambda_1 lambda_2 = det(S)=a c-b^2 >= 0$.

In fact, the last of these works more generally. One can think of isolating a principal submatrix $M$ as choosing a binary matrix $P$ with exactly one non-zero entry per row so that $M=P S P^top$. Clearly, if $S succ.eq 0$ then $M succ.eq 0$, and the determinant of $M$, that is, the product of its eigenvalues, must be non-negative.

The more magical bit is the opposite direction: if the determinant of $P S P^top$ is non-negative for arbitrary $P$, that would be enough to get $S succ.eq 0$. (Why?) The principal submatrix condition says that checking this for binary $P$'s with one non-zero entry per row is enough.

= Second-Order Conditions for Local Minimization
We are back to the unconstrained land. Our aim this time is to get both necessary and sufficient conditions for local minima, although not both at the same time.

#proposition[
	$x$ is a local min for $f => nabla f(x)=0$ and $nabla^2 f(x) succ.eq 0$.
]

#proof[
	We already know $nabla f(x)=0$ is necessary. Say $nabla f(x)=0$ and $nabla^2 f(x) succ.eq.not 0$. Then there exists $u$ for which $u^T nabla^2 f(x) u < 0$ and $norm(u)=1$. Consider $y=x+alpha u$.
	$ f(x+alpha u) = f(x)+angle(nabla f(x),u) alpha + 1/2 alpha^2 u^T nabla^2 f(x) u + r(alpha,x,u) text(" where ") r(alpha,x,u)=o(alpha^2) $
	Hence, there exists $overline(alpha)>0$ such that for all $0 < alpha < overline(alpha)$, $|r(alpha,x,u)| <= 1/4 |alpha^2 u^T nabla^2 f(x) u|$. Now
	$ f(x+alpha u) <= f(x)+1/4 u^T nabla^2 f(x) u alpha^2 < f(x). $
	Yet, $x+alpha u$ lives in an arbitrarily small neighborhood of $x$ since this holds for all $0<alpha<overline(alpha)$. Thus, $x$ cannot be a local min.
]

#proposition[
	$nabla f(x)=0$ and $nabla^2 f(x) succ 0 => x$ is a strict local min at $x$.
]

#proof[
	Consider any $y=x+d epsilon$ with $epsilon != 0$ and unit vector $d$.
	$ f(y)=f(x)+epsilon angle(nabla f(x),d)+epsilon^2/2 d^T nabla^2 f(x)d+r(x,epsilon,d) text(" where ") r=o(epsilon^2) $

	For any fixed $d$, there exists $overline(epsilon)>0$ such that for all $0<epsilon<overline(epsilon)$ we have $|r(x,epsilon,d)| <= epsilon^2/4 |d^T nabla^2 f(x)d|$. So $f(y)>=f(x)+epsilon^2/4 d^T nabla^2 f(x)d > f(x)$. The one complication is that this choice of $overline(epsilon)$ depends on $d$.

	To fix this, define $g(d)=sup{overline(epsilon): forall 0<epsilon< overline(epsilon), f(x+epsilon d)>f(x)}$, which is continuous. Clearly $g(d)>0$ for any $d$; compactness of the set of all unit-norm $d$'s implies that $min_(norm(d)=1) g(d)$ exists, and hence, is positive. Thus, $f(y)>f(x)$ for all $y in B(x,min_(norm(d)=1) g(d)) - {x}$.
]

== Why not the constrained case?
It is also possible to define second-order conditions (sufficient and necessary) for constrained local minima. However, their practical utility is limited by concerns about computational complexity.

For example, is $x=0$ a local min for some $f(x)$ on the non-negative orthant $RR_(>=0) = {x : x>=0}$?

To establish necessary conditions, one would hope for:
+ $angle(nabla f(0),x)>=0$ for all $x in RR_(>=0)$ which is the same as $nabla f(0)>=0$.
+ And for all $d>=0$ such that $angle(nabla f(0), d) = 0$, we have $d^top nabla^2 f(0)d >=0$.

The intuition behind (2) is that for all directions one could move from $x=0$ without changing the first-order expansion, the second-order term is non-negative.

But, also, checking given a matrix $A$ if for all $x>=0$, $x^top A x>=0$ is NP-hard, which is something you would need to do with $nabla f(0)=0$. This is called the copositivity problem. Copositive matrices are of course a superset of PSD matrices, but they only attest to non-negativity for non-negative vectors, and checking membership within this set is computationally hard.

= A Few Topological Definitions

#definition[
	Closure of a set $upright("cl")(A)$ is the smallest closed set containing $A$.
]

Here _smallest_ is defined in the sense of set inclusion. So two sets may not be comparable, yet the smallest in this sense is well defined.

#definition[
	$x$ is an interior point of $A$ if there exists $epsilon>0$ such that $B(x,epsilon) subset.eq A$.

	Interior of a set $upright("int")(A)$ is the set of all interior points of $A$.
]

#definition[
	Boundary of a set $upright("bd")(A)=upright("cl")(A)-upright("int")(A)$.
]

As an example, consider $QQ$, the set of rational numbers. $upright("cl")(QQ)=RR$, $upright("int")(QQ)=emptyset$ since any small ball around a rational contains an irrational and $upright("bd")(QQ)=RR$.

As a second example, consider the open ball $A={x:norm(x)<1}$. Then $upright("cl")(A)={x:norm(x)<=1}$ is the closed ball, $upright("int")(A)=A$ and $upright("bd")(A)={x:norm(x)=1}$ is the sphere.

The annoying bit about these definitions is that the (ambient) space you are in matters. Consider the same set $A={x:x_1^2+x_2^2<1, x_3=0}$ in $RR^3$. Now, $upright("cl")(A)={x:x_1^2+x_2^2<=1, x_3=0}$ same as before, but $upright("int")(A)=emptyset$ since any 3d ball cannot be in the flat set $A$; $upright("bd")(A)=upright("cl")(A)$. Thus, every point in $A$ is now on the boundary. To avoid such nuisances, we will consider notions that are robust to such changes in description.

Recall that $A$ is an affine space if it can be written as $A=x_0+V$ for some point $x_0$ and subspace $V$.

#definition[
	$x$ is a relative interior point of $A$ if there exists $epsilon>0$, $B(x,epsilon) inter upright("aff")(A) subset.eq A$, where $upright("aff")(A)$ is the smallest affine space containing $A$.

	The relative interior of $A$, $upright("relint")(A)$, is the set of all relative interior points of $A$.
]

#definition[
	The relative boundary of a set $A$ is $upright("rbd")(A)=upright("cl")(A)-upright("relint")(A)$.
]

Thus for $A={x:x_1^2+x_2^2<1}$, $upright("relint")(A)=A$, $upright("rbd")(A)={x:x_1^2+x_2^2=1}$ regardless of the space in which $A$ is embedded.

= The Separation Theorem
We will now prove one of the most fundamental properties of convex sets.

#definition[
	$C$ and $D$ are separable if there exist $a != 0 in RR^n$, $b in RR$ such that $a^T x <= b$ for all $x in C$ and $a^T x >= b$ for all $x in D$.

	Similarly, $C$ and $D$ are strictly separable if there exist $a != 0 in RR^n$, $b in RR$ such that $a^T x < b$ for all $x in C$ and $a^T x > b$ for all $x in D$.
]

#theorem[Non-empty convex sets $C$ and $D$ can be separated if $upright("relint")(C) inter upright("relint")(D)=emptyset$.]

This theorem in fact holds with "iff" if we strengthen the definition of separability so that at least one of $C,D$ is not contained in the set ${x:a^top x=b}$. Why do we need this additional condition? Note $upright("relint")$ of a point is itself. Yet, as per the original definition, a point is separable from itself.

#align(center)[#image("assets/no-sep.svg", width: 8cm)]

As seen above, convexity of _both_ sets is necessary for this result. (Imagine B as an ellipse.)

This is a pretty fundamental result about convex sets. We will state a couple of corollaries that are convenient and then prove a weakened version.

#corollary[
	Non-empty closed convex sets $C$ and $D$ such that $C inter D=emptyset$ are separable.
]

#corollary[
	Any non-empty convex set is separable from a point on its boundary.
]

#theorem[If $C$ and $D$ are closed convex non-empty sets such that $C inter D=emptyset$ and at least one of $C,D$ is bounded, then there exist $a != 0 in RR^n$, $b in RR$ such that $a^T x < b$ for all $x in C$ and $a^T x > b$ for all $x in D$.]

We will prove this result. This is weaker than the separation theorem in many ways, but also implies strict separation.

Let us examine why each condition is required.

#align(center)[#image("assets/sep-example.svg", width: 12cm)]

#example[
	Take $C$ to be the open unit diesk centered at $(0,0)$ and $D$ to be the closed unit disk centered at $(2,0)$. These are convex, bounded, and $C inter D=emptyset$. Yet they are not strictly separable.
]

You might say that maybe being closed is not required if we weaken strict separability to strict inequality only on one side. That does not work either.

#example[
	In two dimensions, take $C=({x:x_2<0} union {x:x_2=0,x_1<0}) inter {x:norm(x)_infinity<=1}$, $D=({x:x_2>0} union {x:x_2=0,x_1>0}) inter {x:norm(x)_infinity<=1}$, which are convex and bounded, with $C inter D=emptyset$. Yet they are not strictly separable (even under the one-sided definition).
]

Similarly, one of the sets being bounded is necessary.

#example[
	In two dimensions, take $C={x:x_1<=0}$, $D={x:x_1 x_2>=1,x_1>=0}$, which are closed and convex, with $C inter D=emptyset$. Yet they are not strictly separable.
]

#proof[
	Generically, say $C$ is bounded. Consider $f(x)=inf_(d in D) norm(x-d)$ over $x in C$. Now $C$ is bounded, $f$ is continuous (Why?). Hence, $min_(x in C) f(x)$ is achieved at some $c in C$; let's fix this $c$.

	Now, consider $f(x)=norm(x-c)$ over $x in D$. Now $f$ is continuous with bounded sub-level sets, so $min$ of $f(x)$ over $x in D$ is attained. Hence let $(c,d)=upright("argmin")_(c in C,d in D) norm(c-d)$.

	The idea hereafter is simple. The bisector of the two points is the hyperplane we are looking for. Take $a=d-c !=0$ since else $C inter D != emptyset$ and $b=(norm(d)^2-norm(c)^2)/2$. By construction, $a^top ((c+d)/2)=b$, or equivalently, $(c+d)/2$ lies on the hyperplane.

	We claim $a^top x<b$ for all $x in C$ and $a^top x>b$ for all $x in D$. This might seem obvious, but I encourage you to try to complete the proof before reading further.

	#align(center)[#image("assets/strict-separation.svg")]

	Choose any $d' in D$. Since $f(x)=1/2 norm(x-c)^2$ is minimized at $d$, from @necessary-first-order, we know that $ angle(nabla f(d),d'-d)= angle(d-c,d'-d)>=0, quad forall d' in D. $
	Using this, we get
	$ a^top d' &= (d-c)^top d' >= (d-c)^top d \
		&= (d-c)^top ((d+c)/2 + (d-c)/2) = b + norm(c-d)^2/2 > b. $
	This works (anti)symmetrically for any $c' in C$.
]

= Equivalent Conditions for Convexity

#proposition[
	A differentiable function $f$ is convex on a convex set $K$ if and only if, for all $x,y in K$, we have
	$ f(y)>=f(x)+angle(nabla f(x),y-x). $
] <convex-first-order>

#proof[
	$arrow.l.double$ Consider any $x,y in K$, and $lambda x+(1-lambda)y in K$ for some $lambda in [0,1]$. Now
	$ f(x)>=f(lambda x+(1-lambda)y)+(1-lambda)angle(nabla f(lambda x+(1-lambda)y),x-y), $
	$ f(y)>=f(lambda x+(1-lambda)y)+lambda angle(nabla f(lambda x+(1-lambda)y),y-x). $
	Add (1)×$lambda$ + (2)×$(1-lambda)$ to get $lambda f(x)+(1-lambda)f(y)>=f(lambda x+(1-lambda)y)$.

	$=>$ We know $f(lambda x+(1-lambda)y)<=lambda f(x)+(1-lambda)f(y)$, and hence
	$ angle(nabla f(y),x-y) = lim_(lambda -> 0^+) (f(y+lambda(x-y))-f(y))/lambda <= f(x)-f(y). $
]

#proposition[
	A differentiable function $f$ on a convex set $K$ is convex if and only if, for all $x,y in K$, we have $ angle(nabla f(y)-nabla f(x),y-x)>=0. $
]

#proof[
	$=>$ We know
	$ f(y)>=f(x)+angle(nabla f(x),y-x), $
	$ f(x)>=f(y)+angle(nabla f(y),x-y). $
	Add these.

	$arrow.l.double$ Using the fundamental theorem of calculus, we have
	$ f(y)&=f(x)+integral_0^1 angle(nabla f(x+lambda(y-x)),y-x) dif lambda \
		&= f(x)+angle(nabla f(x),y-x)+integral_0^1 angle(nabla f(x+lambda(y-x))-nabla f(x),y-x) dif lambda\
		&>= f(x)+angle(nabla f(x),y-x). $
]

#proposition[
	A twice-differentiable function $f$ is convex on a convex set $K$ with non-empty interior if and only if, $forall x in upright("int")(K)$, we have $ nabla^2 f(x) succ.eq 0. $
]

#proof[
	$arrow.l.double$  Note that
	$ f(y)&=f(x)+integral_0^1 angle(nabla f(x+lambda(y-x)),y-x) dif lambda \
	 &= f(x)+angle(nabla f(x),y-x)
	 + integral_0^1 integral_0^lambda (y-x)^top nabla^2 f(x+mu(y-x))(y-x) dif mu dif lambda \
	 &>= f(x)+angle(nabla f(x),y-x). $

	$=>$ Consider any $x in upright("int") K$; then there exists $overline(epsilon)>0$ such that $forall 0<epsilon<overline(epsilon), x+epsilon d in K$ for any unit vector $d$. Now, let $y=x+epsilon d$.
	$ 0 &<= (angle(nabla f(y)-nabla f(x),y-x))/epsilon^2\
	&= 1/epsilon angle(nabla f(y)-nabla f(x),d)\
	&= 1/epsilon integral_0^epsilon d^top nabla^2 f(x+mu d)d dif mu. $
	Take $lim_(epsilon->0^+)$. Since this holds for all $d$, we are done. (Why?)
]

== First-Order Local-to-Global for Convex Functions
Finally, our efforts have paid off. We have the following local condition for convexity that can be cheaply checked and characterizes global optimality.

#corollary[
	Let $f$ be a differentiable convex function on $K$. Then, $overline(x)$ is a global minimum on $K$ if and only if $ angle(nabla f(overline(x)),x-overline(x))>=0, quad forall x in K. $
]

This follows by combining @necessary-first-order and @convex-first-order.

Also, conveniently, for unconstrained optimization, or more generally when $overline(x) in upright("int")(K)$, this condition is equivalent to $nabla f(overline(x)) = 0$.

== Subgradients and Their Existence

Convex univariate functions are differentiable everywhere except at countably many points. However, functions of the latter form might not be convex at all, e.g. a step function.

In convex analysis, there is an alternative way to make this concrete with a global definition, that generalizes gradients.

#definition[
	$g$ is a subgradient of $f$ at $x$ on $K$ if for all $y in K$, $f(y)>=f(x)+g^T(y-x)$.

	The set of all subgradients of $f$ at $x$ on $K$ is called the subdifferential $partial f(x)$.
]

#align(center)[#image("assets/supporting-lines.svg")]

#theorem[Let $f$ be a function on a convex set $K$. It is true that:
+ If for all $x in K$ $partial f(x)!=emptyset$, then $f$ is convex, and
+ If $f$ is convex, then, for all $x in upright("int") K$, $partial f(x)!=emptyset$.
]

#proof[
	(1) Same proof as $arrow.l.double$ for @convex-first-order.

	(2) Take any $x in upright("int")(K)$. Then $(f(x),x)$ is on the boundary of $upright("epi")(f)$. Thus, $exists (alpha,beta)!=0$ such that
	$ alpha t+beta^top y >= alpha f(x)+beta^top x quad forall (t,y) in upright("epi")(f). $
	We claim $alpha>0$. If so, taking $t=f(y)$, $f(y)>=f(x)+(-beta/alpha)^top (y-x)$, and $partial f(x)!=emptyset$.

	To finish up: First, $alpha>=0$ follows from taking $t -> infinity$, which is compatible with the epigraph. Finally, there exists $epsilon>0$ such that $y=x+epsilon d in K$ for all unit vectors $d$. Then, if $alpha=0$, note that $beta^top (y-x)=beta^top d epsilon>=0$ for all unit vectors $d$ which can only happen if $beta=0$, but both cannot be zero simultaneously.
]

#exercise[
Prove that:
	+ Non-negative sums of convex functions are convex.
	+ $f(A x+b)$ is convex for any convex $f$.
	+ $max_y f(x,y)$ is convex if for all $y$, $x mapsto f(x,y)$ is convex.
	+ $g(f(x))$ is convex if $f$ is convex and $g$ is convex nondecreasing.
	+ $min_y f(x,y)$ is convex if $f$ is jointly convex in $(x,y)$.
	+ $g(x,y)=y f(x/y)$ on ${(x,y):y>0}$ is convex if $f$ is convex on $RR^n$.
]
