#import "47851.typ": *
#show: lecture.with([Conic Programming and Duality], 3)

This marks the beginning of the second unit of the course. So far, we have treated convexity at an abstract (generic) level. In this unit, we will look at structured optimization problems, where we can peek at the internals and try to make more fine-grained statements. We will begin by trying to generalize the best-understood family of convex programs: linear programs.

#outline()

= From Linear Programs to Conic Programs

Let's recall linear programming duality. We'll build up to conic duality in the next two lectures, and this will serve as a good comparison.

$
	text("Dual: ", weight: "bold") min &c^top x #h(6em) text("Primal: ", weight: "bold") max &b^top y\
	text("subject to") &A x>=b. #h(6em) text("subject to") &A^top y=c,\
	& & y>=0.
$

*Observation:* For any primal feasible $y$ and dual feasible $x$, $c^top x>=b^top y$. This is called weak duality.

#theorem(title: "Strong LP Duality")[
	These conditions are equivalent:
		+ Dual is feasible & bounded below.
		+ Primal is feasible & bounded above.
		+ Dual is solvable.
		+ Primal is solvable.
		+ Both primal & dual are feasible. \
	Whenever they hold, both primal and dual are equal in value.
]

Let us try to abstract what $>=$ in the dual means. Of course, we know that for LPs, this corresponds to a cooridnate-wise comparison of vectors, but our goal is to isolate the essential structure.

For all $a,b,c,d in RR^n$ and $lambda in RR_(>=0)$, it is true that
+ $a>=a$,
+ $a>=b$ & $b>=c$ $=>$ $a>=c$,
+ $a>=b$ & $b>=a$ $=>$ $a=b$,
+ $a>=b$ & $c>=d$ $=>$ $a+c>=b+d$,
+ $a>=b$ $=>$ $lambda a>=lambda b$.

#definition[
	A set $K$ in a real vector space is a cone if for all $lambda >=0$ and  $x in K$, $lambda x in K$.
]

#exercise[
	Show that a cone is convex if and only if it is closed under addition.
]

Similarly, in the primal, $y>=0 equiv y in K=RR_(>=0)^n$, where
+ $K$ is a cone,
+ $K$ is pointed, i.e., $a in K$ & $-a in K$ $=>$ $a=0$,
+ $K$ is closed under addition.

Thankfully, with the interpretation that $a>=b equiv a-b in K$, these two generalizations (or abstractions) can be seen to be equivalent. For example:
+ $a>=a$ follows from $K$ being a cone.
+ $a>=b$ & $b>=c$ $=>$ $a>=c$ follows from $K$ being a cone and closed under addition.
+ $a>=b$ & $b>=a$ $=>$ $a=b$ follows from $K$ being pointed.
+ $a>=b$ & $c>=d$ $=>$ $a+c>=b+d$ follows from $K$ being closed under addition.
+ $a>=b$ $=>$ $lambda a>=lambda b$ follows from $K$ being a cone.

Conic programs generalize LPs from $RR_(>=0)^n$ to arbitrary, closed convex pointed cones with non-empty interiors. (More on the non-empty interior condition in a couple of lectures.) A general conic program involves linear optimization over a feasible set, which is an intersection of an affine space and such a cone.

#align(center)[
  $
  	max &c^top x \
  	text("subject to") &A x=b,\
  	&x in K.
  $
  *A general conic program.*
]

= A Hierarchy of Conic Programs

Conic programs can in fact represent any convex program. Let's do this in two steps. First, we can always make the objective linear.
$
	min_(x in K) f(x) equiv min_((t,x) in upright("epi")(f)) t.
$
To finish, we need to make the set a cone. Consider the cone $C(K) = {lambda (1, x): lambda>=0, x in K}$.
$
	min_(x in K) c^top x equiv &min c^top x\
	&text("subject to") t=1 \
	&text("and ") (t, x) in C(K).
$

So to specialize, we will look at _specific_ cones. Examples of admissible cones:
+ $RR_(>=0)^n$ gives us linear programs.
+ The Lorentz cone $cal(L)^(n+1)={(t,x):t>=norm(x)_2}$ gives us SOCPs (second-order conic programs).
+ $S_+^(n times n)={S:S succ.eq 0}$ gives us SDPs (semi-definite programs).

Let us define a bunch of optimization problems.#align(center)[
  *QPs:* $min x^top Q x+q^top x+d$ with $A x<=b$.

  *QCQPs:* $min x^top Q x+q^top x+d$ with $x^top Q_i x+q_i^top x+d_i<=0$ for all $i in [m]$.

  *SOCPs:* $min f^top x$ with $norm(A_i x+b_i)_2<=c_i^top x+d_i$ for all $i in [m]$.

  *SDPs:* $min tr(C X)$ with $tr(A_i X)=b_i$ for all $i in [m]$, $X⪰0$.
]

SOCPs and SDPs are always convex programs. But that is not true for Quadratic Programs (QPs) and Quadratically Constrained Quadratic Programs (QCQPs). More carefully, Convex QPs are those where $Q succ.eq 0$, and similarly Convex QCQPs are those with $Q, Q_1, dots Q_m succ.eq 0$.

Note that below $angle(A,B)=tr (A^top B) = sum_(i,j=1)^n A_(i j) B_(i j)$ and $angle(A,A)=norm(A)_F^2=sum_(i,j=1)^n A_(i j)^2$.

*On SOCPs:* Let's take a second to appreciate why SOCPs as described can be written as the intersection of an affine space and the Lorentz cone. Clearly:
$ norm(A x+b)_2<=c x+d equiv (c x+d,A x+b) in cal(L)^(n+1). $

This is an affine map. Recall any (linear) subspace can be expressed as the column-space of a matrix (the generator view) and as the null-space of another (the constraint view). Thus, by the fundamental theorem of linear algebra, there exist $E,f$ such that
$ {(c x+d,A x+b):x in RR^n}={y:E y=f}. $

As a concrete example, consider that
$
	{mat(1; 2)x: x in RR} equiv {y: mat(2; -1)^top y = 0}.
$
Essentially, the null-space vectors are all the orthogonal vectors the generators leave out. Thus, this equivalence can also be algorithmically carried out by Gram-Schmidt. Geometrically, we are appealing to the fact that all hyperplanes can be described by an anchoring point and a bunch of vectors to create the span about the anchoring point, and equivalently as the anchoring point together with a normal vector to the hyperplane.

*On SDPs:* Let's talk about SDPs. SDPs are tricky.

Why not impose linear constraints as $A X=B$ instead? Not all linear constraints can be written this way, e.g., $X_(11)=0$, which is same as $0=e_1^top X e_1=tr(X e_1 e_1^top)$.

#example[
  SDPs can have infinitely many extreme points.
  $ S={(x,y):x^2+y^2<=1}={(x,y):mat(1+x,y;y,1-x) succ.eq 0}. $
]

#example[
  For SDPs, feasibility & bounded below $⇏$ solvable.
  $ min x_1 text(" subject to ") mat(x_1,1;1,x_2) succ.eq 0. $
]

#theorem[
	$upright("LPs") subset.eq upright("Convex QPs") subset.eq upright("Convex QCQPs") subset.eq upright("SOCPs") subset.eq upright("SDPs")$.
]

#exercise[
  Write $max c^top x$ subject to $A x>=b$ as directly an SDP.
]

A useful result in proving the above theorem is next.

#lemma(title: "Schur Complements")[
	Let $A succ 0$. Then $ X=mat(A,B;B^top,C) succ.eq 0 text(" if and only if ") C-B^top A^(-1)B succ.eq 0. $
  Here $C-B^top A^(-1)B$ is called the Schur complement of $A$ in $X$.
]

#proof[
	Verify the following factorization:
  $ X=mat(I,0;B^top A^(-1),I) mat(A,0;0,C-B^top A^(-1)B) mat(I,A^(-1)B;0,I) = P^top S P. $
  The first and last matrices are transposes of the other and invertible. In fact, for any triangular matrix, its eigenvalues can be read off the diagonal. (Why?)

  Note that $x^top X x = (P x)^top S (P x)$ and $y^top S y = (P^(-1) y)^top X (P^(-1) y)$. Hence, $X succ.eq 0$ if and only if the middle matrix is PSD if and only if $C-B^top A^(-1)B succ.eq 0$.
]

In fact, such matrices $S$ and $X$ are said to be _congruent_. A well-know fact, although not common in a first linear-algebra course, is that $P$ and $X$ have the same _inertia_, where inertia of a matrix is a 3-tuple composed of the number of positive, negative, and zero eigenvalues. Thus
$ upright("inertia")(X) = upright("inertia")(A) + upright("inertia")(C-B^top A^(-1) B), $
of which we are proving a special case.

A related notion that I sometimes confuse this with is _similar_ matrices. If $S$ and $X$ were similar, then it would hold that $X=P^(-1) S P$, and a much stronger conclusion would follow that all eigenvalues are preserved in similarity transformation. (Correct me if this happens!)

We are now ready to prove the main theorem.

#proof[
	The first two inclusions are straightforward. We can take $Q=0$ and $Q_i=0$, respectively.

	*(Convex QCQPs $subset.eq$ SOCPs)* For any $Q=sum_i lambda_i u_i u_i^top⪰0$, we can define $Q^(1/2)=sum_i sqrt(lambda_i)u_i u_i^top⪰0$. Thus completing the square
	$ x^top Q x+2q^top x+d <= 0 equiv norm(Q^(1/2)x+Q^(-1/2)q)_2 <= sqrt(q^top Q q-d), $
	but this only works if $Q succ 0$ (strictly).

	*(SOCPs $subset.eq$ SDPs)* Using the Schur Complement Lemma, $norm(A x+b)<=c x+d$ if and only if
	$ c^top x+d - norm(A x+b)^2 / (c^top x+d) >= 0 equiv mat((c x+d)I,A x+b;(A x+b)^top,c x+d) succ.eq 0. $
	This is clear if $c^top x+d>0$. Incidentally the end-to-end conclusion also works when $c^top x+d=0$. (Why?) Thus, SOCPs $subset.eq$ SDPs.
]

#exercise[
	Fix the above proof so that QCCPs can be written as SOCPs even when $Q succ.eq 0$.

	_Hint:_ There are many ways to do this. For example, via SVDs, though not directly; not clear why we can tolerate a $sqrt(text("linear function"))$ term on the RHS. One idea is to make use of the following observation:
	$ norm(mat(x; y-1))_2 <= y+1 equiv x^2 <= 4y $
]

= Deriving SDP Duals

Here we will find a way to derive SDP duals. The proof of strong duality will come next week. All dualities start with the following elementary lemma. In words, it says without assuming anything about the function $f$, the second mover has the ultimate say. If min/max do not exist, interpret these as inf/sup for this section.

#lemma[ $min_(x in X) max_(y in Y) f(x,y) >= max_(y in Y) min_(x in X) f(x, y)$. ]

#proof[
  Clearly, for all $x in X$, $y in Y$, $f(x,y)>=min_(x in X) f(x,y)$. Also, $max_(y in y) g(y)>=max_(y in y) h(x)$ is a valid conclusion if $g(y)>=h(y)$ for all $y in y$. Using this, we get for all $x in X$ that
  $ max_(y in Y) f(x,y)>=max_(y in Y) min_(x in X) f(x,y). $

  Because the last inequality holds for all $x in X$, we also have that
  $ min_(x in X) max_(y in Y) f(x,y)>=max_(y in Y) min_(x in X) f(x,y). $
]

Let's try to apply this to SDPs. To begin, we have
$
	&min_X angle(C, X) text(" subject to ") angle(A_i, X)=b_i forall i in [m] text(" and ") X succ.eq 0 \
	= &min_(X succ.eq 0) max_(y in RR^m) angle(C, X)+sum_(i=1)^m y_i (b_i-angle(A_i, X))
$
This works becasue if there exists some $i in [m]$ such that $b_i > angle(A_i, X)$, then $y_i->infinity$ making a minimization program approach $+infinity$ in value, which by convention is same as infeasibility. Similarly, if $b_i < angle(A_i, X)$, then $y_i-> -infinity$ making a minimization program again approach $+infinity$ in value. Now, we use by minimax result as
$
	&min_X angle(C, X) text(" subject to ") angle(A_i, X)=b_i forall i in [m] text(" and ") X succ.eq 0 \
	= &min_(X succ.eq 0) max_(y in RR^m) angle(C, X)+sum_(i=1)^m y_i (b_i-angle(A_i, X))\
	= &min_(X succ.eq 0) max_(y in RR^m) b^top y + angle(C-sum_(i=1)^m y_i A_i, X)\
	>=&max_(y in RR^m) min_(X succ.eq 0)  b^top y + angle(C-sum_(i=1)^m y_i A_i, X).
$
To simplify the inner optimization problem, note that if $C succ.not.eq sum_(i=1)^m y_i A_i$, there exists an eigenvectors $u$ with eigenvector $lambda$ such that $u^top (C - sum_(i=1)^m y_i A_i) u = lambda < 0$. So, we can take $X=alpha u u^top$ with $alpha -> +infinity$, so that $angle(C - sum_(i=1)^m y_i A_i, X) = alpha u^top (C - sum_(i=1)^m y_i A_i) u -> -infinity$, which is the same as infeasibility for a maximization program. Also if $C succ.eq sum_(i=1)^m y_i A_i$ then for all $X succ.eq 0$, we have that $angle(C - sum_(i=1)^m y_i A_i,X)>=0$, in which case the inner minimizer might as well choose $X=0$.

#lemma[$A⪰0$, $B⪰0$ $=>$ $angle(A,B)>=0$.]

#proof[
	Let the spectral decoposition of $A$ be $A=sum_i lambda_i u_i u_i^top$. Here $A succ.eq 0 => lambda_i >= 0$. Hence
	$
		angle(A,B)&=sum_i lambda_i angle(u_i u_i^top, B) \
		&=sum_i lambda_i u_i^top B u_i>=0.
	$
]

Using all this, we get
$
	&min_X angle(C, X) text(" subject to ") angle(A_i, X)=b_i forall i in [m] text(" and ") X succ.eq 0 \
	>=&max_(y in RR^m) min_(X succ.eq 0)  b^top y + angle(C-sum_(i=1)^m y_i A_i, X)\
	= &max_(y in RR^m) cases(b^top y & "if " C succ.eq sum_(i=1)^m y_i A_i, -infinity & "otherwise")\
	= &max_(y in RR^m) b^top y text(" subject to ") C succ.eq sum_(i=1)^m y_i A_i.
$
This is the dual we sought. An inequality of the form $A_0 + sum_(i=1)^m y_i A_i succ.eq 0$ with variables $y_i$ is called a _linear matrix inequality_. Again, it might not obvious at first sight if this dual is a SDP too. However, the LMI constraint is a PSD constraint on an affine map, and hence representable as intersection of an affine space and the PSD cone. Hence, the dual is a SDP too.

By the minimax claim, the only used inequality was a valid one. Hence, we have weak duality by construction.

#lemma[
	For any primal feasible $X$, dual feasible $y$, we have $angle(C, X)>=b^top y$.
]

= Application: Approximation Algorithms for MaxCut (WIP)

Maximum cut on a graph asks for $max_(S subset.eq V) |(S times (V-S)) inter E|$. $upright("Cut")(S)$ = Edges crossing between $S$ and $V-S$.

For a graph $G=(V,E)$, the Laplacian $L in RR^(n times n)$ is such that $L_(i i)=deg(v_i)$ for all $i in [n]$ and for all $i!=j$, $L_(i j)=-1$ if $(v_i,v_j) in E$. Let $x_S in RR^n$ such $x_i=1$ if $i in S$ else $-1$.

Clearly,
$ 1/4 x_S^top L x_S=1/2 sum_((i,j) in E)(1-x_i x_j)=upright("Cut")(S). $

#claim[
  $L⪰0$.
]

#proof[
  $ x^top L x=sum_(i=1)^n deg(i)x_i^2-2 sum_((i,j) in E)x_i x_j=sum_((i,j) in E)(x_i-x_j)^2. $
]

In fact, consider any binary quadratic program with PSD cost matrix; MaxCut is a special case.

$ BQP &equiv max_(x_i^2=1 " for all " i in [n]) x^top A x \
&= max angle(A,x x^top) \
&= max angle(A,X) $
with $X_(i i)=1$ for all $i in [n]$, $X⪰0$, $X$ is rank-1.

Since any rank-1 PSD matrix can be expressed as $v v^top$.

Define $SDP equiv max angle(A,X)$ with $X_(i i)=1$ for all $i in [n]$, $X⪰0$. Clearly $SDP>=BQP$.

#theorem[(Nesterov) $BQP >= (2/pi) SDP$ for BQP with $A⪰0$.]

In fact, the proof of the theorem provides a $2/pi approx 0.636$ approximation algorithm for such BQPs. For MaxCut, this SDP can also be used to produce a $approx 0.878$ approximation (which is unimprovable upto UGC). But the above result is more general; we don't use the fact that $A$ is a Laplacian. Also, random cuts, greedy methods, LP relaxations all produce $1/2$ approximation. So, this is not unimpressive.

#proof[
  Consider any SDP-feasible $X$. There exist $v_1,dots,v_n in RR^n$ such that $X_(i j)=v_i^top v_j$. Since $X_(i i)=1$, $v_i$'s are unit vectors. Take a random $v ∼ upright("Unif")(S^(n-1))$ unit vector. Construct a vector $x$ with $x_i=upright("sign")(angle(v,v_i))$,
  where $upright("sign")(z)=bold(1)_(z>=0)-bold(1)_(z<0)$.
  Clearly, $x$ is BQP feasible.

  $ EE x_i x_j &= Pr(upright("sign")(v,v_i)=upright("sign")(v,v_j))
  -Pr(upright("sign")(v,v_i)!=upright("sign")(v,v_j)) \
  &=1-2 Pr(upright("sign")(v,v_i)!=upright("sign")(v,v_j)) \
  &=1-(4 cos^(-1)(angle(v_i,v_j)))/(2 pi) \
  &=2/pi (pi/2-cos^(-1)(angle(v_i,v_j))) \
  &=2/pi sin^(-1)(angle(v_i,v_j))=2/pi sin^(-1)(X_(i j)). $

  //#image("assets/maxcut_rounding.png", width: 27%)

  Thus,
  $ BQP>=EE x^top A x=sum_(i,j) A_(i j) EE x_i x_j
  =2/pi angle(A,sin^(-1)X), $
  applied elementwise.

  #lemma[$A⪰0$, $B⪰0$ $=>$ $A ⊙ B⪰0$. $⊙$ is element-wise product.]

  #proof[
    There exist $v_1,dots,v_n$ such for all $i,j$, $A_(i j)=v_i^top v_j$. There exist $u_1,dots,u_n$ such for all $i,j$, $B_(i j)=u_i^top u_j$. But now
    $ (A ⊙ B)_(i j)&=(v_i^top v_j)(u_i^top u_j) \
    &=sum_(k,ell) v_(i k)v_(j k)u_(i ell)u_(j ell) \
    &=angle(u_i v_i^top,u_j v_j^top), $
    vectorizing a matrix. Thus $A ⊙ B⪰0$.
  ]

  #corollary[$A⪰0$ $=>$ $A^(⊙ n)⪰0$ for all $n in ZZ_(>=0)$, elementwise powers.]

  #lemma[$sin^(-1)(X)⪰X$ for all $X⪰0$.]

  #proof[
    $ sin^(-1)(x)&=x+x^3/(2 times 3)+(1 times 3 x 5)x^5/(2 times 4 times 5) \
    &quad +(1 times 3 times 5)x^7/(2 times 4 times 6 times 7)+dots \
    &=x+sum_(k>=2) a_k x^k $
    for some $a_k>=0$.

    This power series converges to $sin^(-1)(x)$ uniformly on $[-1,1]$. Using the last corollary, $sin^(-1)(X)⪰X$.
  ]

  To complete the proof,

  $ BQP&>=2/pi angle(A,sin^(-1)(X)) \
  &=2/pi (angle(A,sin^(-1)(X)-X)+angle(A,X)) \
  &>=2/pi SDP. $
]

= Application: Stabilizability of Linear Dynamical Systems

Let's do another application, before we move onto SDP duality, and then finally conic duality.

A fundamental question about dynamical system $x_(t+1)=f(x_t)$ is if their long-term behaviors can be neatly categorized. As a concrete question, consider a linear dynamical system $x_(t+1)=A x_t$. When does such a system end up at 0 (i.e., $x_t->0$ as $t->infinity$) regardless of initialization $(x_0)$? Such systems are called stable.

#definition[A linear dynamical system with transition matrix $A$ is stable if for all $x_0 in RR^n$, $x_t=A^t x_0->0$ as $t->infinity$.]

Consider a diagonalizable matrix $A=U Sigma U^(-1)$. Here, $U$ is the matrix of eigenvectors concatenated together. Note that the eigenvectors may not be orthogonal, neither might eigenvalues be real. Any $x$ can be written as $x=sum_i c_i u_i$. $A^t x=sum_i lambda_i^t c_i u_i$. From this we get:

*Observation:* $A$ is stable iff $lim_(t->infinity) A^t=0$ iff $rho(A)<1$.

Here, $rho(A)=max_i |lambda_i|$, modulus. Note, $rho(dot)$ is non-convex.

#example[
  Consider $A=mat(0,1;0,0)$, $B=mat(0,0;1,0)$. $rho(A)=rho(B)=0$, $norm(A)_2=norm(B)_2=1$. Yet, $rho((A+B)/2)=1/2$, $norm((A+B)/2)_2=1/2$.
]

Operator/spectral norm: $norm(A)_2=max_(norm(x)_2<=1) norm(A x)_2$.

Stability is a property that holds over all trajectories (each with unique starting state). Lyapunov gave a criteria to certify the stability of a system without analyzing all its trajectories. The idea was to find a function $V(x)$ such $V(x)>0$ for all $x!=0$, $V(0)=0$, $V(f(x))<=V(x)$ for all $x$, with strict inequality for $x!=0$. Under such conditions, it is clear $V(x_t)$ decreases asymptotically. The nice bit is that he showed such a certificate always exists for “nice-enough” systems. Today, such analysis is called potential-based.

#lemma[$A$ is stable iff there exists $P⪰I$ such $P succ A^top P A$.]

#proof[
  $<=$ Take $V(x)=x^top P x$ where $P⪰I$, $c P⪰A^top P A$ for $0<=c<1$. $V(0)=0$, $V(x)>=norm(x)^2>0$ for all $x!=0$ and $V(A x)=x^top A^top P A x<=c V(x)$. Therefore, $V(A^t x)<=c^t V(x)->0$
  as $t->infinity$, which also implies $x_t->0$ as $t->infinity$.

  $=>$ Take $P=sum_(t>=0) (A^t)^top A^t$. Clearly, $P⪰I$ as $x_0^top P x_0=sum_(t>=0) norm(x_t)^2$. Also
  $ A^top P A=sum_(t>=0)(A^(t+1))^top A^(t+1)=P-I, $
  hence $P succ A^top P A$. The only bit left out is that $P$ as defined has to exist. We use the following fact that follows from Jordan decomposition:

  *Fact:* $norm(A^t)<=O(t^n rho(A)^t)$.

  Hence $P$ is well-defined, as $rho(A)<1$.
]

Note that ${P:P⪰I,c P⪰A^top P A}$ is a SDP feasibility program; and one can do bisection on guess of $c$ in $[0,1)$.

A more involved question is that of stabilizability. Consider a dynamical system with inputs, $x_(t+1)=A x_t+B u_t$. Many systems are not intrinsically stable, but can be stable under correct choice of control inputs. This is what stabilizability captures.

#definition[A linear dynamical system $(A,B)$ with inputs is stabilizable if for all $x_0 in RR^n$, there exist $u_1,u_2,dots in RR^m$ such $x_t->0$ as $t->infinity$.]

The dynamical systems we have considered so far are Markovian; all information at time $t$ is summarized in the state $x_t$. Using this and linear dynamics, one can gain an alternate characterization.

*Observation:* $(A,B)$ is stabilizable iff there exists $K$ such $rho(A+B K)<1$.

The basic idea is to use inputs of the form $u_t=K x_t$ regardless of $x_t$. If so,
$ x_(t+1)=A x_t+B u_t=(A+B K)x_t=(A+B K)^t x_0. $
Thus the dynamical system looks like one without inputs where $A_("eff")=A+B K$. Since such a $K$ exists with $rho(A+B K)<1$ whenever $(A,B)$ is stabilizable, can one find it efficiently?

Note that $min_K rho(A+B K)$ is a non-convex program.

Instead, let's try to use our lemma.

#align(center)[
  $(A,B)$ is stabilizable

  $equiv$ there exists $K$ such $rho(A+B K)<1$

  $equiv$ there exists $K$ such there exists $P⪰I$,
  $ P succ (A+B K)^top P(A+B K) $

  $equiv$ there exists $K$ such there exists $P⪰I$,
  $ mat(P,((A+B K)P)^top;(A+B K)P,P) succ 0. quad "(1)" $
]

Why? We already know $P succ 0$, thus by Schur complement lemma, satisfaction of the above constraint is same as
$ P-(A+B K)P P^(-1)P(A+B K)^top succ 0. $

Note this does not match what we had in the last step. But that's fine since what we have as same as $rho((A+B K)^top)<1$, but $(A+B K)^top$ and $A+B K$ have the same characteristic equation and hence eigenvalues. However we still don't have a SDP at this stage due to non-linearity produced by multiplication of $K$ & $P$; both are variables.

#align(center)[
  $equiv$ there exists $K in RR^(n times m)$, $P⪰I$,
  $ mat(P,(A P+B L)^top;A P+B L,P) succ 0. quad "(2)" $
]

At first glance, this might seem like a relaxation, we used $L=K P$ to get a (near) SDP feasibility program. But this is in fact a reformulation since one can go back (i.e. construct a feasible solution for the original program) by setting $K=L P^(-1)$.

For this program, one does once again perform bisection with $succ$ constraints.

What we were able to accomplish in the last step (an exact reformulation) was a lucky coincidence. As an example, if we want to impose a box constraint on $K$, we can do so in (1) which was (nominally) non-convex anyways, but now an exact reformulation such as (2) as an SDP is no longer possible. In fact, this problem is NP-hard.

= Conic Duality

Consider $min c^top x$ subject to $A x >=_K b$. Before we get to strong duality, it is worth noting that it's not clear how we even define a dual in the first place. For LPs, our derivation was based on the fact that non-negative combinations of vector inequalities always collapse to valid scalar inequalities. But this is not true for CPs.

#example[
  $cal(L)^(2+1)={(t,x):t>=norm(x)_2}$. Clearly, $(1,0,-1)>=_(cal(L)^(2+1))(0,0,0)$. But
  $ (0,0,1)^top(1,0,-1)=-1 ≱0=(0,0,1)^top(0,0,0). $
]

You might guess that the multiplier should come from $K$ itself. This is also false.

#example[
  $K={(x,y):|y|<=100x}$. Then $(1,-100)>=_K(0,0)$. But $(1,100) in K$, yet
  $ angle((1,100),(1,-100))=1-100^2 ≱0=angle((1,100),(0,0)). $
]

//#image("assets/dual_cone.png", width: 15%)

The true answer follows from a tautology. For all $y in K_*$, $a>=_K b$ $=>$ $y^top a>=y^top b$ for all $a,b in K$.

#definition[The dual cone $K_*$ of a set $K subset.eq RR^n$ is $K_*={y:$ for all $x in K$, $y^top x>=0}$.]

Now, $min_x c^top x$ subject to $A x>=_K b$ is primal. Take $max b^top y$ subject to $A^top y=c$, $y in K_*$ to be the dual.

For any feasible pair $(x,y)$,
$ c^top x=y^top A x=y^top(A x-b)+b^top y>=b^top y. $

#lemma[For primal feasible $x$, dual feasible $y$, we have $c^top x>=b^top y$.]

#lemma[Let $K$ be a non-empty set in $RR^n$. Then:
  + $K_*$ is a closed convex cone.
  + $upright("int")(K)!=emptyset$ $=>$ $K_*$ is pointed.
  + $K$ is a closed convex pointed cone $=>$ $upright("int")(K_*)!=emptyset$.
  + $K$ is a closed convex cone $=>$ $K_(**) = K$.
]

#proof[
  (1)
  $ K_* = ⋂_(x in K) {y:y^top x>=0} $
  is intersection of halfspaces (hence convex) and intersection of closed sets (hence closed). Cone is true from definition, because homogeneous inequalities are involved.

  (2) Say $K_*$ is not pointed. Then $K_*$ contains an entire line passing through origin, say ${lambda a:lambda in RR}$. Then $a^top x>=0$ for all $x in K$ and $(-a)^top x>=0$ for all $x in K$. Thus $a^top x=0$. But then $K$ is flat, it's entirely contained in the maximal subspace orthogonal to $a$. Hence, $upright("int")(K)=emptyset$.

  #exercise[
    Prove part (3) of the above lemma.
  ]

  (4) First $K subset.eq K_(**)$. Why? Pick any $x in K$, clearly for all $y in K_*$, $y^top x>=0$. Thus $x in K_(**)$. For the other direction, we'll invoke separation. Let's do the contrapositive. Say $x ∉ K$. There exists $y !=0$ such $y^top x<y^top z$ for all $z in K$. But then $z^top y$ is lower bounded on $K$, hence $min_(z in K) z^top y=0$. Hence $y in K_*$, and from $y^top x<0$, $x ∉ K_(**)$.
]

#corollary[$K$ is a closed convex pointed cone with non-empty interior iff $K_*$ is a closed convex pointed cone with non-empty interior.]

#align(center)[
  $min c^top x$ subject to $A x>=_K b$, where $K$ is closed convex pointed cone, $upright("int")(K)!=emptyset$: primal (P).
]
#align(center)[
  $max b^top y$ subject to $A^top y=c$, $y in K_*$: dual (D).
]

#theorem[If P is bounded below and strictly feasible (i.e. $A x>_K b$), then D is solvable and optimal values match.]

By symmetry, if D is bounded above and strictly feasible ($y>_(K_*)0$) then P is solvable and optimal values match.

#corollary[If one of P/D is bounded and strictly feasible, and $(x,y)$ is P/D feasible pair, then $(x,y)$ is P,D optimal iff $b^top y=c^top x$ iff $y^top(A x-b)=0$.]

#proof[
  Under given conditions, we know optimal values are equal. We also know from weak duality that for any pair of feasible solutions
  $ b^top y<=D^*<=P^*<=c^top x. $
  Thus attainment of optimal values at $(x,y)$ for $(P,D)$ happens iff $b^top y=c^top x$. From the proof of weak duality this is equivalent to $y^top(A x-b)=0$, the only inequality being tight.
]

#proof(title: "Proof of strong duality")[
  Let $P^*$ be the optimal value for the primal problem. We will construct a dual feasible $y^*$ such $b^top y^*>=P^*$, via the separation theorem.

  Let $V={A x-b:c^top x<=P^*}$.
  It is non-empty, because any affine/linear map of a halfspace can not be empty. Also $V inter upright("int")(K)=emptyset$, otherwise given $x$ such $c^top x<=P^*$ and $A x>_K b$, one can strictly decrease the objective for a close enough feasible point, assuming $c!=0$.

  Aside note: $c!=0$ holds generically. If not, $P^*=0$. Take $y=0$, $A^top y=c$ holds and we have $0=b^top y>=P^*$.

  Now, by (actual) separation theorem, there exists $y!=0 in RR^m$ such that
  $ sup_(x in V) y^top x<=inf_(x in K) y^top x. $
  So, $y^top x$ is bounded below on $K$. But $inf_(x in K) y^top x$ generically only take values ${0,-infinity}$. Thus
  $ sup_(x in V) y^top x<=0<=inf_(x in K) y^top x, $
  and $y in K_*$ follows from the second inequality. But also the first of these is equivalent to $c^top x<=P^* => (A^top y)^top x<=y^top b$. This can only happen if $A^top y=alpha c$ for some $alpha>=0$. If $alpha>0$, then take $y^*=y/alpha in K_*$ for which $A^top y^*=c$ and
  $ b^top y^*>=(A^top y)^top x/alpha=c^top x $
  for all $c^top x<=P^*$. We are done. If $alpha=0$, $A^top y=0$ whereas $b^top y>0$ for there exists $x$, $A x>_K b$. $y in K_*$ is not null, and hence $y^top(A x-b)>0$ $=>$ $c^top x<=P^*$ yet $y^top A x=0>b^top y$. This is a contradiction because then $c^top x<=P^*$ has to empty.
]

= The S-Lemma

Consider a non-convex QCQP $min_x q_a(x)=x^top Q_a x+q_a^top x+c_a$ subject to $q_b(x)=x^top Q_b x+q_b^top x+c_b>=0$. $Q_a,Q_b$ are arbitrary, need not be $Q_a⪰0$, $Q_b⪰0$.

Such programs are efficiently solvable, despite apparent non-convexity.

#theorem[Let there exist $x$ such $q_b(x)>0$ strictly. Then for all $x$, $q_b(x)>=0$ $=>$ $q_a(x)>=0$ implies there exists $lambda>=0$ such that for all $x$, $q_a(x)>=lambda q_b(x)$.]

$<=$ above is evident. Such $lambda$ “certifies” the former clause. This has the same spirit as Farkas' lemma or, theorem of alternatives you might have seen before:

#align(center)[
  $min_x q_a(x)$ subject to $q_b(x)>=0$
]

#align(center)[
  $=max_gamma gamma$ subject to for all $x$, $q_b(x)>=0$ $=>$ $q_a(x)-gamma>=0$
]

#align(center)[
  $=max_(gamma,lambda>=0) gamma$ subject to for all $x$, $q_a(x)-gamma>=lambda q_b(x)$ (using S-lemma)
]

#align(center)[
  $=max_(gamma,lambda>=0) gamma$ subject to for all $x$,
]
$ (x,1)^top mat(Q_a-lambda Q_b,(q_a-lambda q_b)/2;(q_a-lambda q_b)^top/2,c_a-gamma-lambda c_b)(x,1)>=0 $

#align(center)[
  $=max_(gamma,lambda>=0) gamma$ subject to
  $ mat(Q_a-lambda Q_b,(q_a-lambda q_b)/2;(q_a-lambda q_b)^top/2,c_a-gamma-lambda c_b)⪰0. $
]
(This is an SDP.)

#proof[
  This needs a brief justification. $=>$ is clear. You might reasonably suspect we are asking for more than is necessary, e.g., $(x,1)^top M(x,1)>=0$ also, in the latter formulation. To complete the argument, say there exists $y$, $y^top M y<0$. Then,
  (1) last coordinate of $y!=0$, then by rescaling there exists $x$ such that $(x,1)^top M(x,1)<0$.
  (2) last coordinate of $y=0$, since $y^top M y$ is continuous, there exists $z$ with non-zero last coordinate with $z^top M z<0$.
]

Finally, we'll prove S-lemma. Let's start with something simpler.

#lemma[For symmetric $A,B$, $tr(A)>=0>tr(B)$ $=>$ there exists $x$ such $x^top A x>=0$ and $x^top B x<0$.]

#proof[
  Two observations: for any $x in {-1,1}^n$, $x^top D x=tr(D)$ for any diagonal $D$, since $x_i^2=1$. Secondly, for any matrix $M$,
  $ EE_(x ∼ upright("Unif")({-1,1}^n)) x^top M x=tr(M). $

  Now choose $xi ∼ upright("Unif")({-1,1}^n)$, $x=U xi$, where $B=U Lambda U^top$ diagonal.
  $ x^top B x&=xi^top U^top U Lambda U^top U xi \
  &=tr(Lambda)=tr(Lambda U^top U)=tr(B)<0, $
  yet
  $ EE x^top A x&=EE xi^top(U^top A U)xi \
  &=tr(U^top A U)=tr(A U U^top)=tr(A)>=0. $

  Therefore, there exists $x$ such $x^top A x>=0$ and $x^top B x<0$.
]

#theorem[Say there exists $macron(x)$ such $macron(x)^top A macron(x)>0$. Then for all $x$, $x^top A x>=0$ $=>$ $x^top B x>=0$ implies there exists $lambda>=0$ such that $B⪰lambda A$.]

#proof[
  Consider $OPT=min_(x^top A x>=0) x^top B x$. Say for all $x$, $x^top A x>=0$ $=>$ $x^top B x>=0$ and there exists $macron(x)$, $macron(x)^top A macron(x)>0$. Then $OPT=0$. Consider $OPT'=min_(x^top A x>=0, x^top x=n) x^top B x$. $OPT'$ is also strictly feasible; $OPT'>=0$.
  #align(center)[
    $equiv min_X tr(B X)$ subject to $tr(A X)>=0$, $tr(X)=n$, $X⪰0$, $X$ is rank 1.
  ]

  #align(center)[
    $>= min_X tr(B X)$ subject to $tr(A X)>=0$, $tr(X)=n$, $X⪰0$ (SDP).
  ]

  #align(center)[
    SDP dual is $max_(u,lambda>=0) n u$ subject to $B⪰u I+lambda A$.
  ]

  #exercise[
    Show that SDP dual is indeed the dual of the SDP.
  ]

  Let's say we guarantee that SDP $>=0$ and that strong duality holds. Then, we are done as then there exist $u,lambda>=0$ such $B⪰u I+lambda A$ and hence there exists $lambda>=0$, $B⪰lambda A$.

  Take
  $ macron(X)=(macron(x)macron(x)^top+alpha I)/tr(macron(x)macron(x)^top+alpha I) n, $
  for small enough $alpha>0$, is strictly feasible for SDP. And for all $lambda>0$, there exists $u$ such $B succ u I+lambda A$. Thus the SDP dual is also strictly feasible and strong duality holds.

  #exercise[
    In words, describe why adding $x^top x=n$ as a constraint was essential for this argument.
  ]

  Now to argue SDP $>=0$, observe: Compactness ensures that it is solvable; let $X^*$ be the solution. There exists $L$, $X=L L^top$. Consider that
  $ tr(A X^*)=tr(L^top A L)>=0, $
  and
  $ tr(B X^*)=tr(L^top B L)=n u^*. $
  If $u^*<0$, then by our sign lemma, there exists $x$ such $(L x)^top A(L x)>=0$ and $(L x)^top B(L x)<0$. But this contradicts for all $x$, $x^top A x>=0$ $=>$ $x^top B x>=0$.
]

#proof(title: "Proof of the original theorem")[
  The idea is to use homogeneous polynomials.

  $ macron(q)_a(x,t)=x^top Q_a x+q_a^top x t+c_a t^2. $

  $ macron(q)_b(x,t)=x^top Q_b x+q_b^top x t+c_b t^2. $

  Given there exists $x$ with $q_b(x)>0$, we have $macron(q)_b(x,1)>0$.

  We'll prove for all $(x,t)$, $macron(q)_b(x,t)>=0$ $=>$ $macron(q)_a(x,t)>=0$. If so, using our simplified S-lemma, there exists $lambda>=0$, $macron(q)_a(x,t)>=lambda macron(q)_b(x,t)$ for all $(x,t)$. Set $t=1$.

  Let's say our sub-claim is untrue. Then, there exists $(x,t)$ such $macron(q)_b(x,t)>=0$ but $macron(q)_a(x,t)<0$. Consider
  (1) $t!=0$. Then substituting $(x/t,1)$ results in contradiction.
  (2) $t=0$, $macron(q)_b(x,t)>0$, by continuity there exists $(y,s)$ with $s!=0$ with $macron(q)_b(y,s)>=0$. Case (1) applies.

  #exercise[
    Argue why the remaining case $t=0$, $macron(q)_b(x,t)=0$ results in a contradiction too, completing the proof.
  ]
]
