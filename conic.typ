#import "47851.typ": *
#show: lecture.with([Conic Programming and Duality (WIP)], 3)

#outline()

Let's recall linear programming duality. We'll build up to conic duality in the next two lectures, and this will serve as a good comparison.

*Dual:* $min c^T x$, $A x>=b$. *Primal:* $max b^T y$, $A^T y=c$, $y>=0$.

*Observation:* For any primal feasible $y$ and dual feasible $x$, $c^T x>=b^T y$. This is called weak duality.

#theorem[These conditions are equivalent:
(A) dual is feasible & bounded below,
(B) primal is feasible & bounded above,
(C) dual is solvable,
(D) primal is solvable,
(E) both primal & dual are feasible.
Whenever they hold, both primal & dual are equal in value.]

$=>$ in the dual satisfies $y>=0 equiv y in K=RR_(>=0)^n$.

For all $a,b,c,d in RR^n$ which satisfies that and $lambda in RR_(>=0)$, (1) $a>=a$, (2) $a>=b$ & $b>=c$ $=>$ $a>=c$, (3) $a>=b$ & $b>=a$ $=>$ $a=b$, (4) $a>=b$ & $c>=d$ $=>$ $a+c>=b+d$, (5) $a>=b$ $=>$ $lambda a>=lambda b$.

(1) $K$ is a cone. (2) $K$ is pointed, i.e., $a in K$ & $-a in K$ $=>$ $a=0$. (3) $K$ is closed under addition. These are related as $a>=b equiv a-b in K$.

Lines underscore equivalence of the stated conditions. Conic programs generalize LPs from $RR_(>=0)^n$ to arbitrary, closed convex pointed cones with non-empty interiors. A general conic program involves linear optimization over a feasible set, which is an intersection of an affine space and such a cone.

*General conic program:* $max c^T x$ subject to $A x=b$, $x in K$.

*Examples of admissible cones:*
+ $RR_(>=0)^n$ gives linear programs.
+ The Lorentz cone $cal(L)^(n+1)={(t,x):t>=norm(x)_2}$ gives SOCPs (second-order conic programs).
+ $S_+^(n times n)={S:S⪰0}$ gives rise to SDPs (semi-definite programs).

Note $angle(A,B)=upright("Tr")(A^T B)$ below, and $angle(A,A)=norm(A)_F^2=sum_(i,j=1)^n A_(i j)^2$.

$upright("LPs") subset.eq upright("QPs") subset.eq upright("QCQPs") subset.eq upright("SOCPs") subset.eq upright("SDPs").$

*Convex QPs:* $min x^T Q x+q^T x+d$ with $A x<=b$.

*Convex QCQPs:* $min x^T Q x+q^T x+d$ with $x^T Q_i x+q_i^T x+d_i<=0$ for all $i in [m]$.

*SOCPs:* $min f^T x$ with $norm(A_i x+b_i)_2<=c_i^T x+d_i$ for all $i in [m]$.

*SDPs:* $min upright("Tr")(C X)$ with $upright("Tr")(A_i X)=b_i$ for all $i in [m]$, $X⪰0$.

Let's talk about SDPs. SDPs are tricky.
+ SDPs can have infinitely many extreme points. $S={(x,y):x^2+y^2<=1}={(x,y):mat(1+x,y;y,1-x)⪰0}$.
+ For SDPs, feasibility & bounded below $⇏$ solvable. $min x_1$, $mat(x_1,1;1,x_2)⪰0$.
+ Why not impose linear constraints as $A x=B$? Not all linear constraints can be written this way, e.g., $X_(11)=0$, which is same as $0=e_1^T X e_1=upright("Tr")(X e_1 e_1^T)$.

//#image("assets/sdp_nonattainment.png", width: 17%)

*Exercise:* Write $max c^T x$ subject to $A x>=b$ as a SDP.

*Exercise:* Write a convex generic QCQP as a SOCP. For any $Q=sum_i lambda_i u_i u_i^T⪰0$, $Q^(1/2)=sum_i sqrt(lambda_i)u_i u_i^T⪰0$. You might be tempted with $x^T Q x+2q^T x+d<=0 equiv norm(Q^(1/2)x+Q^(-1/2)q)_2<=sqrt(q^T Q q-d)$, but this only works if $Q succ 0$ (strictly). How to handle $Q⪰0$?

Also, let's take a second to understand why SOCPs as described can indeed be written as the intersection of an affine space and the Lorentz cone. Clearly, $norm(A x+b)_2<=c x+d equiv (c x+d,A x+b) in cal(L)^(n+1)$.

This is an affine map. Thus, by the fundamental theorem of linear algebra, there exist $E,f$ such that ${(c x+d,A x+b):x in RR^n}={y:E y=f}$. Recall any (linear) subspace can be expressed as the columnspace of a matrix and as the null-space of another.

*Exercise:* Show that any convex optimization problem can be recast as a conic program (for some cone).

*Schur complement lemma.* Let $A succ 0$. Then $X=mat(A,B;B^T,C)⪰0$ iff $C-B^T A^(-1)B⪰0$. Schur complement of $A$ in $X$.

*Proof:* $X=mat(I,0;B^T A^(-1),I) mat(A,0;0,C-B^T A^(-1)B) mat(I,A^(-1)B;0,I)$. (Verify!) These are transposes of the other and invertible. Hence, $X⪰0$ iff the middle matrix is PSD iff $C-B^T A^(-1)B⪰0$. $square$

Using this, $norm(A x+b)<=c x+d$ iff $mat((c x+d)I,A x+b;(A x+b)^T,c x+b)⪰0$. (Handle the case when $c x+d=0$ separately.) Thus, SOCPs $subset.eq$ SDPs.

= Approximation algorithms for MaxCut

Maximum cut on a graph asks for $max_(S subset.eq V) |(S times (V-S)) inter E|$. $upright("Cut")(S)$ = Edges crossing between $S$ and $V-S$.

For a graph $G=(V,E)$, the Laplacian $L in RR^(n times n)$ is such that $L_(i i)=deg(v_i)$ for all $i in [n]$ and for all $i!=j$, $L_(i j)=-1$ if $(v_i,v_j) in E$. Let $x_S in RR^n$ such $x_i=1$ if $i in S$ else $-1$.

Clearly, $1/4 x_S^T L x_S=1/2 sum_((i,j) in E)(1-x_i x_j)=upright("Cut")(S)$.

*Claim:* $L⪰0$.

*Proof:* $x^T L x=sum_(i=1)^n deg(i)x_i^2-2 sum_((i,j) in E)x_i x_j=sum_((i,j) in E)(x_i-x_j)^2$. $square$

In fact, consider any binary quadratic program with PSD cost matrix; MaxCut is a special case.

$BQP equiv max_(x_i^2=1 " for all " i in [n]) x^T A x
= max angle(A,x x^T)
= max angle(A,X)$ with $X_(i i)=1$ for all $i in [n]$, $X⪰0$, $X$ is rank-1.

Since any rank-1 PSD matrix can be expressed as $v v^T$.

Define $SDP equiv max angle(A,X)$ with $X_(i i)=1$ for all $i in [n]$, $X⪰0$. Clearly $SDP>=BQP$.

#theorem[(Nesterov) $BQP >= (2/pi) SDP$ for BQP with $A⪰0$.]

In fact, the proof of the theorem provides a $2/pi approx 0.636$ approximation algorithm for such BQPs. For MaxCut, this SDP can also be used to produce a $approx 0.878$ approximation (which is unimprovable upto UGC). But the above result is more general; we don't use the fact that $A$ is a Laplacian. Also, random cuts, greedy methods, LP relaxations all produce $1/2$ approximation. So, this is not unimpressive.

*Proof:* Consider any SDP-feasible $X$. There exist $v_1,dots,v_n in RR^n$ such that $X_(i j)=v_i^T v_j$. Since $X_(i i)=1$, $v_i$'s are unit vectors. Take a random $v ∼ upright("Unif")(S^(n-1))$ unit vector. Construct a vector $x$ with $x_i=upright("sign")(angle(v,v_i))$, where $upright("sign")(z)=bold(1)_(z>=0)-bold(1)_(z<0)$. Clearly, $x$ is BQP feasible.

$EE x_i x_j = Pr(upright("sign")(v,v_i)=upright("sign")(v,v_j))
-Pr(upright("sign")(v,v_i)!=upright("sign")(v,v_j))$
$=1-2 Pr(upright("sign")(v,v_i)!=upright("sign")(v,v_j))
=1-(4 cos^(-1)(angle(v_i,v_j)))/(2 pi)$
$=2/pi (pi/2-cos^(-1)(angle(v_i,v_j)))
=2/pi sin^(-1)(angle(v_i,v_j))=2/pi sin^(-1)(X_(i j))$.

//#image("assets/maxcut_rounding.png", width: 27%)

Thus, $BQP>=EE x^T A x=sum_(i,j) A_(i j) EE x_i x_j
=2/pi angle(A,sin^(-1)X)$, applied elementwise.

#lemma[$A⪰0$, $B⪰0$ $=>$ $angle(A,B)>=0$.]

*Proof:* Let $A=sum_i lambda_i u_i u_i^T$. $angle(A,B)=sum_i lambda_i angle(u_i u_i^T,B)=sum_i lambda_i u_i^T B u_i>=0$. $square$

#lemma[$A⪰0$, $B⪰0$ $=>$ $A ⊙ B⪰0$. $⊙$ is element-wise product.]

*Proof:* There exist $v_1,dots,v_n$ such for all $i,j$, $A_(i j)=v_i^T v_j$. There exist $u_1,dots,u_n$ such for all $i,j$, $B_(i j)=u_i^T u_j$. But now $(A ⊙ B)_(i j)=(v_i^T v_j)(u_i^T u_j)=sum_(k,ell) v_(i k)v_(j k)u_(i ell)u_(j ell)
=angle(u_i v_i^T,u_j v_j^T)$, vectorizing a matrix. Thus $A ⊙ B⪰0$. $square$

#corollary[$A⪰0$ $=>$ $A^(⊙ n)⪰0$ for all $n in ZZ_(>=0)$, elementwise powers.]

#lemma[$sin^(-1)(X)⪰X$ for all $X⪰0$.]

*Proof:* $sin^(-1)(x)=x+x^3/(2 times 3)+(1 times 3 x 5)x^5/(2 times 4 times 5)+(1 times 3 times 5)x^7/(2 times 4 times 6 times 7)+dots
=x+sum_(k>=2) a_k x^k$ for some $a_k>=0$.

This power series converges to $sin^(-1)(x)$ uniformly on $[-1,1]$. Using the last corollary, $sin^(-1)(X)⪰X$. $square$

To complete the proof,
$BQP>=2/pi angle(A,sin^(-1)(X))
=2/pi (angle(A,sin^(-1)(X)-X)+angle(A,X))
>=2/pi SDP$. $square$

= Stabilizability of linear systems

Let's do another application, before we move onto SDP duality, and then finally conic duality.

A fundamental question about dynamical system $x_(t+1)=f(x_t)$ is if their long-term behaviors can be neatly categorized. As a concrete question, consider a linear dynamical system $x_(t+1)=A x_t$. When does such a system end up at 0 (i.e., $x_t->0$ as $t->infinity$) regardless of initialization $(x_0)$? Such systems are called stable.

#definition[A linear dynamical system with transition matrix $A$ is stable if for all $x_0 in RR^n$, $x_t=A^t x_0->0$ as $t->infinity$.]

Consider a diagonalizable matrix $A=U Sigma U^(-1)$. Here, $U$ is the matrix of eigenvectors concatenated together. Note that the eigenvectors may not be orthogonal, neither might eigenvalues be real. Any $x$ can be written as $x=sum_i c_i u_i$. $A^t x=sum_i lambda_i^t c_i u_i$. From this we get:

*Observation:* $A$ is stable iff $lim_(t->infinity) A^t=0$ iff $rho(A)<1$.

Here, $rho(A)=max_i |lambda_i|$, modulus. Note, $rho(dot)$ is non-convex. Consider $A=mat(0,1;0,0)$, $B=mat(0,0;1,0)$. $rho(A)=rho(B)=0$, $norm(A)_2=norm(B)_2=1$. Yet, $rho((A+B)/2)=1/2$, $norm((A+B)/2)_2=1/2$. Operator/spectral norm: $norm(A)_2=max_(norm(x)_2<=1) norm(A x)_2$.

Stability is a property that holds over all trajectories (each with unique starting state). Lyapunov gave a criteria to certify the stability of a system without analyzing all its trajectories. The idea was to find a function $V(x)$ such $V(x)>0$ for all $x!=0$, $V(0)=0$, $V(f(x))<=V(x)$ for all $x$, with strict inequality for $x!=0$. Under such conditions, it is clear $V(x_t)$ decreases asymptotically. The nice bit is that he showed such a certificate always exists for “nice-enough” systems. Today, such analysis is called potential-based.

#lemma[$A$ is stable iff there exists $P⪰I$ such $P succ A^T P A$.]

*Proof:* $<=$ Take $V(x)=x^T P x$ where $P⪰I$, $c P⪰A^T P A$ for $0<=c<1$. $V(0)=0$, $V(x)>=norm(x)^2>0$ for all $x!=0$ and $V(A x)=x^T A^T P A x<=c V(x)$. Therefore, $V(A^t x)<=c^t V(x)->0$ as $t->infinity$, which also implies $x_t->0$ as $t->infinity$.

$=>$ Take $P=sum_(t>=0) (A^t)^T A^t$. Clearly, $P⪰I$ as $x_0^T P x_0=sum_(t>=0) norm(x_t)^2$. Also $A^T P A=sum_(t>=0)(A^(t+1))^T A^(t+1)=P-I$, hence $P succ A^T P A$. The only bit left out is that $P$ as defined has to exist. We use the following fact that follows from Jordan decomposition:

*Fact:* $norm(A^t)<=O(t^n rho(A)^t)$.

Hence $P$ is well-defined, as $rho(A)<1$. $square$

Note that ${P:P⪰I,c P⪰A^T P A}$ is a SDP feasibility program; and one can do bisection on guess of $c$ in $[0,1)$.

A more involved question is that of stabilizability. Consider a dynamical system with inputs, $x_(t+1)=A x_t+B u_t$. Many systems are not intrinsically stable, but can be stable under correct choice of control inputs. This is what stabilizability captures.

#definition[A linear dynamical system $(A,B)$ with inputs is stabilizable if for all $x_0 in RR^n$, there exist $u_1,u_2,dots in RR^m$ such $x_t->0$ as $t->infinity$.]

The dynamical systems we have considered so far are Markovian; all information at time $t$ is summarized in the state $x_t$. Using this and linear dynamics, one can gain an alternate characterization.

*Observation:* $(A,B)$ is stabilizable iff there exists $K$ such $rho(A+B K)<1$.

The basic idea is to use inputs of the form $u_t=K x_t$ regardless of $x_t$. If so, $x_(t+1)=A x_t+B u_t=(A+B K)x_t=(A+B K)^t x_0$. Thus the dynamical system looks like one without inputs where $A_("eff")=A+B K$. Since such a $K$ exists with $rho(A+B K)<1$ whenever $(A,B)$ is stabilizable, can one find it efficiently?

Note that $min_K rho(A+B K)$ is a non-convex program.

Instead, let's try to use our lemma.

$(A,B)$ is stabilizable
$equiv$ there exists $K$ such $rho(A+B K)<1$
$equiv$ there exists $K$ such there exists $P⪰I$, $P succ (A+B K)^T P(A+B K)$
$equiv$ there exists $K$ such there exists $P⪰I$, $mat(P,((A+B K)P)^T;(A+B K)P,P) succ 0$. (1)

Why? We already know $P succ 0$, thus by Schur complement lemma, satisfaction of the above constraint is same as $P-(A+B K)P P^(-1)P(A+B K)^T succ 0$.

Note this does not match what we had in the last step. But that's fine since what we have as same as $rho((A+B K)^T)<1$, but $(A+B K)^T$ and $A+B K$ have the same characteristic equation and hence eigenvalues. However we still don't have a SDP at this stage due to non-linearity produced by multiplication of $K$ & $P$; both are variables.

$equiv$ there exists $K in RR^(n times m)$, $P⪰I$, $mat(P,(A P+B L)^T;A P+B L,P) succ 0$. (2)

At first glance, this might seem like a relaxation, we used $L=K P$ to get a (near) SDP feasibility program. But this is in fact a reformulation since one can go back (i.e. construct a feasible solution for the original program) by setting $K=L P^(-1)$.

For this program, one does once again perform bisection with $succ$ constraints.

What we were able to accomplish in the last step (an exact reformulation) was a lucky coincidence. As an example, if we want to impose a box constraint on $K$, we can do so in (1) which was (nominally) non-convex anyways, but now an exact reformulation such as (2) as an SDP is no longer possible. In fact, this problem is NP-hard.

= Deriving SDP duals

*Claim:* $min_(x in X) max_(y in Y) f(x,y) >= max_(y in Y) min_(x in X) f(x,y)$.

If min/max do not exist, interpret these as inf/sup for this section and the next.

*Proof:* Clearly, for all $x in X$, $y in Y$, $f(x,y)>=min_(x in X) f(x,y)$. Also, $max_(x in X) g(x)>=max_(x in X) h(x)$ if $g(x)>=h(x)$ for all $x in X$. Using this, we get
for all $x in X$, $max_(y in Y) f(x,y)>=max_(y in Y) min_(x in X) f(x,y)$;
hence $min_(x in X) max_(y in Y) f(x,y)>=max_(y in Y) min_(x in X) f(x,y)$. $square$

Applying this to SDPs:
$min_X upright("Tr")(C X)$ subject to $upright("Tr")(A_i X)=b_i$ for all $i in [m]$, $X⪰0$ (primal)
$=min_(X⪰0) max_(lambda in RR^m) upright("Tr")(C X)+sum_(i=1)^m lambda_i(b_i-upright("Tr")(A_i X))$
$=min_(X⪰0) max_(lambda in RR^m) b^T lambda+upright("Tr")((C-sum_i lambda_i A_i)X)$.
Why? If $b_i=upright("Tr")(A_i X)$, fine. If $b_i>upright("Tr")(A_i X)$, $lambda_i->infinity$. If $b_i<upright("Tr")(A_i X)$, $lambda_i->-infinity$.

$>=max_(lambda in RR^m) min_(X⪰0) b^T lambda+upright("Tr")((C-sum_i lambda_i A_i)X)$
$=max_(lambda in RR^m) cases(b^T lambda & "if " C⪰sum_i lambda_i A_i, -infinity & "otherwise")$
since $upright("Tr")(Y Z)>=0$ for all $Y,Z⪰0$, and one can choose $X=0$. Otherwise there exists $u$ such that $u^T(C-sum_i lambda_i A_i)u<0$ and one can choose $X=alpha u u^T$, $alpha->infinity$.

$=max_(lambda in RR^m) b^T lambda$ subject to $C⪰sum_i lambda_i A_i$ (dual). This is a PSD constraint on an affine map, hence representable as intersection of an affine space and the PSD cone, hence is a SDP.

By our claim, we always used the correct $>=$ relations. Hence, we have weak duality by construction.

#lemma[For any primal feasible $X$, dual feasible $lambda$, we have $upright("Tr")(C X)>=b^T lambda$.]

= Conic duality

Consider $min c^T x$ subject to $A x >=_K b$. Before we get to strong duality, it is worth noting that it's not clear how we even define a dual in the first place. For LPs, our derivation was based on the fact that non-negative combinations of vector inequalities always collapse to valid scalar inequalities. But this is not true for CPs.

For example, $cal(L)^(2+1)={(t,x):t>=norm(x)_2}$. Clearly, $(1,0,-1)>=_(cal(L)^(2+1))(0,0,0)$. But $(0,0,1)^T(1,0,-1)=-1 ≱0=(0,0,1)^T(0,0,0)$. You might guess that the multiplier should come from $K$ itself. This is also false. For example, $K={(x,y):|y|<=100x}$. Then $(1,-100)>=_K(0,0)$. But $(1,100) in K$, yet $angle((1,100),(1,-100))=1-100^2 ≱0=angle((1,100),(0,0))$.

//#image("assets/dual_cone.png", width: 15%)

The true answer follows from a tautology. For all $y in K_*$, $a>=_K b$ $=>$ $y^T a>=y^T b$ for all $a,b in K$.

#definition[The dual cone $K_*$ of a set $K subset.eq RR^n$ is $K_*={y:$ for all $x in K$, $y^T x>=0}$.]

Now, $min_x c^T x$ subject to $A x>=_K b$ is primal. Take $max b^T y$ subject to $A^T y=c$, $y in K_*$ to be the dual.

For any feasible pair $(x,y)$, $c^T x=y^T A x=y^T(A x-b)+b^T y>=b^T y$.

#lemma[For primal feasible $x$, dual feasible $y$, we have $c^T x>=b^T y$.]

#lemma[Let $K$ be a non-empty set in $RR^n$. Then:
+ $K_*$ is a closed convex cone.
+ $upright("int")(K)!=emptyset$ $=>$ $K_*$ is pointed.
+ $K$ is a closed convex pointed cone $=>$ $upright("int")(K_*)!=emptyset$.
+ $K$ is a closed convex cone $=>$ $K_(**) = K$.
]

*Proof:* (1) $K_* = ⋂_(x in K) {y:y^T x>=0}$ is intersection of halfspaces (hence convex) and intersection of closed sets (hence closed). Cone is true from definition, because homogeneous inequalities are involved.

(2) Say $K_*$ is not pointed. Then $K_*$ contains an entire line passing through origin, say ${lambda a:lambda in RR}$. Then $a^T x>=0$ for all $x in K$ and $(-a)^T x>=0$ for all $x in K$. Thus $a^T x=0$. But then $K$ is flat, it's entirely contained in the maximal subspace orthogonal to $a$. Hence, $upright("int")(K)=emptyset$.

*Exercise:* Prove part (3) of the above lemma.

(4) First $K subset.eq K_(**)$. Why? Pick any $x in K$, clearly for all $y in K_*$, $y^T x>=0$. Thus $x in K_(**)$. For the other direction, we'll invoke separation. Let's do the contrapositive. Say $x ∉ K$. There exists $y !=0$ such $y^T x<y^T z$ for all $z in K$. But then $z^T y$ is lower bounded on $K$, hence $min_(z in K) z^T y=0$. Hence $y in K_*$, and from $y^T x<0$, $x ∉ K_(**)$.

#corollary[$K$ is a closed convex pointed cone with non-empty interior iff $K_*$ is a closed convex pointed cone with non-empty interior.]

$min c^T x$ subject to $A x>=_K b$, where $K$ is closed convex pointed cone, $upright("int")(K)!=emptyset$: primal (P).
$max b^T y$ subject to $A^T y=c$, $y in K_*$: dual (D).

#theorem[If P is bounded below and strictly feasible (i.e. $A x>_K b$), then D is solvable and optimal values match.]

By symmetry, if D is bounded above and strictly feasible ($y>_(K_*)0$) then P is solvable and optimal values match.

#corollary[If one of P/D is bounded and strictly feasible, and $(x,y)$ is P/D feasible pair, then $(x,y)$ is P,D optimal iff $b^T y=c^T x$ iff $y^T(A x-b)=0$.]

Under given conditions, we know optimal values are equal. We also know from weak duality that for any pair of feasible solutions $b^T y<=D^*<=P^*<=c^T x$. Thus attainment of optimal values at $(x,y)$ for $(P,D)$ happens iff $b^T y=c^T x$. From the proof of weak duality this is equivalent to $y^T(A x-b)=0$, the only inequality being tight.

*Proof:* Let $P^*$ be the optimal value for the primal problem. We will construct a dual feasible $y^*$ such $b^T y^*>=P^*$, via the separation theorem.

Let $V={A x-b:c^T x<=P^*}$. It is non-empty, because any affine/linear map of a halfspace can not be empty. Also $V inter upright("int")(K)=emptyset$, otherwise given $x$ such $c^T x<=P^*$ and $A x>_K b$, one can strictly decrease the objective for a close enough feasible point, assuming $c!=0$.

Aside note: $c!=0$ holds generically. If not, $P^*=0$. Take $y=0$, $A^T y=c$ holds and we have $0=b^T y>=P^*$.

Now, by (actual) separation theorem, there exists $y!=0 in RR^m$ such that $sup_(x in V) y^T x<=inf_(x in K) y^T x$. So, $y^T x$ is bounded below on $K$. But $inf_(x in K) y^T x$ generically only take values ${0,-infinity}$. Thus $sup_(x in V) y^T x<=0<=inf_(x in K) y^T x$, and $y in K_*$ follows from the second inequality. But also the first of these is equivalent to $c^T x<=P^* => (A^T y)^T x<=y^T b$. This can only happen if $A^T y=alpha c$ for some $alpha>=0$. If $alpha>0$, then take $y^*=y/alpha in K_*$ for which $A^T y^*=c$ and $b^T y^*>=(A^T y)^T x/alpha=c^T x$ for all $c^T x<=P^*$. We are done. If $alpha=0$, $A^T y=0$ whereas $b^T y>0$ for there exists $x$, $A x>_K b$. $y in K_*$ is not null, and hence $y^T(A x-b)>0$ $=>$ $c^T x<=P^*$ yet $y^T A x=0>b^T y$. This is a contradiction because then $c^T x<=P^*$ has to empty. $square$

= S-lemma

Consider a non-convex QCQP $min_x q_a(x)=x^T Q_a x+q_a^T x+c_a$ subject to $q_b(x)=x^T Q_b x+q_b^T x+c_b>=0$. $Q_a,Q_b$ are arbitrary, need not be $Q_a⪰0$, $Q_b⪰0$.

Such programs are efficiently solvable, despite apparent non-convexity.

#theorem[Let there exist $x$ such $q_b(x)>0$ strictly. Then for all $x$, $q_b(x)>=0$ $=>$ $q_a(x)>=0$ implies there exists $lambda>=0$ such that for all $x$, $q_a(x)>=lambda q_b(x)$.]

$<=$ above is evident. Such $lambda$ “certifies” the former clause. This has the same spirit as Farkas' lemma or, theorem of alternatives you might have seen before:
$min_x q_a(x)$ subject to $q_b(x)>=0$
$=max_gamma gamma$ subject to for all $x$, $q_b(x)>=0$ $=>$ $q_a(x)-gamma>=0$
$=max_(gamma,lambda>=0) gamma$ subject to for all $x$, $q_a(x)-gamma>=lambda q_b(x)$ (using S-lemma)
$=max_(gamma,lambda>=0) gamma$ subject to for all $x$, $(x,1)^T mat(Q_a-lambda Q_b,(q_a-lambda q_b)/2;(q_a-lambda q_b)^T/2,c_a-gamma-lambda c_b)(x,1)>=0$
$=max_(gamma,lambda>=0) gamma$ subject to $mat(Q_a-lambda Q_b,(q_a-lambda q_b)/2;(q_a-lambda q_b)^T/2,c_a-gamma-lambda c_b)⪰0$. (This is an SDP.)

This needs a brief justification. $=>$ is clear. You might reasonably suspect we are asking for more than is necessary, e.g., $(x,1)^T M(x,1)>=0$ also, in the latter formulation. To complete the argument, say there exists $y$, $y^T M y<0$. Then,
(1) last coordinate of $y!=0$, then by rescaling there exists $x$ such that $(x,1)^T M(x,1)<0$.
(2) last coordinate of $y=0$, since $y^T M y$ is continuous, there exists $z$ with non-zero last coordinate with $z^T M z<0$.

Finally, we'll prove S-lemma. Let's start with something simpler.

#lemma[For symmetric $A,B$, $upright("Tr")(A)>=0>upright("Tr")(B)$ $=>$ there exists $x$ such $x^T A x>=0$ and $x^T B x<0$.]

*Proof:* Two observations: for any $x in {-1,1}^n$, $x^T D x=upright("Tr")(D)$ for any diagonal $D$, since $x_i^2=1$. Secondly, for any matrix $M$, $EE_(x ∼ upright("Unif")({-1,1}^n)) x^T M x=upright("Tr")(M)$.

Now choose $xi ∼ upright("Unif")({-1,1}^n)$, $x=U xi$, where $B=U Lambda U^T$ diagonal. $x^T B x=xi^T U^T U Lambda U^T U xi=upright("Tr")(Lambda)=upright("Tr")(Lambda U^T U)=upright("Tr")(B)<0$, yet $EE x^T A x=EE xi^T(U^T A U)xi=upright("Tr")(U^T A U)=upright("Tr")(A U U^T)=upright("Tr")(A)>=0$.

Therefore, there exists $x$ such $x^T A x>=0$ and $x^T B x<0$. $square$

#theorem[Say there exists $macron(x)$ such $macron(x)^T A macron(x)>0$. Then for all $x$, $x^T A x>=0$ $=>$ $x^T B x>=0$ implies there exists $lambda>=0$ such that $B⪰lambda A$.]

*Proof:* Consider $OPT=min_(x^T A x>=0) x^T B x$. Say for all $x$, $x^T A x>=0$ $=>$ $x^T B x>=0$ and there exists $macron(x)$, $macron(x)^T A macron(x)>0$. Then $OPT=0$. Consider $OPT'=min_(x^T A x>=0, x^T x=n) x^T B x$. $OPT'$ is also strictly feasible; $OPT'>=0$.
$equiv min_X upright("Tr")(B X)$ subject to $upright("Tr")(A X)>=0$, $upright("Tr")(X)=n$, $X⪰0$, $X$ is rank 1.
$>= min_X upright("Tr")(B X)$ subject to $upright("Tr")(A X)>=0$, $upright("Tr")(X)=n$, $X⪰0$ (SDP).

SDP dual is $max_(u,lambda>=0) n u$ subject to $B⪰u I+lambda A$.

*Exercise:* Show that SDP dual is indeed the dual of the SDP.

Let's say we guarantee that SDP $>=0$ and that strong duality holds. Then, we are done as then there exist $u,lambda>=0$ such $B⪰u I+lambda A$ and hence there exists $lambda>=0$, $B⪰lambda A$.

Take $macron(X)=(macron(x)macron(x)^T+alpha I)/upright("Tr")(macron(x)macron(x)^T+alpha I) n$, for small enough $alpha>0$, is strictly feasible for SDP. And for all $lambda>0$, there exists $u$ such $B succ u I+lambda A$. Thus the SDP dual is also strictly feasible and strong duality holds.

*Exercise:* In words, describe why adding $x^T x=n$ as a constraint was essential for this argument.

Now to argue SDP $>=0$, observe: Compactness ensures that it is solvable; let $X^*$ be the solution. There exists $L$, $X=L L^T$. Consider that $upright("Tr")(A X^*)=upright("Tr")(L^T A L)>=0$, and $upright("Tr")(B X^*)=upright("Tr")(L^T B L)=n u^*$. If $u^*<0$, then by our sign lemma, there exists $x$ such $(L x)^T A(L x)>=0$ and $(L x)^T B(L x)<0$. But this contradicts for all $x$, $x^T A x>=0$ $=>$ $x^T B x>=0$. $square$

*Proof of the original theorem.* The idea is to use homogeneous polynomials.
$macron(q)_a(x,t)=x^T Q_a x+q_a^T x t+c_a t^2$.
$macron(q)_b(x,t)=x^T Q_b x+q_b^T x t+c_b t^2$.

Given there exists $x$ with $q_b(x)>0$, we have $macron(q)_b(x,1)>0$.

We'll prove for all $(x,t)$, $macron(q)_b(x,t)>=0$ $=>$ $macron(q)_a(x,t)>=0$. If so, using our simplified S-lemma, there exists $lambda>=0$, $macron(q)_a(x,t)>=lambda macron(q)_b(x,t)$ for all $(x,t)$. Set $t=1$.

Let's say our sub-claim is untrue. Then, there exists $(x,t)$ such $macron(q)_b(x,t)>=0$ but $macron(q)_a(x,t)<0$. Consider
(1) $t!=0$. Then substituting $(x/t,1)$ results in contradiction.
(2) $t=0$, $macron(q)_b(x,t)>0$, by continuity there exists $(y,s)$ with $s!=0$ with $macron(q)_b(y,s)>=0$. Case (1) applies.

*Exercise:* Argue why the remaining case $t=0$, $macron(q)_b(x,t)=0$ results in a contradiction too, completing the proof. $square$
