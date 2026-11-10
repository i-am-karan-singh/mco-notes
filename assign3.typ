#import "47851.typ": *
#show: lecture.with([Assignment 3], 8)
#set heading(numbering: "Q1.")

= Grads and Dists (10 points)

*Part A:* Consider a function $f : RR^d -> RR$ that is $mu$-strongly convex and $L$-smooth. Let $x^* in arg min_(x) f(x)$. Prove for all $x in RR^d$ that
$ 1/(2L) norm(nabla f(x))^2<= f(x) - f(x^*) <= 1/(2mu) norm(nabla f(x))^2. $

*Part B.* Prove that $f$ is $mu$-strongly convex if and only if for all $x, y in RR^d$, we have
$ (nabla f(x) - nabla f(y))^top (x-y) >= mu norm(x-y)^2. $

= Just Enough Convexity (10 points)

Consider a $L$-smooth function $f : RR^d -> RR$ that satisfies for all $x$ that
$ f(x) - min_(x^*) f(x^*) <= 1/(2mu) norm(nabla f(x))^2. $

*Part A:* Prove that such a function may not even be convex (much less strongly convex).

*Part B:* Nevertheless, prove that gradient descent on such a function, with an appropriate step size, produces iterates so that
$ f(x_t) - min_(x^*) f(x^*) <= O(e^(-mu t/L)). $

= Small Grad Bounds (10 points)

*Part A:* Fix a convex set $K subset RR^d$ of diameter $D$. Recall that online gradient guarantees the following regret bound against a sequence of adversarially chosen convex $G$-Lipschitz functions $f_t$.
$ upright("Regret") = sum_(t=1)^T f_t (x_t) - min_(x^* in K) sum_(t=1)^T f_t (x^*) <= O(G D sqrt(T)), $
With suitable dynamically chosen step sizes, show that one can guarantee that
$ upright("Regret") <= O(D sqrt(sum_(t=1)^T norm(nabla f_t (x_t))^2)). $

_Hint:_ Set $eta_t prop 1/sqrt(sum_(s=1)^t norm(nabla f_s (x_s))^2)$.

*Part B:* Assume that all loss functions $f_t:RR^d -> RR_(>= 0)$ are in addition non-negative, well-defined and $L$-smooth over all of $RR^d$. Our online algorithm must still play points in $K$ and the regret is still defined with respect to the best point in $K$. Prove that the same algorithm guarantees
$ upright("Regret") <= O( sqrt(L D^2 (min_(x^* in K)sum_(t=1)^T f_t (x^*))) + L D^2). $
This means whenever there is a point with zero net loss, we get constant regret, independent of $T$.

_Hint:_ Use part A from the first question.
