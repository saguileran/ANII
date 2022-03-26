### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ 054454a0-ad17-11ec-3385-29ec65f21c33
begin
	using PlutoUI
	TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)
end

# ╔═╡ 72ea7f70-ad12-11ec-32b3-7fcfb7de94a9
md"
# Homework 2 

"


# ╔═╡ b52c7910-ad12-11ec-3836-73322a4d09bf
md"
## 1.T
Use the Galerkin method with the functions ${\phi_1, \phi_2} = \{x^2-2x, x(x-1)^2\}$ to get an aproximated solution (two parameters) of the BVP solution

```math
\frac{d^2 u}{dx^2} + u = 1; \qquad 0<x<1\qquad \text{with} \quad u(0)=0, \qquad \frac{du}{dx}=0 \quad \text{in} \quad x=1
```
Find the residual error $\epsilon(x)= d^2u_n/dx^2+u_n-1$. Plot in a single plane the aproximate solution $u_n$, the exact solution $u=-\cos(x)-\tan(1)\sin(x)+1$ and the residual error $\epsilon(x)$
"

# ╔═╡ fe3a24d0-ad13-11ec-127e-b99341ce2d42
md"
## 2.T 
Prove that $\ddot{u}$ is a solution of the variational problem. Find $u\in H_0^1(I)$ such that

```math
\int_I u'v' dx = \in_I fv, \qquad \forall v\in H_0^1(I)
```
if and only if $\ddot{u}$ is a functional minimizer.

```math
F(w) = \frac{1}{2} \int_I (w')^2 - \int_I fw dx
```
over the space $H_0^1(I)$
"

# ╔═╡ fec1f270-ad13-11ec-1dc6-47c977b59d4d
md"""
## 3.TP

Consider the problem with Robin conditions

```math
-u''+2u=4+2x^3, \qquad x\in I = (-1,1)
```
```math
u'(-1)=\frac{\alpha}{\kappa}(u(-1)+1)
```
```math
u'(1) = -\frac{\alpha}{\kappa} {\kappa}(u(1)+1)
```

a) Write the weak formulation correspondingt to the boundary value problem $(e)-(3)$ 

b) Prove tue existence and uniquiness of the weak solution.

c) To discretize the weak problem consider the picewise lineal polynomial spaces and write exactly the rigid matrix, load mass and vector.

d) Using the result find vefore, part c., do the modifications to the Matlab codes and solve the problem taking $\alpha =2$ and $\kappa=1$. Plot the numerical solution and compare it with the exact solution $u(x)=1-x^2$.

e) Solve the equation with Dirichlet homogeneous boundary condtions. Plot and compare it with the exact solution $u(x)=1-x^2$.

f) Plot the error made in each mesh point $\left( |u(x_i)-u_h(x_i)|,i=0,1,...,m  \right)$ of the finded solutions in the items d) and 
"""


# ╔═╡ 0c0cd9d0-ad15-11ec-38b8-596132c52440
md"
## 4.TP 

**Bar Deflection**

Consider a beam fix in a wall in both extrems like in the figure. The strcutre, of lenght $L$ is hold in a distributed load $P [Kg \; m^{-1}$ that vary along the $x$ coordinate. Furthermore, supposse that the beam has a rectangular uniform section of width $r$ and depth $s$, inertial moment $J = rs^3/12 [m^4]$ and Young module $E[kgm^{-2}]$.

The beam deflection, suposing small displacements, is describe by the 4th order differential equation
```math
(EJu''(x))'' = P(x), \qquad \qquad u'(0)=u'(L)=0,
```
where the $u=u(x)$ denote the vertical displacement. The following boundary conditions


model the fix beam behaviour in both extrems.

a. Write the weak formulation corresponding to the boundary value problem $(4)-(5)
b. Prove the existence and uniqueness of the weak solution.

To discretize the weak problem we shuold aproximate the space $H_0^2([0,1])$ by a picewise finite dimensional cubic polynomial space

```math
\mathcal{H}_h := \{ \phi \in C^1([0,L]):\quad \phi\Big|_{I_i} \in \mathbb{P}_3(I_i) \quad i=1,2,...,N; 
```
```math
\phi(0)= \phi(1)=0 \quad \text{and} \quad  \phi_x(0)=\phi_x(1)=0 \}
```

where $0 <c_0<x_1<\cdots < x_N =L, \; I_i = (x_{i-1},x_), \; h_i=x_i-x_{i-1}, \; h:=\text{max}_{i=1,...,N}$ and $\mathbb{P}_3(I_i)$ denote the set of all the polynomials over $I_i$ of order less or equal to $3$.

With the goal to build a base for $\mathcal{H}_h$, let's asociate a each inner node $x_k,\; k= 1,...,N-1$ a support $I_k \cup I_{k+1}$ and two functions $\zeta_i$ and $\eta_i$ such that

```math
\zeta_i(x_k) = \delta_{ik}, \qquad \zeta_i'(x_k)=0
```
```math
\eta_i(x_k) = 0, \qquad \eta_i'(x_k)=\delta_{ik}
```

c. Use the previous conditions to determine the $\zeta_i$ and $\eta_i$ explicit formulas and plot them.

d. Write $u_h$ as a linear combination of the base $\zeta_i$ and $\eta_i$ where $i=1,...,N-1$.

e. Consider the discrete weak formulation. Indicate wich are the unknowns.

f. Find the matrix components of the lineal equation system that gives the finite elemebt discretization.

g. Solve the lineal equation systems useing for $N=25,\;50,\;100$ and $200$ elements. Plot the solution. use the folowwing dates: $EJ =1$, $P(x)=1$, and $L=1$.
"

# ╔═╡ feda3562-ad13-11ec-3996-df8b9c71b329
md"

"

# ╔═╡ fef02e60-ad13-11ec-01c9-176af939a84e
md"
```math
```

"

# ╔═╡ ff04eee0-ad13-11ec-2231-99e0a6f74daf


# ╔═╡ Cell order:
# ╠═054454a0-ad17-11ec-3385-29ec65f21c33
# ╟─72ea7f70-ad12-11ec-32b3-7fcfb7de94a9
# ╟─b52c7910-ad12-11ec-3836-73322a4d09bf
# ╟─fe3a24d0-ad13-11ec-127e-b99341ce2d42
# ╟─fec1f270-ad13-11ec-1dc6-47c977b59d4d
# ╟─0c0cd9d0-ad15-11ec-38b8-596132c52440
# ╠═feda3562-ad13-11ec-3996-df8b9c71b329
# ╠═fef02e60-ad13-11ec-01c9-176af939a84e
# ╠═ff04eee0-ad13-11ec-2231-99e0a6f74daf
