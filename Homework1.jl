### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 326d71c0-a9b9-11ec-3986-c9776d176014
begin
	using PlutoUI, QuadGK, Calculus, Roots, LaTeXStrings, PlotlyJS, ImageShow, ImageIO, Images
	TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)
end

# ╔═╡ 2b55a600-a9b9-11ec-3eb1-1f25b963a10f
html"
<center> <h1> Homework #1 <br> 2022</h1> </center>
"

# ╔═╡ 4e476b30-a9b9-11ec-2d51-2fa7754c2747
md"
## 1.
asdasd
"

# ╔═╡ f7dc9db0-a9b8-11ec-15a1-ab2e399d50f1
md"
## 2.  

Let $u(x)$ and $v(x)$ two functions defined as follows

```math 
u(x)=\left\{\begin{array}{ll}
x, & \text { if } 0 \leq x \leq 1 \\
2-x & \text { if } 1 \leq x \leq 2
\end{array} \quad y \quad v(x)=\sin \pi x\right.
```

with their derivatives

```math
u'(x)=\left\{\begin{array}{ll}
1, & \text { if } 0 \leq x \leq 1 \\
-1 & \text { if } 1 \leq x \leq 2
\end{array} \quad y \quad v'(x) = \pi\cos \pi x\right.
```

"

# ╔═╡ 06249f30-a9b9-11ec-1093-a59349277d1a
begin
	function u(x)
		if 0 <= x <=1
			return x
		elseif 1 < x <=2
			return 2 - x
		else 
			return 0
		end
	end

	function up(x)
		if 0 <= x <=1
			return 1
		elseif 1 < x <=2
			return -1
		else 
			return 0
		end
	end

	v(x)  = sin(pi * x)
	vp(x) = pi * cos(pi * x)

	x = range(0, 2, step=0.1) |> collect;
end

# ╔═╡ 74ae8100-a9b9-11ec-1701-1fc4d58c62ad
begin
	plot_u  = scatter(;x=x, y=u.(x), mode="lines+markers", name="u(x)")
	plot_v  = scatter(;x=x, y=v.(x), mode="lines+markers", name="v(x)")
	plot_uv = scatter(;x=x, y=v.(x).*u.(x), mode="lines+markers", fill="tozeroy", name="u(x) * v(x)")

	plot_up = scatter(;x=x, y=up.(x), mode="lines+markers", name="u'(x)")
	plot_vp = scatter(;x=x, y=vp.(x), mode="lines+markers", name="v'(x)")
	plot_uvp = scatter(;x=x, y=vp.(x).*up.(x), mode="lines+markers", fill="tozeroy", name="u'(x) * v'(x)")

	p1 = plot([plot_u, plot_v,plot_uv], Layout(title="Functions",  xaxis_title="x"))
	p2 = plot([plot_up, plot_vp,plot_uvp], Layout(title="Function Derivative",  xaxis_title="x"))

	p_1 = [p1 p2]
	#relayout!(p, height=300, width=700, title_text="Functions", legend_title_text="Legend")
	#p
end

# ╔═╡ 7c7f3c80-a9b9-11ec-1732-4d03f3b17dd6
md"
### Orthogonality 


#### - $L^2(0,2)$ Space

Norm of $L^2$ and its norm $||u||_{L^2}$

```math 
(u,v) =  \int_0^2 u(x) v(x) dx  = \int_0^1 x \; \sin(\pi x) dx + \int_1^2 (2-x)\; \sin(\pi x) dx 
```

Integranting by steps with the sustitution $w = x$ and $dz = \sin (\pi x) \; dx$

```math
(u,v)_{L^2} = \frac{1}{\pi}x\cos(\pi x)\Big|_0^1 - \frac{1}{\pi}\int_0^1 \; \cos(\pi x) dx + \frac{1}{\pi}(2-x)\cos(\pi x)\Big|_1^2 - \frac{1}{\pi} \int_1^2 \; \cos(\pi x) dx \\
(u,v)_{L^2} = \frac{1}{\pi}\cos(\pi) - \frac{1}{\pi^2}\sin(\pi x)\Big|_0^1 - \frac{1}{\pi}\cos(\pi) - \frac{1}{\pi^2}  \sin(\pi x)\Big|_0^1  \\
(u,v)_{L^2} = \frac{1}{\pi}\cos(\pi) - \frac{1}{\pi}\cos(\pi) -  \frac{1}{\pi^2} ( \sin(\pi) - \sin(0) ) - \frac{1}{\pi^2} ( \sin(2\pi) - \sin(\pi) ) = \frac{1}{\pi} - \frac{1}{\pi} -(0-0) -(0-0) \\
\Longrightarrow  \quad \boxed{(u,v)_{L^2} =  0}
```
The function $u$ and $v$ are orthogonals in $L^2(0,2)$.  
"

# ╔═╡ ad4e5a80-a9b9-11ec-2b9f-8b076b06cb32
integral_1, err_1 = quadgk(x -> u(x)*v(x), 0, 2, rtol=1e-8)

# ╔═╡ b0f1c000-a9b9-11ec-32a0-81846014b48d
md"
#### - $H^1(0,2)$ Space

Now, let's check if $u$ and $v$ are orthogonal in $H^1(0,2)$ usign the its norm and their derivatives

```math 
(u,v)_{H^1} =  \int_0^2 u(x) v(x) dx  + \int_0^2 u'(x) v'(x) dx = 0 + \int_0^1 u'(x) v'(x) dx  + \int_1^2 u'(x) v'(x) dx  \\
(u,v)_{H^1} = \pi\int_0^1 1 \; \cos(\pi x) dx + \pi\int_1^2 (-1)\; \cos(\pi x) dx = \sin(\pi x) \Big|_0^1 - \sin(\pi x)\Big|_1^2\\ 
\Longrightarrow \quad \boxed{(u,v)_{H^1} = 0}
```

Then, the functions $u$ and $v$ are orthogonal on $H^1(0,2)$.
"

# ╔═╡ c50f4940-a9b9-11ec-07df-dd60bccf4f04
integral_2, err_2 = quadgk(x -> up(x)*vp(x), 0, 2, rtol=1e-8)

# ╔═╡ c5f4ee50-a9b9-11ec-07ad-8ff5b133082b
md"
### $|| u-v ||$ Difference

Let's calculate the distance between $u$ and $v$ in both spaces $L^2(0,2)$ and $H^1(0,2)$.
"

# ╔═╡ d5f6adc0-a9b9-11ec-3e41-e7cedf88a606
md"
#### -  $L^2(0,2)$ Space

```math 
||u-v||_{L^2}^2 =  \int_0^2 |u-v|^2 dx =  \int_0^2 (u-v)^2 dx\\
||u-v||_{L^2}^2 =  \int_0^1 (x - \sin(\pi x))^2 dx + \int_1^2 (2-x-  \sin(\pi x))^2 dx  \\
||u-v||_{L^2}^2 =  \int_0^1 \left(x^2 - 2x\sin(\pi x) + \sin^2(\pi x)\right) dx + \int_1^2 \left(4 + x^2 + \sin^2(\pi x) -4x -4 \sin(\pi x) + 2x\sin(\pi x) \right) dx  
```

Using integral tables

```math
||u-v||_{L^2}^2 =   \left(\frac{1}{3}x^3 - \frac{2}{\pi^2}\sin(\pi x) + \frac{2}{\pi} x\cos(\pi x) + \frac{x}{2}  - \frac{1}{4\pi} \sin(2\pi x)\right)\Big|_0^1 + \left(4x + \frac{1}{3}x^3 + \frac{x}{2}  - \frac{1}{4\pi} \sin(2\pi x)  -2x^2 + \frac{4}{\pi} \cos(\pi x) + \frac{2}{\pi^2}\sin(\pi x) - \frac{2}{\pi} x\cos(\pi x)  \right)\Big|_1^2  
```

```math
||u-v||_{L^2}^2 =   \left(\frac{1^3}{3} + \frac{2}{\pi} \cos(\pi) + \frac{1}{2} \right)+ \left(4(2) + \frac{2^3}{3} + \frac{2}{2} -(2)2^2 +\frac{4}{\pi} \cos(2\pi)  - \frac{4}{\pi} \cos(2\pi)  \right) -  \left(4 + \frac{1^3}{3} + \frac{1}{2}  -(2)1^2 +\frac{4}{\pi} \cos(\pi ) - \frac{2}{\pi} \cos(\pi)  \right) \\
||u-v||_{L^2}^2 =   \left(\frac{1}{3} - \frac{2}{\pi}  + \frac{1}{2} \right) + \left(8 + \frac{2^3}{3} + 1 - 8 \right) -  \left(4 + \frac{1}{3} + \frac{1}{2} -2 - \frac{4}{\pi} + \frac{2}{\pi}  \right) \\
||u-v||_{L^2}^2 =  \left(\frac{5}{6} - \frac{2}{\pi} \right) + \left( \frac{5}{6} - \frac{2}{\pi}\right) = \frac{5}{3} \qquad 
\Longrightarrow \qquad \boxed{||u-v||_{L^2}^2 =  1.67}
```
"

# ╔═╡ 0b8b9130-a9ba-11ec-2207-774f008ffed7
md"
#### - $H^1(0,2)$ Space

```math
||u-v||_{H^1}^2 =  ||u-v||_{L^2}^2 + ||(u-v)'||_{L^2}^2 =  \int_0^2 |u-v|^2 dx  + \int_0^2 |(u-v)'|^2 dx 
```

The first integral is already done, the second is

```math
||(u-v)'||_{L^2}^2 = \int_0^2 |(u-v)'|^2 dx = \int_0^1 (1 - \pi\cos(\pi x))^2 dx + \int_1^2 (-1 -\pi \cos(\pi x))^2 dx \\
\int_0^2 |(u-v)'|^2 dx = \int_0^1 \left(1 - 2\pi\cos(\pi x) + \pi^2\cos^2(\pi x)\right) dx + \int_1^2 \left(1 +2\pi \cos(\pi x)+ \pi^2 \cos^2(\pi x)\right) dx \\
||(u-v)'||_{L^2}^2 = \int_0^2 |(u-v)'|^2 dx =  \left(x - 2\sin(\pi x) + \frac{\pi^2}{2}x + \frac{\pi}{4}\sin(2\pi x) \right)\Big|_0^1  + \left(x +2\sin(\pi x)+  \frac{\pi^2}{2}x + \frac{\pi}{4}\sin(2\pi x)\right)\Big|_1^2 \\
||(u-v)'||_{L^2}^2 = \int_0^2 |(u-v)'|^2 dx =  \left(1 + \frac{\pi^2}{2} \right)  + \left(1 + \frac{\pi^2}{2}\right) =  \pi^2 +2 \qquad \Longrightarrow \qquad \boxed{||(u-v)'||_{L^2}^2  = 11.8696}
```
"

# ╔═╡ 200e1dd0-a9ba-11ec-13b5-9f9c9bad7ea8
md"
Finally, the distance between $u$ and $v$ in $H^1(0,2)$ is 
```math
||u-v||_{H^1}^2 =  \frac{5}{3}+ \pi^2 +2  = \frac{11}{3} + \pi^2 \qquad \Longrightarrow \qquad \boxed{||u-v||_{H^1} = 3.6792}
```
"

# ╔═╡ 4210ef1e-a9ba-11ec-3e6b-9314377e22aa
umvH = sqrt(integral1 + integral2)

# ╔═╡ 45897410-a9ba-11ec-0ebb-e5c07ed91ae8
md"
## 3.
"

# ╔═╡ 5795f7a0-a9ba-11ec-015f-834a762cdcdf
md"
Now, let's verify the Cauchy-Schwarz inequality in the space $L^2(0,2)$ for the previous functions, $u$ and $u-v$ 

```math
||(u, u-v)||_{L^2} \leq ||u||_{L^2} \cdot |u-v||_{L^2}  
```

Notice that only the last term is calculated. First let's calculate the norm of the inner product between $u$ and $u-v$

```math
||(u, u-v)||_{L^2}^2 = \int_0^2 |u(u-v)|^2 dx = \int_0^1 (x(x-\sin(\pi x)))^2 dx + \int_1^2 ((2-x)(2-x-\sin(\pi x)))^2 dx \\
||(u, u-v)||_{L^2}^2 = \int_0^1 \left(x^4+x^2\sin^2(\pi x) - 2x^3\sin(\pi x)\right) dx + \int_1^2 \left((2-x)^4+(2-x)^2 \sin^2(\pi x)-2(2-x)^3\sin(\pi x) \right) dx \\
||(u, u-v)||_{L^2}^2 = \int_0^1 x^4 dx + \frac{1}{2}\int_0^1 x^2 \left(1-\cos(2\pi x) \right) dx -2 \int_0^1 x^3\sin(\pi x)dx + \int_1^2 (2-x)^4 dx  + \frac{1}{2}\int_1^2 (2-x)^2 \left(1-\cos(2\pi x) \right) dx -2 \int_1^2 (2-x)^3\sin(\pi x) dx \\
||(u, u-v)||_{L^2}^2 = \int_0^1 x^4 dx + \frac{1}{2}\int_0^1 x^2 dx - \frac{1}{2}\int_0^1 x^2 \cos(2\pi x) dx -2 \int_0^1 x^3\sin(\pi x)dx + \int_1^2 (2-x)^4 dx  + \frac{1}{2}\int_1^2 (2-x)^2  dx - \frac{1}{2}\int_1^2 (2-x)^2 \cos(2\pi x)  dx -2 \int_1^2 (2-x)^3\sin(\pi x) dx \\
||(u, u-v)||_{L^2}^2 = \frac{1}{5} + \frac{1}{6} - \frac{1}{2} \frac{1}{2\pi^2} - \frac{2}{\pi} + \frac{12}{\pi^3} + \frac{1}{5} + \frac{1}{6} - \frac{1}{2} \frac{1}{2\pi^2} - \frac{12}{\pi^3} + \frac{2}{\pi} = \frac{2}{5} + \frac{1}{3} - \frac{1}{2\pi^2} = \frac{11}{15}  - \frac{1}{2\pi^2} \qquad \Longrightarrow \qquad \boxed{||(u, u-v)||_{L^2}^2 = 0.683}
```

Now, let's calculate the norm of $u$

```math
||u||_{L^2}^2 = \int_0^2 |u|^2 dx =  \int_0^1 x^2 dx  + \int_1^2 (2-x)^2 dx = \frac{1}{3}x^3\Big|_0^1 - \frac{1}{3} (2-x)^3\Big|_1^2 = \frac{2}{3} 
```


To check the Cauchy-Schwarz inequality let's use it to the power of two 

```math
||(u, u-v)||_{L^2}^2 \leq ||u||_{L^2}^2 \cdot ||u-v||_{L^2}^2  \qquad \Longleftrightarrow \qquad  \boxed{0.683 \leq \frac{2}{3}(1.67) = 1.1111}
```
"

# ╔═╡ 7ca852e0-a9ba-11ec-0eef-c385a17b5f74
sqrt(integraluumv) <= sqrt(integralu2) *  sqrt(umvH)

# ╔═╡ 823f23f0-a9ba-11ec-3747-0b794a9a3151
md"
Now, let's verify the Cauchy-Schwarz inequality in the space $H^1(0,2)$ for the previous functions, $u$ and $u-v$ 

```math
||(u, u-v)||_{H^1} \leq ||u||_{H^1} \cdot ||u-v||_{H^1} = \left( ||u||_{L^2}^2  + ||u'||_{L^2}^2 \right)^{1/2} \left( ||u-v||_{L^2}^2  + ||(u-v)'||_{L^2}^2 \right)^{1/2}
```

Notice that only the last term is calculated. Let's start calculating the norm of the inner product between $u$ and $u-v$

```math
||(u, u-v)||_{H^1}^2 = ||(u, u-v)||_{L^2}^2 +||(u', u'-v')||_{L^2}^2 = \int_0^2 |u (u-v)|^2 dx + \int_0^2 |u' (u'-v')|^2 dx \\
```

Let's start calculating the inner product between the derivatives of $u$ and $u-v$

```math
\int_0^2 |u' (u'-v')|^2 dx = \int_0^1(1- \pi\cos(\pi x))^2 dx + \int_1^2 (1+\pi\cos(\pi x)) ^2dx =  \int_0^1(1- 2\pi\cos(\pi x) + \pi^2\cos^2(\pi x)) dx + \int_1^2 (1+2\pi\cos(\pi x)+\pi^2\cos^2(\pi x))dx 
```
```math
\int_0^2 |u' (u'-v')|^2 dx = \int_0^1 dx -2\pi \int_0^1 \cos(\pi x) dx + \frac{\pi^2}{2}\int_0^1(1+\cos(2\pi x)) dx + \int_1^2 dx + 2\pi\int_1^2 \cos(\pi x) dx + \frac{\pi^2}{2}\int_1^2 (1+\cos(2\pi x)) dx 
```
```math
\int_0^2 |u' (u'-v')|^2 dx = \int_0^1 dx -2\pi \int_0^1 \cos(\pi x) dx + \frac{\pi^2}{2}\int_0^1 dx + \frac{\pi^2}{2}\int_0^1\cos(2\pi x) dx + \int_1^2 dx + 2\pi\int_1^2 \cos(\pi x) dx + \frac{\pi^2}{2}\int_1^2 dx + \frac{\pi^2}{2}\int_1^2 \cos(2\pi x) dx 
```
```math
\int_0^2 |u' (u'-v')|^2 dx = 1 + \frac{\pi^2}{2} + 1 + \frac{\pi^2}{2} = 2+\pi^2 \qquad \Longrightarrow \qquad \boxed{\int_0^2 |u' (u'-v')|^2 dx = 11.87}
```
Then,

```math
||(u, u-v)||_{H^1}^2 =  \frac{11}{15}  - \frac{1}{2\pi^2} + 2 + \pi^2 = \frac{41}{15} - \frac{1}{2\pi^2} + \pi^2 \qquad \Longrightarrow \qquad \boxed{||(u, u-v)||_{H^1}^2 = 12.552 }
```
Now let's calculate the $u$ norm in $H^1$ 
```math
||u||_{H^|}^2 =  ||u||_{L^2}^2  + ||u'||_{L^2}^2 = \frac{2}{3} + \int_0^2 u'^2 dx = \frac{2}{3} + \int_0^1 dx + \int_1^2 dx =  \frac{2}{3} +2 = \frac{8}{3}
```
Finally, checking the Cauchy-Schwarz inequality to the power of two 
```math
||(u, u-v)||_{H^1}^2 \leq ||u||_{H^1}^2 \cdot ||u-v||_{H^1}^2  \qquad \Longleftrightarrow \qquad \boxed{12.552 \leq \frac{8}{3}\cdot (13.536) = 36.0.97}
```

"

# ╔═╡ 0e849070-a9bb-11ec-35c7-99a985571bbb
sqrt(integraluumv + integralupupmp) <= sqrt(integralu2 + 2) * umvH

# ╔═╡ 11431430-a9bb-11ec-1174-6df0edc48380
begin
	plot_u2  = scatter(;x=x, y=broadcast(abs, u.(x).*u.(x)), mode="lines+markers", fill="tozeroy", name="|u|^2")
	plot_up2 = scatter(;x=x, y=broadcast(abs, up.(x).*up.(x)), mode="lines+markers", fill="tozeroy", name="|u'|^2")

	plot([plot_u2, plot_up2], Layout(title="Functions",  xaxis_title="x"))
end

# ╔═╡ 2dedf230-a9bb-11ec-1474-ef4bd3035175
md"
## 4. 

Let $f(x,y)$ a continous function over $\mathbb{R}^2$ such that
```math
\iint_{\Omega} f(x,y) dx fy = 0
```
for all rectangle $\Omega \subseteq \mathbb{R}^2$, then $f(x,y)=0$ for all $(x,y)$.

Since $f$ is continous in $\Omega$, it satisfyies that
```math
\lim_{(x,y) \to (x_0, y_0)} f(x,y) = f(x_0,y_0) = L
```
for all $(x_0,y_0) \in \Omega$. It can be rewriting using  the epsilo-delta definition.

```math
\forall \epsilon > 0, \exists \delta>0  \text{, whenever }  0 < |(x,y)-(x_{0},y_0)| < \delta,  \; \Rightarrow \; |f(x_0,y_0)-L|<\epsilon
```

"

# ╔═╡ 59f981d0-ade6-11ec-01e1-cb295e707273
begin
	philip = load("limit2d.jpg")
	imresize(philip, (450, 500));
end

# ╔═╡ 25be652e-adea-11ec-1667-d9bcf73815e9
md"
Let's start from the limit epsilon delta definition 

```math
|f(x_0,y_0)-L|<\epsilon \qquad \Longleftrightarrow \qquad -\epsilon <  f(x_0,y_0)-L<\epsilon -\epsilon +L<  f(x,y)<\epsilon+L
```
and the following hypothesis
```math
\iint_{\Omega} f(x,y) dx fy = 0
```
that it is valid for any $\epsilon >0$, $\delta >0$, and for any $(x_0,y_0)\in \Omega$, In particular, taking $\epsilon=|L|/2$
```math
\frac{L}{2} <  f(x,y)< \frac{3L}{2}
```
Now, integrating the last inequality over $\Omega$

```math
\iint_{\Omega} \frac{L}{2} \;dxdy< \iint_{\Omega}  f(x,y) \;dxdy< \iint_{\Omega} \frac{3L}{2} \; dxdy
```
then
```math
 \iint_{\Omega}  f(x,y) \;dxdy > \iint_{\Omega} \frac{L}{2} \;dxdy \geq \int_{y_0-\delta}^{y_0+\delta} \int_{x_0-\delta}^{x_0+\delta} \frac{L}{2} \; dx dy = \frac{L}{2} (2\delta) (2\delta) = 2 L \; \cdot \; \delta^2 
```
Thus 
```math
 \iint_{\Omega}  f(x_0,y_0) \;dxdy > 2 L \; \cdot \; \delta^2 
```
Where this integral is zero if and only if $L=0$ but it is contradictory since $\epsilon > 0$, then $L=f(x_0,y_0)$ must be $0$ for any $(x_0,y_0)\in \Omega$. 

Thereby, since $f(x,y)$ is continuos in $\Omega$ it satisfies that
```math
 \boxed{\iint_{\Omega}  f(x,y) \;dxdy = 0}
```
"

# ╔═╡ 91fb9ae0-0f8c-4d2e-9220-70701da28913


# ╔═╡ cdf1c2e0-a9b9-11ec-059a-8b69956003ba
begin
	plot_umv2  = scatter(;x=x, y=broadcast(abs, u.(x) - v.(x)).^2, mode="lines+markers", fill="tozeroy", name="|u(x)-v(x)|^2")
	plot_umv  = scatter(;x=x, y=u.(x) - v.(x), mode="lines+markers", name="u(x)-v(x)")

	plot_upmvp2 = scatter(;x=x, y=broadcast(abs, up.(x) - vp.(x)).^2, mode="lines+markers", fill="tozeroy", name="|u'(x) -v'(x)|^2")
	plot_upmvp  = scatter(;x=x, y= up.(x) - vp.(x), mode="lines+markers", name="u'(x)-v'(x)")

	p3 = plot([plot_umv2, plot_umv], Layout(title="Distance Between Functions",  xaxis_title="x"))
	p4 = plot([plot_upmvp2, plot_upmvp], Layout(title="Distance Between Derivative Functions",  xaxis_title="x"))

	p = [p3 p4]
	#relayout!(p, height=300, width=700, title_text="Functions", legend_title_text="Legend")
	#p
end

# ╔═╡ 1b947c8e-a9ba-11ec-3426-81427ab164ec
integral2, err = quadgk(x -> broadcast(abs, up(x)-vp(x))^2, 0, 2, rtol=1e-8)

# ╔═╡ 78ab34a0-a9ba-11ec-0276-1f5df54ae99e
integraluumv, err = quadgk(x -> broadcast(abs, u.(x).*(u.(x) - v.(x))).^2, 0, 2, rtol=1e-8)

# ╔═╡ 08041820-a9ba-11ec-0caf-0dac024f7add
integral1, err = quadgk(x -> broadcast(abs, u(x)-v(x))^2, 0, 2, rtol=1e-8)

# ╔═╡ 7aaf1af0-a9ba-11ec-2a0a-f3cfb773fe66
integralu2, err = quadgk(x -> broadcast(abs, u.(x)).^2, 0, 2, rtol=1e-8)

# ╔═╡ c0cca110-a9ba-11ec-3eac-07ae73f67832
integralupupmp, err = quadgk(x -> broadcast(abs, up.(x).*(up.(x) - vp.(x))).^2, 0, 2, rtol=1e-8)

# ╔═╡ 4caa1970-a9ba-11ec-2a41-57bbb9c0cd9e
begin
	plot_umumv2    = scatter(;x=x, y=broadcast(abs, u.(x).*(u.(x) - v.(x))).^2, mode="lines+markers", fill="tozeroy", name="|u(u-v)|^2")
	plot_upmupmvp2 = scatter(;x=x, y=broadcast(abs, up.(x).*(up.(x) - vp.(x))).^2, mode="lines+markers", fill="tozeroy", name="|u'(u'-v')|^2")

	p5 = plot([plot_u, plot_umv, plot_umumv2], Layout(title="Functions",  xaxis_title="x"))
	p6 = plot([plot_up, plot_upmvp, plot_upmupmvp2], Layout(title="Derivative Functions",  xaxis_title="x"))

	p = [p5 p6]
	#relayout!(p, height=300, width=700, title_text="Functions", legend_title_text="Legend")
	#p
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Calculus = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"

[compat]
Calculus = "~0.5.1"
ImageIO = "~0.6.1"
ImageShow = "~0.3.3"
Images = "~0.25.1"
LaTeXStrings = "~1.3.0"
PlotlyJS = "~0.18.8"
PlutoUI = "~0.7.37"
QuadGK = "~2.4.2"
Roots = "~1.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "d49f55ff9c7ee06930b0f65b1df2bfa811418475"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.4"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BinDeps]]
deps = ["Libdl", "Pkg", "SHA", "URIParser", "Unicode"]
git-tree-sha1 = "1289b57e8cf019aede076edab0587eb9644175bd"
uuid = "9e28174c-4ba2-5203-b857-d8d62c4213ee"
version = "1.0.2"

[[deps.Blink]]
deps = ["Base64", "BinDeps", "Distributed", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Reexport", "Sockets", "WebIO", "WebSockets"]
git-tree-sha1 = "08d0b679fd7caa49e2bca9214b131289e19808c0"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.5"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "7a20463713d239a19cbad3f6991e404aca876bda"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.15"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "7668b123ecfd39a6ae3fc31c532b588999bdc166"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.1"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1d2d73b14198d10f7f12bf7f8481fd4b3ff5cd61"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.0"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "d0ac64c9bee0aed6fdbb2bc0e5dfa9a3a78e3acc"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.3"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "42fe8de1fe1f80dab37a39d391b6301f7aeaa7b8"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.4"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "11d268adba1869067620659e7cdf07f5e54b6c76"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "cf737764159c66b95cdbf5c10484929b247fecfe"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "bd6c034156b1e7295450a219c4340e32e50b08b1"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.3"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2ef87eeaa28713cb010f9fb0be288b6c1a4ecd53"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.1.0+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "2af69ff3c024d13bde52b34a2a7d6887d4e7b438"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "bfbd6fb946d967794498790aa7a0e6cdf1120f41"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.13"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "Pkg", "Sockets", "WebSockets"]
git-tree-sha1 = "82dfb2cead9895e10ee1b0ca37a01088456c4364"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "0.7.6"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "180d744848ba316a3d0fdf4dbd34b77c7242963a"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.18"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "53d6325e14d3bdb85fd387a085075f36082f35a3"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.8"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "522770af103809e8346aefa4b25c31fbec377ccf"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "6085b8ac184add45b586ed8d74468310948dcfe8"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.4.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "a167638e2cbd8ac41f9cd57282cab9b042fa26e6"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "65068e4b4d10f3c31aaae2e6cb92b6c6cedca610"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "aaa19086bc282630d82f818456bc40b4d314307d"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.4"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "c9529be473e97fa0b3b2642cdafcd0896b4c9494"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.17"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "f91a602e25fe6b89afc93cf02a4ae18ee9384ce3"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.5.9"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "505c31f585405fc375d99d02588f6ceaba791241"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═2b55a600-a9b9-11ec-3eb1-1f25b963a10f
# ╠═326d71c0-a9b9-11ec-3986-c9776d176014
# ╟─4e476b30-a9b9-11ec-2d51-2fa7754c2747
# ╟─f7dc9db0-a9b8-11ec-15a1-ab2e399d50f1
# ╠═06249f30-a9b9-11ec-1093-a59349277d1a
# ╠═74ae8100-a9b9-11ec-1701-1fc4d58c62ad
# ╟─7c7f3c80-a9b9-11ec-1732-4d03f3b17dd6
# ╠═ad4e5a80-a9b9-11ec-2b9f-8b076b06cb32
# ╟─b0f1c000-a9b9-11ec-32a0-81846014b48d
# ╠═c50f4940-a9b9-11ec-07df-dd60bccf4f04
# ╟─c5f4ee50-a9b9-11ec-07ad-8ff5b133082b
# ╠═cdf1c2e0-a9b9-11ec-059a-8b69956003ba
# ╟─d5f6adc0-a9b9-11ec-3e41-e7cedf88a606
# ╠═08041820-a9ba-11ec-0caf-0dac024f7add
# ╟─0b8b9130-a9ba-11ec-2207-774f008ffed7
# ╠═1b947c8e-a9ba-11ec-3426-81427ab164ec
# ╟─200e1dd0-a9ba-11ec-13b5-9f9c9bad7ea8
# ╠═4210ef1e-a9ba-11ec-3e6b-9314377e22aa
# ╟─45897410-a9ba-11ec-0ebb-e5c07ed91ae8
# ╠═4caa1970-a9ba-11ec-2a41-57bbb9c0cd9e
# ╟─5795f7a0-a9ba-11ec-015f-834a762cdcdf
# ╠═78ab34a0-a9ba-11ec-0276-1f5df54ae99e
# ╠═7aaf1af0-a9ba-11ec-2a0a-f3cfb773fe66
# ╠═7ca852e0-a9ba-11ec-0eef-c385a17b5f74
# ╠═823f23f0-a9ba-11ec-3747-0b794a9a3151
# ╠═c0cca110-a9ba-11ec-3eac-07ae73f67832
# ╠═0e849070-a9bb-11ec-35c7-99a985571bbb
# ╠═11431430-a9bb-11ec-1174-6df0edc48380
# ╠═2dedf230-a9bb-11ec-1474-ef4bd3035175
# ╟─59f981d0-ade6-11ec-01e1-cb295e707273
# ╠═25be652e-adea-11ec-1667-d9bcf73815e9
# ╠═91fb9ae0-0f8c-4d2e-9220-70701da28913
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
