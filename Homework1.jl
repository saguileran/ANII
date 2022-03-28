### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ 326d71c0-a9b9-11ec-3986-c9776d176014
begin
	using PlutoUI, QuadGK, Calculus, Roots, LaTeXStrings, PlotlyJS, ImageShow, ImageIO, Images
	TableOfContents(title="📚 Table of Contents Homework 1", indent=true, depth=4, aside=true)
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
	url = "https://files.mtstatic.com/site_4425/3594/0?Expires=1648398822&Signature=XiIbvMeRCn6a6DmhKE6ympokXM1slrsUKh2qWfEeYqoq4lqUjE27RXL~HfOiOoDQUszhX6sZB1T1-AXv2O9UNxDZSWZa7NXxDbyC3xQ1BwpZLriKln82AD8EFP1saE2PLD~9D2bIAGOQYJggPt37KYYVfE6VA6nr0LsTI4q03uI_&Key-Pair-Id=APKAJ5Y6AV4GI7A555NA"
	philip_filename = download(url) # download to a local file. The filename is returned
	philip = load(philip_filename)
	imresize(philip, (450, 500));
end

# ╔═╡ 25be652e-adea-11ec-1667-d9bcf73815e9


# ╔═╡ 08041820-a9ba-11ec-0caf-0dac024f7add
integral1, err = quadgk(x -> broadcast(abs, u(x)-v(x))^2, 0, 2, rtol=1e-8)

# ╔═╡ c0cca110-a9ba-11ec-3eac-07ae73f67832
integralupupmp, err = quadgk(x -> broadcast(abs, up.(x).*(up.(x) - vp.(x))).^2, 0, 2, rtol=1e-8)

# ╔═╡ 1b947c8e-a9ba-11ec-3426-81427ab164ec
integral2, err = quadgk(x -> broadcast(abs, up(x)-vp(x))^2, 0, 2, rtol=1e-8)

# ╔═╡ 7aaf1af0-a9ba-11ec-2a0a-f3cfb773fe66
integralu2, err = quadgk(x -> broadcast(abs, u.(x)).^2, 0, 2, rtol=1e-8)

# ╔═╡ 78ab34a0-a9ba-11ec-0276-1f5df54ae99e
integraluumv, err = quadgk(x -> broadcast(abs, u.(x).*(u.(x) - v.(x))).^2, 0, 2, rtol=1e-8)

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

# ╔═╡ Cell order:
# ╠═2b55a600-a9b9-11ec-3eb1-1f25b963a10f
# ╠═326d71c0-a9b9-11ec-3986-c9776d176014
# ╟─4e476b30-a9b9-11ec-2d51-2fa7754c2747
# ╟─f7dc9db0-a9b8-11ec-15a1-ab2e399d50f1
# ╟─06249f30-a9b9-11ec-1093-a59349277d1a
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
# ╟─823f23f0-a9ba-11ec-3747-0b794a9a3151
# ╠═c0cca110-a9ba-11ec-3eac-07ae73f67832
# ╠═0e849070-a9bb-11ec-35c7-99a985571bbb
# ╠═11431430-a9bb-11ec-1174-6df0edc48380
# ╟─2dedf230-a9bb-11ec-1474-ef4bd3035175
# ╟─59f981d0-ade6-11ec-01e1-cb295e707273
# ╠═25be652e-adea-11ec-1667-d9bcf73815e9
