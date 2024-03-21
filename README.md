# Álgebra Computacional
 Prácticas de la asignatura Álgebra Computacional.

# Práctica 5
Implementación de la función `sqrt_mod(a, p, n)`, que calcula una solución de $x^2 \equiv a \mod p^n$, donde $p$ es un primo impar, $n \geq 1$ y $0 \leq a < p^n$.

Fichero `5_raiz_cuadrada_modular.py`.

# Práctica 6
Implementación de las funciones `mult_pol_mod(f, g, p)` y `mult_ss_mod(f, g, k, p)`.
* `mult_pol_mod(f, g, p)` calcula el producto de los polinomios $f$ y $g$ en el anillo $(\mathbb{Z}/p\mathbb{Z})[x]$.
* `mult_ss_mod(f, g, k, p)` calcula el producto de los polinomios $f$ y $g$ en el anillo $(\mathbb{Z}/p\mathbb{Z})[x]/⟨x^{2^k}+1⟩$ mediante el algoritmo de Schönhage–Strassen.

Fichero `6_schonhage_strassen.py`.

# Práctica 7
Implementación de las operaciones de la ley de grupo definida en $\mathcal{C}_F$ como $P+Q=\mathcal{O}∗(P∗Q)$, al fijar el punto $\mathcal{O} = [1:0:1]$ como neutro, con $F(x, y, z) = x^3 + y^3 - z^3$.

Fichero `7_curva_eliptica_fermat.py`.
