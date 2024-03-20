# %% [markdown]
# ### ENTREGA 6

# %% [markdown]
# Funciones auxiliares:

# %% [markdown]
# Operaciones vectoriales sobre polinomios como listas de coeficientes:

# %%
# Función que suma dos polinomios como listas de enteros módulo p componente a componente.
def sumar_polinomios_mod_p(f, g, p):
    n = len(f)
    return [(f[i] + g[i]) % p for i in range(n)]

# %%
# Función que resta dos polinomios como listas de enteros módulo p componente a componente.
def restar_polinomios_mod_p(f, g, p):
    n = len(f)
    return [(f[i] - g[i]) % p for i in range(n)]

# %%
# Función que cambia de signo un polinomio como lista de enteros módulo p componente a componente.
def cambiar_signo_pol_mod_p(f, p):
    return [(-x) % p for x in f]

# %%
# Función que multiplica el polinomio f como lista de coeficientes enteros por el escalar n, módulo p.
def multiplicar_escalar_mod_p(f, n, p):
    return [(a * n) % p for a in f]

# %% [markdown]
# Operaciones necesarias para llevar a cabo la negaconvolución:

# %%
# Función que "rota" un polinomio como lista de enteros un número de posiciones a la derecha.
# Los elementos que "salen de la lista" y "vuelven a entrar" por la izquierda cambian de signo módulo p.
# Ejemplo: rotar_pol([2, 3, 1, 4], 1, 7) = [-4, 2, 3, 1] mód 7 = [3, 2, 3, 1] mód 7.
def rotar_pol(f, posiciones, p):
    return cambiar_signo_pol_mod_p(f[-posiciones:], p) + f[:-posiciones]

# %%
# Función que "rota negativamente" un polinomio como lista de enteros un número de posiciones a la derecha.
# Los elementos que "no llegan a salir de la lista" cambian de signo módulo p.
# Ejemplo: rotar_pol([2, 3, 1, 4], 1, 7) = [4, -2, -3, -1] mód 7 = [4, 5, 4, 6] mód 7.
def rotar_pol_neg(f, posiciones, p):
    return f[-posiciones:] + cambiar_signo_pol_mod_p(f[:-posiciones], p)

# %%
# Función que multiplica los vectores a = [1_A, xi_2n2, (xi_2n2)^2, ..., (xi_2n2)^(n2-1)] y f = [[f_0], [f_1], ..., [f_(n2-1)]] componente a componente en A[y] / <y^(n2) + 1>,
# donde xi_2n2 es una raíz 2(n2)-ésima de la unidad en A = (Z/pZ)[u] / <u^(2n1) + 1>.
def multiplicar_por_a(f, exp, p):
    # xi_2n2 será [u] o [u^2] (exp = 1 o exp = 2).

    n2 = len(f)

    # En el caso caso exp = 1, a = ([1], [u], [u^2], [u^3], ..., [u^(n2-1)]). Los exponentes avanzan de 1 en 1.
    # En este caso exp = 2, a = ([1], [u^2], [u^4], [u^6], ..., [u^(2n2-2)]). Los exponentes avanzan de 2 en 2.
    # Multipicar un polinomio por [u^i] es aplicarle una rotación positiva de i lugares a la derecha.
    # Los coeficientes cambian de signo al "volver a entrar" por la izquierda.

    for i in range(1, n2): # Obviamos el caso i = 0 porque multiplicar por [1] no provoca cambio.
        f[i] = rotar_pol(f[i], i * exp, p)

    return f

# %%
# Función que multiplica los vectores a^(-1) = ([1_A, xi_2n2, (xi_2n2)^2, ..., (xi_2n2)^(n2-1)])^(-1)
# y f = [[f_0], [f_1], ..., [f_(n2-1)]] componente a componente en A[y] / <y^(n2) + 1>,
# donde xi_2n2 es una raíz 2(n2)-ésima de la unidad en A = (Z/pZ)[u] / <u^(2n1) + 1>.
def multiplicar_por_a_inv(f, exp, p):
    # xi_2n2 será [u] o [u^2] (exp = 1 o exp = 2).

    n2 = len(f)
    dosn1 = len(f[0])

    # En el caso caso exp = 1, a = ([1], [u], [u^2], [u^3], ..., [u^(n2-1)]) y 2n1 = n1.
    # En este caso exp = 2, a = ([1], [u^2], [u^4], [u^6], ..., [u^(2n2-2)]) y n1 = n2.
    
    # El elemento neutro de la operación "multiplicar componente a componente" es ([1], [1], ..., [1]),
    # por lo que el inverso de a es el vector de los inversos de los elementos de a.
    # El inverso de [1] es [1], y para el resto se cambia el exponente por 2n1 - exponente y se cambia el signo de la raíz
    # (en este caso siempre a negativo por ser todas positivas). Es decir,

    # En el caso caso exp = 1, a^(-1) = ([1], [-u^(2n1-1)], [-u^(2n1-2)], [-u^(2n1-3)], ..., [u^(2n1-n2+1)] = [u]). Los exponentes disminuyen de 1 en 1.
    # En este caso exp = 2, a = ([1], [-u^(2n1-2)], [-u^(2n1-4)], [-u^(2n1-6)], ..., [-u^(2n1-2n2+2)] = [u^2]). Los exponentes disminuyen de 2 en 2.

    # Multipicar un polinomio por [-u^i] es aplicarle una rotación negativa de i lugares a la derecha.
    # Los coeficientes que cambian de signo son los que "no salen del polinomio".

    for i in range(1, n2): # Obviamos el caso i = 0 porque multiplicar por [1] no provoca cambio.
        f[i] = rotar_pol_neg(f[i], dosn1 - i * exp, p)

    return f

# %%
# Función que calcula xi^exp = ([+/- u^i])^exp módulo <u^n + 1>.
# Por ejemplo:
# (u^2)^3 = u^6 = u^4 * u^2 = -u^2 mód <u^4 + 1>,
# (-u^2)^3 = -u^6 = -u^4 * u^2 = u^2 mód <u^4 + 1>,
# (-u^2)^6 = u^12 = u^8 * u^4 = -u^4 mód <u^8 + 1>.
def potencia_u_n(xi, exp, n):
    exp_nuevo = exp * xi[1]
    cambios_signo = exp_nuevo // n
    exp_nuevo = exp_nuevo % n
    if exp % 2 == 0: # Si la potencia es par, se ignora el signo antiguo:
        cambios_signo = cambios_signo % 2
    else: # Si es impar, se toma el signo antiguo xi[0] y se le suman los cambios de signo:
        cambios_signo = (cambios_signo + xi[0]) % 2
    return (cambios_signo, exp_nuevo)

# %%
# Función que calcula el inverso de xi mód <u^n + 1> 
def inverso_xi(xi, n):
    if xi[1] == 0: # El inverso de 1 es 1, y el inverso de -1 es -1.
        return xi
    else:
        return (1 - xi[0], n - xi[1])

# %%
# Función que calcula el inverso multiplicativo de un número n en Z/pZ.
# Como n es coprimo con p (n es una potencia de 2 y p es primo distinto de 2), por el teorema de Euler tenemos
# n^(-1) mód p = n ^ (phi(p) - 1) mód p, donde phi es la función de Euler.
# Como p es primo, phi(p) = p - 1, phi(p) - 1 = p - 2.
# El cálculo tiene complejidad O(log(p)) gracias a la exponenciación binaria modular de la función pow() de Python.
def inv_mod_p(n, p):
    return pow(n, p - 2, p)

# %% [markdown]
# Transforma discreta de Fourier:

# %%
# Función que calcula la transformada discreta (rápida) de Fourier de un polinomio f como lista de coeficientes en ((Z/pZ)[u] / <u^(2n1) + 1>) ^ (n2).
# La raíz n2-ésima de la unidad xi será una potencia de [u] y está representada internamente por una tupla (signo, exp)
# donde signo es un entero (0 = positivo, 1 = negativo) y exp es el exponente entre 0 y 2n1-1, (xi = [+/- u^exp]).
def fft(f, xi, p):

    n2 = len(f)
    if n2 == 1:
        return f
    
    mitad = n2 // 2
    dosn1 = len(f[0])

    f_even = [[] for _ in range(mitad)]
    f_odd = [[] for _ in range(mitad)]

    for i in range(mitad):
        f_even[i] = f[2*i]
        f_odd[i]  = f[2*i+1]
        
    xi_cuadrado = potencia_u_n(xi, 2, dosn1)

    res_even = fft(f_even, xi_cuadrado, p)
    res_odd = fft(f_odd, xi_cuadrado, p)
    res = [[] for _ in range(n2)]

    # Separamos el caso i = 0 para no perder tiempo en multiplicar por xi^0 = 1:
    res[0] = sumar_polinomios_mod_p(res_even[0], res_odd[0], p)
    res[mitad] = restar_polinomios_mod_p(res_even[0], res_odd[0], p)

    for i in range(1, mitad):
        potencia = potencia_u_n(xi, i, dosn1)
        if potencia[0] == 0:
            res[i] = sumar_polinomios_mod_p(res_even[i], rotar_pol(res_odd[i], potencia[1], p), p)
            res[i + mitad] = sumar_polinomios_mod_p(res_even[i], rotar_pol_neg(res_odd[i], potencia[1], p), p)
        else:
            res[i] = sumar_polinomios_mod_p(res_even[i], rotar_pol_neg(res_odd[i], potencia[1], p), p)
            res[i + mitad] = sumar_polinomios_mod_p(res_even[i], rotar_pol(res_odd[i], potencia[1], p), p)
    return res

# %% [markdown]
# Transformada discreta inversa de Fourier:

# %%
# Función que calcula la transformada discreta (rápida) inversa de Fourier de un polinomio f como lista de coeficientes en ((Z/pZ)[u] / <u^(2n1) + 1>) ^ (n2).
# La raíz n2-ésima de la unidad xi será una potencia de [u] y está representada internamente por una tupla (signo, exp)
# donde signo es un entero (0 = positivo, 1 = negativo) y exp es el exponente entre 0 y 2n1-1, (xi = [+/- u^exp]).
def ifft(f, xi, p):
    n2 = len(f)
    dosn1 = len(f[0])

    # Calculamos la transformada discreta de Fourier de f con el inverso de la raíz n2-ésima de la unidad xi.
    res = fft(f, inverso_xi(xi, dosn1), p)

    # para "dividir" el resultado entre n, calculamos el inverso de n módulo p y multiplicamos por él cada una de sus componentes.
    n2_inv = inv_mod_p(n2, p)
    for i in range(n2):
        res[i] = multiplicar_escalar_mod_p(res[i], n2_inv, p)
    return res

# %% [markdown]
# Negaconvolución:

# %%
# Función que calcula la negaconvolución de f y g en ((Z/pZ)[u] / <u^(2n1) + 1>) ^ (n2) con xi_2n2.
# xi_2n2 es una raíz 2(n2)-ésima de la unidad, que será [u] (exp = 1) o [u^2] (exp = 2).
# El parámetro k es necesario para la recursión (negaconv(f, g, k, exp, p) se invoca desde mult_ss_mod(f, g, k, p)
# e invoca a su vez a mult_ss_mod(f, g, k, p)).
def negaconv(f, g, k, exp, p):

    n2 = len(f)

    # Se calcula la transformada discreta de Fourier de a * f y de a * g.
    # Como xi_2n2 es una raíz 2(n2)-ésima de la unidad, (xi_2n2)^2 es una raíz n2-ésima de la unidad.
    # Si xi_2n2 = [u], (xi_2n2)^2 = [u^2]. Si xi_2n2 = [u^2], (xi_2n2)^2 = [u^4].
    # No es necesario invocar a potencia_u_n(xi, n) porque el exponente siempre va a ser menor que 2n1.
    # Si exp = 1 (2n1 = n2), lo mínimo que puede valer n1 (sin ser un caso base) es 2^1 = 2 (k = 3, k1 = 1, k2 = 2), y 2n1 = 2*2 = 4 > 2 = exp * 2.
    # Si exp = 2 (n1 = n2), lo mínimo que puede valer n1 (sin ser un caso base) es 2^2 = 4 (k = 4, k1 = 2, k2 = 2), y 2n1 = 2*4 = 8 > 4 = exp * 2.
    f_conv = fft(multiplicar_por_a(f, exp, p), (0, exp << 1), p)
    g_conv = fft(multiplicar_por_a(g, exp, p), (0, exp << 1), p)

    # Se multiplican f_conv y g_conv componente a componente.
    # Cada una de sus n2 componentes es una clase de polinomio de (Z/pZ)[u] / <u^(2n1) + 1>, por lo que se multiplican mediante
    # n2 llamadas recursivas a mult_ss_mod(f, g, k, p) con k = k1 + 1
    # (el valor de k ya viene correcto desde la invocación a negaconv(f, g, k, exp, p)).
    rec = [[] for _ in range(n2)]
    for i in range(n2):
        rec[i] = mult_ss_mod(f_conv[i], g_conv[i], k, p)

    # Se calcula la transformada discreta de Fourier inversa de f_conv * g_conv.
    inv = ifft(rec, (0, exp << 1), p)

    # Se calcula y se devuelve a^(-1) * inv.
    return multiplicar_por_a_inv(inv, exp, p)

# %% [markdown]
# Últimas funciones auxiliares necesarias para implementar el algoritmo de Shönhage-Strassen:

# %%
# Función que calcula el polinomio f_til = [f_0] + [f_1]y + [f_2]y^2 + ... + [f_(n2-1)]y^(n2-1) de ((Z/pZ)[u] / <u^(2n1) + 1>)[y] a partir de f.
def f_a_f_til(f, n1, n2):
    # La función debe devolver [f[i * n1 : (i + 1) * n1] + ([0] * n1) for i in range(n2)].
    # Para evitar realizar tantas multiplicaciones por el valor n1, que es una potencia de 2 que puede ser muy grande,
    # llevamos un acumulador. Únicamente realizamos un recorrido lineal de la lista que representa al polinomio f.
    f_til = [[] for i in range(n2)]
    acum = 0
    for i in range(n2):
        aux = acum
        acum += n1
        f_til[i] = f[aux : acum] + ([0] * n1)
    return f_til

# %%
# Función que realiza la operación inversa a f_a_f_til(f, n1, n2).
def f_til_a_f(f_til, n1, n2, p):
    n = n1 * n2
    f = [0] * n
    # Rellenamos los últimos n1 elementos de f (n1 operaciones de resta en Z/pZ).
    f[:n1] = restar_polinomios_mod_p(f_til[0][:n1], f_til[-1][n1:], p)
    acum = n1
    # En cada iteración del bucle rellenamos n1 elementos de f (n1 operaciones de suma en Z/pZ) y hay n2-1 iteraciones.
    for i in range(1, n2):
        aux = acum
        acum += n1
        f[aux : acum] = sumar_polinomios_mod_p(f_til[i][:n1], f_til[i - 1][n1:], p)
    # En total se han realizado n1 + n1 * (n2 - 1) = n1 * n2 = n operaciones de suma o resta en Z/pZ.
    return f

# %% [markdown]
# Función mult_ss_pol_mod(f, g, k, p) que calcula el producto de las clases de polinomios [f] y [g] en (Z/pZ)[x] / <x^(2^k) + 1> mediante el método de Schönhage-Strassen:

# %%
def mult_ss_mod(f, g, k, p):

    # CASOS BASE:
    # k = 0 (f, g son polinomios constantes).

    if k == 0:
        return [(f[0] * g[0]) % p]

    # k = 1
    # (a + bx)*(c + dx) = ac + (bc + ad)x + bdx^2 = (ac - bd) + (bc + ad)x mód <x^2 + 1>.

    if k == 1:
        return [(((f[0] * g[0]) % p) - ((f[1] * g[1]) % p) % p), (((f[1] * g[0]) % p) + ((f[0] * g[1]) % p) % p)]

    # k = 2
    # (a + bx + cx^2 + dx^3)*(A + Bx + Cx^2 + Dx^3) = aA + (aB+bA)x^2 + (aC+bB+cA)x^2 + (aD+bC+cB+dA)x^3 + (bD+cC+dB)x^4 + (cD+dC)x^5 + dDx^6 =
    # = (aA - bD - cC - dB) + (aB + bA - cD - dC)x + (aC + bB + cA - dD)x^2 + (aD + bC + cB + dA)x^3

    if k == 2:
        r0 = ((f[0] * g[0]) % p - (f[1] * g[3]) % p - (f[2] * g[2]) % p - (f[3] * g[1]) % p) % p
        r1 = ((f[0] * g[1]) % p + (f[1] * g[0]) % p - (f[2] * g[3]) % p - (f[3] * g[2]) % p) % p
        r2 = ((f[0] * g[2]) % p + (f[1] * g[1]) % p + (f[2] * g[0]) % p - (f[3] * g[3]) % p) % p
        r3 = ((f[0] * g[3]) % p + (f[1] * g[2]) % p + (f[2] * g[1]) % p + (f[3] * g[0]) % p) % p
        return [r0, r1, r2, r3]

    # CASO RECURSIVO:

    # Elegimos k1, k2:
    k1 = k // 2
    k2 = k - k1

    n1 = 1 << k1
    n2 = 1 << k2

    # Calculamos f_til y g_til:

    f_til = f_a_f_til(f, n1, n2)
    g_til = f_a_f_til(g, n1, n2)

    # Realizamos una negaconvolución con xi_2n2, que será [u] (exp = 1) o [u^2] (exp = 2)
    # según si k es impar (2n1 = n2, 2n1/n2 = 1) o par (n1 = n2, 2n1/n2 = 2).
    # En ella se realizarán n2 llamadas recursivas a mult_ss_mod(f, g, k, p) con listas de tamaño
    # 2n1 en lugar de n. El parámetro k de la llamada debe valer k1 + 1, ya que 2^(k1 + 1) = (2^k1)*2 = 2*n1.

    exp = 2 if k % 2 == 0 else 1

    h_til = negaconv(f_til, g_til, k1 + 1, exp, p)

    # Calculamos h:

    h = f_til_a_f(h_til, n1, n2, p)

    return h

# %% [markdown]
# Función mult_pol_mod(f, g, p) que calcula el producto de los polinomios f y g en (Z/pZ)[x] como caso particular de mult_ss_pol_mod(f, g, k, p):

# %%
def mult_pol_mod(f, g, p):
    deg_f = len(f) - 1
    deg_g = len(g) - 1
    if deg_f < 0 or deg_g < 0: # Si alguno de ellos era nulo (len(f) = 0 o len(g) = 0), el producto es nulo.
        return []
    # Calculamos k, con n = 2^k la menor potencia de 2 tal que deg(f) + deg(g) < n:
    # Si consideramos la representación binaria de deg(f) + deg(g), k es exactamente su número de bits,
    # o lo que es lo mismo, la cantidad de desplazamientos lógicos a la derecha necesarios para convertirlo en 0.
    suma = deg_f + deg_g
    aux = suma
    k = 0
    while aux > 0: # Este bucle es lineal en el número de bits de deg(f) + deg(g), es decir, logarítmico en el valor de deg(f) + deg(g).
        aux = aux >> 1
        k += 1
    # Calculamos la cantidad de coeficientes 0 que hay que añadir a f y a g para que cada lista tenga exactamente n = 2^k elementos:
    n = 1 << k
    ceros_extra_f = n - deg_f - 1
    ceros_extra_g = n - deg_g - 1
    # Además, sabemos que el polinomio resultado tendrá exactamente grado deg(f) + deg(g) por ser p primo:
    # El coeficiente del término de grado deg(f) de f y el de grado deg(g) de g deben ser no nulos (o no tendrían dicho grado),
    # es decir, no múltiplos de p. Llamémoslos a y b respectivamente.
    # Si el coeficiente de grado deg(f) + deg(g) de f*g, que es ab, fuera nulo, tendríamos que ab es múltiplo de p sin serlo ni a ni b,
    # y esto es imposible porque p es primo. Por tanto, no necesitamos un bucle ni comprobaciones para eliminar los ceros sobrantes.
    return mult_ss_mod(f + [0] * ceros_extra_f, g + [0] * ceros_extra_g, k, p) [:suma + 1] # El formato pedido, sin los ceros del final.