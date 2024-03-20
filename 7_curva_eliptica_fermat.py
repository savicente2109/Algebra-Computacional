# %% [markdown]
# ENTREGA 7

# %%
# Función que toma un primo p=/=3 y devuelve el conjunto de puntos distintos de CF(Fp).
# En el caso peor, el coste en tiempo es O(p + p^2) = O(p^2).
# Sí se nos ocurrió una forma de hacerlo en O(p), pero a cambio pasaría a ser O(p)
# en espacio adicional también.
def points(p):

    # Puntos con z=0 (puntos de la curva que están en la recta del infinito):

    res = {(1, p-1, 0)} # [1:-1:0] está siempre en la curva.

    if p % 3 == 1:
        # Recorrido lineal en p que busca un elemento de orden 3 de Fp*.
        for i in range(2, p):
            if pow(i, 3, p) == 1:
                res.add((1, p - i, 0)) # [1:-xi_3:0]
                res.add((1, p - ((i*i) % p), 0)) # [1:-xi_3^2:0]
                break
    
    # Puntos con z=1 (puntos de la curva que están en el plano afín):

    for x in range(p):
        for y in range(p):
            if (x**3 + y**3) % p == 1:
                res.add((x, y, 1))
    
    return res

# %%
# Elemento neutro de la ley de adición definida en la curva:
O = (1, 0, 1)

# %%
# Función que calcula el inverso multiplicativo de un número a en Fp*.
# Como a es coprimo con p (es estrictamente menor), por el teorema de Euler tenemos
# a^(-1) mód p = a ^ (phi(p) - 1) mód p, donde phi es la función de Euler.
# Como p es primo, phi(p) = p - 1, phi(p) - 1 = p - 2.
# El cálculo tiene complejidad O(log(p)) gracias a la exponenciación binaria modular de la función pow() de Python.
def inverso_mod_p(a, p):
    return pow(a, p - 2, p)

# %%
# Función que normaliza un punto P = [x:y:z] de la curva por si su representación no es exactamente la misma
# que la que precisan las demás funciones.
def normalizar(P, p):
    if P[2] == 0:
        if P[0] == 1:
            return P
        # P[0] no puede ser 0 porque no hay ningún punto de la curva con x, z ambas igual a 0.
        return (1, (P[1] * inverso_mod_p(P[0], p)) % p, 0)
    
    if P[2] == 1:
        return P
    
    inv = inverso_mod_p(P[2], p)
    return ((P[0] * inv) % p, (P[1] * inv) % p, 1)

# %%
# Función que avergiua si un punto P de la curva es uno de sus puntos de inflexión.
def es_pto_inflexion(P, p):
    # La distinción de casos no es necesaria, pero evita calcular potencias innecesarias cuando p = 2 mód 3.
    # Si p = 1 mód 3, los puntos de inflexión son los tres que hay con z = 0, O, [1:0:xi_3] = [xi_3^2:0:xi_3^3] = [xi_3^2:0:1],
    # [1:0:xi_3^2] = [xi_3:0:xi_3^3] = [xi_3:0:1], [0:1:1], [0:xi_3:1] y [0:xi_3^2:1].
    if p % 3 == 1:
        return P == O or P[2] == 0 or P[1] == 0 and (P[0]**3) % p == 1 or P[0] == 0 and (P[1]**3) % p == 1
    
    # Si p = 2 mód 3, los puntos de inflexión son O, [1:-1:0] (que es el único con z = 0) y [0:1:1]
    return P == O or P[2] == 0 or P[0] == 0 and P[1] == 1

# %%
# Función que calcula P * Q distinguiendo diversos casos.
def ast(P, Q, p):

    P = normalizar(P, p)
    Q = normalizar(Q, p)

    if P == Q:
        # Si P es punto de inflexión, P * P = P.
        if es_pto_inflexion(P, p):
            return P
        
        # Si P no es punto de inflexión, hay que buscar el punto de corte de la curva con la recta tangente en P.
        # Además, en este caso tenemos asegurado que z = 1, ya que todos los puntos con z = 0 son de inflexión.
        # P = [a:b:1].
        # gradiente(F) = (3x^2, 3y^2, -3z^2), gradiente(F)(P) = (3a^2, 3b^2, -3).
        # La recta tangente en P es 3a^2x + 3b^2y - 3z = 0, es decir, a^2x + b^2y = z.
        # El vector director de esa recta (considerando solo los puntos afines), es el (-b^2, a^2) y, por tanto,
        # la parametrización de la recta es
        # x = a - b^2t
        # y = b + a^2t
        # z = 1 (solo podemos obtener los puntos afines)

        # Sustituyendo en F(x,y,z) = 0, tenemos (a-b^2t)^3 + (b+a^2t)^3 - 1^3 = 0.
        # Desarrollando, obtenemos una ecuación cúbica At^3 + Bt^2 + Ct + D = 0 con
        # A = a^6 - b^6, B = 3(ab^4 + ba^4), C = 0, D = a^3 + b^3 - 1.

        # El valor para t = 0 corresponde al propio punto P, por lo que t = 0 siempre es raíz doble de la ecuación.
        # Por ello, C = D = 0 y At^3 + Bt^2 = 0. Sacando factor común, At + B = 0.

        A = (pow(P[0], 6, p) - pow(P[1], 6, p)) % p

        # Si A =/= 0, t = -B / A define el tercer punto de corte de la recta.
        if A != 0:
            B = (3 * (P[0] * pow(P[1], 4, p) + P[1] * pow(P[0], 4, p))) % p
            t = p - ((B * inverso_mod_p(A, p)) % p)
            return ((P[0] - P[1]**2 * t) % p, (P[1] + P[0]**2 * t) % p, 1)
        
        # Si A = 0, la ecuación no es cúbica y no tiene una tercera raíz. Esto se debe a que el punto de corte está
        # en la recta del infinito.
        # Si p = 2 mód 3, solo hay un punto de la curva con z = 0, el [1:-1:0], por lo que necesariamente debe ser ese.
        if p % 3 == 2:
            return (1, p - 1, 0)
        
        # En el caso p = 1 mód 3, calculamos el punto del infinito de todas las rectas con pendiente la de la recta tangente que pasa por P.
        # Recordamos que su vector director es (-b^2, a^2), así que la pendiente es -a^2/b^2.
        m = p - (((P[0] * P[0]) % p) * inverso_mod_p(P[1] * P[1], p) % p)
        return (1, m, 0)
    
    else: # Caso P =/= Q:
        # Si tanto P como Q están en la recta del infinito, hay que devolver el tercer punto de la recta del infinito,
        # ya que están alineados y son los únicos puntos con z = 0 de la curva (este caso solo se puede dar con p = 1 mód 3).
        if P[2] == 0 and Q[2] == 0:
            if P[1] == p - 1:
                return (1, p - ((Q[1] * Q[1]) % p), 0)
            if Q[1] == p - 1:
                return (1, p - ((P[1] * P[1]) % p), 0)
            return (1, p - 1, 0)
        
        # Si P está en la recta del infinito y Q está en el plano afín, R = P * Q estará también en el plano afín (z = 1),
        # ya que en caso contrario R y P estarían alineados en la recta del infinito y por tanto R * P también estaría en ella,
        # pero R * P = Q (!).
        # Buscamos por tanto el punto de corte (en el plano afín) de la curva con la recta que pasa por Q = [x_1:y_1:1]
        # y tal que su punto del infinito es P = [1:y_2:0]
        # Como su punto del infinito es P, la pendiente debe ser (y_2)/1 = y_2.
        # Por tanto, la recta en paramétricas es
        # x = x_1 + t
        # y = y_1 + y_2t
        # z = 1

        # Sustituyendo en F(x,y,z) = 0, tenemos (x_1 + t)^3 + (y_1 + y_2t)^3 - 1^3 = 0.
        # Desarrollando, obtenemos una ecuación cúbica At^3 + Bt^2 + Ct + D = 0 con
        # A = 1 + y_2^3, B = 3(x_1 + y_1y_2^2), C = 3(x_1^2 + y_1^2y_2), D = x_1^3 + y_1^3 - 1.

        # Como t = 0 corresponde al punto Q, siempre tenemos que t = 0 es raíz. Por tanto, D = 0.
        # Además, como la ecuación solo va a tener dos raíces (pues el punto P está en la recta del infinito),
        # también tendremos A = 0. Es decir, Bt^2 + C^t = 0. Sacando factor común, Bt + C = 0, y t = -C/B.
        if P[2] == 0:
            # Llamamos C a C/3 y D a D/3
            y1y2 = Q[1] * P[1] % p
            C = (Q[0] * Q[0] + Q[1] * y1y2) % p
            B = (Q[0] + y1y2 * P[1]) % p
            t = p - ((C * inverso_mod_p(B, p)) % p)
            return ((Q[0] + t) % p, (Q[1] + P[1] * t) % p, 1)
        
        # Caso análogo al anterior pero con Q en la recta del infinito y P en el plano afín.
        if Q[2] == 0:
            y1y2 = P[1] * Q[1] % p
            C = (P[0] * P[0] + P[1] * y1y2) % p
            B = (P[0] + y1y2 * Q[1]) % p
            t = p - ((C * inverso_mod_p(B, p)) % p)
            return ((P[0] + t) % p, (P[1] + Q[1] * t) % p, 1)

        # Si tanto P = [x_1:y_1:1] como Q = [x_2:y_2:1] están en el plano afín, la recta parametrizada que los une es tP + (1-t)Q = 0.
        # Es decir,
        # x = x_2 + (x_1 - x_2)t
        # y = y_2 + (y_1 - y_2)t
        # z = 1
        # Sustituyendo en F(x,y,z) = 0, tenemos (x_2 + (x_1 - x_2)t)^3 + (y_2 + (y_1 - y_2)t)^3 - 1^3 = 0.
        # Desarrollando, obtenemos una ecuación cúbica At^3 + Bt^2 + Ct + D = 0 con
        # A = (x_1 - x_2)^3 + (y_1 - y_2)^3, B = 3(x_2(x_1 - x_2)^2 + y_2(y_1 - y_2)^2), C = 3(x_2^2(x_1 - x_2) + y_2^2(y_1 - y_2)), D = x_2^3 + y_2^3 - 1.

        # Como t = 0 corresponde al punto P, tenemos que t = 0 siempre es raíz de la ecuación y por tanto D = 0.
        # Tenemos entonces que At^3 + Bt^2 + Ct = 0. Sacando factor común, At^2 + Bt + C = 0.

        # Como t = 1 corresponde al punto Q, tenemos que t = 1 siempre es raíz de la ecuación y por tanto A + B + C = 0.
        # Dividimos entre el factor (t + 1) y obtenemos At + (A + B) = 0.
            
        restax = (P[0] - Q[0]) % p
        restay = (P[1] - Q[1]) % p
        restax_3 = (restax**3) % p
        restay_3 = (restay**3) % p

        A = (restax_3 + restay_3) % p

        # Si A =/= 0, tenemos que t = - (A + B) / A = C / A.
        if A != 0:
            C = (3 * (Q[0] * Q[0] * restax + Q[1] * Q[1] * restay)) % p
            t = (C * inverso_mod_p(A, p)) % p
            return ((t * restax + Q[0]) % p, (t * restay + Q[1]) % p, 1)
        
        # Si A = 0, la ecuación solo tiene dos raíces (t=0, t=1). Esto se debe a que el tercer punto de corte se encuentra
        # en la recta del infinito.
        # En el caso p = 2 mód 3, solo hay un punto en dicha recta, por lo que forzosamente debe ser ese.
        if p % 3 == 2:
            return (1, p - 1, 0)

        # En el caso p = 1 mód 3, calculamos el punto del infinito de todas las rectas con pendiente la de la recta afín que une P y Q.
        # La pendiente es (y_2 - y_1)/(x_2 - x_1).
        m = (restay * inverso_mod_p(restax, p)) % p
        return (1, m, 0)

# %%
# Función que calcula P + Q según la ley de adición definida por O.
def add(P, Q, p):
    # P + Q = O * (P * Q).
    return ast(O, ast(P,Q,p), p)

# Función que calcula el inverso de P, esto es, el punto Q tal que P + Q = 0.
def inv(P, p):
    # El inverso de P es (O * O) * P y, como O = [1:0:1] es punto de inflexión, O * O = O.
    return ast(O, P, p)

# Función que calcula P + P + ... + P (k veces) tomada de "eliptica.py" (vista en clase).
def mul(k, P, p):
    if k < 0:
        return mul(-k, inv(P, p), p)
    if k == 0:
        return O
    if k % 2 == 0:
        Q = mul(k//2, P, p)
        return add(Q, Q, p)
    Q = mul(k-1, P, p)
    return add(P, Q, p)