# %% [markdown]
# #ENTREGA 5

# %% [markdown]
# ## Funciones auxiliares
# 
# Para resolver el problema de encontrar una solución de la ecuación x*i = y mod p, en lugar de buscar la solución recorriendo todo el rango de posibles soluciónes, calculamos el inverso de x y obtenemos i como i = x^-1 * y mod p. El inverso de x lo hayamos con gcd_extendido_binario(x,y), ya que el primer valor que devuelve (a) es el inverso de x.  (1 = x*a + p*b = x*a mod p).

# %%
def gcd_extendido_binario(x,y):
  if x < 0:
    a, b, gcd = gcd_extendido_binario(-x,y)
    return -a, b, gcd
  if y < 0:
    a, b, gcd = gcd_extendido_binario(x,-y)
    return a, -b, gcd
  xespar = x%2 == 0
  yespar = y%2 == 0
  if x == 0: # caso base : gcd (0,y)=y
      return 0, 1, y
  elif y == 0: # caso base : gcd(x ,0)= x
      return 1, 0, x
  elif xespar and yespar :
      a, b , gcd = gcd_extendido_binario(x//2, y//2)
      return a, b, 2*gcd
  elif xespar :
      a, b , gcd = gcd_extendido_binario(x//2, y)
      if a % 2 == 0:
          return a//2, b, gcd
      else:
          return (a+y)//2, (2*b-x)//2, gcd
  elif yespar :
      a, b , gcd = gcd_extendido_binario(x, y//2)
      if b % 2 == 0:
          return a, b//2, gcd
      else:
          return (2*a+y)//2, (b-x)//2, gcd
  elif x > y:
      a, b, gcd = gcd_extendido_binario(y, x-y)
      return b, a-b, gcd
  else:
    a, b, gcd = gcd_extendido_binario(x, y-x)
    return  a-b, b, gcd

# %%
def legendre_symbol(a, p):
   "Calculo del simbolo de legendre visto en clase"
   if(a % p == 0):
    return 0
   ls = pow(a, (p - 1)//2, p)
   if ls == 1:
    return 1
   else:
    return -1

# %% [markdown]
# Calculamos un **no resto cuadrático** módulo p, en lugar de una raíz primitiva. Lo hacemos aleatoriamente porque la probabilidada de hallar que un candidato elegido al azar lo sea es muy alta.

# %%
import random
def resto_no_cuadratico(p):
  a = random.randint(2, p - 1)
  while(legendre_symbol(a, p) != -1):
    a = random.randint(2, p - 1)
  return a

# %%
def tonelli_shanks(a,p):
  "Implementación del algoritmo de tonelli shanks visto en clase"
  if legendre_symbol(a, p) != 1:
      return None
  r = p-1
  k = 0
  while(r % 2 == 0): # p = r*2^k+1
    k+=1
    r=r//2
  z = pow(a, r, p)
  #z = pow(a,r)
  x = pow(a,(r+1)//2,p)
  if (z == 1):
  #if ((z%p) ==1):
    return x
  s = k-1
  potencia = pow(2,s-1)
  potencia2 = pow(2, k -s -1)
  while(s>0):
    if(pow(z,potencia,p)!=1):
    #if(pow(z,pow(2,s-1),p)!=1):
      t = pow(resto_no_cuadratico(p), r * potencia2, p)
      #t = pow(resto_no_cuadratico(p), r * pow(2, k -s -1), p)
      #z = (z*pow(t,2))%p
      z = ((z%p) * pow(t, 2, p))%p
      #x = (x*t)%p
      x = ((x%p) * (t%p))%p
      if (z == 1):
        return x
    potencia2 = potencia2 << 1
    potencia = potencia >> 1
    s-=1
  return x

# %% [markdown]
# ## Función sqrt_mod

# %%
def sqrt_mod(a,p,n):
  if(a == 0):
    return 0
  lim = pow(p,n)
  a = a % lim
  if(n == 1):
    return tonelli_shanks(a,p)

  if(a % p == 0): #p divide a a
    s = 0
    z = a
    while(z % p == 0): # x^2 = a = p^s*z (p^n) && x = p^(s//2) * raíz(z)
      s += 1
      z = z // p
    # a = p^s * z
    if(s % 2 != 0):
      return None

    z0 = sqrt_mod(z,p,n)
    if(z0 == None):
      return None
    #return (pow(p, s // 2) * z0) % lim
    return (pow(p, s // 2, lim) * (z0%lim))%lim

  else: #p no divide a a
    x0 = tonelli_shanks(a,p)
    if(x0 == None):
      return None
    potencia = 1
    for k in range(2,n+1):
        potencia *= p
        t = (pow(x0,2)-a)//potencia
        l = (((gcd_extendido_binario((p-2)*x0,p)[0])%p)*(t%p))%p
        x0 = x0 + potencia * l
        #x0 = x0+pow(p,k-1)*l
    return x0

  return None

