# 30 agosto 2025

1. Se declaró la función `construir_CF_flujo_de_masa` en
`condiciones_de_frontera_MDot.hpp` con el objetivo de construir las condiciones
de frontera para el flujo de masa.

2. Se declaró el `struct` `Parches_Flujo_de_Masa` en
`condiciones_de_frontera_MDot.hpp` para tener una copia de los parches
obtenidos de la `Malla::Mallador` de los archivos `malla_por_bloques.*` y así
poder manipular los nodos frontera a conveniencia.

3. En el `main.cpp` se comenzó a probar lo de las copias hacia el `struct`
`Parches_Flujo_de_Masa`, luego de eso sería conveniente ponerlo todo en una
función.

TODO (para la siguiente sesión):

1. [x] Implementar funciones miembro del `struct` `condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa`.
2. [x] Hacer que el código corra y arroje resultados.


# 31 agosto 2025

1. Se implementaron las funciones miembro en
`condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa` de forma exitosa.

2. El código ya compila y corre sin problemas (por ahora).

TODO (para la siguiente sesión):

1. [x] Implementar la función `condiciones_de_frontera_MDot.hpp::construir_CF_flujo_de_masa`


# 1 septiembre 2025

1. Se implementó la función `asignar_condiciones_de_frontera_MDot.hpp::construir_CF_flujo_de_masa`

2. Se comenzó a realizar la función `condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`


# 2 septiembre 2025

1. Se implementó la función miembro
`condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa::añadir_tipo_de_CF`
para poder almacenar el tipo de condición de frontera en cada objeto de tipo
`Parches_Flujo_de_Masa`. Esto se hizo para hacer más sencilla a la función que
se implementará después:
`condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`


# 3 septiembre 2025

1. Se comenzó a implementar la función
`condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`,
dándome cuenta que era necesario redefinir el constructor de la clase
`condiciones_de_frontera_MDot.hpp::CF_MDot<Dirichlet_MDot>`.

TODO (para la siguiente sesión):

1. [x] Pre calcular el `desfase` en el struct `condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa`
2. [x] Dividir `listas` en `listas_verticales` y `listas_horizontales` en `condiciones_de_frontera_MDot.hpp::construir_CF_flujo_de_masa`


# 5 septiembre 2025

1. Se implementó la función miembro `condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa::calcular_desfase`
2. Se implementó una lista para los parches verticales y otra para los horizontales en `condiciones_de_frontera_MDot.hpp::construir_CF_flujo_de_masa`, con el objetivo de evitar pasar velocidades que no sean normales a las caras

TODO (para la siguiente sesión):

1. [x] Terminar de implementar la función `condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`
2. [x] Obtener las deltas en el struct `condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa::calcular_desfase`


# 7 septiembre 2025

1. Se eliminó al flujo de masa del struct `Campo.hpp::Momentum` porque no era necesario tenerlo, ya que el flujo de masa está implícito en los coeficientes agrupados `A_coef`.
2. Se implementó la función miembro `condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`.
3. Se obtuvieron las deltas en el struct `condiciones_de_frontera_MDot.hpp::Parches_Flujo_de_Masa::calcular_desfase` con éxito.


# 8 septiembre 2025

1. Se corrigió el código de la función `condiciones_de_frontera_MDot.cpp::Parches_Flujo_de_Masa::cortar_nodos_esquina` haciendo un bucle en reversa para evitar saltos de elementos ocasionado por el anterior enfoque. Si tienes dudas de lo que hablo ponte a hacerlo a mano, flojo.
2. Se especificó en `main.cpp` que es necesario robustecer el código relacionado con el struct `Parches_Flujo_de_Masa`, ya que las funciones miembro deben de ponerse en un determinado orden para que se obtenga lo deseado, lo cual hace que el código sea muy frágil.
3. Se eliminaron las variables locales `listas_verticales` y `listas_horizontales` en `condiciones_de_frontera_MDot.cpp::construir_CF_flujo_de_masa` porque ocurrían errores de segmentación de memoria, pues, al salir dejar de usar la función estas se destruían, causando un comportamiento indefinido.
4. Se implementó con éxito la función miembro `aplicar` de `condiciones_de_frontera_MDot.hpp::CF_MDot<Dirichlet_MDot>`.

TODO (para la siguiente sesión):

1. [ ] Implementar la función miembro `aplicar` de la clase `condiciones_de_frontera_MDot.hpp::CF_MDot<Zero_Neumann_MDot>`.
2. [ ] Investigar una herramienta para pruebas unitarias.
