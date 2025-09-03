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

TODO (para la siguiente sesión):

1. [ ] Implementar la función `condiciones_de_frontera_MDot.hpp::asignar_condiciones_de_frontera_MDot`
