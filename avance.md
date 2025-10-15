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

1. [x] Implementar la función miembro `aplicar` de la clase `condiciones_de_frontera_MDot.hpp::CF_MDot<Zero_Neumann_MDot>`.
2. [x] No olvidar la multiplicación por el vector normal en `condiciones_de_frontera_MDot.hpp::CF_MDot<Dirichlet_MDot>`.
3. [ ] Investigar una herramienta para pruebas unitarias.


# 9 septiembre 2025

1. Se implementó la función miembro `aplicar` de la clase `condiciones_de_frontera_MDot.hpp::CF_MDot<Zero_Neumann_MDot>`, añadiendo el vector normal. Se comprobó que el cálculo está bien.


# 10 septiembre 2025

1. Se retrasó lo del cálculo del flujo de masa debido a que primero se deben de calcular los coeficientes agrupados para la ecuación de momentum.
2. Se comenzó el cálculo de la conductancia difusiva, agregando los archivos `conductancia_difusiva.*`.

TODO (para la siguiente sesión):

1. [x] Codificar los constructores en `variables_discretizacion.cpp` y crear las variables en el `main` para uso posterior.
2. [x] Implementar el cálculo de la conductancia difusiva en `conductancia_difusiva.*`


# 11 septiembre 2025

1. Se crearon los archivos `ecuacion_momentum.*`, donde se definió una clase
   para facilitar el cálculo del campo de velocidad.
2. Se refactorizó la función miembro
   `malla_por_bloques.cpp::obtener_factores_de_interpolacion`, quitando a
   `push_back`.
3. Se refactorizó la función `esquemas_de_discretizacion.*::gradiente` para
   hacerla más legible.
4. Se añadieron `namespaces` en `esquemas_de_discretizacion.*` para diferenciar
   entre esquemas implícitos y explícitos.

TODO (para la siguiente sesión):

1. [x] Crear en la struct `Ecuacion_Momentum` de los archivos `ecuacion_momentum.*` una función miembro que calcule la conductancia difusiva.
2. [x] Investigar una herramienta para pruebas unitarias.


# 12 septiembre 2025

1. Se implementó la función miembro `Ecuacion_Momentum::ecuacion_momentum.*::calcular_conductancia_difusiva`.
2. Se investigó la herramienta Catch2.

TODO (para la siguiente sesión):

1. [x] Investigar a `std::variant` para poder elegir en tiempo de compilación y usarlo en la elección de solvers lineales.
2. [ ] Implementar función `esquemas_de_discretizacion.*::construccion_matriz_A_momentum`
3. [ ] Implementar una prueba unitaria con Catch2, [referencia](https://youtu.be/hNaj9AOGFGM?si=i51ppoGVgKTg71PU).


# 13 septiembre 2025

1. Se construyeron las funciones `esquemas_de_discretizacion.cpp::Discretizacion::construccion_matriz_A_momentum` y `esquemas_de_discretizacion.cpp::Discretizacion::construccion_coeficiente_b_momemtum`, aunque aún falta testeo.
2. Se usaron las funciones mencionadas en el punto uno en la clase `ecuacion_momentum.cpp::Ecuacion_Momentum`.
3. Se usó `std::variant` en `ecuacion_momentum.cpp::Ecuacion_Momentum`. Se implementó en `solvers_lineales.hpp::solverVariant`.


TODO (para la siguiente sesión):

1. [x] Revisar cuidadosamente los cálculos para la obtención del coeficiente A de la ecuación de momentum.


# 14 septiembre 2025

1. Se revisó el cálculo de la matriz de coeficientes A. No encontré nada raro, pero bueno, puede ser que esté mal algo que no haya notado y sea la raíz de un debuggeo del infierno.


TODO (para la siguiente sesión):

1. [x] Investigar como se trata la ecuación de corrección en las fronteras para cuando sean tipo Dirichlet y Zero Neumann.


# 15 septiembre 2025

1. Se investigó el tratamiento adecuado para las condiciones de frontera Dirichlet y Zero Neumann para la ecuación de corrección de presión.


# 16 septiembre 2025

1. Se elaboró la función para calcular el coeficiente $d$ en `flujo_de_masa.cpp::Coeficiente_d::calcular`, aunque aún no se ha probado.
2. Se elaboró una función de apoyo para llamar al solver lineal elegido en `ecuacion_momentum.hpp::Ecuacion_Momentum::resolver_con`.


TODO (para la siguiente sesión):

1. [x] Probar la función para el cálculo del coeficiente $d$ (`flujo_de_masa.cpp::Coeficiente_d::calcular`).
2. [x] Corregir la implementación del coeficiente $d$ en `flujo_de_masa.cpp::Coeficiente_d::calcular`, de tal forma que no haya coeficientes $d$ para las caras $w$ y $s$.


# 21 septiembre 2025

1. Se probó la función `flujo_de_masa.cpp::Coeficiente_d::calcular`. No le vi nada raro, pero es uno de los puntos donde revisar cuando venga el debuggeo del infierno.


# 22 septiembre 2025

1. Se comprobó que el coeficiente $d_{n} = 0$ en la frontera correspondiente para la ecuación de corrección de presión, por lo que no será necesario el uso del coeficiente $a_{N}$ en las fronteras para ninguna de las condiciones de frontera.

TODO (para la siguiente sesión):

1. [x] Comprobar que para la condición de frontera tipo Dirichlet para la ecuación de **corrección de presión** el coeficiente $d_{n}$ ($^{V_{n}}/_{a_{n}}$) es cero, dado que si se interpola se requiere el valor de $V_{N}$, y este es cero porque $\Delta_x$ o $\Delta_y$ pueden ser cero. En dado caso de ser cierto entonces será necesario hacer $a_{N} = 0$ en la frontera $N$ (se considera a $N$ como una frontera genérica, entonces no necesariamente corresponde con la frontera norte).


# 24 septiembre 2025

1. Se empezó la construcción de la función `flujo_de_masa.cpp::calcular_flujo_de_masa`. Falta por probar y hacer el cálculo para la dirección $y$.
2. Se llegó a la conclusión de que pasar una instancia de `Campo::velFace` para crear un objeto del tipo `Ecuacion_Momentum` no es necesario, pues se puede obtener a través de este sin mayor inconveniente (eso creo).


# 26 septiembre 2025

1. Se concluyó la implementación del cálculo de flujo de masa en `flujo_de_masa.hpp::calcular_flujo_de_masa`. Esta es una de las partes cruciales que se deben de revisar cuando comience el debuggeo del infierno.


TODO (para la siguiente sesión):

1. [x] Revisar de nuevo la implementación del cálculo del flujo de masa en `flujo_de_masa.hpp::calcular_flujo_de_masa`.
2. [x] Comenzar con la implementación de la ecuación de corrección de presión, con base en lo aprendido del coeficiente "d".


# 27 septiembre 2025

1. Se revisó un poco la implementación de `flujo_de_masa.hpp::calcular_flujo_de_masa` para no perder tiempo. Luego veré si eso da problemas.
2. Se comenzó con la implementación de `ecuacion_presion.*` para calcular la presión de corrección.


# 30 septiembre 2025

1. Se culminó la implementación de la función `esquemas_de_discretizacion.hpp::laplaciano_lineal_presion`
2. Se comenzó con la implementación de la función `esquemas_de_discretizacion.cpp::construccion_matriz_A_presion`


TODO (para la siguiente sesión):

1. [x] Terminar la implementación de `ecuacion_presion.*`.
2. [x] Terminar la implementación de de la función `esquemas_de_discretizacion.cpp::construccion_matriz_A_presion`


# 3 octubre 2025

1. Se implementó la función `esquemas_de_discretizacion.cpp::construccion_matriz_A_presion`.
2. Se implementó a `ecuacion_presion.cpp::Ecuacion_Presion::resolver`.


TODO (para la siguiente sesión):

1. [x] Probar la implementación de `ecuacion_presion.cpp::Ecuacion_Presion::resolver`.
2. [x] Comenzar con la corrección de los campos.


# 6 octubre 2025

1. No se probó la implementación de `ecuacion_presion.cpp::Ecuacion_Presion::resolver`, pero se revisó lo suficiente como para confiar de que está bien (eso espero). Esta es otra de las partes cruciales para revisar cuando comienzo el debuggeo del infierno.
2. Se comenzó con la corrección de los campos en `correccion_campos.*`.
3. Se implementó un `for_each` para el cálculo de la corrección de campos en `correccion_campos.cpp::corregir`.


TODO (para la siguiente sesión):

1. [x] Terminar la implementación de la corrección de los campos en `correccion_campos.*`.


# 6 octubre 2025 (segunda parte)

1. Se terminó de implementar la corrección de campos (`correccion_campos.cpp::corregir`). Además se utilizaron flags de pre procesamiento para elegir entre un bucle en paralelo o serial.


TODO (para la siguiente sesión):

1. [x] Revisar la implementación de la corrección de los campos en `correccion_campos.cpp::corregir` en busca de errores.
2. [x] Comenzar con la función para verificar la convergencia.


# 8 octubre 2025

1. Se modificó a `CMakeLists.txt` con ayuda de chatGPT para definir los siguientes modos de construcción: `Debug`, `Release` y `RelWithDebInfo`.
2. Se revisó la implementación de `correccion_campos.cpp::corregir`, esperando que esté bien. De todos modos, esta es una de las partes que hay que revisar en el debuggeo del infierno.
3. Se terminó la función para verificar la convergencia en `convergencia.*`.


TODO (para la siguiente sesión):

1. [x] Revisar la implementación de `convergencia.*`.
2. [x] Implementar una función para la reasignación de campos para la nueva iteración.


# 13 octubre 2025

1. Se revisó la implementación de `convergencia.*`, pero no lo suficiente como para descartarla cuando comience el debuggeo del infierno.
2. Se implementó la función `reasignacion.*::reasignar`.


TODO (para la siguiente sesión):

1. [x] Preparar todo para la primer corrida del SIMPLE.


# 14 octubre 2025

1. Comenzó el debuggeo del infierno.


TODO (para la siguiente sesión):

1. [x] Revisar la función `reasignacion.*::reasignar`.
2. [ ] Revisar a `ecuacion_presion.*` en busca de bugs.
