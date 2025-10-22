# Solver de flujo de fluidos sencillos

Este solver posee la capacidad de solucionar problemas de flujo de fluidos
newtonianos, laminares, con transferencia de calor en geometrías rectangulares
con un número arbitrario de parches.

## jueves 21 agosto 2025

**NOTE**: el polimorfismo en tiempo de ejecución volvió horrible a mi código,
muy complicado y difícil de mantener, por lo que decidí reescribirlo usando
plantillas (polimorfismo en tiempo de compilación) y usando el paradigma
orientada a datos.

## lunes 8 agosto 2025

**NOTE** (xaman pag. 232 PDF): Para cuando toque el desarrollo de la ecuación
de corrección de presión:
- En Dirichlet se hace todo el calculo normal pero
P_{b}^{'} (presión en la frontera) es cero, haciendo que a_{N} sea cero del
lado derecho de la igualdad pero no en su contribución al coeficiente central
a_{C}.
- En Zero Neumann se sabe cuál es el flujo de masa desde el principio, por lo
que no hay contribución al coeficiente central por parte de a_{N} porque este
nunca se forma.
