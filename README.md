# Sincronización cuántica en la representación de fase

*Código de Mathematica del Trabajo de Fin de Máster*

![TFM](https://img.shields.io/badge/TFM-UCM%202026--27-1f6feb)
![Wolfram Mathematica](https://img.shields.io/badge/Wolfram-Mathematica-dd1100?logo=wolframmathematica&logoColor=white)

**Mario Córdoba Burgueño**, con **Gloria Platero** y **Álvaro Gómez-León**<br>
Máster en Física Teórica, Universidad Complutense de Madrid, 2026–27<br>
Título oficial: *Low Dimensional Topological Lattices for the Transfer of Quantum Information*

---

El objetivo es resolver el oscilador cuántico de Van der Pol (y después varios
osciladores acoplados) con la ecuación de Lindblad **directamente en la representación
de fase**, con el núcleo $K(\theta,\theta') = \langle\theta|\rho|\theta'\rangle$, y
caracterizar desde ahí la sincronización. Lo numérico en la base de Fock sirve de
comprobación.

## Qué hay

| Archivo | Qué es |
|---|---|
| `Mathematica/RepresentacionFase.wl` | Paquete: liouvilliano en la base de Fock (espacio de Fock–Liouville), estado estacionario y paso a la representación de fase: $K(\theta,\theta')$, $P(\theta)$, $D_\theta$. Desde el 2026-11-12, también el diccionario de fase para dos osciladores. |
| `Mathematica/vdp_sin_forzamiento.nb` | VdP sin forzamiento: estacionario y distribución de fase. |
| `Mathematica/vdp_forzado.nb` | VdP forzado en Fock: $P(\theta)$, lengua de Arnold y núcleo $K$. |
| `Mathematica/vdp_funcion_generatriz.nb` | Sector $k=0$ sin truncar Fock: función generatriz $G(z)$, ecuación de Kummer y comparación con Fock. |
| `Mathematica/vdp_forzado_fase.nb` | El caso forzado escrito directamente en la representación de fase. |
| `Mathematica/vdp_corte_k.nb` | Corte en sectores $\vert k\vert \le K$: estacionario en la banda y cómo cambia la distribución de fase con el corte. |
| `Mathematica/vdp_dos_osciladores_fase.nb` | Dos VdP con acoplo disipativo $\gamma_c\mathcal{D}[a_1+a_2]$ solo en la representación de fase, sin Fock: ecuaciones desde el diccionario de fase, límites con solución cerrada y estacionario con $r=0.05$ como serie en $r$; distribución de la fase relativa, núcleo y lenguas de Arnold con frecuencias distintas. |
| `Mathematica/figuras/` | Figuras exportadas de los cuadernos. |

## Cómo se usa

Cada cuaderno carga en su sección 0 el paquete `RepresentacionFase.wl` de su misma
carpeta. Basta con abrirlo y evaluarlo.

