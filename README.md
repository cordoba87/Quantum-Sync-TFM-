# Sincronización cuántica en la representación de fase

*Código y notas del Trabajo de Fin de Máster*

![TFM](https://img.shields.io/badge/TFM-UCM%202026--27-1f6feb)
![Wolfram Mathematica](https://img.shields.io/badge/Wolfram-Mathematica-dd1100?logo=wolframmathematica&logoColor=white)
![LaTeX](https://img.shields.io/badge/LaTeX-biblatex%20%2B%20biber-008080?logo=latex&logoColor=white)

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
| `Mathematica/RepresentacionFase.wl` | Paquete: liouvilliano en la base de Fock (espacio de Fock–Liouville), estado estacionario y paso a la representación de fase: $K(\theta,\theta')$, $P(\theta)$, $D_\theta$. |
| `Mathematica/vdp_sin_forzamiento.nb` | VdP sin forzamiento: estacionario y distribución de fase. |
| `Mathematica/vdp_forzado.nb` | VdP forzado en Fock: $P(\theta)$, lengua de Arnold y núcleo $K$. |
| `Mathematica/vdp_funcion_generatriz.nb` | Sector $k=0$ sin truncar Fock: función generatriz $G(z)$, ecuación de Kummer y comparación con Fock. |
| `Mathematica/vdp_forzado_fase.nb` | El caso forzado escrito directamente en la representación de fase. |
| `Mathematica/vdp_corte_k.nb` | Corte en sectores $\vert k\vert \le K$: estacionario en la banda y cómo cambia la distribución de fase con el corte. |
| `Mathematica/figuras/` | Figuras exportadas de los cuadernos. |
| `Notas/notas-tutores.tex` | Notas para los tutores. |
| `Notas/apuntes-vdp-fase.tex` | Apuntes de trabajo: el VdP en la representación de fase. |
| `Notas/bibliografia.bib` | Bibliografía, exportada de Zotero (no se edita a mano). |

## Cómo se usa

- **Mathematica.** Cada cuaderno carga en su sección 0 el paquete
  `RepresentacionFase.wl` de su misma carpeta. Basta con abrirlo y evaluarlo.
- **LaTeX.** pdflatex → biber → pdflatex ×2. En VS Code (LaTeX Workshop) la receta
  está en `.vscode/settings.json` y la salida va a `build/`.

Los PDF de artículos y libros, los compilados de LaTeX y el material de trabajo
que no es para los tutores no se suben: la lista está en `.gitignore`.

## Guardar cambios

```bash
git add .
git commit -m "Qué he hecho hoy"
git push
```

Si se trabaja en otro ordenador, antes de empezar: `git pull`.

## Citar

GitHub ofrece la cita en el botón **Cite this repository** (sale de `CITATION.cff`).
En LaTeX, la entrada está en `Notas/bibliografia.bib` con la clave
`github2026TFM`.
