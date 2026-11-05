(* ::Package:: *)

(* RepresentacionFase.wl
   Oscilador cuántico abierto: Lindblad en la base de Fock, estado estacionario en el
   espacio de Fock–Liouville (FLS) y representación de fase de las notas de los tutores
   (|θ〉 = Σ e^{inθ}|n〉, K(θ,θ') = 〈θ|ρ|θ'〉, D_θ = i∂θ, D_θ' = -i∂θ').

   Convenios:
   - ħ = 1. Disipador habitual: D[L]ρ = LρL† - (L†Lρ + ρL†L)/2, con la tasa fuera de L.
   - Truncado a d niveles |0〉…|d-1〉. Los productos se hacen con las matrices ya truncadas
     (a a† vale 0 en |d-1〉), así que el modelo truncado conserva la traza.
   - vec(ρ) apila columnas: vec(AρB) = (Bᵀ⊗A) vec(ρ).
   - Los parámetros simbólicos se tratan como reales (daga no conjuga símbolos).
   - Operadores: se escriben con a, ad y ** (NonCommutativeMultiply), p. ej. Δ ad**a.
     Ojo: E es el número e en Mathematica; para la amplitud del forzamiento usa ℰ.
   - Desde el 2026-11-05: centro y diferencia (Θ, φ), sectores k y corte |k| ≤ K
     (nucleoCD, vectorSector, nucleoSector, ordenPorSectores, posicionesBanda,
     estacionarioEnBanda). Lo anterior se revisó y comprobó ese día sin cambios. *)

BeginPackage["RepresentacionFase`"];

a::usage = "a: operador de destrucción. Los productos se escriben con **, p. ej. ad**a.";
ad::usage = "ad: operador de creación a†.";
\[Theta]::usage = "θ: ángulo del ket en el núcleo K(θ, θp).";
\[Theta]p::usage = "θp: ángulo del bra en el núcleo K(θ, θp).";

matrizFock::usage = "matrizFock[op, d] da la matriz d×d del operador op (escrito con a, ad, **) en la base |0〉…|d-1〉.";
daga::usage = "daga[op] da el adjunto de op: invierte los productos, cambia a↔ad y conjuga los números complejos (los símbolos se suponen reales).";
salto::usage = "salto[k] es el operador de salto que cambia el número de fonones en k: salto[1] = ad, salto[-1] = a, salto[-2] = a**a.";

liouvilliano::usage = "liouvilliano[H, saltos, d] da la matriz d²×d² de Lindblad en el espacio de Fock–Liouville. H es un operador (con a, ad, **) o una matriz. saltos = {{γ1, L1}, {γ2, L2}, …}; cada L es un operador, una matriz o un entero k (= salto[k]).";
vecFLS::usage = "vecFLS[ρ] apila las columnas de ρ (vectorización de Fock–Liouville).";
matrizDeVec::usage = "matrizDeVec[v, d] deshace vecFLS.";
etiquetasFLS::usage = "etiquetasFLS[d] da los nombres ρnm en el orden de vecFLS.";
ordenPorK::usage = "ordenPorK[d] da las posiciones de vecFLS ordenadas por bloques k = n - m (k = 0, 1, 2, …, -1, -2, …). ℒ[[o, o]] con o = ordenPorK[d] muestra la estructura de bloques.";
estadoEstacionario::usage = "estadoEstacionario[ℒ, d] resuelve ℒ·vec(ρ) = 0 con Tr ρ = 1 (sustituye la fila de ρ00 por la de la traza) y devuelve ρ.";

nucleo::usage = "nucleo[ρ] = K(θ, θp) = Σ ρnm e^{-inθ} e^{imθp}.";
matrizDeNucleo::usage = "matrizDeNucleo[K, d] recupera ρnm a partir del núcleo K (inversa de nucleo).";
distribucionFase::usage = "distribucionFase[ρ] = P(θ) = K(θ, θ)/(2π), escrita en forma real.";
funcionD\[Theta]::usage = "funcionDθ[f, K] aplica f(D_θ) a K, con D_θ = i∂θ definido espectralmente: e^{-inθ} → f(n) e^{-inθ}.";
funcionD\[Theta]p::usage = "funcionDθp[f, K] aplica f(D_θp) a K, con D_θp = -i∂θp: e^{imθp} → f(m) e^{imθp}.";
accionIzq::usage = "accionIzq[op, K, d] es el núcleo de op·ρ (ec. 41 de las notas: aρ → e^{iθ}√D_θ K, a†ρ → e^{-iθ}√(D_θ+1) K). Con d = Infinity no trunca.";
accionDer::usage = "accionDer[op, K, d] es el núcleo de ρ·op (ec. 42: ρa → e^{iθp}√(D_θp+1) K, ρa† → e^{-iθp}√D_θp K). Con d = Infinity no trunca.";
lindbladFase::usage = "lindbladFase[K, H, saltos, d] = ∂t K: la ecuación de Lindblad escrita en la representación de fase con el diccionario de las notas.";
liouvillianoDesdeFase::usage = "liouvillianoDesdeFase[H, saltos, d] construye la matriz de FLS desde la representación de fase: aplica lindbladFase a cada modo de Fourier e^{-inθ}e^{imθp} y lee sus coeficientes (columna por columna). Tiene que coincidir con liouvilliano[H, saltos, d].";
comprobarDiccionario::usage = "comprobarDiccionario[H, saltos, d] aplica lindbladFase al núcleo de una ρ genérica y lo compara con el Liouvilliano de Fock. Devuelve True si coinciden.";

\[CapitalTheta]::usage = "Θ = (θ + θp)/2: ángulo centro (la fase).";
\[CurlyPhi]::usage = "φ = θ - θp: ángulo diferencia (el conjugado del número).";
nucleoCD::usage = "nucleoCD[ρ] = K(Θ, φ) = Σ ρnm e^{-i(n-m)Θ} e^{-i(n+m)φ/2}: el núcleo con θ = Θ + φ/2 y θp = Θ - φ/2.";
vectorSector::usage = "vectorSector[ρ, k] = x_k: la diagonal k de ρ en columna, (ρ_{k,0}, ρ_{k+1,1}, …) = (c_0, c_1, …) si k ≥ 0, y (ρ_{0,|k|}, ρ_{1,|k|+1}, …) si k < 0.";
nucleoSector::usage = "nucleoSector[ρ, k] = K_k(φ) = Σ_n c_n e^{-i(n+|k|/2)φ}, con c_n = vectorSector[ρ, k]. Cumple K(Θ, φ) = Σ_k e^{-ikΘ} K_k(φ).";
ordenPorSectores::usage = "ordenPorSectores[d] da las posiciones de vecFLS ordenadas como en la §6 de las notas: k = -(d-1), …, d-1 y, dentro de cada sector, n creciente. ℒ[[o, o]] con o = ordenPorSectores[d] es la matriz tridiagonal por bloques.";
posicionesBanda::usage = "posicionesBanda[d, K] da las posiciones de vecFLS con |n - m| ≤ K (el corte en k), en el mismo orden que ordenPorSectores.";
estacionarioEnBanda::usage = "estacionarioEnBanda[M, d, K] resuelve M·x = 0 con Tr ρ = 1, donde M es la matriz restringida a |n - m| ≤ K y escrita en el orden de posicionesBanda[d, K]. Devuelve ρ (d×d, con ceros fuera de la banda). Trabaja con números.";

matrizFock::orden = "Hay un producto de operadores escrito con * (conmutativo). Usa ** para fijar el orden.";
matrizFock::forma = "No sé traducir `1` a matriz.";

Begin["`Private`"];

opQ[x_] := ! FreeQ[x, a | ad];
conjR[x_] := x /. Complex[re_, im_] :> Complex[re, -im];
dagaM[m_] := Transpose[conjR[m]];

(* ---------- Base de Fock ---------- *)

matA[d_] := Table[If[j == i + 1, Sqrt[i], 0], {i, d}, {j, d}];

separarEscalar[x_Times] := {Times @@ Select[List @@ x, ! opQ[#] &], Select[List @@ x, opQ]};

matrizFock[x_, d_Integer] := Which[
  ! opQ[x], x IdentityMatrix[d],
  x === a, matA[d],
  x === ad, Transpose[matA[d]],
  Head[x] === Plus, Total[matrizFock[#, d] & /@ List @@ x],
  Head[x] === NonCommutativeMultiply, Dot @@ (matrizFock[#, d] & /@ List @@ x),
  Head[x] === Times, With[{so = separarEscalar[x]},
    If[Length[so[[2]]] > 1, Message[matrizFock::orden]; $Failed, so[[1]] matrizFock[First[so[[2]]], d]]],
  Head[x] === Power && IntegerQ[x[[2]]] && x[[2]] > 0, MatrixPower[matrizFock[x[[1]], d], x[[2]]],
  True, Message[matrizFock::forma, x]; $Failed];

daga[x_] := Which[
  ! opQ[x], conjR[x],
  x === a, ad,
  x === ad, a,
  Head[x] === Plus || Head[x] === Times, daga /@ x,
  Head[x] === NonCommutativeMultiply, NonCommutativeMultiply @@ Reverse[daga /@ List @@ x],
  Head[x] === Power, daga[x[[1]]]^x[[2]],
  True, $Failed];

salto[1] = ad;
salto[-1] = a;
salto[k_Integer?Positive] := NonCommutativeMultiply @@ ConstantArray[ad, k];
salto[k_Integer?Negative] := NonCommutativeMultiply @@ ConstantArray[a, -k];

operadorDeSalto[l_] := If[IntegerQ[l], salto[l], l];
aMatriz[x_, d_] := If[MatrixQ[x], x, matrizFock[x, d]];

(* ---------- Espacio de Fock–Liouville ---------- *)

vecFLS[r_?MatrixQ] := Flatten[Transpose[r]];
matrizDeVec[v_List, d_Integer] := Transpose[Partition[v, d]];
etiquetasFLS[d_Integer] := Flatten[Transpose[Table["\[Rho]" <> ToString[n] <> ToString[m], {n, 0, d - 1}, {m, 0, d - 1}]]];

(* posición en vecFLS de ρnm: n + d m + 1 *)
ordenPorK[d_Integer] := Module[{pares = Flatten[Table[{n, m}, {m, 0, d - 1}, {n, 0, d - 1}], 1]},
  Ordering[{If[#[[1]] - #[[2]] >= 0, 0, 1], Abs[#[[1]] - #[[2]]], #[[2]]} & /@ pares]];

superDisipador[L_, id_] := With[{LdL = dagaM[L] . L},
  KroneckerProduct[conjR[L], L] - 1/2 KroneckerProduct[id, LdL] - 1/2 KroneckerProduct[Transpose[LdL], id]];

liouvilliano[h_, saltos_List, d_Integer] := Module[{H = aMatriz[h, d], id = IdentityMatrix[d]},
  -I (KroneckerProduct[id, H] - KroneckerProduct[Transpose[H], id]) +
   Total[Table[s[[1]] superDisipador[aMatriz[operadorDeSalto[s[[2]]], d], id], {s, saltos}]]];

estadoEstacionario[L_?MatrixQ, d_Integer] := Module[{A = L, x},
  A[[1]] = vecFLS[IdentityMatrix[d]];          (* fila de ρ00 → Tr ρ = 1 *)
  x = LinearSolve[A, UnitVector[d^2, 1]];
  Map[Together, matrizDeVec[x, d], {2}]];

(* ---------- Representación de fase ----------
   Por dentro, K se escribe como polinomio en u = e^{-iθ} y v = e^{iθp}: el término
   ρnm e^{-inθ} e^{imθp} es ρnm u^n v^m. Así D_θ cuenta la potencia de u y D_θp la de v. *)

aUV[K_] := Expand[(TrigToExp[K] /. {\[Theta] -> I Log[uu], \[Theta]p -> -I Log[vv]}) /.
     Power[E, x_] :> PowerExpand[Power[E, Expand[x]]]];
desdeUV[p_] := p /. {uu -> Exp[-I \[Theta]], vv -> Exp[I \[Theta]p]};

(* f(D) espectral sobre la variable w (uu o vv) *)
fUV[f_, p_, w_] := Module[{q = Expand[p], nmax},
  nmax = Exponent[q, w];
  If[q === 0 || nmax < 0, 0, Sum[f[n] Coefficient[q, w, n] w^n, {n, 0, nmax}]]];

coefsUV[p_, d_Integer] := With[{q = Expand[p]},
  Table[Coefficient[Coefficient[q, uu, n], vv, m], {n, 0, d - 1}, {m, 0, d - 1}]];
desdeCoefs[c_] := Sum[c[[n + 1, m + 1]] uu^n vv^m, {n, 0, Length[c] - 1}, {m, 0, Length[c] - 1}];

(* reglas de un solo operador; con d finito, lo que sube por encima de |d-1〉 se anula *)
debajo[k_, d_] := If[d === Infinity, True, k < d];
izq1[a, p_, d_] := Expand[uu^-1 fUV[Sqrt[#] &, p, uu]];
izq1[ad, p_, d_] := Expand[uu fUV[If[debajo[# + 1, d], Sqrt[# + 1], 0] &, p, uu]];
der1[a, p_, d_] := Expand[vv fUV[If[debajo[# + 1, d], Sqrt[# + 1], 0] &, p, vv]];
der1[ad, p_, d_] := Expand[vv^-1 fUV[Sqrt[#] &, p, vv]];

(* op·ρ : el factor de más a la derecha actúa primero *)
izqUV[x_, p_, d_] := Which[
  MatrixQ[x], desdeCoefs[x . coefsUV[p, Length[x]]],
  ! opQ[x], x p,
  x === a || x === ad, izq1[x, p, d],
  Head[x] === Plus, Total[izqUV[#, p, d] & /@ List @@ x],
  Head[x] === Times, With[{so = separarEscalar[x]}, so[[1]] izqUV[First[so[[2]]], p, d]],
  Head[x] === NonCommutativeMultiply, Fold[izqUV[#2, #1, d] &, p, Reverse[List @@ x]],
  Head[x] === Power, Nest[izqUV[x[[1]], #, d] &, p, x[[2]]],
  True, $Failed];

(* ρ·op : el factor de más a la izquierda actúa primero *)
derUV[x_, p_, d_] := Which[
  MatrixQ[x], desdeCoefs[coefsUV[p, Length[x]] . x],
  ! opQ[x], x p,
  x === a || x === ad, der1[x, p, d],
  Head[x] === Plus, Total[derUV[#, p, d] & /@ List @@ x],
  Head[x] === Times, With[{so = separarEscalar[x]}, so[[1]] derUV[First[so[[2]]], p, d]],
  Head[x] === NonCommutativeMultiply, Fold[derUV[#2, #1, d] &, p, List @@ x],
  Head[x] === Power, Nest[derUV[x[[1]], #, d] &, p, x[[2]]],
  True, $Failed];

nucleo[r_?MatrixQ] := With[{d = Length[r]},
  Sum[r[[n + 1, m + 1]] Exp[-I n \[Theta]] Exp[I m \[Theta]p], {n, 0, d - 1}, {m, 0, d - 1}]];
matrizDeNucleo[K_, d_Integer] := coefsUV[aUV[K], d];
distribucionFase[r_?MatrixQ] := Simplify[ComplexExpand[(nucleo[r] /. \[Theta]p -> \[Theta])/(2 Pi)]];

funcionD\[Theta][f_, K_] := desdeUV[fUV[f, aUV[K], uu]];
funcionD\[Theta]p[f_, K_] := desdeUV[fUV[f, aUV[K], vv]];
accionIzq[op_, K_, d_: Infinity] := desdeUV[izqUV[op, aUV[K], d]];
accionDer[op_, K_, d_: Infinity] := desdeUV[derUV[op, aUV[K], d]];

lindbladFaseUV[p_, h_, saltos_List, d_] := Expand[
  -I (izqUV[h, p, d] - derUV[h, p, d]) +
   Sum[With[{g = s[[1]], L = operadorDeSalto[s[[2]]]},
     With[{Ld = If[MatrixQ[L], dagaM[L], daga[L]]},
      With[{LdL = If[MatrixQ[L], Ld . L, Ld ** L]},
       g (izqUV[L, derUV[Ld, p, d], d] - 1/2 izqUV[LdL, p, d] - 1/2 derUV[LdL, p, d])]]], {s, saltos}]];

lindbladFase[K_, h_, saltos_List, d_: Infinity] := desdeUV[lindbladFaseUV[aUV[K], h, saltos, d]];

(* columna (n,m) = coeficientes de Fourier de ∂t K cuando K = e^{-inθ} e^{imθp}, en el orden de vecFLS *)
liouvillianoDesdeFase[h_, saltos_List, d_Integer] := Transpose[Flatten[
   Table[vecFLS[coefsUV[lindbladFaseUV[uu^n vv^m, h, saltos, d], d]], {m, 0, d - 1}, {n, 0, d - 1}], 1]];

comprobarDiccionario[h_, saltos_List, d_Integer] := Module[{r, rho, fase, fock},
  rho = Table[r[n, m], {n, 0, d - 1}, {m, 0, d - 1}];
  fase = coefsUV[lindbladFaseUV[aUV[nucleo[rho]], h, saltos, d], d];
  fock = matrizDeVec[liouvilliano[h, saltos, d] . vecFLS[rho], d];
  Simplify[fase - fock] === ConstantArray[0, {d, d}]];

(* ---------- Centro y diferencia: Θ = (θ+θp)/2, φ = θ - θp ----------
   ρnm e^{-inθ} e^{imθp} = ρnm e^{-ikΘ} e^{-i(n+m)φ/2}, con k = n - m. El sector k
   es la diagonal k de ρ; dentro de él, \hat N = i∂φ vale n + |k|/2 sobre c_n. *)

nucleoCD[r_?MatrixQ] := With[{d = Length[r]},
  Sum[r[[n + 1, m + 1]] Exp[-I (n - m) \[CapitalTheta]] Exp[-I (n + m) \[CurlyPhi]/2], {n, 0, d - 1}, {m, 0, d - 1}]];

vectorSector[r_?MatrixQ, k_Integer] := With[{d = Length[r]},
  If[k >= 0, Table[r[[n + k + 1, n + 1]], {n, 0, d - 1 - k}],
   Table[r[[n + 1, n - k + 1]], {n, 0, d - 1 + k}]]];

nucleoSector[r_?MatrixQ, k_Integer] := With[{c = vectorSector[r, k]},
  Sum[c[[n + 1]] Exp[-I (n + Abs[k]/2) \[CurlyPhi]], {n, 0, Length[c] - 1}]];

(* posición en vecFLS de ρnm: n + d m + 1. Sector k ≥ 0: ρ_{n+k,n}; k < 0: ρ_{n,n-k} *)
posicionSector[d_, k_, n_] := If[k >= 0, (n + k) + d n + 1, n + d (n - k) + 1];
posicionesBanda[d_Integer, K_Integer] := Flatten[
  Table[posicionSector[d, k, n], {k, -K, K}, {n, 0, d - 1 - Abs[k]}]];
ordenPorSectores[d_Integer] := posicionesBanda[d, d - 1];

estacionarioEnBanda[M_?MatrixQ, d_Integer, K_Integer] := Module[
  {pos = posicionesBanda[d, K], A = N[M], fila00, x, v},
  fila00 = First[FirstPosition[pos, 1]];                         (* ρ00 *)
  A[[fila00]] = If[MemberQ[Table[n + d n + 1, {n, 0, d - 1}], #], 1., 0.] & /@ pos;  (* Tr ρ = 1 *)
  x = LinearSolve[A, UnitVector[Length[pos], fila00]];
  v = ConstantArray[0., d^2]; v[[pos]] = x;
  Chop[matrizDeVec[v, d]]];

End[];
EndPackage[];
