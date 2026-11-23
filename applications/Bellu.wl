(* Bellu.wl - Tick co-infection models from Belluccini et al. *)
(* Two viral strains V1 (resident) and V2 (invasive) in tick population *)

(*
\begin{table}[h!]
\centering
\caption{Models in Bellu.wl}
\begin{tabular}{|c|c|c|l|}
\hline
\textbf{Model} & \textbf{RHS} & \textbf{Variables} & \textbf{Description} \\
\hline\hline
1 & RHS & $m_0, m_1, m_2, m_c$ & Basic co-feeding \\
\hline
2 & RHS2 & $m_0, m_1, m_2, m_c$ & Within-host ($\phi$ factors) \\
\hline
3 & RHS3 & $m_0, m_1, m_2, M_1, M_2, m_c$ & Alizon \\
\hline
4 & RHS4 & $m_0, m_1, m_2, M_1, M_2, m_c$ & Two-slot \\
\hline
\end{tabular}
\end{table}

Two-letter abbreviations:
\begin{itemize}
\item n0, n1, n2, nc, nM1, nM2 $\to$ $\nu_0, \nu_1, \nu_2, \nu_c, \nu_{M_1}, \nu_{M_2}$
\item a1, a2 $\to$ $\alpha_1, \alpha_2$
\item d1, d2 $\to$ $\delta_1, \delta_2$
\item ec $\to$ $\epsilon_c$
\item ph $\to$ $\Phi$
\item La1, La2, Lac $\to$ $\Lambda_1, \Lambda_2, \Lambda_c$
\item l1c, l2c $\to$ $\lambda_{1,c}, \lambda_{2,c}$
\item p21, p12 $\to$ $\phi_{2|1}, \phi_{1|2}$
\end{itemize}
*)

(* ===== Model 1: Basic co-feeding transmission (Eq. 1) ===== *)
(* Variables: m0=susceptible, m1=infected V1, m2=infected V2, mc=co-infected *)

prF[ex_] := ex /. {
  m0 -> Subscript[m, 0], m1 -> Subscript[m, 1],
  m2 -> Subscript[m, 2], mc -> Subscript[m, c],
  n0 -> Subscript[\[Nu], 0], n1 -> Subscript[\[Nu], 1],
  n2 -> Subscript[\[Nu], 2], nc -> Subscript[\[Nu], c],
  a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2],
  d1 -> Subscript[\[Delta], 1], d2 -> Subscript[\[Delta], 2],
  ec -> Subscript[\[Epsilon], c], ph -> \[CapitalPhi],
  la1 -> Subscript[\[Lambda], 1], la2 -> Subscript[\[Lambda], 2],
  lac -> Subscript[\[Lambda], c]
};

var = {m0, m1, m2, mc};

(* Force of infection terms *)
la1 = a1 m1 + d1 (1 - ec) mc;
la2 = a2 m2 + d2 (1 - ec) mc;
lac = (d1 + d2) ec mc;

RHS = {
  ph - n0 m0 - m0 (la1 + la2 + lac),
  -n1 m1 + m0 la1 - m1 (la2 + lac),
  -n2 m2 + m0 la2 - m2 (la1 + lac),
  -nc mc + m0 lac + m1 (la2 + lac) + m2 (la1 + lac)
};

(* Expanded RHS *)
RHS1 = {
  ph - n0 m0 - m0 (a1 m1 + d1 (1-ec) mc + a2 m2 + d2 (1-ec) mc + (d1+d2) ec mc),
  -n1 m1 + m0 (a1 m1 + d1 (1-ec) mc) - m1 (a2 m2 + d2 (1-ec) mc + (d1+d2) ec mc),
  -n2 m2 + m0 (a2 m2 + d2 (1-ec) mc) - m2 (a1 m1 + d1 (1-ec) mc + (d1+d2) ec mc),
  -nc mc + m0 (d1+d2) ec mc + m1 (a2 m2 + d2 (1-ec) mc + (d1+d2) ec mc) + m2 (a1 m1 + d1 (1-ec) mc + (d1+d2) ec mc)
};

(* ===== Model 2: Within-host (Eq. b) - adds phi factors ===== *)
(* phi_{2|1} = susceptibility to V2 given infected with V1 *)
(* phi_{1|2} = susceptibility to V1 given infected with V2 *)

prF2[ex_] := ex /. {
  m0 -> Subscript[m, 0], m1 -> Subscript[m, 1],
  m2 -> Subscript[m, 2], mc -> Subscript[m, c],
  n0 -> Subscript[\[Nu], 0], n1 -> Subscript[\[Nu], 1],
  n2 -> Subscript[\[Nu], 2], nc -> Subscript[\[Nu], c],
  a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2],
  d1 -> Subscript[\[Delta], 1], d2 -> Subscript[\[Delta], 2],
  ec -> Subscript[\[Epsilon], c], ph -> \[CapitalPhi],
  p21 -> Subscript[\[Phi], "2|1"], p12 -> Subscript[\[Phi], "1|2"]
};

var2 = {m0, m1, m2, mc};

RHS2 = {
  ph - n0 m0 - m0 (la1 + la2 + lac),
  -n1 m1 + m0 la1 - m1 p21 (la2 + lac),
  -n2 m2 + m0 la2 - m2 p12 (la1 + lac),
  -nc mc + m0 lac + m1 p21 (la2 + lac) + m2 p12 (la1 + lac)
};

(* ===== Model 3: Alizon's model (Eq. c) ===== *)
(* 6 compartments: m0, m1, m2, M1 (super-V1), M2 (super-V2), mc *)
(* New rates: Λ_1, Λ_2, λ_{1,c}, λ_{2,c} *)

prF3[ex_] := ex /. {
  m0 -> Subscript[m, 0], m1 -> Subscript[m, 1],
  m2 -> Subscript[m, 2], mc -> Subscript[m, c],
  M1 -> Subscript[M, 1], M2 -> Subscript[M, 2],
  n0 -> Subscript[\[Nu], 0], n1 -> Subscript[\[Nu], 1],
  n2 -> Subscript[\[Nu], 2], nc -> Subscript[\[Nu], c],
  nM1 -> Subscript[\[Nu], Subscript[M, 1]], nM2 -> Subscript[\[Nu], Subscript[M, 2]],
  a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2],
  d1 -> Subscript[\[Delta], 1], d2 -> Subscript[\[Delta], 2],
  ec -> Subscript[\[Epsilon], c], ph -> \[CapitalPhi],
  La1 -> Subscript[\[CapitalLambda], 1], La2 -> Subscript[\[CapitalLambda], 2],
  l1c -> Subscript[\[Lambda], "1,c"], l2c -> Subscript[\[Lambda], "2,c"]
};

var3 = {m0, m1, m2, M1, M2, mc};

(* Additional force of infection from super-infected and co-infected *)
(* La1 = Λ_1, La2 = Λ_2 from super-infected M1, M2 *)
(* l1c = λ_{1,c}, l2c = λ_{2,c} from co-infected mc *)

RHS3 = {
  ph - n0 m0 - m0 (la1 + la2 + lac + La1 + La2),
  -n1 m1 + m0 la1 - m1 (la1 + l1c + La1) - m1 (la2 + l2c + La2),
  -n2 m2 + m0 la2 - m2 (la2 + l2c + La2) - m2 (la1 + l1c + La1),
  -nM1 M1 + m0 La1 + m1 (la1 + l1c + La1),
  -nM2 M2 + m0 La2 + m2 (la2 + l2c + La2),
  -nc mc + m0 lac + m1 (la2 + l2c + La2) + m2 (la1 + l1c + La1)
};

(* ===== Model 4: Two-slot model (Eq. d) ===== *)
(* Same 6 compartments, different m0 → M1, M2, mc transitions *)

prF4[ex_] := ex /. {
  m0 -> Subscript[m, 0], m1 -> Subscript[m, 1],
  m2 -> Subscript[m, 2], mc -> Subscript[m, c],
  M1 -> Subscript[M, 1], M2 -> Subscript[M, 2],
  n0 -> Subscript[\[Nu], 0], n1 -> Subscript[\[Nu], 1],
  n2 -> Subscript[\[Nu], 2], nc -> Subscript[\[Nu], c],
  nM1 -> Subscript[\[Nu], Subscript[M, 1]], nM2 -> Subscript[\[Nu], Subscript[M, 2]],
  a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2],
  d1 -> Subscript[\[Delta], 1], d2 -> Subscript[\[Delta], 2],
  ec -> Subscript[\[Epsilon], c], ph -> \[CapitalPhi],
  La1 -> Subscript[\[CapitalLambda], 1], La2 -> Subscript[\[CapitalLambda], 2],
  Lac -> Subscript[\[CapitalLambda], c],
  l1c -> Subscript[\[Lambda], "1,c"], l2c -> Subscript[\[Lambda], "2,c"]
};

var4 = {m0, m1, m2, M1, M2, mc};

(* Two-slot: m0 → M1 with (d1/(d1+d2)) l1c + La1 *)
(* Two-slot: m0 → M2 with (d2/(d1+d2)) l2c + La2 *)
(* Two-slot: m0 → mc with Lac (not lac) *)

RHS4 = {
  ph - n0 m0 - m0 (la1 + la2 + d1/(d1+d2) l1c + La1 + d2/(d1+d2) l2c + La2 + Lac),
  -n1 m1 + m0 la1 - m1 (la1 + l1c + La1) - m1 (la2 + l2c + La2),
  -n2 m2 + m0 la2 - m2 (la2 + l2c + La2) - m2 (la1 + l1c + La1),
  -nM1 M1 + m0 (d1/(d1+d2) l1c + La1) + m1 (la1 + l1c + La1),
  -nM2 M2 + m0 (d2/(d1+d2) l2c + La2) + m2 (la2 + l2c + La2),
  -nc mc + m0 Lac + m1 (la2 + l2c + La2) + m2 (la1 + l1c + La1)
};
