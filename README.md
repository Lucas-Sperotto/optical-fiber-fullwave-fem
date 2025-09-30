# Optical Fiber Full-Wave FEM (Vector / Nédélec)

Solver vetorial 2D (elementos de aresta – Nédélec/Whitney) para modos de fibra óptica.
Resolve o problema generalizado:

\[
(K - k_0^2\,M_\varepsilon)\,e \;=\; \beta^2\,M_\mu\,e,\qquad n_\text{eff} = \beta/k_0
\]

- **Malha**: Gmsh v2 (triângulos).
- **Bases**: Whitney 1-forms (edge) – anti-espúrio para Maxwell.
- **Saída**: `modes_summary.csv` (β², β, n_eff) e `mode0_edges.csv` (coeficientes por aresta).
- **Roadmap**: integrar ARPACK/Spectra (shift-invert), PML anular, export VTK.

## Requisitos
- **C++17** (ou superior)
- **Eigen3** (`sudo apt install libeigen3-dev`)
- **Gmsh** para gerar a malha (`gmsh` CLI/GUI)

## Estrutura

```text
src/fiber_fem_nedelec.cpp   # montagem K, M_eps, M_mu + solver denso demo
meshes/fibra.geo            # exemplo de geometria (núcleo + casca)
results/                    # saídas .csv / .vtk
scripts/build.sh            # build helper
````

## Build rápido

```bash
# Ubuntu/Debian (Eigen3):
sudo apt update && sudo apt install -y build-essential libeigen3-dev gmsh

# Compilar:
g++ -O3 -march=native -DNDEBUG src/fiber_fem_nedelec.cpp -o fiber_fem -I/usr/include/eigen3
````

> Se preferir: `scripts/build.sh` já chama o comando acima.

## Gerar malha

Edite `meshes/fibra.geo` (raios/etiquetas) e rode:

```bash
gmsh -2 meshes/fibra.geo -o meshes/fibra.msh
```

## Executar

```bash
./fiber_fem \
  --msh meshes/fibra.msh \
  --lambda 1.55e-6 \
  --ncore 1.450 \
  --nclad 1.444 \
  --nev 6
```

### Saídas

* `results/modes_summary.csv` → `mode, beta2, beta, neff`
* `results/mode0_edges.csv` → `edge_id, node_a, node_b, coef`

> **Modos guiados** devem satisfazer: (k_0 n_\text{clad} < \beta < k_0 n_\text{core}).

## Validação

* Refinar malha e verificar convergência de (n_\text{eff}).
* Aumentar raio externo quando não houver PML.
* Comparar HE(*{11}) com aproximação LP(*{01}) para (V \lesssim 2.405).

## Roadmap

* [ ] Integrar **Spectra/ARPACK** (shift-invert) para malhas grandes
* [ ] **PML** anular (estiramento complexo)
* [ ] Export **VTK (.vtu)** com (E_x,E_y) reconstruídos
* [ ] Perfis GRIN e perdas (ε complexo)

## Licença

Defina conforme o seu objetivo (privado/fechado ou MIT/BSD). *Placeholder*.
