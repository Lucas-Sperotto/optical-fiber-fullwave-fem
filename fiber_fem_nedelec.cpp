/*
 * fiber_fem_nedelec.cpp
 * ---------------------------------------------------------------
 * Full‑wave vectorial 2D FEM (edge/Nédélec, 1ª ordem) para modos
 * de uma fibra óptica (núcleo + casca, com opção de PML anular).
 * Resolve o problema generalizado:
 *     (K - k0^2 M_eps) e = beta^2 M_mu e
 * onde e agrega os graus de liberdade tangenciais por aresta.
 *
 * STATUS: protótipo didático com montagem local TRI3 (1ª ordem).
 * - Monta K (curl‑curl) e M_eps (massa dielétrica) com bases de Whitney 1‑formas.
 * - Lê malha Gmsh v2 ASCII ($Nodes/$Elements) com triângulos (tipo=2).
 * - Constrói arestas globais, orientações locais e conectividade.
 * - Resolve ingenuamente com fallback denso (Eigen) para malhas pequenas.
 *   Para malhas reais, conecte Spectra/ARPACK (ver TODO no fim).
 *
 * COMPILAR (ex.: Ubuntu):
 *   g++ -O3 -march=native -DNDEBUG fiber_fem_nedelec.cpp -o fiber_fem \
 *       -I/usr/include/eigen3
 *
 * RODAR:
 *   ./fiber_fem --msh fibra.msh --lambda 1.55e-6 --ncore 1.450 --nclad 1.444 --pml 0
 *   (defina --pml 1 se quiser marcar região de PML por Physical Group)
 *
 * MALHA (Gmsh):
 *   Crie um .geo com dois discos concêntricos (raios a e R). Atribua
 *   Physical Surface 1 = núcleo, 2 = casca, (opcional) 3 = PML.
 *   Exemplo rápido (cole num arquivo fibra.geo e malhe no Gmsh GUI/CLI):
 *
 *   a = 4e-6; R = 20e-6; // raio núcleo e domínio externo
 *   Point(1) = {0,0,0, a/3};
 *   Circle(1) = {0,0,0, a};
 *   Circle(2) = {0,0,0, R};
 *   Curve Loop(1) = {1}; Plane Surface(1) = {1}; // disco núcleo
 *   Curve Loop(2) = {2}; Plane Surface(2) = {2}; // disco externo
 *   Surface(3) = Surface{2}; Surface{3} -= {1}; // casca = anel R\a
 *   Physical Surface(1) = {1}; // núcleo
 *   Physical Surface(2) = {3}; // casca
 *   Mesh.Algorithm = 6; Mesh.CharacteristicLengthMin = a/6; Mesh.CharacteristicLengthMax = a/4;
 *
 *   Salve e gere: gmsh -2 fibra.geo -o fibra.msh
 *
 * VALIDAÇÃO: Modos guiados devem ter k0*nclad < beta < k0*ncore.
 * ---------------------------------------------------------------
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

// ------------------------------- Util -------------------------------
struct Cmd {
    std::string msh = "fibra.msh";
    double lambda = 1.55e-6;
    double ncore = 1.450;
    double nclad = 1.444;
    int pml_on = 0; // 1 se Physical Surface 3 for PML
    int nev = 6;    // quantos modos (apenas para o solver denso demonstrativo)
};

Cmd parse_cmd(int argc, char** argv){
    Cmd c; for(int i=1;i<argc;i++){
        std::string t = argv[i];
        auto need = [&](const char* flag){ if(i+1>=argc){fprintf(stderr,"Falta arg de %s\n",flag); exit(1);} return std::string(argv[++i]); };
        if(t=="--msh") c.msh = need("--msh");
        else if(t=="--lambda") c.lambda = std::stod(need("--lambda"));
        else if(t=="--ncore") c.ncore = std::stod(need("--ncore"));
        else if(t=="--nclad") c.nclad = std::stod(need("--nclad"));
        else if(t=="--pml") c.pml_on = std::stoi(need("--pml"));
        else if(t=="--nev") c.nev = std::stoi(need("--nev"));
        else { fprintf(stderr,"Arg desconhecido: %s\n", t.c_str()); exit(1);} }
    return c;
}

// ------------------------------- Malha -------------------------------
struct Node { double x,y; };
struct Tri { int tag; int v[3]; }; // tag = Physical; v[] = índices de nós (1-based do .msh será convertido p/ 0-based)

struct Mesh {
    std::vector<Node> nodes;      // tamanho NN
    std::vector<Tri>  tris;       // tamanho NT, só triângulos
    // Arestas globais (min,max) -> id
    std::vector<std::array<int,2>> edge2verts; // id -> {a,b} com a<b
    std::vector<std::vector<int>>  tri2edges;  // NT x 3, mapeia para arestas globais
    std::vector<std::array<int,3>> tri_edge_sign; // +1 se orientação local coincide, -1 caso contrário
};

// Leitor mínimo para Gmsh v2 ASCII ($Nodes/$Elements) com triângulos tipo=2
Mesh read_msh_v2_ascii(const std::string& path){
    std::ifstream in(path); if(!in){ throw std::runtime_error("Não abriu .msh"); }
    std::string line; Mesh M; int numNodes=0, numElem=0;
    auto expect = [&](const char* token){ if(line!=token) throw std::runtime_error(std::string("Esperado ")+token+", veio "+line); };
    while(std::getline(in,line)){
        if(line=="$Nodes"){
            std::getline(in,line); numNodes = std::stoi(line); M.nodes.resize(numNodes);
            for(int i=0;i<numNodes;i++){
                std::getline(in,line); std::istringstream ss(line);
                int id; double x,y,z; ss>>id>>x>>y>>z; M.nodes[id-1]={x,y};
            }
            std::getline(in,line); expect("$EndNodes");
        } else if(line=="$Elements"){
            std::getline(in,line); numElem = std::stoi(line);
            for(int i=0;i<numElem;i++){
                std::getline(in,line); std::istringstream ss(line);
                int id, type, ntags; ss>>id>>type>>ntags;
                std::vector<int> tags(ntags); for(int k=0;k<ntags;k++) ss>>tags[k];
                if(type==2){ // triangle
                    int n1,n2,n3; ss>>n1>>n2>>n3; Tri t; t.tag = (ntags>=1?tags[0]:0); t.v[0]=n1-1; t.v[1]=n2-1; t.v[2]=n3-1; M.tris.push_back(t);
                } else {
                    // ignorar outros elementos
                }
            }
            std::getline(in,line); expect("$EndElements");
        }
    }
    // Construir arestas globais e conectividade
    std::unordered_map<uint64_t,int> mapEdge; // chave compacta (min,max) -> id
    auto key = [&](int a,int b){ if(a>b) std::swap(a,b); return ( (uint64_t)a<<32 ) | (uint64_t)b; };
    for(auto &t : M.tris){
        int vs[3] = {t.v[0],t.v[1],t.v[2]};
        int einds[3]; int esgn[3];
        int loc[3][2] = {{1,2},{2,0},{0,1}}; // arestas locais: (v1,v2),(v2,v0),(v0,v1)
        for(int e=0;e<3;e++){
            int a=vs[loc[e][0]], b=vs[loc[e][1]];
            uint64_t k = key(a,b);
            auto it = mapEdge.find(k);
            int eid;
            if(it==mapEdge.end()){
                eid = (int)M.edge2verts.size();
                mapEdge[k]=eid; M.edge2verts.push_back({{std::min(a,b), std::max(a,b)}});
            } else eid=it->second;
            einds[e]=eid;
            esgn[e] = (a<b)? +1 : -1; // sinal local: +1 se orientação local (a->b) coincide com a<b
        }
        M.tri2edges.push_back({einds[0],einds[1],einds[2]});
        M.tri_edge_sign.push_back({esgn[0],esgn[1],esgn[2]});
    }
    return M;
}

// ---------------------- Bases de Whitney (triângulo) ----------------------
// Para triângulo com vértices (r1,r2,r3), funções barycêntricas λi têm gradientes
// constantes no elemento: ∇λi = ni / (2A), com ni = (y_j - y_k, x_k - x_j), (i,j,k) cíclico.
// Bases de aresta (Whitney 1‑forms) associadas às arestas (i,j): w_ij = λ_i ∇λ_j - λ_j ∇λ_i.
// Arestas locais em TRI: e0=(v1,v2) => w_12, e1=(v2,v0) => w_23, e2=(v0,v1) => w_31.

struct TriGeom {
    double A;            // área
    Vector2d gradL[3];   // ∇λ1, ∇λ2, ∇λ3
};

TriGeom tri_geom(const Node& n0, const Node& n1, const Node& n2){
    TriGeom g{};
    double x1=n0.x, y1=n0.y;
    double x2=n1.x, y2=n1.y;
    double x3=n2.x, y3=n2.y;
    double det = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
    g.A = 0.5*det;
    // normais "rotacionadas" para gradientes barycêntricos
    Vector2d n1v(y2 - y3, x3 - x2);
    Vector2d n2v(y3 - y1, x1 - x3);
    Vector2d n3v(y1 - y2, x2 - x1);
    g.gradL[0] = n1v / (2.0*g.A);
    g.gradL[1] = n2v / (2.0*g.A);
    g.gradL[2] = n3v / (2.0*g.A);
    return g;
}

// Integrais de momentos barycêntricos em TRI linear:
// ∫ λ_i^2 dΩ = A/6 ; ∫ λ_i λ_j dΩ = A/12 (i≠j).
inline double I_li_lj(double A, int i, int j){ return (i==j) ? (A/6.0) : (A/12.0); }

// Curl 2D (escalar fora do plano) de w_ij é constante: ∇×w_ij = 2 (∇λ_i × ∇λ_j)
inline double cross2(const Vector2d& a, const Vector2d& b){ return a.x()*b.y() - a.y()*b.x(); }

// ------------------ Propriedades materiais por elemento ------------------
struct MatProps { double eps_r; double mu_r; int is_pml; };

// Exemplo simples: define eps_r por Physical tag
MatProps element_material(const Cmd& cmd, int physTag){
    MatProps m{}; m.mu_r = 1.0; m.is_pml = 0;
    if(physTag==1)      m.eps_r = cmd.ncore*cmd.ncore; // núcleo
    else if(physTag==2) m.eps_r = cmd.nclad*cmd.nclad; // casca
    else if(physTag==3) { m.eps_r = cmd.nclad*cmd.nclad; m.is_pml = cmd.pml_on; }
    else m.eps_r = cmd.nclad*cmd.nclad; // default
    return m;
}

// -------------------------- Montagem de matrizes --------------------------
struct Assembled {
    SparseMatrix<double> K;     // curl‑curl
    SparseMatrix<double> M_eps; // massa dielétrica
    SparseMatrix<double> M_mu;  // massa magnética (aqui = identidade ponderada, mu_r)
};

Assembled assemble_system(const Cmd& cmd, const Mesh& M){
    const int NE = (int)M.edge2verts.size();
    std::vector<Triplet<double>> tk, tmep, tmu;
    tk.reserve(M.tris.size()*9);
    tmep.reserve(M.tris.size()*9);
    tmu.reserve(M.tris.size()*9);

    for(size_t it=0; it<M.tris.size(); ++it){
        const Tri& T = M.tris[it];
        const Node& n0 = M.nodes[T.v[0]];
        const Node& n1 = M.nodes[T.v[1]];
        const Node& n2 = M.nodes[T.v[2]];
        TriGeom g = tri_geom(n0,n1,n2);
        if(g.A<=0){ fprintf(stderr,"Tri degenerado em %zu\n", it); continue; }

        // Materiais
        MatProps mp = element_material(cmd, T.tag);
        double epsr = mp.eps_r; double mur = mp.mu_r; double inv_mur = 1.0/mur;
        // (PML real não implementado – aqui apenas etiqueta)

        // Mapas locais
        int EID[3] = { M.tri2edges[it][0], M.tri2edges[it][1], M.tri2edges[it][2] };
        int SGN[3] = { M.tri_edge_sign[it][0], M.tri_edge_sign[it][1], M.tri_edge_sign[it][2] };

        // Gradientes barycêntricos
        const Vector2d& dL1 = g.gradL[0];
        const Vector2d& dL2 = g.gradL[1];
        const Vector2d& dL3 = g.gradL[2];

        // Arestas locais (seguindo convenção acima):
        // e0 ~ (v1,v2) => w_12 = λ1 dλ2 − λ2 dλ1
        // e1 ~ (v2,v0) => w_23 = λ2 dλ3 − λ3 dλ2
        // e2 ~ (v0,v1) => w_31 = λ3 dλ1 − λ1 dλ3
        // Guardar pares (i,j) para cada aresta local – indices em {0,1,2} (λ1,λ2,λ3)
        int IJ[3][2] = {{0,1},{1,2},{2,0}}; // (1,2),(2,3),(3,1) em 0‑based
        const Vector2d dL[3] = {dL1,dL2,dL3};

        // Matrizes locais 3x3
        double Ke[3][3]={0}, Meps[3][3]={0}, Mmu[3][3]={0};

        // --- K local: ∫ mu^{-1} (curl w_a)(curl w_b) dΩ ---
        // curl w_ij = 2 (dLi × dLj)
        double curlLoc[3];
        for(int a=0;a<3;a++){
            int i = IJ[a][0], j = IJ[a][1];
            curlLoc[a] = 2.0 * cross2(dL[i], dL[j]);
        }
        for(int a=0;a<3;a++) for(int b=0;b<3;b++){
            Ke[a][b] = inv_mur * (curlLoc[a]*curlLoc[b]) * g.A; // constantes no elemento
        }

        // --- M_eps local: ∫ eps_r w_a · w_b dΩ ---
        // w_ij = λ_i dλ_j − λ_j dλ_i
        auto dot = [&](const Vector2d& u, const Vector2d& v){ return u.dot(v); };
        for(int a=0;a<3;a++){
            int ia = IJ[a][0], ja = IJ[a][1];
            for(int b=0;b<3;b++){
                int ib = IJ[b][0], jb = IJ[b][1];
                // w_a · w_b = (λ_ia λ_ib)(dL_ja·dL_jb) − (λ_ia λ_jb)(dL_ja·dL_ib)
                //            − (λ_ja λ_ib)(dL_ia·dL_jb) + (λ_ja λ_jb)(dL_ia·dL_ib)
                double t1 = I_li_lj(g.A, ia, ib) * dot(dL[ja], dL[jb]);
                double t2 = I_li_lj(g.A, ia, jb) * dot(dL[ja], dL[ib]);
                double t3 = I_li_lj(g.A, ja, ib) * dot(dL[ia], dL[jb]);
                double t4 = I_li_lj(g.A, ja, jb) * dot(dL[ia], dL[ib]);
                Meps[a][b] = epsr * (t1 - t2 - t3 + t4);
            }
        }

        // --- M_mu local: ∫ mu^{-1} w_a · w_b dΩ --- (útil como T no autoproblema)
        for(int a=0;a<3;a++){
            int ia = IJ[a][0], ja = IJ[a][1];
            for(int b=0;b<3;b++){
                int ib = IJ[b][0], jb = IJ[b][1];
                double t1 = I_li_lj(g.A, ia, ib) * dot(dL[ja], dL[jb]);
                double t2 = I_li_lj(g.A, ia, jb) * dot(dL[ja], dL[ib]);
                double t3 = I_li_lj(g.A, ja, ib) * dot(dL[ia], dL[jb]);
                double t4 = I_li_lj(g.A, ja, jb) * dot(dL[ia], dL[ib]);
                Mmu[a][b] = (inv_mur) * (t1 - t2 - t3 + t4);
            }
        }

        // Acumular com sinais de orientação (SGN)
        for(int a=0;a<3;a++){
            int Aglob = EID[a]; double sa = (double)SGN[a];
            for(int b=0;b<3;b++){
                int Bglob = EID[b]; double sb = (double)SGN[b];
                tk.emplace_back(Aglob,Bglob, sa*sb*Ke[a][b]);
                tmep.emplace_back(Aglob,Bglob, sa*sb*Meps[a][b]);
                tmu.emplace_back(Aglob,Bglob, sa*sb*Mmu[a][b]);
            }
        }
    }

    Assembled A;
    const int NE = (int)M.edge2verts.size();
    A.K.resize(NE,NE); A.M_eps.resize(NE,NE); A.M_mu.resize(NE,NE);
    A.K.setFromTriplets(tk.begin(), tk.end());
    A.M_eps.setFromTriplets(tmep.begin(), tmep.end());
    A.M_mu.setFromTriplets(tmu.begin(), tmu.end());
    return A;
}

// ------------------------ Solver (denso de demonstração) ------------------------
// Para malhas pequenas: converte as matrizes esparsas em densas e usa Eigen denso
// no problema generalizado simétrico (auto-adjunto): S e = beta^2 T e
// com S = K - k0^2 M_eps, T = M_mu.

struct Modes { std::vector<double> beta2; MatrixXd V; };

Modes solve_dense_demo(const SparseMatrix<double>& K,
                       const SparseMatrix<double>& M_eps,
                       const SparseMatrix<double>& M_mu,
                       double k0, int nev)
{
    MatrixXd Sd = MatrixXd(K) - (k0*k0)*MatrixXd(M_eps);
    MatrixXd Td = MatrixXd(M_mu);
    // Generalized self-adjoint eigendecomposition: Sd v = lambda Td v
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> ges(Sd, Td);
    if(ges.info()!=Eigen::Success) throw std::runtime_error("Falha no eigensolver denso");
    VectorXd lam = ges.eigenvalues(); // crescente
    MatrixXd V   = ges.eigenvectors();

    // Selecionar os nev maiores lambdas fisicamente relevantes (beta^2)
    // OBS: queremos beta^2 entre (k0^2 n_clad^2, k0^2 n_core^2). Aqui apenas
    // ordenamos decrescente para inspecionar; usuário filtra por faixa depois.
    std::vector<int> idx(lam.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a,int b){ return lam[a] > lam[b]; });

    Modes m; m.V.resize(V.rows(), std::min<int>(nev, V.cols()));
    for(int k=0;k<std::min<int>(nev, (int)idx.size()); ++k){
        m.beta2.push_back(lam[idx[k]]);
        m.V.col(k) = V.col(idx[k]]);
    }
    return m;
}

// ------------------------------- Main -------------------------------
int main(int argc, char** argv){
    try{
        Cmd cmd = parse_cmd(argc, argv);
        double k0 = 2.0*M_PI / cmd.lambda;
        std::cerr << "Lendo malha: "<<cmd.msh<<"\n";
        Mesh mesh = read_msh_v2_ascii(cmd.msh);
        std::cerr << "Nós: "<<mesh.nodes.size()<<", Triângulos: "<<mesh.tris.size()
                  <<", Arestas: "<<mesh.edge2verts.size()<<"\n";
        auto A = assemble_system(cmd, mesh);
        std::cerr << "Montagem ok. K nnz="<<A.K.nonZeros()<<", M_eps nnz="<<A.M_eps.nonZeros()
                  <<", M_mu nnz="<<A.M_mu.nonZeros()<<"\n";

        // Resolver (demo denso)
        Modes modes = solve_dense_demo(A.K, A.M_eps, A.M_mu, k0, cmd.nev);

        // Reportar beta, neff e salvar primeiro modo em CSV
        std::ofstream out("modes_summary.csv");
        out << "mode,beta2,beta,neff\n";
        for(size_t i=0;i<modes.beta2.size();++i){
            double b2 = modes.beta2[i];
            if(b2 <= 0) continue; // ignorar não físicos (ou contínuo)
            double beta = std::sqrt(b2);
            double neff = beta / k0;
            out << (i+1) << "," << std::setprecision(16) << b2 << "," << beta << "," << neff << "\n";
        }
        out.close();
        std::cerr << "Salvo modes_summary.csv\n";

        // Exportar autovetor do 1º modo (coluna 0) em CSV por aresta
        if(modes.V.cols()>0){
            std::ofstream f("mode0_edges.csv");
            f << "edge_id,node_a,node_b,coef\n";
            for(int e=0;e<(int)mesh.edge2verts.size();++e){
                int a = mesh.edge2verts[e][0];
                int b = mesh.edge2verts[e][1];
                f << e << "," << a << "," << b << "," << std::setprecision(16) << modes.V(e,0) << "\n";
            }
            f.close();
            std::cerr << "Salvo mode0_edges.csv\n";
        }

        std::cerr << "Pronto.\n";
        return 0;
    } catch(const std::exception& e){
        std::cerr << "Erro: "<<e.what()<<"\n";
        return 1;
    }
}

/* -------------------------------------------------------------------------
 * NOTAS TÉCNICAS
 * 1) Orientação de arestas: O sinal SGN assegura consistência tangencial entre
 *    orientação local da base e a orientação global (min->max). Isso evita
 *    cancelamentos errados no assembly global.
 *
 * 2) Faixa física: Após obter beta^2, filtre por k0^2*n_clad^2 < beta^2 < k0^2*n_core^2.
 *    Use isso para identificar modos guiados. O summary CSV ajuda a inspecionar.
 *
 * 3) PML: Aqui não implementamos o estiramento complexo. Uma rota simples é
 *    multiplicar os operadores por métricas de PML (s_x,s_y) no elemento. Para
 *    uma primeira versão, amplie o domínio (R grande) e aceite erro de truncamento.
 *
 * 4) Estabilização anti‑espúrios: Bases de Nédélec já eliminam boa parte de modos
 *    espúrios. Se necessário, use formul. mistas (adicionando escalar nodal) ou
 *    um pequeno termo de penalidade de grad‑div em um espaço compatível.
 *
 * 5) Solver grande escala: Para malhas reais, use ARPACK/Spectra com shift‑invert
 *    no problema (K - k0^2 M_eps) e = beta^2 M_mu e. Apontar o alvo próximo de
 *    k0^2*n_core^2 (ou usar recíproca) acelera a convergência dos primeiros modos.
 *    Ex.: Spectra::SymGEigsShiftSolver com OP = (S - sigma T)^{-1} T.
 *
 * 6) Campo (Ex,Ey) no pós‑processamento: O autovetor retorna coeficientes por
 *    aresta. Para recuperar E_t em pontos (ex.: centróide), avalie combinações
 *    das bases w_ij (respeitando o sinal local de cada tri) e acumule por elemento.
 *    Para visualização VTK, amostre em nós ou quadratura e exporte como cell data.
 *
 * 7) Checagem com analítico (step‑index circular): Compare o modo fundamental
 *    HE11 com aproximação LP01 quando V <= 2.405. Isso dá sanity check dos valores
 *    de neff (erro típico < 1e‑3 com malha moderada).
 *
 * 8) Extensões: dispersão (n(lambda)), perdas (eps complexo), curvatura, fibras GDF,
 *    perfis GRIN (eps_r(r)). Basta trocar element_material() e, para PML, métrica.
 * ------------------------------------------------------------------------- */
