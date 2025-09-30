a = 4e-6;       // raio do núcleo
R = 20e-6;      // raio do domínio externo

// Pontos "fakes" para controlar tamanho de malha via lc ~ a/3
Point(1) = {0, 0, 0, a/3};

// Circunferências de núcleo e domínio
Circle(1) = {0, 0, 0, a};
Circle(2) = {0, 0, 0, R};

// Superfícies
Curve Loop(1) = {1};  Plane Surface(1) = {1};        // disco núcleo
Curve Loop(2) = {2};  Plane Surface(2) = {2};        // disco total
Surface(3) = Surface{2}; Surface{3} -= {1};          // casca = anel R\a

// Physical tags: 1 = núcleo, 2 = casca (3 reservado para PML, se quiser)
Physical Surface(1) = {1};
Physical Surface(2) = {3};

// Malha
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = a/6;
Mesh.CharacteristicLengthMax = a/4;
