# Python Toolbox calculo de distância entre dois pontos

---

Tool de geoprocessamento que realiza o cálculo do menor caminho, utilizando o algoritmo A*, entre dois pontos do mapa com base em uma malha de linhas.
<h3>Recebe três entradas</h3>

Estrutura de grafo gerada pela malha de linhas, sendo os vértices suas intersecções e arestas seus caminhos (criada por ferramenta auxiliar).

Camada de pontos de origem.

Camada de pontos de destino.

---

<h3>Saída</h3>

Armazena o resultado da distância percorrida de cada registro relacionado na camada de origem.

Gera uma camada com o caminho de cada registro.

---

<h3>Exemplo de execução</h3>


<h4>Os dois pontos de entrada:</h4>
<img src='https://github.com/NickolasCrema/imagens_readmes/blob/main/projeto_peixes/input_points.PNG?raw=true' alt='entrada'/>


<h4>Saída:</h4>
<img src='https://github.com/NickolasCrema/imagens_readmes/blob/main/projeto_peixes/output_path.png?raw=true', alt='saida'/>
<img src='https://github.com/NickolasCrema/imagens_readmes/blob/main/projeto_peixes/output_db.PNG?raw=true' width=700px, alt='saida-db'/>
