# E-Digraph / S-Digraph — implementação em Python

Implementação em Python do procedimento geométrico descrito no artigo para os E-Digraph e S-Digraph.

## Arquivos

```text
edigraph-sdigraph/
├── README.md
└── edigraph_sdigraph.py
```

O projeto não possui dependências externas.

## Requisitos

Python 3.9 ou superior.

## Uso

Coloque `edigraph_sdigraph.py` no diretório desejado e importe as funções diretamente:

```python
from edigraph_sdigraph import derivar

resultado = derivar("energia", "U", "V")
print(resultado)
```

O resultado corresponde a:

```text
(∂U/∂V)_(S,N) = -P
```

## Execução direta

Para executar toda a implementação diretamente pelo terminal:

```bash
python edigraph_sdigraph.py
```

A execução é dividida em três partes:

1. **E-Digraph (representação da Energia):** percorre as 7 funções não nulas, reconstrói suas variáveis naturais, executa as 3 derivadas parciais de primeira ordem de cada função e analisa os 3 pares possíveis de derivadas sucessivas nas duas ordens.
2. **S-Digraph (representação da Entropia):** executa o mesmo procedimento para as 7 funções não nulas dessa representação.
3. **Validação específica do S-Digraph:** verifica passo a passo as 36 relações geométricas usadas como conjunto de referência, mostrando o resultado determinado geometricamente e o resultado esperado.

O Potencial Nulo permanece fora da rotina geral de derivação. Nas derivadas sucessivas, uma diferença entre os rótulos geométricos obtidos nas duas ordens é apresentada como informação da representação geométrica e não, por si só, como falha da igualdade analítica entre derivadas cruzadas.

## Validação das relações do S-Digraph

```python
from edigraph_sdigraph import validar_relacoes_sdigraph

resultado = validar_relacoes_sdigraph()
print(resultado)
```

Para o conjunto de referência utilizado no artigo, são verificadas 36 relações geométricas: 24 no primeiro bloco e 12 no segundo.

## Notação usada no código

Alguns símbolos do artigo são escritos com nomes textuais no Python:

```text
μ   -> mu
ψ'  -> psi_linha
ℱ   -> F_script
ℋ   -> H_script
```

As explicações específicas da implementação estão preservadas nas docstrings e comentários do próprio arquivo `edigraph_sdigraph.py`.
