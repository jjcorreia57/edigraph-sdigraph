# -*- coding: utf-8 -*-

"""Implementação computacional dos E-Digraph e S-Digraph.


O módulo reproduz o procedimento geométrico usado no artigo para derivadas
parciais de primeira ordem e para a preparação/análise de derivadas sucessivas.
Também contém a verificação das 36 relações geométricas do S-Digraph usadas
como conjunto de referência.

Convenções de nomes usadas no código:
    mu         -> μ
    psi_linha  -> ψ'
    F_script   -> ℱ
    H_script   -> ℋ

O Potencial Nulo é tratado separadamente, pois as três grandezas intensivas
associadas a esse vértice não constituem três variáveis naturais independentes.
"""
vertices = {'energia': {'U': (0, 0, 0), 'F': (1, 0, 0), 'H': (0, 1, 0), 'I': (0, 0, 1), 'G': (1, 1, 0), 'Omega': (1, 0, 1), 'Gamma': (0, 1, 1), 'zero': (1, 1, 1)}, 'entropia': {'S': (0, 0, 0), 'psi': (1, 0, 0), 'varphi': (0, 1, 0), 'phi': (0, 0, 1), 'psi_linha': (1, 1, 0), 'F_script': (1, 0, 1), 'H_script': (0, 1, 1), 'zero': (1, 1, 1)}}
eixos = {'energia': [{0: ('S', '+', 'T'), 1: ('T', '-', 'S')}, {0: ('V', '-', 'P'), 1: ('P', '+', 'V')}, {0: ('N', '+', 'mu'), 1: ('mu', '-', 'N')}], 'entropia': [{0: ('U', '+', '1/T'), 1: ('1/T', '-', 'U')}, {0: ('V', '+', 'P/T'), 1: ('P/T', '-', 'V')}, {0: ('N', '-', 'mu/T'), 1: ('mu/T', '+', 'N')}]}

def montar_regras(representacao, funcao):
    """Reconstrói as variáveis naturais e as regras de derivação de uma função.

A coordenada binária do vértice seleciona, em cada uma das três direções,
a extremidade do eixo que atua como variável natural. Para cada variável,
é mantido o par (sinal, variável conjugada)."""
    if representacao not in vertices:
        raise ValueError('Representacao desconhecida.')
    if funcao not in vertices[representacao]:
        raise ValueError('Funcao desconhecida.')
    coordenada = vertices[representacao][funcao]
    naturais = []
    regras = {}
    for numero_do_eixo, bit in enumerate(coordenada):
        dados = eixos[representacao][numero_do_eixo][bit]
        natural, sinal, conjugada = dados
        naturais.append(natural)
        regras[natural] = (sinal, conjugada)
    return (naturais, regras)

def derivar(representacao, funcao, variavel):
    """Calcula uma derivada parcial de primeira ordem.

Recebe a representação, a função termodinâmica e a variável de derivação.
Reconstrói as variáveis naturais e identifica a grandeza conjugada e seu
sinal. A saída mantém separadamente o sinal algébrico e a variável-resultado,
além da forma já formatada para exibição; isso permite reutilizar o resultado
como variável intermediária em uma derivação sucessiva."""
    if funcao == 'zero':
        raise ValueError('O Potencial Nulo requer tratamento separado.')
    naturais, regras = montar_regras(representacao, funcao)
    if variavel not in regras:
        raise ValueError('A variavel escolhida nao e natural da funcao.')
    sinal, conjugada = regras[variavel]
    constantes = [x for x in naturais if x != variavel]
    sinal_numerico = 1 if sinal == '+' else -1
    if sinal_numerico == 1:
        resultado = conjugada
    else:
        resultado = '-' + conjugada
    return {'funcao': funcao, 'variavel': variavel, 'constantes': constantes, 'sinal_resultado': sinal_numerico, 'variavel_resultado': conjugada, 'resultado': resultado}

def preparar_derivacao_sucessiva(representacao, funcao, primeira_variavel, segunda_variavel):
    """Organiza os dados necessários para uma derivação sucessiva.

Esta função ainda não determina a variável-resultado da segunda derivação.
Ela registra a ordem das derivações, a variável intermediária, seu sinal
algébrico e a única variável natural mantida constante, preparando os dados
para a etapa geométrica baseada na orientação tridimensional e na Regra da
Mão Direita."""
    naturais, _ = montar_regras(representacao, funcao)
    if primeira_variavel not in naturais:
        raise ValueError('A primeira variavel escolhida nao e natural da funcao.')
    if segunda_variavel not in naturais:
        raise ValueError('A segunda variavel escolhida nao e natural da funcao.')
    if primeira_variavel == segunda_variavel:
        raise ValueError('As duas variaveis de derivacao devem ser distintas.')
    primeira_derivacao = derivar(representacao, funcao, primeira_variavel)
    constantes = [x for x in naturais if x != primeira_variavel and x != segunda_variavel]
    if len(constantes) != 1:
        raise RuntimeError('Nao foi possivel identificar a variavel mantida constante.')
    return {'representacao': representacao, 'funcao': funcao, 'primeira_variavel': primeira_variavel, 'segunda_variavel': segunda_variavel, 'variavel_constante': constantes[0], 'sinal_intermediaria': primeira_derivacao['sinal_resultado'], 'variavel_intermediaria': primeira_derivacao['variavel_resultado'], 'intermediaria_formatada': primeira_derivacao['resultado']}
base_cartesiana = {0: (1, 0, 0), 1: (0, 1, 0), 2: (0, 0, 1)}

def produto_vetorial(vetor_a, vetor_b):
    """Calcula o produto vetorial entre dois vetores tridimensionais."""
    ax, ay, az = vetor_a
    bx, by, bz = vetor_b
    return (ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)

def orientar_par_de_eixos(eixo_primeiro, eixo_segundo):
    """Determina o eixo remanescente e a orientação associada à ordem de dois eixos."""
    eixos_validos = (0, 1, 2)
    if eixo_primeiro not in eixos_validos:
        raise ValueError('O primeiro eixo deve ser 0, 1 ou 2.')
    if eixo_segundo not in eixos_validos:
        raise ValueError('O segundo eixo deve ser 0, 1 ou 2.')
    if eixo_primeiro == eixo_segundo:
        raise ValueError('O produto vetorial exige dois eixos distintos.')
    vetor_primeiro = base_cartesiana[eixo_primeiro]
    vetor_segundo = base_cartesiana[eixo_segundo]
    vetor_resultante = produto_vetorial(vetor_primeiro, vetor_segundo)
    componentes_nao_nulas = [indice for indice, componente in enumerate(vetor_resultante) if componente != 0]
    if len(componentes_nao_nulas) != 1:
        raise RuntimeError('O produto vetorial nao produziu um unico eixo remanescente.')
    eixo_remanescente = componentes_nao_nulas[0]
    sinal_geometrico = vetor_resultante[eixo_remanescente]
    return {'eixo_primeiro': eixo_primeiro, 'eixo_segundo': eixo_segundo, 'eixo_remanescente': eixo_remanescente, 'vetor_resultante': vetor_resultante, 'sinal_geometrico': sinal_geometrico}

def localizar_eixo(representacao, variavel):
    """Localiza em qual das três direções coordenadas uma variável está representada."""
    if representacao not in eixos:
        raise ValueError('Representacao desconhecida.')
    for numero_do_eixo, extremidades in enumerate(eixos[representacao]):
        variaveis_do_eixo = {dados[0] for dados in extremidades.values()}
        if variavel in variaveis_do_eixo:
            return numero_do_eixo
    raise ValueError('A variavel nao pertence a nenhum eixo coordenado.')

def identificar_eixo_remanescente(preparacao):
    """Identifica o eixo não usado pelo par de variáveis de uma derivação sucessiva."""
    representacao = preparacao['representacao']
    eixo_intermediaria = localizar_eixo(representacao, preparacao['variavel_intermediaria'])
    eixo_segunda_variavel = localizar_eixo(representacao, preparacao['segunda_variavel'])
    eixos_utilizados = {eixo_intermediaria, eixo_segunda_variavel}
    eixos_remanescentes = [eixo for eixo in (0, 1, 2) if eixo not in eixos_utilizados]
    if len(eixos_remanescentes) != 1:
        raise RuntimeError('Nao foi possivel identificar um unico eixo remanescente.')
    return {'eixo_intermediaria': eixo_intermediaria, 'eixo_segunda_variavel': eixo_segunda_variavel, 'eixo_remanescente': eixos_remanescentes[0]}

def diagnosticar_orientacao_dos_eixos(preparacao):
    """Reúne informações de orientação dos eixos envolvidos em uma derivação sucessiva."""
    representacao = preparacao['representacao']
    dados_dos_eixos = identificar_eixo_remanescente(preparacao)
    eixo_intermediaria = dados_dos_eixos['eixo_intermediaria']
    eixo_segunda_variavel = dados_dos_eixos['eixo_segunda_variavel']
    eixo_por_exclusao = dados_dos_eixos['eixo_remanescente']
    orientacao_cartesiana = orientar_par_de_eixos(eixo_intermediaria, eixo_segunda_variavel)
    eixo_por_produto_cartesiano = orientacao_cartesiana['eixo_remanescente']
    if eixo_por_exclusao != eixo_por_produto_cartesiano:
        raise RuntimeError('A exclusao e o produto vetorial produziram eixos diferentes.')
    seta_intermediaria = obter_vetor_seta_coordenada(representacao, eixo_intermediaria)
    seta_segunda_variavel = obter_vetor_seta_coordenada(representacao, eixo_segunda_variavel)
    seta_eixo_remanescente = obter_vetor_seta_coordenada(representacao, eixo_por_exclusao)
    posicao_intermediaria = obter_vetor_posicional_da_variavel(representacao, preparacao['variavel_intermediaria'])
    posicao_segunda_variavel = obter_vetor_posicional_da_variavel(representacao, preparacao['segunda_variavel'])
    if posicao_intermediaria['eixo'] != eixo_intermediaria:
        raise RuntimeError('O vetor posicional da variavel intermediaria pertence ao eixo errado.')
    if posicao_segunda_variavel['eixo'] != eixo_segunda_variavel:
        raise RuntimeError('O vetor posicional da segunda variavel pertence ao eixo errado.')
    vetor_resultante_posicional = produto_vetorial(posicao_intermediaria['vetor_posicional'], posicao_segunda_variavel['vetor_posicional'])
    componentes_nao_nulas = [indice for indice, componente in enumerate(vetor_resultante_posicional) if componente != 0]
    if len(componentes_nao_nulas) != 1:
        raise RuntimeError('O produto dos vetores posicionais nao produziu um unico eixo.')
    eixo_resultante_posicional = componentes_nao_nulas[0]
    if eixo_resultante_posicional != eixo_por_exclusao:
        raise RuntimeError('O produto dos vetores posicionais nao coincide com o eixo remanescente.')
    vetor_seta_remanescente = seta_eixo_remanescente['vetor_orientado']
    vetor_oposto = tuple((-componente for componente in vetor_seta_remanescente))
    if vetor_resultante_posicional == vetor_seta_remanescente:
        sinal_geometrico_digraph = 1
    elif vetor_resultante_posicional == vetor_oposto:
        sinal_geometrico_digraph = -1
    else:
        raise RuntimeError('Nao foi possivel comparar as orientacoes dos vetores.')
    return {'eixo_intermediaria': eixo_intermediaria, 'eixo_segunda_variavel': eixo_segunda_variavel, 'eixo_remanescente': eixo_por_exclusao, 'vetor_resultante': orientacao_cartesiana['vetor_resultante'], 'sinal_geometrico_cartesiano': orientacao_cartesiana['sinal_geometrico'], 'vetor_seta_intermediaria': seta_intermediaria['vetor_orientado'], 'vetor_seta_segunda_variavel': seta_segunda_variavel['vetor_orientado'], 'vetor_posicional_intermediaria': posicao_intermediaria['vetor_posicional'], 'vetor_posicional_segunda_variavel': posicao_segunda_variavel['vetor_posicional'], 'vetor_resultante_orientado': vetor_resultante_posicional, 'vetor_seta_eixo_remanescente': vetor_seta_remanescente, 'sinal_geometrico_digraph': sinal_geometrico_digraph}

def identificar_lado_indicado_pela_regra(sinal_geometrico_digraph):
    """Identifica a extremidade inicialmente indicada pela orientação geométrica."""
    if sinal_geometrico_digraph not in (-1, 1):
        raise ValueError('O sinal geometrico do Digraph deve ser -1 ou 1.')
    if sinal_geometrico_digraph == 1:
        lado_indicado = 'ponto'
    else:
        lado_indicado = 'cruz'
    return {'sinal_geometrico_digraph': sinal_geometrico_digraph, 'lado_indicado_pela_regra': lado_indicado}

def determinar_variavel_indicada_geometricamente(preparacao):
    """Determina a variável indicada na extremidade final do eixo remanescente.

A Regra da Mão Direita fornece a extremidade inicialmente indicada. Depois,
o sinal da variável intermediária atua sobre essa escolha: sinal positivo
conserva a extremidade; sinal negativo seleciona a extremidade oposta."""
    orientacao_dos_eixos = diagnosticar_orientacao_dos_eixos(preparacao)
    identificacao_do_lado = identificar_lado_indicado_pela_regra(orientacao_dos_eixos['sinal_geometrico_digraph'])
    extremidades = identificar_extremidades_eixo_remanescente(preparacao)
    lado_inicial = identificacao_do_lado['lado_indicado_pela_regra']
    sinal_intermediaria = preparacao['sinal_intermediaria']
    if sinal_intermediaria == 1:
        lado_resultante = lado_inicial
        houve_inversao = False
    elif sinal_intermediaria == -1:
        if lado_inicial == 'cruz':
            lado_resultante = 'ponto'
        elif lado_inicial == 'ponto':
            lado_resultante = 'cruz'
        else:
            raise RuntimeError("O lado inicial deve ser 'cruz' ou 'ponto'.")
        houve_inversao = True
    else:
        raise ValueError('O sinal da variavel intermediaria deve ser -1 ou 1.')
    if lado_resultante == 'cruz':
        variavel_indicada = extremidades['extremidade_cruz']
    elif lado_resultante == 'ponto':
        variavel_indicada = extremidades['extremidade_ponto']
    else:
        raise RuntimeError("O lado resultante deve ser 'cruz' ou 'ponto'.")
    return {'eixo_remanescente': extremidades['eixo_remanescente'], 'lado_indicado_pela_regra': lado_inicial, 'sinal_intermediaria': sinal_intermediaria, 'houve_inversao': houve_inversao, 'lado_resultante': lado_resultante, 'variavel_indicada_geometricamente': variavel_indicada, 'sinal_geometrico_digraph': orientacao_dos_eixos['sinal_geometrico_digraph'], 'resultado_formatado': variavel_indicada}

def localizar_extremidade(representacao, variavel):
    """Localiza se uma variável ocupa o lado cruz ou o lado ponto de um eixo."""
    if representacao not in eixos:
        raise ValueError('Representacao desconhecida.')
    for numero_do_eixo, extremidades in enumerate(eixos[representacao]):
        for valor_binario, dados in extremidades.items():
            variavel_do_eixo = dados[0]
            if variavel == variavel_do_eixo:
                return {'eixo': numero_do_eixo, 'extremidade_binaria': valor_binario}
    raise ValueError('A variavel nao pertence a nenhum eixo coordenado.')
orientacoes = {'energia': [{'cruz': 0, 'ponto': 1}, {'cruz': 0, 'ponto': 1}, {'cruz': 1, 'ponto': 0}], 'entropia': [{'cruz': 1, 'ponto': 0}, {'cruz': 1, 'ponto': 0}, {'cruz': 0, 'ponto': 1}]}

def obter_vetor_seta_coordenada(representacao, numero_do_eixo):
    """Obtém o vetor associado à orientação de uma Seta Coordenada."""
    if representacao not in orientacoes:
        raise ValueError('Representacao desconhecida.')
    if numero_do_eixo not in (0, 1, 2):
        raise ValueError('O numero do eixo deve ser 0, 1 ou 2.')
    orientacao = orientacoes[representacao][numero_do_eixo]
    extremidade_cruz = orientacao['cruz']
    extremidade_ponto = orientacao['ponto']
    sentido_binario = extremidade_ponto - extremidade_cruz
    if sentido_binario not in (-1, 1):
        raise RuntimeError('A orientacao do eixo deve ligar extremidades binárias distintas.')
    vetor_base = base_cartesiana[numero_do_eixo]
    vetor_orientado = tuple((sentido_binario * componente for componente in vetor_base))
    return {'representacao': representacao, 'eixo': numero_do_eixo, 'extremidade_cruz': extremidade_cruz, 'extremidade_ponto': extremidade_ponto, 'sentido_binario': sentido_binario, 'vetor_orientado': vetor_orientado}

def identificar_lado(representacao, variavel):
    """Identifica o lado cruz ou ponto a partir da orientação de um vetor."""
    posicao = localizar_extremidade(representacao, variavel)
    numero_do_eixo = posicao['eixo']
    extremidade_binaria = posicao['extremidade_binaria']
    orientacao = orientacoes[representacao][numero_do_eixo]
    if extremidade_binaria == orientacao['cruz']:
        lado = 'cruz'
    elif extremidade_binaria == orientacao['ponto']:
        lado = 'ponto'
    else:
        raise RuntimeError('A extremidade nao corresponde a orientacao registrada.')
    return {'eixo': numero_do_eixo, 'extremidade_binaria': extremidade_binaria, 'lado': lado}

def obter_vetor_posicional_da_variavel(representacao, variavel):
    """Obtém o vetor posicional associado à extremidade ocupada por uma variável."""
    dados_do_lado = identificar_lado(representacao, variavel)
    numero_do_eixo = dados_do_lado['eixo']
    lado = dados_do_lado['lado']
    dados_da_seta = obter_vetor_seta_coordenada(representacao, numero_do_eixo)
    vetor_da_seta = dados_da_seta['vetor_orientado']
    if lado == 'ponto':
        fator_posicional = 1
    elif lado == 'cruz':
        fator_posicional = -1
    else:
        raise RuntimeError("O lado da variavel deve ser 'cruz' ou 'ponto'.")
    vetor_posicional = tuple((fator_posicional * componente for componente in vetor_da_seta))
    return {'representacao': representacao, 'variavel': variavel, 'eixo': numero_do_eixo, 'lado': lado, 'vetor_da_seta': vetor_da_seta, 'vetor_posicional': vetor_posicional}

def obter_extremidades_orientadas(representacao, numero_do_eixo):
    """Retorna as variáveis nas extremidades cruz e ponto de um eixo."""
    orientacao = orientacoes[representacao][numero_do_eixo]
    dados_do_eixo = eixos[representacao][numero_do_eixo]
    extremidade_cruz = orientacao['cruz']
    extremidade_ponto = orientacao['ponto']
    variavel_cruz = dados_do_eixo[extremidade_cruz][0]
    variavel_ponto = dados_do_eixo[extremidade_ponto][0]
    return {'variavel_cruz': variavel_cruz, 'variavel_ponto': variavel_ponto}

def identificar_extremidades_eixo_remanescente(preparacao):
    """Identifica as duas extremidades orientadas do eixo remanescente."""
    dados_dos_eixos = identificar_eixo_remanescente(preparacao)
    eixo_remanescente = dados_dos_eixos['eixo_remanescente']
    extremidades = obter_extremidades_orientadas(preparacao['representacao'], eixo_remanescente)
    return {'eixo_remanescente': eixo_remanescente, 'extremidade_cruz': extremidades['variavel_cruz'], 'extremidade_ponto': extremidades['variavel_ponto']}

def diagnosticar_derivacao_sucessiva(preparacao):
    """Executa o diagnóstico geométrico completo de uma derivação sucessiva."""
    representacao = preparacao['representacao']
    lado_intermediaria = identificar_lado(representacao, preparacao['variavel_intermediaria'])
    lado_segunda_variavel = identificar_lado(representacao, preparacao['segunda_variavel'])
    extremidades_remanescentes = identificar_extremidades_eixo_remanescente(preparacao)
    orientacao_dos_eixos = diagnosticar_orientacao_dos_eixos(preparacao)
    return {'representacao': representacao, 'funcao': preparacao['funcao'], 'primeira_variavel': preparacao['primeira_variavel'], 'segunda_variavel': preparacao['segunda_variavel'], 'variavel_constante': preparacao['variavel_constante'], 'sinal_intermediaria': preparacao['sinal_intermediaria'], 'variavel_intermediaria': preparacao['variavel_intermediaria'], 'intermediaria_formatada': preparacao['intermediaria_formatada'], 'eixo_intermediaria': lado_intermediaria['eixo'], 'lado_intermediaria': lado_intermediaria['lado'], 'eixo_segunda_variavel': lado_segunda_variavel['eixo'], 'lado_segunda_variavel': lado_segunda_variavel['lado'], 'eixo_remanescente': extremidades_remanescentes['eixo_remanescente'], 'vetor_resultante': orientacao_dos_eixos['vetor_resultante'], 'sinal_geometrico_cartesiano': orientacao_dos_eixos['sinal_geometrico_cartesiano'], 'vetor_seta_intermediaria': orientacao_dos_eixos['vetor_seta_intermediaria'], 'vetor_seta_segunda_variavel': orientacao_dos_eixos['vetor_seta_segunda_variavel'], 'vetor_posicional_intermediaria': orientacao_dos_eixos['vetor_posicional_intermediaria'], 'vetor_posicional_segunda_variavel': orientacao_dos_eixos['vetor_posicional_segunda_variavel'], 'vetor_resultante_orientado': orientacao_dos_eixos['vetor_resultante_orientado'], 'vetor_seta_eixo_remanescente': orientacao_dos_eixos['vetor_seta_eixo_remanescente'], 'sinal_geometrico_digraph': orientacao_dos_eixos['sinal_geometrico_digraph'], 'extremidade_cruz_eixo_remanescente': extremidades_remanescentes['extremidade_cruz'], 'extremidade_ponto_eixo_remanescente': extremidades_remanescentes['extremidade_ponto']}

def preparar_duas_ordens(representacao, funcao, variavel_a, variavel_b):
    """Prepara as duas ordens possíveis de derivação para um mesmo par de variáveis."""
    if variavel_a == variavel_b:
        raise ValueError('As duas variaveis devem pertencer a eixos distintos.')
    ordem_a_entao_b = preparar_derivacao_sucessiva(representacao, funcao, variavel_a, variavel_b)
    ordem_b_entao_a = preparar_derivacao_sucessiva(representacao, funcao, variavel_b, variavel_a)
    return {'ordem_a_entao_b': ordem_a_entao_b, 'ordem_b_entao_a': ordem_b_entao_a}

def diagnosticar_duas_ordens(representacao, funcao, variavel_a, variavel_b):
    """Diagnostica geometricamente as duas ordens de uma derivação cruzada."""
    preparacoes = preparar_duas_ordens(representacao, funcao, variavel_a, variavel_b)
    return {'ordem_a_entao_b': diagnosticar_derivacao_sucessiva(preparacoes['ordem_a_entao_b']), 'ordem_b_entao_a': diagnosticar_derivacao_sucessiva(preparacoes['ordem_b_entao_a'])}

def validar_variaveis_duas_ordens(representacao, funcao, variavel_a, variavel_b, esperado_a_entao_b, esperado_b_entao_a):
    """Compara os rótulos geométricos produzidos pelas duas ordens de derivação."""
    preparacoes = preparar_duas_ordens(representacao, funcao, variavel_a, variavel_b)
    resultado_a_entao_b = determinar_variavel_indicada_geometricamente(preparacoes['ordem_a_entao_b'])
    resultado_b_entao_a = determinar_variavel_indicada_geometricamente(preparacoes['ordem_b_entao_a'])
    obtido_a_entao_b = resultado_a_entao_b['variavel_indicada_geometricamente']
    obtido_b_entao_a = resultado_b_entao_a['variavel_indicada_geometricamente']
    confere_a_entao_b = obtido_a_entao_b == esperado_a_entao_b
    confere_b_entao_a = obtido_b_entao_a == esperado_b_entao_a
    return {'ordem_a_entao_b': {'esperado': esperado_a_entao_b, 'obtido': obtido_a_entao_b, 'confere': confere_a_entao_b}, 'ordem_b_entao_a': {'esperado': esperado_b_entao_a, 'obtido': obtido_b_entao_a, 'confere': confere_b_entao_a}, 'validacao_global': confere_a_entao_b and confere_b_entao_a}

def diagnosticar_pares_naturais(representacao, funcao):
    """Executa o diagnóstico para os pares de variáveis naturais de uma função."""
    naturais, _ = montar_regras(representacao, funcao)
    if len(naturais) != 3:
        raise RuntimeError('A funcao deve possuir exatamente tres variaveis naturais.')
    pares = [(naturais[0], naturais[1]), (naturais[0], naturais[2]), (naturais[1], naturais[2])]
    diagnosticos = {}
    for numero_do_par, par in enumerate(pares, start=1):
        variavel_a, variavel_b = par
        diagnosticos['par_' + str(numero_do_par)] = {'variavel_a': variavel_a, 'variavel_b': variavel_b, 'duas_ordens': diagnosticar_duas_ordens(representacao, funcao, variavel_a, variavel_b)}
    return {'representacao': representacao, 'funcao': funcao, 'variaveis_naturais': naturais, 'numero_de_pares': len(pares), 'diagnosticos': diagnosticos}

def diagnosticar_representacao(representacao):
    """Executa diagnósticos para as funções de uma representação termodinâmica."""
    if representacao not in vertices:
        raise ValueError('Representacao desconhecida.')
    funcoes = [funcao for funcao in vertices[representacao] if funcao != 'zero']
    diagnosticos_das_funcoes = {}
    for funcao in funcoes:
        diagnosticos_das_funcoes[funcao] = diagnosticar_pares_naturais(representacao, funcao)
    numero_total_de_pares = sum((diagnostico['numero_de_pares'] for diagnostico in diagnosticos_das_funcoes.values()))
    numero_total_de_sequencias = 2 * numero_total_de_pares
    return {'representacao': representacao, 'numero_de_funcoes': len(funcoes), 'funcoes': funcoes, 'funcoes_excluidas': ['zero'], 'numero_total_de_pares': numero_total_de_pares, 'numero_total_de_sequencias': numero_total_de_sequencias, 'diagnosticos': diagnosticos_das_funcoes}

def comparar_variaveis_duas_ordens(representacao, funcao, variavel_a, variavel_b):
    """Compara as variáveis indicadas geometricamente nas duas ordens de derivação."""
    preparacoes = preparar_duas_ordens(representacao, funcao, variavel_a, variavel_b)
    resultado_a_entao_b = determinar_variavel_indicada_geometricamente(preparacoes['ordem_a_entao_b'])
    resultado_b_entao_a = determinar_variavel_indicada_geometricamente(preparacoes['ordem_b_entao_a'])
    variavel_indicada_a_entao_b = resultado_a_entao_b['variavel_indicada_geometricamente']
    variavel_indicada_b_entao_a = resultado_b_entao_a['variavel_indicada_geometricamente']
    return {'representacao': representacao, 'funcao': funcao, 'variavel_a': variavel_a, 'variavel_b': variavel_b, 'ordem_a_entao_b': {'primeira_variavel': variavel_a, 'segunda_variavel': variavel_b, 'variavel_indicada': variavel_indicada_a_entao_b}, 'ordem_b_entao_a': {'primeira_variavel': variavel_b, 'segunda_variavel': variavel_a, 'variavel_indicada': variavel_indicada_b_entao_a}, 'coincidencia_geometrica': variavel_indicada_a_entao_b == variavel_indicada_b_entao_a}

def comparar_pares_da_funcao(representacao, funcao):
    """Compara as duas ordens para todos os pares naturais de uma função."""
    naturais, _ = montar_regras(representacao, funcao)
    if len(naturais) != 3:
        raise RuntimeError('A funcao deve possuir exatamente tres variaveis naturais.')
    pares = [(naturais[0], naturais[1]), (naturais[0], naturais[2]), (naturais[1], naturais[2])]
    comparacoes = {}
    numero_de_coincidencias = 0
    for numero_do_par, par in enumerate(pares, start=1):
        variavel_a, variavel_b = par
        comparacao = comparar_variaveis_duas_ordens(representacao, funcao, variavel_a, variavel_b)
        nome_do_par = 'par_' + str(numero_do_par)
        comparacoes[nome_do_par] = comparacao
        if comparacao['coincidencia_geometrica']:
            numero_de_coincidencias += 1
    numero_de_divergencias = len(pares) - numero_de_coincidencias
    return {'representacao': representacao, 'funcao': funcao, 'variaveis_naturais': naturais, 'numero_de_pares': len(pares), 'numero_de_coincidencias': numero_de_coincidencias, 'numero_de_divergencias': numero_de_divergencias, 'coincidencia_global': numero_de_divergencias == 0, 'comparacoes': comparacoes}
relacoes_sdigraph = {'positivas': [{'resultado': 'N', 'operacoes': [('1/T', 'P/T'), ('V', '1/T'), ('U', 'V'), ('P/T', 'U')]}, {'resultado': 'mu/T', 'operacoes': [('P/T', '1/T'), ('U', 'P/T'), ('V', 'U'), ('1/T', 'V')]}, {'resultado': '1/T', 'operacoes': [('P/T', 'N'), ('mu/T', 'P/T'), ('V', 'mu/T'), ('N', 'V')]}, {'resultado': 'U', 'operacoes': [('N', 'P/T'), ('V', 'N'), ('mu/T', 'V'), ('P/T', 'mu/T')]}, {'resultado': 'V', 'operacoes': [('N', 'U'), ('1/T', 'N'), ('mu/T', '1/T'), ('U', 'mu/T')]}, {'resultado': 'P/T', 'operacoes': [('U', 'N'), ('mu/T', 'U'), ('1/T', 'mu/T'), ('N', '1/T')]}], 'negativas': [{'resultado': 'mu/T', 'operacoes': [('1/T', '-P/T'), ('P/T', '-U')]}, {'resultado': 'N', 'operacoes': [('U', '-P/T'), ('V', '-U')]}, {'resultado': 'U', 'operacoes': [('P/T', '-N'), ('mu/T', '-P/T')]}, {'resultado': '1/T', 'operacoes': [('V', '-N'), ('N', '-P/T')]}, {'resultado': 'P/T', 'operacoes': [('N', '-U'), ('1/T', '-N')]}, {'resultado': 'V', 'operacoes': [('mu/T', '-U'), ('U', '-N')]}]}

def contar_relacoes(bloco):
    """Conta o número de operações/relações armazenadas em um bloco de referência."""
    return sum((len(cadeia['operacoes']) for cadeia in bloco))

def formatar_operador(variavel):
    """Formata a variável de um operador de derivação para exibição."""
    if '/' in variavel:
        return f'({variavel})'
    return variavel

def formatar_operacao(operacao):
    """Formata uma operação geométrica como uma expressão de derivação."""
    variavel, grandeza = operacao
    return '\\partial/\\partial ' + formatar_operador(variavel) + '(' + grandeza + ')'

def formatar_cadeia(cadeia):
    """Formata uma cadeia de operações do conjunto de referência do S-Digraph."""
    operacoes_formatadas = [formatar_operacao(operacao) for operacao in cadeia['operacoes']]
    return ' \\sim_{S} '.join(operacoes_formatadas) + ' \\sim_{S} ' + cadeia['resultado']

def imprimir_bloco(nome_do_bloco):
    """Imprime as cadeias de um bloco de relações de referência."""
    print(nome_do_bloco)
    for numero, cadeia in enumerate(relacoes_sdigraph[nome_do_bloco], start=1):
        print(numero, formatar_cadeia(cadeia))

def preparar_operacao_geometrica_sdigraph(operacao):
    """Converte uma operação de referência em dados compatíveis com a rotina geométrica.

Cada tupla é interpretada como (variável do operador, grandeza sobre a qual
o operador atua). O sinal da grandeza é separado de sua identidade antes da
avaliação geométrica."""
    variavel_operador, grandeza = operacao
    if grandeza.startswith('-'):
        sinal_intermediaria = -1
        variavel_intermediaria = grandeza[1:]
    else:
        sinal_intermediaria = 1
        variavel_intermediaria = grandeza
    return {'representacao': 'entropia', 'funcao': 'relacao_geometrica', 'primeira_variavel': None, 'segunda_variavel': variavel_operador, 'variavel_constante': None, 'sinal_intermediaria': sinal_intermediaria, 'variavel_intermediaria': variavel_intermediaria, 'intermediaria_formatada': grandeza}

def avaliar_operacao_sdigraph(operacao, resultado_esperado):
    """Avalia geometricamente uma operação do S-Digraph e compara com a referência."""
    preparacao = preparar_operacao_geometrica_sdigraph(operacao)
    indicacao = determinar_variavel_indicada_geometricamente(preparacao)
    resultado_geometrico = indicacao['variavel_indicada_geometricamente']
    return {'operador': operacao[0], 'grandeza': operacao[1], 'resultado_geometrico': resultado_geometrico, 'resultado_esperado': resultado_esperado, 'coincide': resultado_geometrico == resultado_esperado}

def validar_relacoes_sdigraph():
    """Verifica automaticamente as 36 relações geométricas de referência do S-Digraph.

A variável esperada não é fornecida como resposta à rotina geométrica. O
algoritmo determina primeiro o resultado a partir dos eixos, orientações,
Regra da Mão Direita e sinal intermediário; somente depois o compara com o
resultado armazenado no conjunto de referência."""
    numero_total = 0
    numero_de_coincidencias = 0
    divergencias = []
    for nome_do_bloco, cadeias in relacoes_sdigraph.items():
        for numero_da_cadeia, cadeia in enumerate(cadeias, start=1):
            resultado_esperado = cadeia['resultado']
            for numero_da_operacao, operacao in enumerate(cadeia['operacoes'], start=1):
                numero_total += 1
                teste = avaliar_operacao_sdigraph(operacao, resultado_esperado)
                if teste['coincide']:
                    numero_de_coincidencias += 1
                else:
                    divergencias.append({'bloco': nome_do_bloco, 'cadeia': numero_da_cadeia, 'operacao': numero_da_operacao, 'operador': teste['operador'], 'grandeza': teste['grandeza'], 'resultado_geometrico': teste['resultado_geometrico'], 'resultado_esperado': teste['resultado_esperado']})
    numero_de_divergencias = numero_total - numero_de_coincidencias
    return {'numero_total_de_relacoes': numero_total, 'numero_de_coincidencias': numero_de_coincidencias, 'numero_de_divergencias': numero_de_divergencias, 'validacao_global': numero_de_divergencias == 0, 'divergencias': divergencias}


def executar_validacao_verbose():
    """Executa a validação das relações do S-Digraph mostrando cada etapa.

    Esta função é usada quando o arquivo é executado diretamente. Ela não
    altera a lógica de validação: apenas apresenta, de forma legível, cada
    operação geométrica, o resultado determinado pelo algoritmo e o valor
    usado posteriormente como referência para comparação.
    """
    print("=" * 72)
    print("VALIDAÇÃO GEOMÉTRICA DO S-DIGRAPH")
    print("=" * 72)
    print()
    print("O algoritmo determina primeiro a variável indicada geometricamente")
    print("e somente depois compara esse resultado com o conjunto de referência.")
    print()

    total = 0
    coincidencias = 0
    divergencias = []

    for nome_do_bloco, cadeias in relacoes_sdigraph.items():
        quantidade = contar_relacoes(cadeias)
        print("-" * 72)
        print(f"Bloco: {nome_do_bloco.upper()} ({quantidade} relações)")
        print("-" * 72)

        for numero_da_cadeia, cadeia in enumerate(cadeias, start=1):
            resultado_esperado = cadeia["resultado"]
            print()
            print(
                f"Cadeia {numero_da_cadeia}: "
                f"resultado de referência = {resultado_esperado}"
            )

            for numero_da_operacao, operacao in enumerate(
                cadeia["operacoes"], start=1
            ):
                total += 1
                teste = avaliar_operacao_sdigraph(
                    operacao, resultado_esperado
                )

                if teste["coincide"]:
                    coincidencias += 1
                    estado = "OK"
                else:
                    estado = "DIVERGÊNCIA"
                    divergencias.append(
                        {
                            "bloco": nome_do_bloco,
                            "cadeia": numero_da_cadeia,
                            "operacao": numero_da_operacao,
                            "operador": teste["operador"],
                            "grandeza": teste["grandeza"],
                            "resultado_geometrico":
                                teste["resultado_geometrico"],
                            "resultado_esperado":
                                teste["resultado_esperado"],
                        }
                    )

                expressao = formatar_operacao(operacao)
                print(
                    f"  [{estado:10}] operação {numero_da_operacao}: "
                    f"{expressao} -> {teste['resultado_geometrico']} "
                    f"(referência: {teste['resultado_esperado']})"
                )

        print()

    numero_de_divergencias = total - coincidencias

    print("=" * 72)
    print("RESUMO DA VALIDAÇÃO")
    print("=" * 72)
    print(f"Relações verificadas : {total}")
    print(f"Coincidências         : {coincidencias}")
    print(f"Divergências          : {numero_de_divergencias}")
    print(
        "Validação global      : "
        + ("SUCESSO" if numero_de_divergencias == 0 else "FALHA")
    )

    if divergencias:
        print()
        print("Divergências encontradas:")
        for divergencia in divergencias:
            print("  ", divergencia)


def _nome_representacao(representacao):
    """Retorna um nome legível para a representação usada na saída."""
    if representacao == "energia":
        return "E-DIGRAPH — representação da Energia"
    if representacao == "entropia":
        return "S-DIGRAPH — representação da Entropia"
    return representacao


def _formatar_sinal(sinal):
    """Formata um sinal numérico para exibição."""
    return "+" if sinal == 1 else "-"


def executar_representacao_verbose(representacao):
    """Executa todas as funções não nulas de uma representação.

    Para cada função, o programa:
      1. reconstrói as três variáveis naturais a partir do vértice;
      2. executa as três derivadas parciais de primeira ordem;
      3. forma os três pares possíveis de variáveis naturais;
      4. executa cada par nas duas ordens;
      5. mostra a variável intermediária, os eixos envolvidos e a variável
         indicada geometricamente em cada ordem.

    A comparação entre as duas ordens é apenas uma comparação de rótulos
    geométricos. Rótulos distintos não significam, por si só, desigualdade
    analítica entre derivadas cruzadas.
    """
    diagnostico = diagnosticar_representacao(representacao)

    print("=" * 78)
    print(_nome_representacao(representacao))
    print("=" * 78)
    print(
        f"Funções analisadas: {diagnostico['numero_de_funcoes']} | "
        f"pares de variáveis: {diagnostico['numero_total_de_pares']} | "
        f"ordens sucessivas: {diagnostico['numero_total_de_sequencias']}"
    )
    print("Potencial Nulo: excluído da rotina geral de derivação.")
    print()

    total_primeira_ordem = 0
    total_ordens_sucessivas = 0

    for funcao in diagnostico["funcoes"]:
        coordenada = vertices[representacao][funcao]
        naturais, _ = montar_regras(representacao, funcao)

        print("-" * 78)
        print(f"Função: {funcao} | vértice: {coordenada}")
        print("Variáveis naturais: " + ", ".join(naturais))
        print()
        print("Derivadas parciais de primeira ordem:")

        for variavel in naturais:
            resultado = derivar(representacao, funcao, variavel)
            total_primeira_ordem += 1
            constantes = ", ".join(resultado["constantes"])
            print(
                f"  ∂{funcao}/∂{variavel} | constantes: {constantes} "
                f"-> {resultado['resultado']}"
            )

        pares = [
            (naturais[0], naturais[1]),
            (naturais[0], naturais[2]),
            (naturais[1], naturais[2]),
        ]

        print()
        print("Derivadas sucessivas — três pares, nas duas ordens:")

        for numero_do_par, (variavel_a, variavel_b) in enumerate(pares, start=1):
            preparacoes = preparar_duas_ordens(
                representacao, funcao, variavel_a, variavel_b
            )

            resultados = []
            print(f"  Par {numero_do_par}: ({variavel_a}, {variavel_b})")

            for chave in ("ordem_a_entao_b", "ordem_b_entao_a"):
                preparacao = preparacoes[chave]
                diagnostico_ordem = diagnosticar_derivacao_sucessiva(preparacao)
                indicacao = determinar_variavel_indicada_geometricamente(preparacao)
                total_ordens_sucessivas += 1

                primeira = preparacao["primeira_variavel"]
                segunda = preparacao["segunda_variavel"]
                intermediaria = preparacao["intermediaria_formatada"]
                indicada = indicacao["variavel_indicada_geometricamente"]
                resultados.append(indicada)

                print(
                    f"    ordem {primeira} -> {segunda}: "
                    f"intermediária = {intermediaria}; "
                    f"eixos = {diagnostico_ordem['eixo_intermediaria']} x "
                    f"{diagnostico_ordem['eixo_segunda_variavel']} -> "
                    f"remanescente {diagnostico_ordem['eixo_remanescente']}; "
                    f"sinal geométrico = "
                    f"{_formatar_sinal(diagnostico_ordem['sinal_geometrico_digraph'])}; "
                    f"variável indicada = {indicada}"
                )

            if resultados[0] == resultados[1]:
                print(f"    comparação geométrica: mesmo rótulo ({resultados[0]})")
            else:
                print(
                    "    comparação geométrica: rótulos distintos "
                    f"({resultados[0]} e {resultados[1]})"
                )
                print(
                    "    observação: isso descreve a indicação geométrica; "
                    "não é um teste da igualdade analítica de Schwarz."
                )

        print()

    print("-" * 78)
    print(f"RESUMO — {_nome_representacao(representacao)}")
    print(f"Derivadas de primeira ordem executadas : {total_primeira_ordem}")
    print(f"Ordens de derivação sucessiva executadas: {total_ordens_sucessivas}")
    print()


def executar_tudo_verbose():
    """Executa E-Digraph, S-Digraph e a validação das 36 relações do S-Digraph."""
    print()
    print("IMPLEMENTAÇÃO COMPUTACIONAL DOS E-DIGRAPH E S-DIGRAPH")
    print("Execução completa das duas representações")
    print()

    executar_representacao_verbose("energia")
    executar_representacao_verbose("entropia")

    print()
    print("A seguir é executada uma validação adicional específica do S-Digraph:")
    print("as 36 relações geométricas usadas como conjunto de referência.")
    print()
    executar_validacao_verbose()


if __name__ == "__main__":
    executar_tudo_verbose()
