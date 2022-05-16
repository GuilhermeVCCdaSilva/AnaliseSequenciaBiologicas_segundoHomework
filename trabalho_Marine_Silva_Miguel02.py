#!/usr/bin/env python3
import sys


def listaParaString(listaDeSequencias):
     '''
    Converte uma lista em string, acrescentando cada elementos 
    da "listaDeSequencias" a uma string "listaDeSequencias".


    arguments: list - listaDaSequencias
    return: string - sequencia
    '''
    sequencia = "" 
    for i in listaDeSequencias: 
        sequencia += i  
    return sequencia


def obterNomeMaiorTamanho(dicionarioNomeSequencia):
    '''
    Devolve o tamanho da chave da sequencia que tem maior tamanho.


    argumento: dict - dicionarioNomeSequencia 
    return: int - maiorTamanho
    '''
    listaDeChaves = dicionarioNomeSequencia.keys()
    maiorTamanho = max(map(len,listaDeChaves))
    return maiorTamanho


def obterFicheiroFasta(argumento):
    '''
    Devolve as linhas do ficheiro "fasta" 
    escrito como argumento


    argumento: string - argumento
    return: list - listaDeLinhas
    '''
    listaDeLinhas = []
    argumento = sys.argv[1] 
    with open(argumento, 'r') as fasta:
        for linha in fasta:
            listaDeLinhas.append(linha)
    return listaDeLinhas


def obterDicionarioNomeSequenciaNtax(listaDeLinhas,dicionarioNomeSequencia):
    '''
    Percorre linha a linha o ficheiro "FASTA", começando por tirar os "\n". 
    Em seguida, verifica se têm o sinal de ">" e se for verdade, 
    cria uma nova chave no dicionário com o nome da sequência  se este nome ainda não existir no dicionário.
    Asociado ao nome vai uma lista vazia que vai ser acrescentada com as sequências linha a linha até surgir o novo sinal ">" e o processo volta a iterar.
    No final retorna o dicionário Nome;chave/Sequência;valor

    argumento: list - listaDeLinhas; dict - dicionarioNomeSequencia 
    return: dict - dicionarioNomeSequencia
    '''
    dicionarioNomeSequencia = {}
    for linha in listaDeLinhas:
        linha = linha.strip("\n") 
        if linha.startswith(">"):
            nome = linha[1:]
            if nome not in dicionarioNomeSequencia:
                dicionarioNomeSequencia[nome] = []
            continue    
        sequencia = linha
        if len(sequencia) > 99:
            sequencia = "A sequência têm mais de 99 caracteres!"
            break
        dicionarioNomeSequencia[nome].append(sequencia)
    return dicionarioNomeSequencia


def escreverNexusAlinhado(dicionarioNomeSequencia):
    '''
    Acrescenta ao dicionário os nomes ajustados à direta com o espaçamento do nome com maior tamanho.
    No Final faz o print com os nomes ajustados e a lista de sequência associada, passada em string única.

    argumento: dict - dicionarioNomeSequencia  
    '''
    dicionarioNomeSequenciaAlinhado = {nome.rjust(obterNomeMaiorTamanho(dicionarioNomeSequencia)+1) : sequencia for nome, sequencia in dicionarioNomeSequencia.items()}
    for nome, sequencia in dicionarioNomeSequenciaAlinhado.items():
            print(nome,listaParaString(sequencia))    


def formatarNexusMrBayes(dicionarioNomeSequencia):
    '''
    Formata corretamente o "dicionarioNomeSequencia" para Nexus com o "ntax" e "nchar".
    É acrescentado também o MrBayes associado com os argumentos "outgroup" e "ngen" escritos pelo utilizador. 
    argumento: dict - dicionarioNomeSequencia  
    '''
    ntax = len(dicionarioNomeSequencia.keys())
    nchar = sum(len(i) for i in list(dicionarioNomeSequencia.items())[0][1])
    print("""#NEXUS

BEGIN DATA;
DIMENSIONS NTAX={} NCHAR={};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX""".format(ntax, nchar))
    escreverNexusAlinhado(dicionarioNomeSequencia)
    print("""  ;
END;""")
    print("""
begin mrbayes;
    set autoclose=yes;
    outgroup={};
    mcmcp ngen={} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=ChangeMyNameMrBayes;
    mcmc;
    sumt filename=ChangeMyNameMrBayes;
end;""".format(sys.argv[2],sys.argv[3]))


if __name__ == '__main__':
    argumento = ""
    dicionarioNomeSequencia = ""
    listaDeLinhas = obterFicheiroFasta(argumento)
    dicionarioNomeSequencia = obterDicionarioNomeSequenciaNtax(listaDeLinhas,dicionarioNomeSequencia)
    formatarNexusMrBayes(dicionarioNomeSequencia)
