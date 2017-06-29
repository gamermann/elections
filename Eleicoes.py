#!/usr/bin/python
#encoding=latin-1

import os
from math import log, exp
from Gnuplot import Gnuplot as gplot
from Gnuplot import Data as gdata
from time import time
from ExGUtils.uts import drand, stats
from scipy.stats import chi2 as chi2_distr

g = gplot(persist=1)


#
# This file contains classes to read and structure the data contained
# in the TSE files concerning electoral campaigns donations. 
# It also contains functions and stuff used in the analysis to produce
# the plots and tables (in latex). 
#

# tex stuff
hline = "\\hline"
cline = "\\cline{%i-%i}"
table_start = "\\btab[H]\n\\small\\caption{%s}\\label{%s}\n\\bc\n\\bt{%s}"
table_start_tiny = "\\btab[H]\n\\caption{%s}\\label{%s}\n\\tiny\n\\bc\n\\bt{%s}"
table_end = "\\et\n\\ec\n\\etab"



# Benford
digis = range(1, 10)
bf = [log(1.+1./ele)/log(10) for ele in digis]
dbf = gdata(digis, bf, with_="boxes", title="Benford")



# functions

def acumulator(lista):
    xs = list(set(lista))
    xs.sort()
    ys = [lista.count(ele) for ele in xs]
    tot = sum(ys)*1.
    ys = [ele/tot for ele in ys]
    y = [sum(ys[:ii]) for ii in xrange(len(ys))]
    return xs, y



def get_digit(num):
    """
    given a number (num) return its first significant digit
    """
    if not num:
        return 0
    if num<0.:
        num = -num
    while num<1.:
        num *= 10
    return int(str(num)[0])
        

def benford_stuff(numbers, name="", wtot=False):
    """
    Given a list of numbers, returns a string containing in latex format a
    line for the table and a gnuplot object with the histogram of digits
    compared with benford distribution.
    """
    digs = [get_digit(ele) for ele in numbers]
    N = len(digs)
    distr = [digs.count(ele)*1./N for ele in digis]
    errs = [((1.-ele)/ N)**.5  if ele!=0 else 0. for ele in distr]
    dd = gdata(digis, distr, errs, with_="err", title="%s"%name)
    chi2 = sum([N*(distr[ii]-ele)**2 / ele for ii, ele in enumerate(bf)])
    p_value = (1 - chi2_distr.cdf(chi2, 8))
    if wtot:
        stri = " & ".join(["%4.3f"%ele for ele in distr])
        line = "%s & %s & %i & %4.3f & %4.3f & %8.2f & %8.2f & %8.2f"%(name, stri, N, chi2, p_value, min(numbers), max(numbers), sum(numbers))
    else:
        stri = " & ".join(["%4.3f"%ele for ele in distr])
        line = "%s & %s & %i & %4.3f & %4.3f"%(name, stri, N, chi2, p_value)
    return line, dd





##############
# Distribution stuff
##############


def drand_distr(gamma, csi0, csim=100., delt=0.):
    F = drand()
    bla = csi0/(csim-delt)
    blabla = F/(bla**gamma+1-F)
    return exp(csi0*blabla**(1./gamma)+delt)

def func(x, gamma, csi0, csim=100., delt=0.):
    # cumulative distribution
    num1 = log(x)-delt
    num2 = (csi0/num1)**gamma
    num3 = (csi0/(csim-delt))**gamma
    term2 = 1.+num3
    term3 = 1./(num2+1.)
    return term2*term3


def func2(x, gamma, csi0, csim=100., delt=0.):
    # distribution
    num1 = log(x)-delt
    num2 = (csi0/num1)**gamma
    num3 = (csi0/(csim-delt))**gamma
    term1 = gamma/x
    term2 = 1.+num3
    term3 = num2/( (num2+1.)**2 )
    term3 /= num1
    return term1*term2*term3

    


def lkhd_distr(nums, gamma, csi0, csim, delt):
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum5 = 0.
    sum6 = 0.
    N = len(nums)
    num1 = csi0/(csim-delt)
    num1g = num1**gamma
    termb = num1g/(1.+num1g)
    for ele in nums:
        num2 = log(ele)
        num3 = num2-delt
        num4 = csi0/(num3)
        num4g = num4**gamma
        term1 = num4g/(num4g+1.)
        # for lkhd
        sum1 += num2
        sum2 += log(num3)
        sum3 += log(1.+num4g)
        # for dlnldcsi0
        sum4 += term1
        # for dlnldgamma
        sum5 += log(num4)*term1
    lkhd = N*(log(gamma)+log(1.+num1g)+gamma*log(csi0))-sum1-(gamma+1.)*sum2-2.*sum3
    g1 = N*(gamma*termb/csi0 + gamma/csi0)-2.*gamma*sum4/csi0 #csi0
    g2 = -N*gamma*termb/(csim-delt) #csim
    g3 = N*(1./gamma + termb*log(num1)+log(csi0)) - sum2-2.*sum5 #gamma
    return lkhd, (g1, g2, g3)


### This algorithm is very unstable for some small sets of numbers
def maxLKHD_distr(nums, gamma=3., csi0=3., csim=100., delt=0, lamb=1., eps=1.e-4):
    lkhd, grad = lkhd_distr(nums, gamma, csi0, csim, delt)
    norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
    ncsi0 = csi0+lamb*grad[0]/norm
    ncsim = csim#+lamb*grad[1]/norm
    ngamma = gamma+lamb*grad[2]/norm
    nlkhd, ngrad = lkhd_distr(nums, ngamma, ncsi0, ncsim, delt)
    while lamb*norm > eps:
        #print lamb, norm
        if nlkhd>lkhd: # accepted
            grad = ngrad
            csi0 = ncsi0
            csim = ncsim
            gamma = ngamma
            lkhd = nlkhd
            lamb *= 1.2
        else:
            lamb *= .01
        norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
        ncsi0 = csi0+lamb*grad[0]/norm
        ncsim = csim#+lamb*grad[1]/norm
        ngamma = gamma+lamb*grad[2]/norm
        nlkhd, ngrad = lkhd_distr(nums, ngamma, ncsi0, ncsim, delt)
    return gamma, csi0, csim






#################
# classes
##################

class Doacao_C:
    """
    Structures the information about a single donation made to a candidate. 
    The list are the elements contained in a line of a doacoes_candidatos file.
    """    
    def __init__(self, ele):
        self.is_candidato = True
        self.is_partido = False
        self.is_comite = False
        self.data_hora = ele[2]
        self.CNPJprestador = ele[3]
        self.UF = ele[5]
        self.partido = ele[6]
        self.numero_cand = ele[7]
        self.cargo = ele[8]
        self.nome_cand = ele[9]
        self.CPF_cand = ele[10]
        self.numero_recibo = ele[11]
        self.numero_doc = ele[12]
        cp = ele[13]
        self.cp = cp
        if len(cp)==11:
            self.is_CPF = True
            self.is_CNPJ = False
            self.CPF = cp
        elif len(cp)==14:
            self.is_CNPJ = True
            self.is_CPF = False
            self.CNPJ = cp
        else:
            self.is_CNPJ = False
            self.is_CPF = False
        self.Nome = ele[14]
        self.Nome_Receita = ele[15]
        self.UE_doador = ele[16]
        self.Numero_partido = ele[17]
        self.numero_cand_doador = ele[18]
        self.Cod_Setor_Economico = ele[19]
        self.Setor_Economico = ele[20]
        self.data_hora_receita = ele[21]
        self.valor = float(ele[22])
        self.tipo = ele[23]
        self.fonte = ele[24]
        self.especie = ele[25]
        self.descricao = ele[26]
        cp = ele[27]
        self.cp_ori = cp
        if len(cp)==11:
            self.is_CPF_ori = True
            self.is_CNPJ_ori = False
            self.CPF_ori = cp
            self.is_ori = False
        elif len(cp)==14:
            self.is_CNPJ_ori = True
            self.is_CPF_ori = False
            self.CNPJ_ori = cp
            self.is_ori = False
        else:
            self.is_CNPJ_ori = False
            self.is_CPF_ori = False
            self.is_ori = True
        self.Nome_ori = ele[28]        
        self.tipo_ori = ele[29]
        self.Setor_Economico_ori = ele[30]
        self.Nome_Receita_ori = ele[31]
    
    def __repr__(self):
        stri = "<%s recebeu de %s %f R$.>"%(self.nome_cand, self.Nome, self.valor)
        return stri


class Doacao_P:
    """
    Structures the information about a single donation made to a party. 
    The list are the elements contained in a line of a doacoes_partidoss file.
    """        
    def __init__(self, ele):
        self.is_candidato = False
        self.is_partido = True
        self.is_comite = False
        self.data_hora = ele[2]
        self.CNPJprestador = ele[3]
        self.UF = ele[5]
        self.Tipo_Diretorio = ele[6]
        self.partido = ele[7]
        self.Tipo_doc = ele[8]
        self.Numero_doc = ele[9]
        cp = ele[10]
        if len(cp)==11:
            self.is_CPF = True
            self.is_CNPJ = False
            self.CPF = cp
        elif len(cp)==14:
            self.is_CNPJ = True
            self.is_CPF = False
            self.CNPJ = cp
        else:
            self.is_CNPJ = False
            self.is_CPF = False
        self.Nome = ele[11]
        self.Nome_Receita = ele[12]
        self.UE_doador = ele[13]
        self.Numero_partido = ele[14]
        self.Numero_cand_doador = ele[15]
        self.Cod_Setor_Economico = ele[16]
        self.Setor_Economico = ele[17]
        self.data_hora_receita = ele[18]
        self.valor = float(ele[19])
        self.tipo = ele[20]
        self.fonte = ele[21]
        self.especie = ele[22]
        self.descricao = ele[23]
        cp = ele[24]
        self.cp_ori = cp
        if len(cp)==11:
            self.is_CPF_ori = True
            self.is_CNPJ_ori = False
            self.CPF_ori = cp
            self.is_ori = False
        elif len(cp)==14:
            self.is_CNPJ_ori = True
            self.is_CPF_ori = False
            self.CNPJ_ori = cp
            self.is_ori = False
        else:
            self.is_CNPJ_ori = False
            self.is_CPF_ori = False
            self.is_ori = True
        self.Nome_ori = ele[25]        
        self.tipo_ori = ele[26]
        self.Setor_Economico_ori = ele[27]
        self.Nome_Receita_ori = ele[28]
    
    def __repr__(self):
        stri = "<%s recebeu de %s %f R$.>"%(self.CNPJprestador, self.Nome, self.valor)
        return stri
        

class Doacao_Comite:
    """
    Structures the information about a single donation made to a committee. 
    The list are the elements contained in a line of a doacoes_comites file.
    """    
    def __init__(self, ele):
        self.is_candidato = False
        self.is_partido = False
        self.is_comite = True
        self.data_hora = ele[2]
        self.CNPJprestador = ele[3]
        self.UF = ele[5]
        self.Tipo_Comite = ele[6]
        self.partido = ele[7]
        self.Tipo_doc = ele[8]
        self.Numero_doc = ele[9]
        cp = ele[10]
        if len(cp)==11:
            self.is_CPF = True
            self.is_CNPJ = False
            self.CPF = cp
        elif len(cp)==14:
            self.is_CNPJ = True
            self.is_CPF = False
            self.CNPJ = cp
        else:
            self.is_CNPJ = False
            self.is_CPF = False
        self.Nome = ele[11]
        self.Nome_Receita = ele[12]
        self.UE_doador = ele[13]
        self.Numero_partido = ele[14]
        self.Numero_cand_doador = ele[15]
        self.Cod_Setor_Economico = ele[16]
        self.Setor_Economico = ele[17]
        self.data_hora_receita = ele[18]
        self.valor = float(ele[19])
        self.tipo = ele[20]
        self.fonte = ele[21]
        self.especie = ele[22]
        self.descricao = ele[23]
        cp = ele[24]
        self.cp_ori = cp
        if len(cp)==11:
            self.is_CPF_ori = True
            self.is_CNPJ_ori = False
            self.CPF_ori = cp
            self.is_ori = False
        elif len(cp)==14:
            self.is_CNPJ_ori = True
            self.is_CPF_ori = False
            self.CNPJ_ori = cp
            self.is_ori = False
        else:
            self.is_CNPJ_ori = False
            self.is_CPF_ori = False
            self.is_ori = True
        self.Nome_ori = ele[25]        
        self.tipo_ori = ele[26]
        self.Setor_Economico_ori = ele[27]
        self.Nome_Receita_ori = ele[28]
    
    def __repr__(self):
        stri = "<%s recebeu de %s %f R$.>"%(self.CNPJprestador, self.Nome, self.valor)
        return stri
        

class Doacoes_Candidatos:
    """
    Structures the information for all donations made to the candidates. 
    The list contains the elements contained in all lines of a doacoes_candidatos file.
    """        
    def __init__(self, data, cargo=""):
        self.data = [Doacao_C(ele) for ele in data[1:]]
        if cargo:
            self.data = [ele for ele in self.data if ele.cargo==cargo]
        else:
            self.data = [ele for ele in self.data]
        self.N = len(self.data)
        self.cargo = ""
        self.candidatos = list(set([ele.nome_cand for ele in self.data]))
        self.doadores = list(set([ele.Nome for ele in self.data]))
        self.doadores_ori = list(set([ele.Nome_ori for ele in self.data]))
        self.Ncands = len(self.candidatos)
        self.Ndoad = len(self.doadores)
        self.Ndoad_ori = len(self.doadores_ori)
        self.iCPFs = [ele for ele in xrange(self.N) if self.data[ele].is_CPF]
        self.iCNPJs = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ]
        self.iCPFs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CPF_ori]
        self.iCNPJs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ_ori]
        self.moneys = [ele.valor for ele in self.data]
        self.partidos = [ele.partido for ele in self.data]
        self.money_tot = round(sum(self.moneys), 2)
        self.totals_cand = [self.total_candidato(ele) for ele in self.candidatos]
        self.totals_doad = [self.total_doador(ele) for ele in self.doadores]
        self.totals_doad_ori = [self.total_doador_ori(ele) for ele in self.doadores_ori]
    
    def __repr__(self):
        if self.cargo:
            stri = "< Candidatos a %s : %i \n  Doadores : %i \n  Doadores_ori : %i \n  $ total : %f >"%(self.cargo, self.Ncands, self.Ndoad, self.Ndoad_ori, self.money_tot)
        else:
            stri = "< Todos candidatos : %i \n  Doadores : %i \n  Doadores_ori : %i \n  $ total : %f >"%(self.Ncands, self.Ndoad, self.Ndoad_ori, self.money_tot)
        return stri
    
    def total_candidato(self, cand):
        soma = sum([ele.valor for ele in self.data if ele.nome_cand==cand])
        return soma
        
    def total_cpf(self, cpf):
        soma = sum([ele.valor for ele in self.data if ele.CPF_cand==cpf])
        return soma
        
    def doacoes_por_cpf(self, cpf):
        soma = [ele.valor for ele in self.data if ele.CPF_cand==cpf]
        return soma
        
    def total_doador(self, doad):
        soma = sum([ele.valor for ele in self.data if ele.Nome==doad])
        return soma
        
    def total_doador_ori(self, doad_ori):
        soma = sum([ele.valor for ele in self.data if ele.Nome_ori==doad_ori])
        return soma
    
    def calc_dics(self):
        CNPJs = set([self.data[ii].CNPJ for ii in self.iCNPJs])
        CNPJs_ori = set([self.data[ii].CNPJ_ori for ii in self.iCNPJs_ori])
        CPFs = set([self.data[ii].CPF for ii in self.iCPFs])
        CPFs_ori = set([self.data[ii].CPF_ori for ii in self.iCPFs_ori])
        CNPJs_t = list(CNPJs|CNPJs_ori)
        CPFs_t = list(CPFs|CPFs_ori)
        Ncnpj = len(CNPJs_t)
        Ncpf = len(CPFs_t)
        self.CNPJs_t = CNPJs_t
        self.CPFs_t = CPFs_t
        dic_CNPJ = dict(zip(CNPJs_t, [[] for ii in xrange(Ncnpj)]))
        dic_CPF = dict(zip(CPFs_t, [[] for ii in xrange(Ncpf)]))
        for ii in list(set(self.iCNPJs)):
            cnpj = self.data[ii].CNPJ
            nome = self.data[ii].Nome
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCNPJs_ori)):
            cnpj = self.data[ii].CNPJ_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCPFs)):
            cpf = self.data[ii].CPF
            nome = self.data[ii].Nome
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        for ii in list(set(self.iCPFs_ori)):
            cpf = self.data[ii].CPF_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        self.dic_CPF = dic_CPF
        self.dic_CNPJ = dic_CNPJ






class Doacoes_Partidos:
    """
    Structures the information for all donations made to the parties. 
    The list contains the elements contained in all lines of a doacoes_partidos file.
    """        
    def __init__(self, data):
        self.data = [Doacao_P(ele) for ele in data[1:]]
        self.N = len(self.data)
        self.Diretorios = list(set([ele.Tipo_Diretorio for ele in self.data]))
        self.Docs = list(set([ele.Tipo_doc for ele in self.data]))
        self.Numeros_doc = [ele.Numero_doc for ele in self.data]
        self.doadores = list(set([ele.Nome for ele in self.data]))
        self.doadores_ori = list(set([ele.Nome_ori for ele in self.data]))
        self.NDiretorios = len(self.Diretorios)
        self.NDocs = len(self.Docs)
        self.Ndoad = len(self.doadores)
        self.Ndoad_ori = len(self.doadores_ori)
        self.iCPFs = [ele for ele in xrange(self.N) if self.data[ele].is_CPF]
        self.iCNPJs = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ]
        self.iCPFs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CPF_ori]
        self.iCNPJs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ_ori]
        self.moneys = [ele.valor for ele in self.data]
        self.partidos = [ele.partido for ele in self.data]
        self.money_tot = round(sum(self.moneys), 2)
        self.totals_doad = [self.total_doador(ele) for ele in self.doadores]
        self.totals_doad_ori = [self.total_doador_ori(ele) for ele in self.doadores_ori]
    
    def __repr__(self):
        stri = "< Doadores : %i \n  Doadores_ori : %i \n  $ total : %f >"%(self.Ndoad, self.Ndoad_ori, self.money_tot)
        return stri
    
    def total_doador(self, doad):
        soma = sum([ele.valor for ele in self.data if ele.Nome==doad])
        return soma
        
    def total_doador_ori(self, doad_ori):
        soma = sum([ele.valor for ele in self.data if ele.Nome_ori==doad_ori])
        return soma
    
    def calc_dics(self):
        CNPJs = set([self.data[ii].CNPJ for ii in self.iCNPJs])
        CNPJs_ori = set([self.data[ii].CNPJ_ori for ii in self.iCNPJs_ori])
        CPFs = set([self.data[ii].CPF for ii in self.iCPFs])
        CPFs_ori = set([self.data[ii].CPF_ori for ii in self.iCPFs_ori])
        CNPJs_t = list(CNPJs|CNPJs_ori)
        CPFs_t = list(CPFs|CPFs_ori)
        Ncnpj = len(CNPJs_t)
        Ncpf = len(CPFs_t)
        self.CNPJs_t = CNPJs_t
        self.CPFs_t = CPFs_t
        dic_CNPJ = dict(zip(CNPJs_t, [[] for ii in xrange(Ncnpj)]))
        dic_CPF = dict(zip(CPFs_t, [[] for ii in xrange(Ncpf)]))
        for ii in list(set(self.iCNPJs)):
            cnpj = self.data[ii].CNPJ
            nome = self.data[ii].Nome
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCNPJs_ori)):
            cnpj = self.data[ii].CNPJ_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCPFs)):
            cpf = self.data[ii].CPF
            nome = self.data[ii].Nome
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        for ii in list(set(self.iCPFs_ori)):
            cpf = self.data[ii].CPF_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        self.dic_CPF = dic_CPF
        self.dic_CNPJ = dic_CNPJ
            



class Doacoes_Comites:
    """
    Structures the information for all donations made to the comitteess. 
    The list contains the elements contained in all lines of a doacoes_comites file.
    """        
    def __init__(self, data):
        self.data = [Doacao_Comite(ele) for ele in data[1:]]
        self.N = len(self.data)
        self.comites = list(set([ele.Tipo_Comite for ele in self.data]))
        self.Docs = list(set([ele.Tipo_doc for ele in self.data]))
        self.Numeros_doc = [ele.Numero_doc for ele in self.data]
        self.doadores = list(set([ele.Nome for ele in self.data]))
        self.doadores_ori = list(set([ele.Nome_ori for ele in self.data]))
        self.Ncomites = len(self.comites)
        self.NDocs = len(self.Docs)
        self.Ndoad = len(self.doadores)
        self.Ndoad_ori = len(self.doadores_ori)
        self.iCPFs = [ele for ele in xrange(self.N) if self.data[ele].is_CPF]
        self.iCNPJs = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ]
        self.iCPFs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CPF_ori]
        self.iCNPJs_ori = [ele for ele in xrange(self.N) if self.data[ele].is_CNPJ_ori]
        self.moneys = [ele.valor for ele in self.data]
        self.partidos = [ele.partido for ele in self.data]
        self.money_tot = round(sum(self.moneys), 2)
        self.totals_doad = [self.total_doador(ele) for ele in self.doadores]
        self.totals_doad_ori = [self.total_doador_ori(ele) for ele in self.doadores_ori]
    
    def __repr__(self):
        stri = "< Doadores : %i \n  Doadores_ori : %i \n  $ total : %f >"%(self.Ndoad, self.Ndoad_ori, self.money_tot)
        return stri
    
    def total_doador(self, doad):
        soma = sum([ele.valor for ele in self.data if ele.Nome==doad])
        return soma
        
    def total_doador_ori(self, doad_ori):
        soma = sum([ele.valor for ele in self.data if ele.Nome_ori==doad_ori])
        return soma
    
    def calc_dics(self):
        CNPJs = set([self.data[ii].CNPJ for ii in self.iCNPJs])
        CNPJs_ori = set([self.data[ii].CNPJ_ori for ii in self.iCNPJs_ori])
        CPFs = set([self.data[ii].CPF for ii in self.iCPFs])
        CPFs_ori = set([self.data[ii].CPF_ori for ii in self.iCPFs_ori])
        CNPJs_t = list(CNPJs|CNPJs_ori)
        CPFs_t = list(CPFs|CPFs_ori)
        Ncnpj = len(CNPJs_t)
        Ncpf = len(CPFs_t)
        self.CNPJs_t = CNPJs_t
        self.CPFs_t = CPFs_t
        dic_CNPJ = dict(zip(CNPJs_t, [[] for ii in xrange(Ncnpj)]))
        dic_CPF = dict(zip(CPFs_t, [[] for ii in xrange(Ncpf)]))
        for ii in list(set(self.iCNPJs)):
            cnpj = self.data[ii].CNPJ
            nome = self.data[ii].Nome
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCNPJs_ori)):
            cnpj = self.data[ii].CNPJ_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CNPJ[cnpj]:
                dic_CNPJ[cnpj].append(nome)
        for ii in list(set(self.iCPFs)):
            cpf = self.data[ii].CPF
            nome = self.data[ii].Nome
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        for ii in list(set(self.iCPFs_ori)):
            cpf = self.data[ii].CPF_ori
            nome = self.data[ii].Nome_ori
            if nome not in dic_CPF[cpf]:
                dic_CPF[cpf].append(nome)
        self.dic_CPF = dic_CPF
        self.dic_CNPJ = dic_CNPJ



class Eleicao_UF:
    """
    Reads all donnations for a given election.
    """    
    def __init__(self, candidatos="", partidos="", comites="", cargo=""):
        if candidatos:
            bla = open(candidatos).read()
            bla = bla.replace(',', ".")
            bla = bla.split("\n")
            data = [[ele2.replace('"', '') for ele2 in ele.split('";"')] for ele in bla]
            data.pop()
            cargos = list(set([ele[8] for ele in data[1:]]))
            if cargo:
                self.candidatos = Doacoes_Candidatos(data, cargo=cargo)
            else:
                self.candidatos = Doacoes_Candidatos(data)
            self.cargos = list(set([ele.cargo for ele in self.candidatos.data]))
        else:
            self.candidatos = Doacoes_Candidatos([])
        #
        if partidos:
            bla = open(partidos).read()
            bla = bla.replace(',', ".")
            bla = bla.split("\n")
            data = [[ele2.replace('"', '') for ele2 in ele.split('";"')] for ele in bla]
            data.pop()
            self.partidos = Doacoes_Partidos(data)
        else:
            self.partidos = Doacoes_Partidos([])
        #
        if comites:
            bla = open(comites).read()
            bla = bla.replace(',', ".")
            bla = bla.split("\n")
            data = [[ele2.replace('"', '') for ele2 in ele.split('";"')] for ele in bla]
            data.pop()
            self.comites = Doacoes_Comites(data)
        else:
            self.comites = Doacoes_Comites([])
        





class Candidato:
    """
    Creates a candidato structure
    """    
    def __init__(self, nome, num, cpf, cargo, total, votos, partido, coligacao, situacao, validos):
        self.nome = nome
        self.numero = num
        self.cpf = cpf
        self.cargo = cargo
        self.money = total
        self.votos = votos
        self.partido = partido
        self.coligacao = coligacao
        self.situacao = situacao
        self.validos = validos
    
    def __repr__(self):
        stri = "< %s candidato a %s recebeu %f R$ e foi %s>"%(self.nome, self.cargo, self.money, self.situacao)
        return stri






class Candidatos:
    
    def __init__(self, donations="", csv="", cargo=""):
        # reads donnations
        bla = open(donations).read()
        bla = bla.replace(',', ".")
        bla = bla.split("\n")
        data = [[ele2.replace('"', '') for ele2 in ele.split('";"')] for ele in bla]
        data.pop()
        cargos = list(set([ele[8] for ele in data[1:]]))
        if cargo:
            donations = Doacoes_Candidatos(data, cargo=cargo)
        else:
            donations = Doacoes_Candidatos(data)
        self.donations = donations
        numeros = list(set([ele.numero_cand for ele in donations.data]))
        totals = [sum([ele.valor for ele in donations.data if ele.numero_cand==num]) for num in numeros]
        dic_nomes = {}
        nomes = ["" for ii in xrange(len(numeros))]
        cpfs = ["" for ii in xrange(len(numeros))]
        for ele in donations.data:
            num = ele.numero_cand
            nome = ele.nome_cand
            cpf = ele.CPF_cand
            if nome not in nomes:
                ii = numeros.index(num)
                nomes[ii] = nome
                cpfs[ii] = cpf
                dic_nomes[num] = nome
            #else:
            #    if nome != dic_nomes[num]:
            #        print "SOMETHNG IS WRONG 1!!!!!!!!!!!!!", nome, dic_nomes[num]
        # reads votes and creates data structure
        bla = open(csv).read()
        bla = bla.replace(',', ".")
        bla = bla.split("\n")
        data = [[ele2.replace('"', '') for ele2 in ele.split('";"')] for ele in bla]
        data.pop(0)
        data.pop()
        data.pop()
        data.pop()
        data.pop()
        data.pop()
        struc = []
        for ele in data:
            num = ele[2]
            name = ele[3]
            party = ele[4]
            coalition = ele[5]
            situation = ele[6]
            votes = int(ele[7].replace(".", ""))
            percentage = ele[8].replace(",", ".")
            if num in numeros:
                ii = numeros.index(num)
                cpf = cpfs[ii]
                total = totals[ii]
                struc.append(Candidato(name, num, cpf, cargo, total, votes, party, coalition, situation, percentage))
                # check
                nome = nomes[ii].replace(" ", "").replace("#", "").upper()
                namer = name.replace(" ", "").replace("#", "").upper()
                #if nome!=namer:
                #    print "SOMETHNG IS WRONG 2!!!!!!!!!!!!!", nome, namer, dic_nomes[num]
            else:
                cpf = "?"
                total = 0.
                struc.append(Candidato(name, num, cpf, cargo, total, votes, party, coalition, situation, percentage))
                dic_nomes[num] = name
                numeros.append(num)
                totals.append(0.)
                nomes.append(name)
                cpfs.append(cpf)
        self.data = struc[:]
        self.numbers = numeros[:]
        self.names = nomes[:]
        self.moneys = totals[:]
        self.CPFs = cpfs[:]
        self.dic_names = dic_nomes
                







