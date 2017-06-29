#!/usr/bin/python
#encoding=latin-1

#
# script that
# Produces tables for the logistic regression fits
#


import os
from Eleicoes import *
from math import log
from Logit import *





##############################################################
# Logistic regression
######################################################
# filters
#uf = "RS"
ufs = ["RS", "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR", "SC", "SE", "SP", "TO"]
#cargo = "Senador"
cargo = "Deputado Federal"
#cargo = "Deputado Estadual"

# tags in the csv files
if cargo == "Deputado Estadual":
    cargo_fil = "dep_est"
elif cargo == "Deputado Federal":
    cargo_fil = "dep_fed"
elif cargo == "Senador":
    cargo_fil = "senador"


sit = "Eleito" # success means the string sit is in the candidate situation attribute
# note that in the csv files the difference between beeing elected or not apears as
# Eleito, Eleito por QP, Eleito por media or Nao eleito, so in this python script, the only difference 
# is the upper or lower case E in the "Eleito" or "nao eleito" string.


print table_start
print "\\small"
print "\\bt{c|cc|cccc|c}"
print "UF & $\\beta\\pm\\sigma$ & p-value (Wald) & $N$ & $n$ & Deviance & p-value & Total Money [R\\$]\\\\"
print hline
ndeps = 0
betas1 = []
betas0 = []
totmoneys =[]
errbetas1 = []
errbetas0 = []
for uf in ufs:
    donation_fil = "../Eleicoes/prestacao_final_2014/receitas_candidatos_2014_%s.txt"%uf # file with donnations to candidates in a given UF
    csv_fil = "../Eleicoes/resultados/%s_%s.csv"%(uf, cargo_fil) # CSV file with election results for a given office (cargo) in a given UF
    candidatos = Candidatos(donation_fil, csv_fil, cargo)
    candidatos = [ele for ele in candidatos.data]
    # preps
    moneys = [ele.money for ele in candidatos]
    eleitos = [int(sit in ele.situacao) for ele in candidatos]
    N = len(moneys)
    pops = [moneys[ii] for ii in xrange(N)]
    tot = sum(pops)
    pops = [moneys[ii]*1./tot for ii in xrange(N)]
    [xm, ss] = stats(pops)
    ys = [eleitos[ii] for ii in xrange(N)]
    ndeps += sum(ys)
    #
    try:
        betas = get_betas(pops, ys)
        #
        xx = pops[:]
        xx.sort()
        yy = [logit(betas, [ele]) for ele in xx]
        d1 = gdata(xx, yy, with_="lines")
        d2 = gdata(pops, ys, with_="points lc 0 pt 3")
        #g.plot(d1, d2)
        pva, incs = pvals(pops, betas)
        dev, pva2 = Deviance(pops, ys, betas)
        #print uf, betas, len(candidatos), N, pva, Deviance(pops, ys, betas)
        print "%s & %3.6f $\\pm$ %3.6f & %3.6f & %i & %i & %3.6f & %3.6f & %12.2f\\\\"%(uf, betas[0], incs[0]**.5, pva[0], N, sum(ys), dev, pva2, sum(moneys))
        print "   & %3.6f $\\pm$ %3.6f & %3.6f &    &    &       &       & \\\\"%(betas[1], incs[1]**.5, pva[1])
        #print cline%(2, 7)
        print hline
        betas1.append(betas[1])
        errbetas1.append(incs[1]**.5)
        betas0.append(betas[0])
        errbetas0.append(incs[0]**.5)
        totmoneys.append(sum(moneys))
    except:
        print "error for %s"%uf

print "\\et"
print table_end

print 

print ndeps

print


