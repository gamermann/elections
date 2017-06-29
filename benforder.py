#!/usr/bin/python
#encoding=latin-1

#
# This script generates tables, fits and plots for testing the benford 
# distribution agains the donnations made to candidates, parties, etc...
#

import os
from Eleicoes import *
from ExGUtils.uts import int_points_gauss, intsum, histogram # For integrating 


# Directories for output
plot_dir = "plots/"
tex_dir = "texs/"
table_dir = "pdfs/"


####################################
# Reads
####################################
UF = "BR"
fil_cands = "../Eleicoes/prestacao_final_2014/receitas_candidatos_2014_%s.txt"%UF # file from TSE with donations to candidates
fil_parts = "../Eleicoes/prestacao_final_2014/receitas_partidos_2014_%s.txt"%UF # file from TSE with donations to parties
fil_comts = "../Eleicoes/prestacao_final_2014/receitas_comites_2014_%s.txt"%UF # file from TSE with donations to committes
tttti = time()
ti = tttti
# Removing declarations from parties, for Br, One is left with presidentional campaign only. 
eleic = Eleicao_UF(candidatos=fil_cands, partidos="", comites=fil_comts, cargo="") 
tf = time()
print "Donations read in %f s"%(tf-ti)


alldons = [ele for ele in eleic.candidatos.data]
alldons += [ele for ele in eleic.partidos.data if ele.is_ori]
alldons += [ele for ele in eleic.comites.data]


parties = list(set([ele.partido for ele in alldons]))

#
# Analysis done for the following sets of data:
# original Donations from CNPJ
# original donnations From CPF
# Non-original donnations
# Only for presidentional campaign 
#


# number of each kind of donnation per party:
print table_start%("Donations per party", "tab:dpp", "c|cccc|c")
print hline
print "Party & CNPJ & CPF & Non-ori & Unknown & Total \\\\"
for party in parties:
    dons = [ele for ele in alldons if ele.partido==party]
    cnpjs = [ele for ele in dons if ele.is_ori and ele.is_CNPJ]
    cpfs = [ele for ele in dons if ele.is_ori and ele.is_CPF]
    nonori = [ele for ele in dons if not ele.is_ori]
    unknown = [ele for ele in dons if ele.is_ori and not ele.is_CNPJ and not ele.is_CPF]
    #
    tot = len(dons)
    ncnpj = len(cnpjs)
    ncpf = len(cpfs)
    nnori = len(nonori)
    nunk = len(unknown)
    print "%8s & %4i  & %4i  & %4i  & %4i  & %4i \\\\"%(party, ncnpj, ncpf, nnori, nunk, tot)

print table_end

###############
### Plot stuff
g("set term post eps enhanced color font 'Helvetica, 25'")


# stuff
ii = 1
base = exp(1.) # Use base=10. to change to logarithm base 10 in the plots
l10 = log(base) 
delt = log(.01)



ttti = time()
##################
# Prepare big Table (tex stuff)
##################
lines = []
lines.append(table_start_tiny%("Benfor - All ", "tab%i"%ii, "c|ccccccccc|ccc|ccccc"))
lines.append("Partido & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & N & $\\chi^2$ & p-val & Min & Max & Sum & $\\gamma$ & $\\xi_0$\\\\")
lines.append(hline)
lines.append(hline)
stri = " & ".join(["%4.3f"%ele for ele in bf])
lines.append("Benford & %s & & & & &\\\\"%stri)
lines.append(hline)
##################
# Starts
##################
for party in parties:
    print "________________________"
    print "Party: %s:"%party
    all_nums = []
    #
    # Fit All
    #
    print "All"
    vals = [ele.valor for ele in alldons if ele.partido==party]
    N = len(vals)
    Nt = N
    if N>20:
        csim = log(max(vals)+1.)
        ti = time()
        gamma, csi0, csim = maxLKHD_distr(vals, csim=csim, delt=delt, eps=1.e-8)
        tf = time()
        gammat = gamma
        csi0t = csi0
        csimt = csim
        lkhd, gg = lkhd_distr(vals, gamma, csi0, csim, delt)
        print "Fit time: %f s"%(tf-ti)
        print "  gamma=%1.10f x0=%8.2f  lkhd=%f"%(gamma, round(exp(csi0), 2), lkhd)
        # Plots
        [x, y] = histogram([log(ele)/l10 for ele in vals])
        dx = .5*(x[1]-x[0])
        dd = gdata(x, y, with_="boxes", title="Data") # data distribution with x in log scale
        y2 = []
        for xp in x:
            xi = base**(xp-dx)
            xf = base**(xp+dx)
            xss = int_points_gauss(xi, xf, 200)
            yss = [func2(ele, gamma, csi0, csim, delt) for ele in xss]
            y2.append(N*intsum(xi, xf, yss))
        dt = gdata(x, y2, with_="lines lw 3", title="Fit")
        g("set output '%shist_%s_all.eps'"%(plot_dir, party))
        g("set title '%s - Distribution for all Donations'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Counts'")
        g("set key")
        g.plot(dd, dt)
        # Chi2
        chi2 = sum([(ele-y[ii])**2/ele for ii, ele in enumerate(y2)])
        p_value = (1. - chi2_distr.cdf(chi2, len(x)-1))
        if p_value>0.: print "%s: p_val=%f"%(party, p_value)
        # Cumulative
        x, y = acumulator(vals)
        y2 = [func(ele, gamma, csi0, csim, delt) for ele in x]
        x = [log(ele)/l10 for ele in x]
        dd = gdata(x, y, title="Data")
        dt = gdata(x, y2, with_="lines lw 4", title="Fit")
        g("set output '%shist_accu_%s_all.eps'"%(plot_dir, party))
        g("set title '%s - All Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd, dt)
        g("set output '%shist_accu_%s_all_nf.eps'"%(plot_dir, party))
        g("set title '%s - All Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd)
        # Benford stuff
        line, dbenf = benford_stuff(vals, "%s - All"%party, True)
        g("set output '%sbenford_%s_all.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for all donnations'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        # random
        nums = [drand_distr(gamma, csi0, csim, delt) for ii in xrange(N)]
        line, dbenf = benford_stuff(nums, "%s - All Rand"%party, True)
        g("set output '%sbenford_%s_all_rand.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for all donnations random generated'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
    #
    # Fit CNPJ
    #
    print "CNPJ"
    vals = [ele.valor for ele in alldons if ele.partido==party and ele.is_ori and ele.is_CNPJ]
    N = len(vals)
    if N>20:
        csim = log(max(vals)+1.)
        ti = time()
        gamma, csi0, csim = maxLKHD_distr(vals, csim=csim, delt=delt, eps=1.e-8)
        tf = time()
        lkhd, gg = lkhd_distr(vals, gamma, csi0, csim, delt)
        print "Fit time: %f s"%(tf-ti)
        print "  gamma=%1.10f x0=%8.2f  lkhd=%f"%(gamma, round(exp(csi0), 2), lkhd)
        # Plots
        [x, y] = histogram([log(ele)/l10 for ele in vals])
        dx = .5*(x[1]-x[0])
        dd = gdata(x, y, with_="boxes", title="Data") # data distribution with x in log scale
        y2 = []
        for xp in x:
            xi = base**(xp-dx)
            xf = base**(xp+dx)
            xss = int_points_gauss(xi, xf, 200)
            yss = [func2(ele, gamma, csi0, csim, delt) for ele in xss]
            y2.append(N*intsum(xi, xf, yss))
        dt = gdata(x, y2, with_="lines lw 3", title="Fit")
        g("set output '%shist_%s_cnpj.eps'"%(plot_dir, party))
        g("set title '%s - Distribution for donations from CNPJ'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Counts'")
        g("set key")
        g.plot(dd, dt)
        # Chi2
        chi2 = sum([(ele-y[ii])**2/ele for ii, ele in enumerate(y2)])
        p_value = (1. - chi2_distr.cdf(chi2, len(x)-1))
        if p_value>0.: print "%s: p_val=%f"%(party, p_value)
        # Cumulative
        x, y = acumulator(vals)
        y2 = [func(ele, gamma, csi0, csim, delt) for ele in x]
        x = [log(ele)/l10 for ele in x]
        dd = gdata(x, y, title="Data")
        dt = gdata(x, y2, with_="lines lw 4", title="Fit")
        g("set output '%shist_accu_%s_cnpj.eps'"%(plot_dir, party))
        g("set title '%s - CNPJ Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd, dt)
        g("set output '%shist_accu_%s_cnpj_nf.eps'"%(plot_dir, party))
        g("set title '%s - CNPJ Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd)
        # Benford stuff
        line, dbenf = benford_stuff(vals, "%s - CNPJ"%party, True)
        g("set output '%sbenford_%s_cnpj.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for CNPJ donnations'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        # random
        nums = [drand_distr(gamma, csi0, csim, delt) for ii in xrange(N)]
        line, dbenf = benford_stuff(nums, "%s - CNPJ Rand"%party, True)
        g("set output '%sbenford_%s_cnpj_rand.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for CNPJ donnations randomly generated'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        all_nums += nums
    #
    # Fit CPF
    #
    print "CPF"
    vals = [ele.valor for ele in alldons if ele.partido==party and ele.is_ori and ele.is_CPF]
    N = len(vals)
    if N>20:
        csim = log(max(vals)+1.)
        ti = time()
        gamma, csi0, csim = maxLKHD_distr(vals, csim=csim, delt=delt, eps=1.e-8)
        tf = time()
        lkhd, gg = lkhd_distr(vals, gamma, csi0, csim, delt)
        print "Fit time: %f s"%(tf-ti)
        print "  gamma=%1.10f x0=%8.2f  lkhd=%f"%(gamma, round(exp(csi0), 2), lkhd)
        # Plots
        [x, y] = histogram([log(ele)/l10 for ele in vals])
        dx = .5*(x[1]-x[0])
        dd = gdata(x, y, with_="boxes", title="Data") # data distribution with x in log scale
        y2 = []
        for xp in x:
            xi = base**(xp-dx)
            xf = base**(xp+dx)
            xss = int_points_gauss(xi, xf, 200)
            yss = [func2(ele, gamma, csi0, csim, delt) for ele in xss]
            y2.append(N*intsum(xi, xf, yss))
        dt = gdata(x, y2, with_="lines lw 3", title="Fit")
        g("set output '%shist_%s_cpf.eps'"%(plot_dir, party))
        g("set title '%s - Distribution for donations from CPF'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Counts'")
        g("set key")
        g.plot(dd, dt)
        # Chi2
        chi2 = sum([(ele-y[ii])**2/ele for ii, ele in enumerate(y2)])
        p_value = (1. - chi2_distr.cdf(chi2, len(x)-1))
        if p_value>0.: print "%s: p_val=%f"%(party, p_value)
        # Cumulative
        x, y = acumulator(vals)
        y2 = [func(ele, gamma, csi0, csim, delt) for ele in x]
        x = [log(ele)/l10 for ele in x]
        dd = gdata(x, y, title="Data")
        dt = gdata(x, y2, with_="lines lw 4", title="Fit")
        g("set output '%shist_accu_%s_cpf.eps'"%(plot_dir, party))
        g("set title '%s - CPF Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd, dt)
        g("set output '%shist_accu_%s_cpf_nf.eps'"%(plot_dir, party))
        g("set title '%s - CPF Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd)
        # Benford stuff
        line, dbenf = benford_stuff(vals, "%s - CPF"%party, True)
        g("set output '%sbenford_%s_cpf.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for CPF donnations'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        # random
        nums = [drand_distr(gamma, csi0, csim, delt) for ii in xrange(N)]
        line, dbenf = benford_stuff(nums, "%s - CPF Rand"%party, True)
        g("set output '%sbenford_%s_cpf_rand.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for CPF donnations randomly generated'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        all_nums += nums
    #
    # Fit nonori
    #
    print "Non-ori"
    vals = [ele.valor for ele in alldons if ele.partido==party and not ele.is_ori]
    N = len(vals)
    if N>20:
        csim = log(max(vals)+1.)
        ti = time()
        gamma, csi0, csim = maxLKHD_distr(vals, csim=csim, delt=delt, eps=1.e-8)
        tf = time()
        lkhd, gg = lkhd_distr(vals, gamma, csi0, csim, delt)
        print "Fit time: %f s"%(tf-ti)
        print "  gamma=%1.10f x0=%8.2f  lkhd=%f"%(gamma, round(exp(csi0), 2), lkhd)
        # Plots
        [x, y] = histogram([log(ele)/l10 for ele in vals])
        dx = .5*(x[1]-x[0])
        dd = gdata(x, y, with_="boxes", title="Data") # data distribution with x in log scale
        y2 = []
        for xp in x:
            xi = base**(xp-dx)
            xf = base**(xp+dx)
            xss = int_points_gauss(xi, xf, 200)
            yss = [func2(ele, gamma, csi0, csim, delt) for ele in xss]
            y2.append(N*intsum(xi, xf, yss))
        dt = gdata(x, y2, with_="lines lw 3", title="Fit")
        g("set output '%shist_%s_nonori.eps'"%(plot_dir, party))
        g("set title '%s - Distribution for donations from non-original'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Counts'")
        g("set key")
        g.plot(dd, dt)
        # Chi2
        chi2 = sum([(ele-y[ii])**2/ele for ii, ele in enumerate(y2)])
        p_value = (1. - chi2_distr.cdf(chi2, len(x)-1))
        if p_value>0.: print "%s: p_val=%f"%(party, p_value)
        # Cumulative
        x, y = acumulator(vals)
        y2 = [func(ele, gamma, csi0, csim, delt) for ele in x]
        x = [log(ele)/l10 for ele in x]
        dd = gdata(x, y, title="Data")
        dt = gdata(x, y2, with_="lines lw 4", title="Fit")
        g("set output '%shist_accu_%s_nonori.eps'"%(plot_dir, party))
        g("set title '%s - Non-original Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd, dt)
        g("set output '%shist_accu_%s_nonori_nf.eps'"%(plot_dir, party))
        g("set title '%s - Non-original Donations - Cumulative Distribution'"%party)
        g("set xlab 'ln(Amount)'")
        g("set ylab 'Sum(Counts)'")
        g("set key")
        g.plot(dd)
        # Benford stuff
        line, dbenf = benford_stuff(vals, "%s - Non-original"%party, True)
        g("set output '%sbenford_%s_nonori.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for Non-original donnations'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        # random
        nums = [drand_distr(gamma, csi0, csim, delt) for ii in xrange(N)]
        line, dbenf = benford_stuff(nums, "%s - Non-original Rand"%party, True)
        g("set output '%sbenford_%s_nonori_rand.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for Non-Original donnations randomly generated'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & %2.4f & %3.4f \\\\"%(gamma, csi0))
        all_nums += nums
    #
    # Fit unknown
    #
    print "Unknown"
    vals = [ele.valor for ele in alldons if ele.partido==party and ele.is_ori and not ele.is_CNPJ and not ele.is_CPF]
    N = len(vals)
    if N>20:
        csim = log(max(vals)+1.)
        ti = time()
        gamma, csi0, csim = maxLKHD_distr(vals, csim=csim, delt=delt, eps=1.e-8)
        tf = time()
        lkhd, gg = lkhd_distr(vals, gamma, csi0, csim, delt)
        print "Fit time: %f s"%(tf-ti)
        print "  gamma=%1.10f x0=%8.2f  lkhd=%f"%(gamma, round(exp(csi0), 2), lkhd)
        nums = [drand_distr(gamma, csi0, csim, delt) for ii in xrange(N)]
        all_nums += nums
    #
    # Last model
    #
    if all_nums:
        if len(all_nums)!=Nt:
            all_nums += [drand_distr(gammat, csi0t, csimt) for ii in xrange(Nt-len(all_nums))]
        line, dbenf = benford_stuff(all_nums, "%s - Model"%party, True)
        g("set output '%sbenford_%s_model.eps'"%(plot_dir, party))
        g("set title '%s - First Digit Distribution for modeled donnations'"%party)
        g("set xlab 'Digit'")
        g("set ylab 'frequency'")
        g("unset key")
        g.plot(dbf, dbenf)
        lines.append(line + " & - & - \\\\")
    #
    # Finishes
    #
    if Nt>20:
        lines.append(hline)
    


#
#
# Finalization
#
#
lines.append(table_end)
tab = [ele.decode("utf-8").strip() for ele in open("%sbegin.tex"%tex_dir).readlines()]
tab += lines
tab += ["\n", "\n","\n", "\n", "\\end{document}"]
tab = [ele.replace("_", "\_").encode("utf-8") for ele in tab]
print >>open("%sfinal_benford_table_%s.tex"%(tex_dir, UF), "w"), "\n".join(tab)
os.system("pdflatex -output-directory=%s %sfinal_benford_table_%s.tex"%(table_dir, tex_dir, UF))

#

tttf = time()
print
print
print "Fitting time: %f s"%(tttf-ttti)
print
print
print "Total run time: %f s"%(tttf-tttti)









g("set term x11")

