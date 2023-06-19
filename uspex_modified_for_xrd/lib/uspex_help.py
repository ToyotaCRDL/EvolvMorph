#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
USPEX 9.4.4 release
2015 Oganov's Lab. All rights reserved.
"""

import os
import sys


# -------------------------------------------------------------------------------
def uspex_help(parm):
    if not parm:
        return

    run_dir = os.getcwd()
    bin_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(bin_dir)
    os.chdir('doc')
    doc_dir = os.getcwd()
    tex_file = doc_dir + '/uspex_keywords.xml'
    # tex_file = 'uspex_keywords.xml'
    os.chdir(run_dir)

    f = open(tex_file, 'rb')
    full_tex_content = f.readlines()
    f.close()

    tex_content = []
    for row in full_tex_content:
#        if row.strip().replace(' ', '') <> '':
        if row.strip().replace(' ', '') != '':   # modified on 210414
            tex_content.append(row. \
                               replace(r'{max\\_coordination\\_number}', '/max_coordination_number'). \
                               replace(r'{\\% ', '{\n% '). \
                               replace(r'\\\\', ''). \
                               replace(r'{\\rm ', '{'). \
                               replace(r'\\rm ', ''). \
                               replace(r'\\rm', ''). \
                               replace('`', "'"). \
                               replace("''", "'"). \
                               replace('[]', '').
                               replace(r'\\begin{itemize}', ''). \
                               replace(r'\\end{itemize}', ''). \
                               replace(r'\\begin{center}', ''). \
                               replace(r'\\begin{tabular}{|l|c|l|}', ''). \
                               replace(r'\\begin{figure}[h]', ''). \
                               replace('\centering', ''). \
                               replace('\includegraphics[scale=0.18]{pic/hardness_example}', ''). \
                               replace(
                '\caption{\footnotesize \textbf{Predictions of the hardest structure of TiO$_2$.}}', ''). \
                               replace('\label{fig:hardness_example}', ''). \
                               replace(r'\\end{figure}', ''). \
                               replace(r'\\vspace{0.5cm}', ''). \
                               replace(r'\\end{tabular}', ''). \
                               replace(r'\\end{center}', ''). \
                               replace('\hline', ''). \
                               replace(r'\\textbf{', ''). \
                               replace(r'\\textshift{', ''). \
                               replace('---', '-'). \
                               replace('--', '-'). \
                               replace(r'\\item', '-'). \
                               replace('\emph{', ''). \
                               replace('\keyword{', ''). \
                               replace(r'\\cite{Chen2011}', ''). \
                               replace('}', ''). \
                               replace('\paramacro{', ''). \
                               replace('$', ''). \
                               replace(r'\\times', '*'). \
                               replace(r'\\frac{', ''). \
                               replace(r'\\_', '_'). \
                               replace(r'\\sim', '~'). \
                               replace(r'\\r{', ''). \
                               replace(r'\\%', '%'). \
                               replace(r'\\v{n', 'n'). \
                               replace(r'\\\'{a', 'a'). \
                               replace(r'\\ldots', '...'). \
                               replace(r'\\file{', ''). \
                               replace(r'\\textshiftleft{', ''). \
                               replace(r'\\keyword{', ''). \
                               replace('\geq', '>='). \
                               replace(r'\\texttt', ''). \
                               replace(r'~\\ref{fig:', ''). \
                               replace('\cite{Oganov:2010', ' (Oganov:2010)'). \
                               replace('\cite{Dubrovinsky:2001', ' (Dubrovinsky:2001)'). \
                               replace('\Delta', 'Delta'). \
                               replace('~', ''). \
                               replace(r'\\ref{appendix17', ''). \
                               replace('\includegraphics[scale=0.3]{pic/Wyckoff_positions', ''). \
                               replace('\caption{\\footnotesize Example of merging atoms onto special Wyckoff', ''). \
                               replace('positions (from Ref.\cite{Lyakhov:2013aa).', ''). \
                               replace('\label{fig:Wyckoff_positions', ''). \
                               replace('\label{fracGene', ''). \
                               replace('\label{fracRand', ''). \
                               replace('\min', 'min'). \
                               replace(r'\\alpha', 'alpha'). \
                               replace(r'\\beta', 'beta'). \
                               replace('\gamma', 'gamma'). \
                               replace(r'\\textcolor{blue{\url', ''). \
                               replace('2\piA^{-1', '2*pi*A^(-1)'). \
                               replace('\leq', '<='). \
                               replace('\delta', 'delta'). \
                               replace(r'\\omega', 'omega'). \
                               replace(r'\\ref{eq:pso2', '8, see online manual'). \
                               replace(r'\\var', ''). \
                               replace(r'\\rightarrow', '->'). \
                               replace('\&', '&'). \
                               replace('{', ''). \
                               replace('\caption\\footnotesize Predictions of the hardest structure of TiO_2.', '')
                               )   # added r'' because of unicode error of python3 on 210414

    keywords_list = []
    for i, row in enumerate(tex_content):
#        if row.find('<tocitem target="name">') >= 0 and row[0] <> '%':
        if row.find('<tocitem target="name">') >= 0 and row[0] != '%':   # modified on 210414
            current_keyword = row.split('<tocitem target="name">')[1].split('</tocitem>')[0].replace('\r\n', '')
            keywords_list.append(current_keyword)
    keywords_list.sort()

#    if parm <> 'all':
    if parm != 'all':   # modified on 210414
        if parm in keywords_list:
            keywords_list = [parm]
        else:
            keywords_list = []
            print ('Keyword does not exist!')   # modified on 210414
            return

#    if parm <> 'all':
    if parm != 'all':   # modified on 210414
        keyword = keywords_list[0]
        brackets_count = 0
        for i, row in enumerate(tex_content):
            if row.find('<tocitem target="name">' + keyword) >= 0:
                keyword_start = i
                brackets_count = tex_content[keyword_start].count('<tocitem target=')
                j = 1
                while brackets_count <= 6:
                    count = tex_content[keyword_start + j].count('<tocitem target=')
                    brackets_count = brackets_count + count
                    j += 1
                keyword_end = keyword_start + j
                break

        keyword_content = ''.join(tex_content[keyword_start:keyword_end])

        print ('Keyword:', keyword_content.split('<tocitem target="name">')[1].split('</tocitem>')[0])   # modified on 210414
        print ('\nDescription:')   # modified on 210414
        print ('------------')   # modified on 210414
        print (keyword_content.split('<tocitem target="description">')[1].split('</tocitem>')[0])   # modified on 210414
        print (keyword_content.split('<tocitem target="possible values">')[1].split('</tocitem>')[0].replace('&', '|'))   # modified on 210414
        print ('\nDefault:', keyword_content.split('<tocitem target="default">')[1].split('</tocitem>')[0])   # modified on 210414
        print ('')   # modified on 210414
        syntax = keyword_content.split('<tocitem target="example">')[1].split('</tocitem>')[0].replace('%\r\n', '% ')
        if syntax.find('%') < 0:
            syntax = syntax.replace('\r\n', ' ')
        print ('\nSyntax:\n', syntax)   # modified on 210414
        notes = keyword_content.split('<tocitem target="notes">')[1].split('</tocitem>')[0]
#        if notes <> '':
        if notes != '':   # modified on 210414
            print ('\n', notes)   # modified on 210414
#        print '\nCalculation method:', \
#            keyword_content.split('<tocitem target="calculationMethod">')[1].split('</tocitem>')[0]
        print ('\nCalculation method:', keyword_content.split('<tocitem target="calculationMethod">')[1].split('</tocitem>')[0])   # modified on 210414
    else:
        print ('')   # modified on 210414
        for key in keywords_list:
            print ('\t' + key + '\n')   # modified on 210414

    if len(keywords_list) > 1:
        print ('\nNumber of keywords:', len(keywords_list))   # modified on 210414


# -------------------------------------------------------------------------------

# Usage:
# uspex_help('trajectoryFile')
# uspex_help('MDrestartFile')

if __name__ == "__main__":
    from optparse import OptionParser

    usage = "usage: %prog OPTIONS"
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--parameter", dest="parm",
                      help="specify parameter to get help. If no value or 'all' value is specified, all INPUT.txt parameters will be shown",
                      metavar="PARM")

    (options, args) = parser.parse_args()

    if options.parm:
        parm = options.parm
        uspex_help(parm)
        sys.exit(0)
