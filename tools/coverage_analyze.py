#!/usr/bin/env python

import os, subprocess, shutil, glob

# Create a clean working dir
if os.path.exists('gcov'):
    shutil.rmtree('gcov')
os.mkdir('gcov')
os.chdir('gcov')

# Run gcov on all relevant files
for root, dirs, fns in os.walk('../'):
    for fn in fns:
        if fn.endswith('.gcno'):
             fn_o = '%s/%s.o' % (root, fn[:-5])
             #print 'Coverage %s' % fn_o
             subprocess.call('gcov %s > /dev/null' % fn_o, shell=True)

# Sometimes gcov signals non-executed lines that are not relevant
def is_relevant(line):
    #print line[16:].strip()
    if line[16:].strip().startswith('class'):
        return False
    if line[16:].split() == ['do', '{']:
        return False
    if 'IGNORE_COVERAGE' in line:
        return False
    return True

# Loop over all gcov files and extract relevant information
coverage_results = {}
for fn in os.listdir('.'):
    info_lines = []
    path = ''
    if fn.endswith('.gcov'):
        with open(fn) as f:
            first = f.next()
            if '/usr' in first:
                continue
            path = first[23:-1]
            for line in f:
                if line.startswith('    #####:') and is_relevant(line):
                    info_lines.append(line[:-1])
    dn, bn = os.path.split(path)
    l = coverage_results.setdefault(dn, [])
    l.append((bn, info_lines, fn))

# Nicely format things on screen
GREEN = '\033[92m'
ENDC = '\033[0m'
RED = '\033[91m'
print 'Number of lines not covered in unit tests:'
for dn, bns in sorted(coverage_results.iteritems()):
    print '   ', dn
    for bn, info_lines, fn_gcov in sorted(bns):
        if len(info_lines) == 0:
            print '%8i  %s' % (len(info_lines), bn)
        else:
            print '%s%8i  %s gcov/%s%s' % (RED, len(info_lines), bn.ljust(30), fn_gcov, ENDC)
