# -*- coding: utf-8 -*-
"""
Created on: Fri Jun 26 21:45:28 2020

@author:    Mark
            Data Scientist
            Neurovine
            m.campbell@neurovine.ai
"""

bib = r'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\bib\thesisBib.bib'
f = open(bib, 'r')
entries = {}
entryString = ''
for line in f:
    if '@' in line and '{' in line:
        key = line[line.index('{')+1 : line.index(',')]
    entryString += line
    
    if '}' in line and line.index('}') == 0:
        entries[key] = entryString
        entryString = ''
        
# alphabetize dictionary
alphEntries = sorted(entries)
f = open(bib, 'w')
for a in alphEntries:
    f.write(entries[a])
f.close()

