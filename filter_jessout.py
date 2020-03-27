#!/usr/bin/env python3
"""Takes the output of Jess in stdin, and for each template-target pair
it keeps only the best aligment (the one with the lowest score)"""

import sys
import os
import re


def run():

    logE_values = dict()
    templates = dict()
    currentTemplate = []
    previous_target = None
    previous_logE = 999.0

    for line in sys.stdin: 
        if re.search("^REMARK", line):
            remarkFields = line.split()
            target = remarkFields[1]
            query = remarkFields[3]
            logE = float(remarkFields[7])
            currentTemplate.append(line)

        if re.search("^ATOM", line):
            currentTemplate.append(line)

        if re.search("^ENDMDL", line):
            currentTemplate.append(line)
            if previous_target == target and logE < previous_logE: #and query != target:
                templates[target] = currentTemplate.copy()
                currentTemplate.clear()
            else:
                currentTemplate.clear()
            previous_logE = logE

        previous_target = target

    for key, template in templates.items():
        for line in template:   
            print(line.strip())
        print()

if __name__ == "__main__":
    run()
