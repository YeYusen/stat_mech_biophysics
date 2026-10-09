"""Extract positive-well Cq values from the authors' ExtremeDiluteR.xlsx (github.com/Azuresky99/quPCR).

usage: python3 extract_exp.py ExtremeDiluteR.xlsx exp_cq.npy      (needs openpyxl)
"""
import sys
import numpy as np
import openpyxl

wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
rows = list(wb.worksheets[0].iter_rows(values_only=True))[1:]
cq = np.array([r[5] for r in rows if isinstance(r[5], (int, float))], float)
np.save(sys.argv[2], cq)
print(f"{len(cq)} positive wells, {np.sum((cq >= 38.07) & (cq <= 39.5))} single-copy in [38.07, 39.5]")
