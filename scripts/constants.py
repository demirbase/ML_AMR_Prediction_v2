#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared constants for the AMR Prediction pipeline.

Single source of truth for values that were previously copy-pasted across
multiple scripts (resolves the duplicated ANTIBIOTIC_CLASSES dictionary that
lived in both 01_data_validation.py and 01b_data_validation.py).

Import as:
    from constants import ANTIBIOTIC_CLASSES

This works because each pipeline script is launched directly
(``python scripts/0X_*.py``), which puts the ``scripts/`` directory on
``sys.path``.
"""

# ============================================================================
# ANTIBIOTIC CLASSIFICATION SYSTEM
# ============================================================================
# Maps each antibiotic class to its member drugs. Used for grouped reporting
# (01_data_validation.py) and class-representation plots (01b_data_validation.py).
ANTIBIOTIC_CLASSES = {
    'Penicillins': ['ampicillin', 'amoxicillin', 'amoxicillin/clavulanic acid', 'piperacillin/tazobactam', 'ampicillin/sulbactam', 'penicillin', 'carbenicillin', 'piperacillin', 'ticarcillin/clavulanic acid'],
    'Cephalosporins': ['ceftazidime', 'cefotaxime', 'cefuroxime', 'ceftriaxone', 'cefepime', 'cefoxitin', 'cephalothin', 'cefazolin', 'ceftiofur', 'cefpodoxime', 'cefotetan', 'ceftazidime/avibactam', 'ceftaroline', 'cephalexin', 'cefpodoxime_clavulanic_acid', 'ceftolozane/tazobactam', 'cefotaxime/clavulanic acid'],
    'Beta-Lactams: Carbapenems & Others': ['meropenem', 'imipenem', 'ertapenem', 'doripenem', 'aztreonam', 'beta-lactam', 'sulbactam'],
    'Aminoglycosides': ['gentamicin', 'amikacin', 'tobramycin', 'streptomycin', 'kanamycin', 'apramycin', 'neomycin', 'netilmicin'],
    'Quinolones': ['ciprofloxacin', 'norfloxacin', 'levofloxacin', 'nalidixic acid', 'moxifloxacin', 'ofloxacin'],
    'Folate Pathway Inhibitors': ['trimethoprim/sulfamethoxazole', 'trimethoprim', 'sulfamethoxazole', 'sulfisoxazole'],
    'Tetracyclines': ['tigecycline', 'tetracycline', 'doxycycline', 'minocycline', 'oxytetracycline'],
    'Others': ['chloramphenicol', 'nitrofurantoin', 'azithromycin', 'colistin', 'fosfomycin', 'erythromycin', 'lincomycin', 'rifampin', 'clindamycin', 'clarithromycin', 'daptomycin', 'linezolid', 'polymyxin B', 'teicoplanin', 'vancomycin']
}
