import re
import sys
import os

def check_red_flags(path, text):
    red_flags = [
        (r"\b0\.0716\b", "0.0716 (superseded lambda_min)"),
        (r"threshold\s*=\s*0\.0239|threshold\s*=\s*\+?0\.0239", "0.0239 as development threshold"),
        (r"###\s*Prediction Error and Calibration", "prediction error and calibration heading"),
        (r"IBS (?:proves|confirms|establishes) calibration", "calibration based only on IBS/Brier"),
        (r"primary prognostic driver", "primary prognostic driver"),
        (r"orthogonal latent factor", "orthogonal latent factors"),
        (r"\bfirst\b|\bfirst-ever\b|\bfirst in HCC\b", "first / first-ever / first in HCC"),
        (r"ElasticNet[- ]Cox was the best[- ]performing|best[- ]performing ElasticNet", "best-performing ElasticNet"),
        (r"external validation for the (?:N\s*=\s*65\s*)?holdout", "external validation for the N=65 holdout"),
        (r"\b0\.6696\b", "leaked holdout 0.6696"),
        (r"drives survival", "drives survival (causal claim)"),
        (r"proves mechanism", "proves mechanism"),
        (r"statistically significant superiority based on repeated-fold tests", "statistically significant superiority based on repeated-fold tests"),
        (r"\bbreakthrough\b", "breakthrough"),
        (r"\bstate-of-the-art\b", "state-of-the-art"),
        (r"\bclinically ready\b|\bready for the clinic\b|\bready for bedside\b", "clinically ready"),
        (r"four[- ]way (?:mechanistic )?synergy", "four-way synergy"),
        (r"mutations are (?:biologically )?unimportant", "mutations are unimportant")
    ]

    violations = []
    for pattern, name in red_flags:
        for m in re.finditer(pattern, text, re.IGNORECASE):
            start = max(0, m.start() - 50)
            end = min(len(text), m.end() + 50)
            snippet = text[start:end].replace("\n", " ")
            lower_snippet = snippet.lower()
            if any(w in lower_snippet for w in ["never", "prohibited", "disallowed", "cannot", "not", "audit", "banned", "none", "zero", "check:", "reconciliation", "resolved", "corrected", "superseded"]):
                continue
            violations.append((name, snippet))

    if violations:
        print(f"FAILED: Found {len(violations)} red-flag violations in {path}:")
        for v, snip in violations:
            print(f"  - [{v}]: ...{snip}...")
        return False
    else:
        print(f"PASSED: Zero red-flag violations in {path}")
        return True

def check_numerical_consistency(path, text):
    required_numbers = [
        ("285", "Development cohort N=285"),
        ("103", "Development deaths=103"),
        ("182", "Development censored=182"),
        ("65", "Holdout N=65"),
        ("19", "Holdout deaths=19"),
        ("46", "Holdout censored=46"),
        ("350", "Primary locked N=350"),
        ("377", "Initial universe N=377"),
        ("11,402", "Total features=11,402"),
        ("0.5875", "ElasticNet CV C-index=0.5875"),
        ("0.5910", "LASSO CV C-index=0.5910"),
        ("0.6197", "Holdout C-index=0.6197"),
        ("0.5974", "Clinical-only C-index=0.5974"),
        ("0.6120", "Combined C-index=0.6120"),
        ("0.0146", "Delta C=0.0146"),
        ("0.0214", "Pooled Delta C=0.0214"),
        ("0.2717", "Paired t-test p=0.2717"),
        ("0.2584", "Paired Wilcoxon p=0.2584"),
        ("0.7751", "Factor 2 HR=0.7751"),
        ("0.000224", "Factor 2 p=0.000224"),
        ("0.1665", "Factor 2 beta=-0.1665"),
        ("0.9784", "Factor 2 stability r=0.9784"),
        ("0.0892", "SNV variance=0.0892%"),
        ("66.27", "CNV variance=66.27%"),
        ("58.61", "Meth variance=58.61%"),
        ("27.55", "RNA variance=27.55%"),
        ("0.1928", "Holdout IBS=0.1928"),
        ("0.0859", "Holdout 1y Brier=0.0859"),
        ("0.2180", "Holdout 3y Brier=0.2180"),
        ("0.2744", "Holdout 5y Brier=0.2744"),
        ("0.08055", "Production lambda_min=0.08055"),
        ("-0.0233", "Development median linear predictor threshold=-0.0233")
    ]

    missing_nums = []
    for num, desc in required_numbers:
        if num not in text:
            missing_nums.append((num, desc))

    if missing_nums:
        print(f"FAILED: Missing {len(missing_nums)} required numbers in {path}:")
        for num, desc in missing_nums:
            print(f"  - Missing {num}: {desc}")
        return False
    else:
        print(f"PASSED: All {len(required_numbers)} required study metrics present in {path}")
        return True

master_path = "./docs/MANUSCRIPT_MASTER.md"
print(f"==================================================")
print(f"Auditing Master Manuscript: {master_path}")
print(f"==================================================")
with open(master_path, "r", encoding="utf-8") as f:
    master_text = f.read()

master_rf_ok = check_red_flags(master_path, master_text)
master_num_ok = check_numerical_consistency(master_path, master_text)

sections_dir = "./docs/manuscript_sections"
print(f"\n==================================================")
print(f"Auditing Individual Manuscript Sections for Red Flags")
print(f"==================================================")
sections_rf_ok = True
for fname in sorted(os.listdir(sections_dir)):
    if fname.endswith(".md"):
        fpath = os.path.join(sections_dir, fname)
        with open(fpath, "r", encoding="utf-8") as f:
            sec_text = f.read()
        if not check_red_flags(fpath, sec_text):
            sections_rf_ok = False

if master_rf_ok and master_num_ok and sections_rf_ok:
    print("\n==================================================")
    print("ALL RED-TEAM VERIFICATIONS AND AUDIT CHECKS PASSED!")
    print("==================================================")
    sys.exit(0)
else:
    print("\n==================================================")
    print("AUDIT FAILED: Please review reported issues above.")
    print("==================================================")
    sys.exit(1)
