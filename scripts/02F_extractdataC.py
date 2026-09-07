#!/usr/bin/env python3
"""
02F_extractdataC.py

Purpose:
  1. Derive location using both 'Location' and 'Food origin' columns:
     - Keep state detail if 'Location' has it (e.g., USA:CA).
     - Check 'Food origin' for state detail if 'Location' is strictly 'USA' or empty.
     - Fall back to 'USA' if both only contain 'USA' without state variables.
  2. Filter for records that remain strictly 'USA' (unresolved state).
  3. Save those records to data6c.csv.
  4. Inspect all remaining columns in data6c.csv for other potential state metadata.
"""

import os
import sys
import pandas as pd

# File Paths
DATA6A_PATH = "/home/samuelajulo/epi/data6a.csv"
DATA6C_PATH = "/home/samuelajulo/epi/data6c.csv"

print("=== SCRIPT 02F: EXTRACTING 'USA'-ONLY RECORDS (UPDATED LOCATION LOGIC) ===")

if not os.path.exists(DATA6A_PATH):
    print(f"[ERROR] Source file not found: {DATA6A_PATH}")
    sys.exit(1)

# Load data6a.csv
df = pd.read_csv(DATA6A_PATH, low_memory=False)

# Ensure required columns exist
for col in ["Location", "Food origin"]:
    if col not in df.columns:
        df[col] = None

def is_empty(val):
    if pd.isna(val):
        return True
    s = str(val).strip()
    return s in ["", "NA", "nan", "NaN", "None"]

def is_just_usa(val):
    if is_empty(val):
        return False
    s_upper = str(val).strip().upper()
    return s_upper in ["USA", "USA:", "USA :", "UNITED STATES", "US"]

def derive_location(row):
    loc = row.get("Location")
    food = row.get("Food origin")

    loc_clean = "" if is_empty(loc) else str(loc).strip()
    food_clean = "" if is_empty(food) else str(food).strip()

    # 1. If Location has detailed state info (not just USA), keep it
    if not is_empty(loc_clean) and not is_just_usa(loc_clean):
        return loc_clean

    # 2. If Location is just "USA" or empty, check Food origin for state info
    if not is_empty(food_clean) and not is_just_usa(food_clean):
        return food_clean

    # 3. If either column was "USA" (without state info), return "USA"
    if is_just_usa(loc_clean) or is_just_usa(food_clean):
        return "USA"

    # 4. Fallback if both are missing/empty
    return "NA"

# Apply updated location logic
df['Derived_Location'] = df.apply(derive_location, axis=1)

# Filter strictly for unresolved 'USA'
usa_df = df[df['Derived_Location'] == 'USA'].copy()

print(f"[INFO] Total records in data6a: {len(df)}")
print(f"[INFO] Found {len(usa_df)} records where derived Location remains strictly 'USA'.")

# Save to data6c.csv
usa_df.to_csv(DATA6C_PATH, index=False)
print(f"[SUCCESS] Saved data6c.csv to: {DATA6C_PATH}")

# Inspect all columns for hidden state/regional information in remaining USA records
print("\n========================================================")
print("   COLUMN METADATA SUMMARY FOR DERIVED LOCATION == 'USA'")
print("========================================================")

for col in usa_df.columns:
    if col == 'Derived_Location':
        continue
    valid_series = usa_df[col].dropna()
    valid_series = valid_series[~valid_series.astype(str).str.strip().isin(["", "NA", "nan", "NaN", "None"])]
    n_count = len(valid_series)
    unique_vals = valid_series.unique()
    
    if n_count > 0:
        sample_str = " | ".join([str(x) for x in unique_vals[:5]])
        print(f"\n[Column: {col}]")
        print(f"  - Non-empty rows : {n_count} / {len(usa_df)}")
        print(f"  - Unique values  : {len(unique_vals)}")
        print(f"  - Sample values  : {sample_str}")

print("\n========================================================")
print(" [DONE] Inspection complete.")
