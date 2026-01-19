"""
Demographic characteristics analysis for clinical data.
Generates summary tables of patient demographics by diagnosis group.
"""

import pandas as pd
import numpy as np
from scipy import stats
from typing import Dict, Tuple


def generate_demographic_table(
    metadata_df: pd.DataFrame,
    output_file: str = None
) -> pd.DataFrame:
    """
    Generate a demographic characteristics table by diagnosis group.
    
    This function creates a summary table of demographic characteristics
    (age, sex, BMI, waist circumference) for each diagnosis group,
    including statistical tests for differences between groups.
    
    Args:
        metadata_df (pd.DataFrame): DataFrame with patient metadata including:
            - 'patient_id': unique patient identifier
            - 'Diagnosis': diagnosis classification
            - 'Age[A2]': patient age (numeric)
            - 'Gender[A1]': patient sex (F or M)
            - 'BMI': body mass index (numeric)
            - 'Height[cm][G1]': height in cm (for BMI calculation if needed)
            - 'Mass[kg][G4]': mass in kg (for BMI calculation if needed)
            - 'Waist Circumference[cm][G7]': waist circumference measurement in cm (numeric)
        output_file (str, optional): Path to save the output table as CSV
    
    Returns:
        pd.DataFrame: Demographic characteristics table with diagnosis groups as columns
    """
    
    # Get one row per patient (in case there are multiple samples per patient)
    patient_df = metadata_df.groupby('patient_id').first().reset_index()
    
    # Convert numeric columns to proper types
    # Replace '?' with NaN first
    patient_df = patient_df.replace('?', np.nan)
    
    numeric_columns = ['Age[A2]', 'BMI', 'Height[cm][G1]', 'Mass[kg][G4]', 'Waist Circumference[cm][G7]', 'Waist Circumfurence(cm)[G7]']
    for col in numeric_columns:
        if col in patient_df.columns:
            patient_df[col] = pd.to_numeric(patient_df[col], errors='coerce')
    
    # Rename diagnoses to match publication format
    patient_df['Diagnosis'] = patient_df['Diagnosis'].replace({
        'Healthy': 'Controls',
        'Metastatic': 'Metastatic OAC'
    })
    
    # Define the order of diagnosis columns
    diagnosis_order = ['Controls', 'GORD', 'BO', 'Dysplasia', 'OAC', 'Metastatic OAC']
    
    # Initialize results dictionary
    results = {}
    
    # 1. Patient counts (N)
    patient_counts = patient_df['Diagnosis'].value_counts()
    results['Patients (N)'] = {diag: int(patient_counts.get(diag, 0)) for diag in diagnosis_order}
    
    # 2. Age (mean, range)
    age_stats = {}
    age_values_for_test = []
    age_group_names = []
    
    # Try different variations of age column name
    age_col = None
    possible_age_cols = ['Age[A2]', 'Age', 'age']
    for col_name in possible_age_cols:
        if col_name in patient_df.columns:
            age_col = col_name
            break
    
    if age_col:
        for diag in diagnosis_order:
            subset = patient_df[patient_df['Diagnosis'] == diag][age_col].dropna()
            if len(subset) > 0:
                mean_age = subset.mean()
                min_age = subset.min()
                max_age = subset.max()
                age_stats[diag] = f"{mean_age:.1f} ({int(min_age)}-{int(max_age)})"
                age_values_for_test.append(subset.values)
                age_group_names.append(diag)
            else:
                age_stats[diag] = "N/A"
    else:
        for diag in diagnosis_order:
            age_stats[diag] = "N/A"
    
    # Perform Kruskal-Wallis test for age
    if len(age_values_for_test) > 1:
        _, age_pval = stats.kruskal(*age_values_for_test)
        age_stats['p value'] = f"{age_pval:.3f}"
    else:
        age_stats['p value'] = "N/A"
    
    results['Age (mean,range)'] = age_stats
    
    # 3. Sex (f/m)
    sex_stats = {}
    sex_contingency = []
    sex_group_names = []
    
    # Try different variations of gender column name
    gender_col = None
    possible_gender_cols = ['Gender[A1]', 'Sex', 'sex', 'Gender', 'gender']
    for col_name in possible_gender_cols:
        if col_name in patient_df.columns:
            gender_col = col_name
            break
    
    if gender_col:
        for diag in diagnosis_order:
            subset = patient_df[patient_df['Diagnosis'] == diag]
            if len(subset) > 0:
                f_count = ((subset[gender_col] == 'F') | (subset[gender_col] == 'Female')).sum()
                m_count = ((subset[gender_col] == 'M') | (subset[gender_col] == 'Male')).sum()
                sex_stats[diag] = f"{int(f_count)}/{int(m_count)}"
                sex_contingency.append([f_count, m_count])
                sex_group_names.append(diag)
            else:
                sex_stats[diag] = "0/0"
                sex_contingency.append([0, 0])
    else:
        for diag in diagnosis_order:
            sex_stats[diag] = "N/A"
            sex_contingency.append([0, 0])
    
    # Chi-square test for sex
    if len(sex_contingency) > 1:
        chi2, sex_pval, _, _ = stats.chi2_contingency(sex_contingency)
        sex_stats['p value'] = f"{sex_pval:.3f}"
    else:
        sex_stats['p value'] = "N/A"
    
    results['Sex (f/m)'] = sex_stats
    
    # 4. BMI (mean, range)
    bmi_stats = {}
    bmi_values_for_test = []
    bmi_group_names = []
    
    # Check if BMI column exists, if not calculate from Height and Mass
    if 'BMI' in patient_df.columns:
        bmi_col = 'BMI'
    else:
        # Calculate BMI from Height[cm][G1] and Mass[kg][G4]
        patient_df['BMI'] = (patient_df['Mass[kg][G4]'] / ((patient_df['Height[cm][G1]'] / 100) ** 2)).round(1)
        bmi_col = 'BMI'
    
    for diag in diagnosis_order:
        subset = patient_df[patient_df['Diagnosis'] == diag][bmi_col].dropna()
        if len(subset) > 0:
            mean_bmi = subset.mean()
            min_bmi = subset.min()
            max_bmi = subset.max()
            bmi_stats[diag] = f"{mean_bmi:.1f} ({min_bmi:.1f}-{max_bmi:.1f})"
            bmi_values_for_test.append(subset.values)
            bmi_group_names.append(diag)
        else:
            bmi_stats[diag] = "N/A"
    
    # Kruskal-Wallis test for BMI
    if len(bmi_values_for_test) > 1:
        _, bmi_pval = stats.kruskal(*bmi_values_for_test)
        bmi_stats['p value'] = f"{bmi_pval:.3f}"
    else:
        bmi_stats['p value'] = "N/A"
    
    results['BMI (mean,range)'] = bmi_stats
    
    # 5. Waist Circumference (mean, range)
    wc_stats = {}
    wc_values_for_test = []
    wc_group_names = []
    
    # Try different variations of the waist circumference column name (note: typo in original "Circumfurence")
    wc_col = None
    possible_wc_cols = ['Waist Circumfurence(cm)[G7]', 'Waist Circumference[cm][G7]', 'Waist_Circumference', 'Waist Circumference', 'Waist_Circumference[cm][G7]']
    for col_name in possible_wc_cols:
        if col_name in patient_df.columns:
            wc_col = col_name
            break
    
    if wc_col:
        for diag in diagnosis_order:
            subset = patient_df[patient_df['Diagnosis'] == diag][wc_col].dropna()
            if len(subset) > 0:
                mean_wc = subset.mean()
                min_wc = subset.min()
                max_wc = subset.max()
                wc_stats[diag] = f"{mean_wc:.1f} ({int(min_wc)}-{int(max_wc)})"
                wc_values_for_test.append(subset.values)
                wc_group_names.append(diag)
            else:
                wc_stats[diag] = "N/A"
    else:
        print("Warning: Waist Circumference column not found. Available columns:")
        print(patient_df.columns.tolist())
        for diag in diagnosis_order:
            wc_stats[diag] = "N/A"
    
    # Kruskal-Wallis test for Waist Circumference
    if len(wc_values_for_test) > 1:
        _, wc_pval = stats.kruskal(*wc_values_for_test)
        wc_stats['p value'] = f"{wc_pval:.3f}"
    else:
        wc_stats['p value'] = "N/A"
    
    results['Waist Circumference (mean,range)'] = wc_stats
    
    # Create DataFrame from results
    table_df = pd.DataFrame(results).T
    
    # Reorder columns
    columns_order = diagnosis_order + ['p value']
    table_df = table_df[columns_order]
    
    # Save to CSV if output file is provided
    if output_file:
        table_df.to_csv(output_file)
        print(f"Demographic characteristics table saved to {output_file}")
    
    return table_df

