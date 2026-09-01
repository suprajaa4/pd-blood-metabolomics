import pandas as pd
import numpy as np
import re
from collections import Counter



INPUT_FILE = "data/raw/pd_met_final.csv"
OUTPUT_FILE = "data/processed/metabolites_grouped_inchikey.xlsx"


def clean_string(x):
   
    if pd.isna(x):
        return np.nan

    x = str(x).strip()

    if x == "" or x.lower() in {"nan", "none", "na", "n/a"}:
        return np.nan

    # Collapse repeated spaces
    x = re.sub(r"\s+", " ", x)

    return x


def clean_inchikey(x):
    
    x = clean_string(x)

    if pd.isna(x):
        return np.nan

    return str(x).strip().upper()


def normalize_name_for_grouping(name):

    name = clean_string(name)

    if pd.isna(name):
        return ""

    x = name.lower()

    # Normalize unicode dashes
    x = x.replace("–", "-").replace("—", "-")

    # Replace punctuation with spaces
    x = re.sub(r"[^a-z0-9]+", " ", x)

    # Collapse spaces
    x = re.sub(r"\s+", " ", x).strip()

    return x


def pretty_name(name):

    name = clean_string(name)

    if pd.isna(name):
        return np.nan


    name = re.sub(r"\s+", " ", name).strip()

    replacements = {
        r"^l-": "L-",
        r"^d-": "D-",
        r"^dl-": "DL-",
    }

    for pattern, repl in replacements.items():
        name = re.sub(pattern, repl, name, flags=re.IGNORECASE)

    return name


def unique_join(series, sep=" | "):

    output = []
    seen = set()

    for value in series:
        value = clean_string(value)

        if pd.isna(value):
            continue

        key = str(value).lower()

        if key not in seen:
            output.append(str(value))
            seen.add(key)

    if not output:
        return np.nan

    return sep.join(output)


def unique_numeric_join(series, decimals=6):

    values = pd.to_numeric(series, errors="coerce").dropna()

    if len(values) == 0:
        return np.nan


    values = values.round(decimals)

    unique_values = pd.unique(values)

    if len(unique_values) == 1:
        return float(unique_values[0])

    return " | ".join(
        f"{v:g}" for v in unique_values
    )


def first_nonmissing(series):

    for x in series:
        if pd.notna(x):
            if str(x).strip() != "":
                return x

    return np.nan


def most_common_nonmissing(series):

    values = []

    for x in series:
        x = clean_string(x)

        if pd.notna(x):
            values.append(x)

    if not values:
        return np.nan

    counts = Counter(values)


    max_count = max(counts.values())

    tied = {
        value for value, count in counts.items()
        if count == max_count
    }

    for value in values:
        if value in tied:
            return value

    return values[0]


def choose_canonical_name(series):

    names = [
        pretty_name(x)
        for x in series
        if pd.notna(clean_string(x))
    ]

    if not names:
        return np.nan

    grouped = {}

    for name in names:
        normalized = normalize_name_for_grouping(name)

        grouped.setdefault(normalized, []).append(name)

    # frequency 
    best_key = sorted(
        grouped.keys(),
        key=lambda k: (
            -len(grouped[k]),
            len(k),
            k
        )
    )[0]

    candidate_names = grouped[best_key]

    # Choose most common capitalization/display variant
    counts = Counter(candidate_names)
    max_count = max(counts.values())

    candidates = [
        x for x, c in counts.items()
        if c == max_count
    ]

    candidates.sort(key=lambda x: (len(x), x))

    return candidates[0]


def log2_safe(x):
    x = pd.to_numeric(x, errors="coerce")

    if pd.isna(x) or x <= 0:
        return np.nan

    return np.log2(x)


def direction_from_fc(fc):
    fc = pd.to_numeric(fc, errors="coerce")

    if pd.isna(fc):
        return np.nan

    if fc > 1:
        return "up"
    elif fc < 1:
        return "down"
    else:
        return "unchanged"




df = pd.read_csv(INPUT_FILE, low_memory=False)


unnamed = [
    c for c in df.columns
    if str(c).startswith("Unnamed:")
]

df = df.drop(columns=unnamed, errors="ignore")


text_columns = [
    "repository",
    "study",
    "phase",
    "mode",
    "metabolite",
    "plasma/serum",
    "PubChem CID",
    "KEGG ID",
    "HMDB ID",
    "Formula",
    "InChIKey",
    "SMILES"
]

for col in text_columns:
    if col in df.columns:
        df[col] = df[col].apply(clean_string)


df["InChIKey"] = df["InChIKey"].apply(clean_inchikey)
df["metabolite"] = df["metabolite"].apply(pretty_name)

# ids
df["study"] = (
    df["study"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
    .str.strip()
)

df.loc[
    df["study"].str.lower().isin(["nan", "none", ""]),
    "study"
] = np.nan




numeric_columns = [
    "m/z ratio",
    "RT(min)",
    "Samples control",
    "Intensity mean control",
    "Intensity median control",
    "Samples",
    "Intensity mean PD",
    "Median Intensity PD",
    "p-value",
    "fold-change"
]

for col in numeric_columns:
    if col in df.columns:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )


# Log2-intensities
df["log_control_mean"] = df["Intensity mean control"].apply(log2_safe)
df["log_control_median"] = df["Intensity median control"].apply(log2_safe)

df["log_PD_mean"] = df["Intensity mean PD"].apply(log2_safe)
df["log_PD_median"] = df["Median Intensity PD"].apply(log2_safe)


# fold change where missing
fc_missing = df["fold-change"].isna()

valid_intensity = (
    df["Intensity mean PD"].notna()
    & df["Intensity mean control"].notna()
    & (df["Intensity mean control"] > 0)
)

df.loc[
    fc_missing & valid_intensity,
    "fold-change"
] = (
    df.loc[fc_missing & valid_intensity, "Intensity mean PD"]
    /
    df.loc[fc_missing & valid_intensity, "Intensity mean control"]
)


# Log2 fold change
df["Log2FC"] = df["fold-change"].apply(log2_safe)

# Direction
df["direction"] = df["fold-change"].apply(direction_from_fc)



df["z_control"] = np.nan
df["z_PD"] = np.nan

for study, idx in df.groupby("study").groups.items():

    control = df.loc[idx, "log_control_mean"]

    if control.notna().sum() >= 2 and control.std(skipna=True) > 0:
        df.loc[idx, "z_control"] = (
            control - control.mean(skipna=True)
        ) / control.std(skipna=True)

    pd_values = df.loc[idx, "log_PD_mean"]

    if pd_values.notna().sum() >= 2 and pd_values.std(skipna=True) > 0:
        df.loc[idx, "z_PD"] = (
            pd_values - pd_values.mean(skipna=True)
        ) / pd_values.std(skipna=True)


df["z_PD_minus_control"] = (
    df["z_PD"] - df["z_control"]
)



df["_clean_name_key"] = df["metabolite"].apply(
    normalize_name_for_grouping
)


def make_group_key(row):


    inchikey = row["InChIKey"]
    name_key = row["_clean_name_key"]
    metabolite = clean_string(row["metabolite"])

    if pd.notna(inchikey):
        return f"IK:{inchikey}"

    if (
        pd.isna(metabolite)
        or name_key == ""
        or name_key.startswith("unknown")
    ):
        return (
            f"UNKNOWN_ROW:"
            f"{row.name}"
        )


    return f"NAME:{name_key}"


df["_group_key"] = df.apply(
    make_group_key,
    axis=1
)



entity_rows = []


for group_key, group in df.groupby(
    "_group_key",
    sort=False
):

    canonical_name = choose_canonical_name(
        group["metabolite"]
    )

    all_names = unique_join(
        group["metabolite"]
    )

    inchikey = most_common_nonmissing(
        group["InChIKey"]
    )

    formula = most_common_nonmissing(
        group["Formula"]
    )

    smiles = most_common_nonmissing(
        group["SMILES"]
    )

    hmdb = most_common_nonmissing(
        group["HMDB ID"]
    )

    kegg = most_common_nonmissing(
        group["KEGG ID"]
    )

    pubchem = most_common_nonmissing(
        group["PubChem CID"]
    )




    if pd.notna(inchikey):

        entity_group = "Inchi"
        entity_id = f"IK_{inchikey}"

    elif group_key.startswith("UNKNOWN_ROW"):

        entity_group = "unknown"

        study = first_nonmissing(
            group["study"]
        )

        name = canonical_name

        if pd.isna(name):
            name = "unknown"

        entity_id = (
            f"UNKNOWN_ROW__"
            f"{study}__"
            f"{name}"
        )

    else:

        entity_group = "No Inchi"

        clean_key = normalize_name_for_grouping(
            canonical_name
        )

        entity_id = (
            "NO_INCHIKEY__"
            + clean_key
        )


    studies = [
        str(x)
        for x in group["study"].dropna().unique()
    ]




    has_formula = pd.notna(formula)
    has_smiles = pd.notna(smiles)
    has_hmdb = pd.notna(hmdb)
    has_kegg = pd.notna(kegg)
    has_pubchem = pd.notna(pubchem)

    annotation_score = sum([
        pd.notna(inchikey),
        has_formula,
        has_smiles,
        has_hmdb,
        has_kegg,
        has_pubchem
    ])


    row = {
        "entity_group": entity_group,
        "entity_id": entity_id,
        "canonical_metabolite": canonical_name,
        "all_metabolite_names": all_names,

        "n_input_rows": len(group),
        "n_studies": len(studies),

        "studies_present":
            " | ".join(studies)
            if studies else np.nan,

        "InChIKey": inchikey,
        "formula": formula,
        "SMILES": smiles,
        "HMDB ID": hmdb,
        "KEGG ID": kegg,
        "PubChem CID": pubchem,

        "annotation_score": annotation_score,

        "has_formula": has_formula,
        "has_SMILES": has_smiles,
        "has_HMDB": has_hmdb,
        "has_KEGG": has_kegg,
        "has_PubChem": has_pubchem
    }



    for study, study_data in group.groupby(
        "study",
        sort=False
    ):

        prefix = f"S_{study}__"


        #metadata
        row[prefix + "repository"] = unique_join(
            study_data["repository"]
        )

        row[prefix + "phase"] = unique_join(
            study_data["phase"]
        )

        row[prefix + "mode"] = unique_join(
            study_data["mode"]
        )

        row[prefix + "matrix"] = unique_join(
            study_data["plasma/serum"]
        )


        # Analytical measurements
        row[prefix + "mz"] = unique_numeric_join(
            study_data["m/z ratio"]
        )

        row[prefix + "RT_min"] = unique_numeric_join(
            study_data["RT(min)"]
        )


        # Sample counts
        row[prefix + "n_control"] = first_nonmissing(
            study_data["Samples control"]
        )

        row[prefix + "n_PD"] = first_nonmissing(
            study_data["Samples"]
        )


        # Intensities
        row[prefix + "int_mean_control"] = first_nonmissing(
            study_data["Intensity mean control"]
        )

        row[prefix + "int_median_control"] = first_nonmissing(
            study_data["Intensity median control"]
        )

        row[prefix + "int_mean_PD"] = first_nonmissing(
            study_data["Intensity mean PD"]
        )

        row[prefix + "int_median_PD"] = first_nonmissing(
            study_data["Median Intensity PD"]
        )


        # Logged intensities
        row[prefix + "log_control_mean"] = first_nonmissing(
            study_data["log_control_mean"]
        )

        row[prefix + "log_control_median"] = first_nonmissing(
            study_data["log_control_median"]
        )

        row[prefix + "log_PD_mean"] = first_nonmissing(
            study_data["log_PD_mean"]
        )

        row[prefix + "log_PD_median"] = first_nonmissing(
            study_data["log_PD_median"]
        )


        # Z-scores
        row[prefix + "z_control"] = first_nonmissing(
            study_data["z_control"]
        )

        row[prefix + "z_PD"] = first_nonmissing(
            study_data["z_PD"]
        )

        row[prefix + "z_PD_minus_control"] = first_nonmissing(
            study_data["z_PD_minus_control"]
        )


        # Statistical information
        row[prefix + "p_value"] = first_nonmissing(
            study_data["p-value"]
        )

        row[prefix + "fold_change"] = first_nonmissing(
            study_data["fold-change"]
        )

        row[prefix + "Log2FC"] = first_nonmissing(
            study_data["Log2FC"]
        )

        row[prefix + "direction"] = first_nonmissing(
            study_data["direction"]
        )


    entity_rows.append(row)




grouped = pd.DataFrame(entity_rows)




group_order = {
    "Inchi": 0,
    "No Inchi": 1,
    "unknown": 2
}

grouped["_order"] = grouped["entity_group"].map(
    group_order
)

grouped = grouped.sort_values(
    [
        "_order",
        "canonical_metabolite"
    ],
    na_position="last"
)

grouped = grouped.drop(
    columns="_order"
).reset_index(drop=True)




core_columns = [
    "entity_group",
    "entity_id",
    "canonical_metabolite",
    "all_metabolite_names",
    "n_input_rows",
    "n_studies",
    "studies_present",
    "InChIKey",
    "formula",
    "SMILES",
    "HMDB ID",
    "KEGG ID",
    "PubChem CID",
    "annotation_score",
    "has_formula",
    "has_SMILES",
    "has_HMDB",
    "has_KEGG",
    "has_PubChem"
]


study_variable_order = [
    "repository",
    "phase",
    "mode",
    "matrix",
    "mz",
    "RT_min",

    "n_control",
    "int_mean_control",
    "int_median_control",
    "log_control_mean",
    "log_control_median",
    "z_control",

    "n_PD",
    "int_mean_PD",
    "int_median_PD",
    "log_PD_mean",
    "log_PD_median",
    "z_PD",

    "z_PD_minus_control",
    "p_value",
    "fold_change",
    "Log2FC",
    "direction"
]


studies = sorted(
    [
        str(x)
        for x in df["study"].dropna().unique()
    ]
)


study_columns = []

for study in studies:

    for variable in study_variable_order:

        column = f"S_{study}__{variable}"

        if column in grouped.columns:
            study_columns.append(column)


other_columns = [
    c for c in grouped.columns
    if c not in core_columns
    and c not in study_columns
]


grouped = grouped[
    core_columns
    + study_columns
    + other_columns
]




summary = pd.DataFrame({
    "metric": [
        "input_rows",
        "output_entity_rows",
        "input_rows_with_InChIKey",
        "unique_input_InChIKeys",
        "grouped_InChIKey_entities",
        "no_InChIKey_entities",
        "unknown_entities",
        "number_of_studies"
    ],

    "value": [
        len(df),
        len(grouped),

        df["InChIKey"].notna().sum(),
        df["InChIKey"].dropna().nunique(),

        (grouped["entity_group"] == "Inchi").sum(),
        (grouped["entity_group"] == "No Inchi").sum(),
        (grouped["entity_group"] == "unknown").sum(),

        df["study"].nunique()
    ]
})




with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:

    grouped.to_excel(
        writer,
        sheet_name="Updated_Grouped",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )


print()
print("DONE!")
print(f"Input rows:        {len(df):,}")
print(f"Grouped entities:  {len(grouped):,}")
print(
    "InChIKey groups:   ",
    (grouped["entity_group"] == "Inchi").sum()
)
print(
    "No-InChI groups:   ",
    (grouped["entity_group"] == "No Inchi").sum()
)
print(
    "Unknown rows:      ",
    (grouped["entity_group"] == "unknown").sum()
)
print()
print(f"saved {OUTPUT_FILE}")
