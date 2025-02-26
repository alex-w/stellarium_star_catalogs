import pathlib
import re
import warnings

import pandas as pd
import tqdm
from astropy.table import Table, vstack
from astroquery.gaia import Gaia
import astropy.units as u

from py.utils import custom_simbad

base_path = pathlib.Path("simbad_query_results")
hip_combined_path = base_path / "hip_combined.dat"
sao_combined_path = base_path / "sao_combined.dat"
hd_combined_path = base_path / "hd_combined.dat"
hr_combined_path = base_path / "hr_combined.dat"
# check if the combined file already exists, if not then raise
for combined_path in [hip_combined_path, sao_combined_path, hd_combined_path, hr_combined_path]:
    if not combined_path.exists():
        raise FileNotFoundError(f"{combined_path} does not exist. Please run simbad_query_hipsaohdhr.py first.")
hip_combined_table = Table.read(hip_combined_path, format="ascii")
sao_combined_table = Table.read(sao_combined_path, format="ascii")
hd_combined_table = Table.read(hd_combined_path, format="ascii")
hr_combined_table = Table.read(hr_combined_path, format="ascii")
hip_to_remove = [1902, 54948]

# create an empty dataframe for cross ID of integers
cross_id_df = pd.DataFrame(
    columns=["hip", "gaia_dr3", "component", "sao", "hd", "hr"], dtype="Int64"
)
cross_id_df["component"] = cross_id_df["component"].astype(pd.StringDtype())


def extract_gaia_number(text, strict_dr3=False):
    """
    This function will extract the Gaia DR number from the text. Try Gaia DR3 first, if not found try Gaia DR2.

    Parameters
    ----------
    text : str
        The text to extract the Gaia DR number from
    strict_dr3 : bool
        If True, only extract Gaia DR3 number. Otherwise, try to extract other Gaia DR number if Gaia DR3 is not found
    """
    if text == text:  # check if not NaN
        match = re.search(r"Gaia DR3\s+(\d+)", text)
        if not match and not strict_dr3:  # if not try to find Gaia DR2
            match = re.search(r"Gaia DR2\s+(\d+)", text)
            if match:
                # make sure the source_id is valid, if not then try to cone search the source_id
                if len(Gaia.launch_job(f"SELECT source FROM gaiadr3.gaia_source WHERE source_id={match.group(1)}").results) < 1:
                    # in case we are upgrading to DR4, this need to be upgraded to DR3
                    tmp = Gaia.cone_search(f"Gaia DR2 {match.group(1)}", radius=15 * u.arcsecond).results
                    if len(tmp) > 0:
                        match = re.search(r"Gaia DR2\s+(\d+)", tmp[0]["designation"])
                    else:
                        print(f"Gaia DR2 {match.group(1)} not found in DR3")
                        return pd.NA
        if match and match.group(1) != "":
            if len(tmp := match.group(1)) > 2:
                return int(tmp)
            else:
                return tmp
    return pd.NA


def parse_result(simbad_table):
    """
    This function will parse the result from the query:
    - Remove rows that are not stars
    - Add additional rows to the table for stars that are potentially in a binary system

    To attempt to resolve potential binary system, the logic is as follows:
    1. If the stars have V and J magnitude dimmer than 3.0, then it should have Gaia source_id
    2. When no Gaia source_id is found, query the children of the star because it might be a binary system resolved by Gaia hence no Gaia source_id
    3. If childen have no Gaia source_id, then assume this is not a binary system. Just for whatever reason Gaia did not observe the parent star
    4. If the star has children with source_id, then add the children to the table and remove the parent star

    Parameters
    ----------
    simbad_table : astropy.table.Table
        Table from the query result

    Returns
    -------
    astropy.table.Table
        Table with additional rows from the query result
    """
    simbad_table = simbad_table.to_pandas()
    simbad_table["otype"] = simbad_table["otype"].fillna("")
    simbad_table = simbad_table[simbad_table["otype"].apply(lambda x: "err" not in x)]

    # more strict on resolving the source_id in case of binary stars
    source_id = simbad_table["ids"].apply(lambda x: extract_gaia_number(x))
    # these stars probably should have source_id, so we want to resolve possible binary stars
    hip_wo_source_id_idx = (pd.isna(source_id) & (simbad_table["V"].fillna(99.99) > 3.0) & (simbad_table["J"].fillna(99.99) > 3.0))
    simbad_table_temp = simbad_table[hip_wo_source_id_idx]
    result = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in tqdm.tqdm(simbad_table_temp.iterrows(), total=len(simbad_table_temp)):
            temp = custom_simbad.query_hierarchy(i[1]["main_id"], hierarchy="children")
            temp = temp.to_pandas()
            temp["user_specified_id"] = i[1]["user_specified_id"].strip()
            temp["source_id"] = temp["ids"].apply(lambda x: extract_gaia_number(x))
            # add letter to the end of the user_specified_id, so first row in temp will end with "A", second row in temp will end with "B", etc
            for idx, row in temp.iterrows():
                row["user_specified_id"] = f"{row['user_specified_id']}{chr(65 + idx)}"
                temp.loc[idx, "user_specified_id"] = row["user_specified_id"]
                # add user_specified_id to the table ids column
                temp.loc[idx, "ids"] = f"{row['user_specified_id']}|{row['ids']}"
            if (
                len(temp) > 1
            ):  # in case some stars only have one component which is itself
                # the second row and beyond should have source_id, if not then assume it is not a binary system
                if pd.isna(temp["source_id"].iloc[1:]).any():
                    continue
                # add temp_source_id for each row to the ids column if it is not NaN
                temp["ids"] = temp.apply(lambda x: f"{x['ids']}|Gaia DR3 {x['source_id']}" if not pd.isna(x["source_id"]) else x["ids"], axis=1)
                temp.drop(columns=["source_id"], inplace=True)
                result.append(Table.from_pandas(temp))
                simbad_table.drop(i[0], inplace=True)
    if len(result) > 0:
        result = vstack(result, metadata_conflicts="silent")
        simbad_table = pd.concat([simbad_table, result.to_pandas()])
    # turn back to astropy table
    simbad_table = Table.from_pandas(simbad_table)
    return simbad_table

def parse_cross_id(result):
    """
    Parse the cross id from the query result

    Parameters
    ----------
    result : pandas.DataFrame
        The result from the query
    """
    str_id = ["HIP", "SAO", "HD", "HR"]

    def extract_id(text, sid, group_id=1):
        if text == text:  # check if not NaN
            match = re.search(rf"{sid}\s+(\d+)([A-Z]?)", text)
            if match and match.group(group_id) != "":
                return match.group(group_id).strip()
        return pd.NA

    def parse_component(text):
        if pd.isna(text) or text == "" or text == 0:
            return pd.NA
        else:
            return text
        
    df = pd.DataFrame(
    columns=["hip", "gaia_dr3", "component", "sao", "hd", "hr"], dtype="Int64"
    )
    df["component"] = df["component"].astype(pd.StringDtype())
    hip, sao, hd, hr = [result["ids"].apply(lambda x: extract_id(x, sid)).astype("Int64") for sid in str_id]
    gaia = result["ids"].apply(lambda x: extract_gaia_number(x)).astype("Int64")
    component = result["ids"].apply(lambda x: parse_component(extract_id(x, "HIP", group_id=2))).astype(pd.StringDtype())

    # fill the dataframe
    df["hip"] = hip
    df["sao"] = sao
    df["hd"] = hd
    df["hr"] = hr
    df["gaia_dr3"] = gaia
    df["component"] = component

    return df


hip_with_binary = parse_result(hip_combined_table)
# deal with Sirius A and Sirius B, because SIMBAD said Sirius does not have children
hip_with_binary["ids"][hip_with_binary["main_id"] == "* alf CMa"] = "HIP 32349A|" + hip_with_binary[hip_with_binary["main_id"] == "* alf CMa"]["ids"]
result = custom_simbad.query_objects(["Sirius B"])
result["ids"] = "HIP 32349B|" + result["ids"]
# cast to pandas and back to prevent potential dtype incompatibility issue
result_df = result.to_pandas()
result = Table.from_pandas(result_df)
simbad_table = vstack([hip_with_binary, result])

cross_id_df = pd.concat([cross_id_df, parse_cross_id(simbad_table.to_pandas())])
# remove rows with just NaN
cross_id_df = cross_id_df.dropna(how="all")
# find which SAO id is missing, should be SAO 1 to SAO max_sao_id. Only query the missing SAO id
missing_sao_idx = [
    i for i in range(1, len(sao_combined_table) + 1) if i not in cross_id_df["sao"].values
]
# get the missing SAO id with corresponding row
missing_sao = sao_combined_table[[i - 1 for i in missing_sao_idx]]
cross_id_df = pd.concat([cross_id_df, parse_cross_id(missing_sao.to_pandas())])
# find which HD id is missing, should be HD 1 to HD max_hd_id. Only query the missing HD id
missing_hd_idx = [
    f"HD {i}" for i in range(1, len(hd_combined_table) + 1) if i not in cross_id_df["hd"].values
]
# get the missing HD id with corresponding row
missing_hd = custom_simbad.query_objects(missing_hd_idx)
cross_id_df = pd.concat([cross_id_df, parse_cross_id(missing_hd.to_pandas())])

cross_id_df.to_csv(base_path / "cross_id_with_gaiaid.dat", index=False, sep="\t")
# if hip is missing, then set hip to the corresponding gaia id and then delete the gaia column
cross_id_df.loc[pd.isna(cross_id_df["hip"]), "hip"] = cross_id_df.loc[pd.isna(cross_id_df["hip"]), "gaia_dr3"].values
cross_id_df = cross_id_df.drop(columns=["gaia_dr3"])
# if all sao, hd, hr is missing, then delete those rows. No point of keeping them
cross_id_df = cross_id_df[
    (~pd.isna(cross_id_df["sao"]) | ~pd.isna(cross_id_df["hd"]) | ~pd.isna(cross_id_df["hr"])) & ~pd.isna(cross_id_df["hip"])
]
cross_id_df = cross_id_df.sort_values(["hip", "component"]).astype({"component": "object"})
# replace 0 to empty string
# cross_id_df = cross_id_df.replace(0, "")
cross_id_df.to_csv(base_path / "cross_id.dat", index=False, sep="\t")
